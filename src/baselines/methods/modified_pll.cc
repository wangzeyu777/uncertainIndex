#include "modified_pll.h"

MPLL::MPLL(Graph &graph, IndexOrder order): order(order), graph(graph), label_backward(NULL), label_forward(NULL) {}

MPLL::~MPLL() {
    delete_index();
}

void MPLL::delete_index() {
    delete []label_forward;
    delete []label_backward;
}

void MPLL::index() {
    delete_index();
    const int v_n = graph.vertex_number();
    label_forward = new MLabels[v_n];
    label_backward = new MLabels[v_n];

    DijkstraBreadthFirstSearch dijk_search(graph);
    Vertex *v = NULL;
    if (order == BY_DEGREE) {
        v = new Vertex[v_n];
        auto v_sort = new std::pair<int, int>[v_n];
        for (int i = 0; i < v_n; ++i) {
            VertexLoad load = graph.get_vertex_load({i});
            v_sort[i] = {- load.in_degree - load.out_degree, i};
        }
        std::sort(v_sort, v_sort+v_n);
        for (int i = 0; i < v_n; ++i) {
            // logger_info << v_sort[i].first << " " << v_sort[i].second << " id " << graph.id({v_sort[i].second}) << std::endl; 
            v[i] = {v_sort[i].second};
        }
        delete []v_sort;
    }

    long long total_size = 0;
    for (int rank = 0; rank < v_n; ++rank) {
        int root;
        if (order == BY_DEGREE) {
            root = v[rank].v_pos;
            logger_info << "round " << rank << " select root " << graph.id({root}) << std::endl;
        } else if (order == BY_INPUT_ORDER) {
            root = {rank};
        }
        Distance d;
        label_backward[root].emplace_back(root, -1, 0, rank, 1);
        label_forward[root].emplace_back(root, -1, 0, rank, 1);
        for (int direction = 0; direction < 2; ++direction) {
            logger_info << "round " << rank << " index " << root << " " << (direction ? "backward" : "forward") << std::endl;
            dijk_search.reset({{root}}, direction);
            dijk_search.next();
            MLabels *labels = direction ? label_forward : label_backward;
            while (graph.valid((d = dijk_search.front()).node)) {
                double queried = query({root}, d.node, direction);
                logger_info << " " << graph.id(d.node) << " (" << graph.id({root}) << " " << d.distance << ")" << " queried " << queried << std::endl;
                if (d.distance < queried + EPS) {
                    dijk_search.next(true);
                } else {
                    Vertex p = d.parent;
                    labels[d.node.v_pos].emplace_back(root, p.v_pos, labels[p.v_pos].size()-1, rank, d.distance);
                    logger_info << "push node " << graph.id(d.node) << " parent " << graph.id(p) << " p_pos "<< labels[p.v_pos].size()-1 << " dis " << d.distance << std::endl;
                    dijk_search.next();
                    total_size += 1;
                }
            }
            // logger_info << "index " << (direction ? "forward" : "backward") << std::endl;
        }
        // logger_info << "total size " << total_size <<  " avg size " << total_size / (double) (rank+1) << std::endl;
    }
    logger_info << "total size " << total_size <<  " avg size " << total_size / (double) graph.vertex_number() << std::endl;

    delete []v;
}

double MPLL::query(const DistanceQuery &q) {
    return MPLL::query(q.s, q.t);
}

MLabel::MLabel() {}
MLabel::MLabel(int node, int parent, int parent_pos, int rank, double distance): node(node),
    parent(parent), parent_pos(parent_pos), rank(rank), distance(distance) {}

double MPLL::path_length(MLabel *s_label, Vertex s, MLabel *t_label, Vertex t, int direction) {
    // logger_info << "query path length " << graph.id(s) << " -> " << graph.id(t) << " direction " << direction << std::endl;
    double res = 1;
    if (direction) {
        // backward
        while (t_label->parent != -1 && t_label != s_label) {
            int find = 0;
            for (auto e=graph.first_in(t); graph.valid(e); e=graph.next_in(e)) {
                if (graph.get_source(e).v_pos == t_label->parent) {
                    res *= graph.get_edge_load(e);
                    t = {t_label->parent};
                    t_label = &label_backward[t_label->parent][t_label->parent_pos];
                    find = 1;
                    break;
                }
            }
            if (!find) {
                logger_error << "backward parent of " << graph.id({t_label->node}) << " not found" << std::endl;
                return 0;
            }
        }
        if (t_label != s_label) return 0;
    } else {
        // forward
        while (s_label->parent != -1 && t_label != s_label) {
            int find = 0;
            for (auto e=graph.first_out(s); graph.valid(e); e=graph.next_out(e)) {
                if (graph.get_target(e).v_pos == s_label->parent) {
                    res *= graph.get_edge_load(e);
                    s = {s_label->parent};
                    s_label = &label_forward[s_label->parent][s_label->parent_pos];
                    find = 1;
                    break;
                }
            }
            if (!find) {
                logger_error << "forward parent of " << graph.id({s_label->node}) << " not found" << std::endl;
                return 0;
            }
        }
        if (t_label != s_label) return 0;
    }
    return res;
}

double MPLL::query(Vertex s, Vertex t, int direction) {
    if (direction) std::swap(s, t);
    MLabels & forward_s = label_forward[s.v_pos];
    MLabels & forward_t = label_forward[t.v_pos];
    MLabels & backward_t = label_backward[t.v_pos];
    MLabels & backward_s = label_backward[s.v_pos];

    auto i = forward_s.begin(), j = backward_t.begin(), k = forward_t.begin();
    double res = 0;
    for (; i != forward_s.end() && j != backward_t.end(); ++i) {
        while (j != backward_t.end() && j->rank < i->rank) {
            ++j;
        }
        if (j != backward_t.end() && j->rank == i->rank) {
            res = std::max(res, j->distance * i->distance);
        }
        while (k != forward_t.end() && k->rank < i->rank) {
            ++k;
        }
        if (k != forward_t.end() && k->rank == i->rank) {
            res = std::max(res, path_length(&*i, s, &*k, t, 0));
        }
    }
    i = backward_s.begin(), j = backward_t.begin();
    for (; i != backward_s.end(); ++i) {
        while (j != backward_t.end() && j->rank < i->rank) {
            ++j;
        }
        if (j != backward_t.end() && j->rank == i->rank) {
            res = std::max(res, path_length(&*i, s, &*j, t, 1));
        }
    }
    return res;
}

void MPLL::load_index(const char *filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        logger_error << "Cannot Open File " << filename << std::endl;
    }
    delete_index();
    label_forward = new MLabels[graph.vertex_number()];
    label_backward = new MLabels[graph.vertex_number()];

    int i, s;
    char c;
    MLabel label;
    MLabels *labels;
    while (file.peek() != EOF) {
        file.read((char *)&i, sizeof(i));
        file.read(&c, sizeof(c));
        if (c == 'f') {
            labels = label_forward;
        } else if (c == 'b') {
            labels = label_backward;
        } else {
            // TODO: wrong format
        }
        file.read((char *)&s, sizeof(s));
        MLabels &l = labels[i];
        l.reserve(s);
        while (s--) {
            file.read((char *)&label, sizeof(label));
            l.emplace_back(label);
        }
    }
    file.close();
}


void MPLL::save_index(const char* filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        logger_error << "Cannot Open File " << filename << std::endl;
        return ;
    }
    for (int i = 0; i < graph.vertex_number(); ++i) {
        auto x = label_forward[i];
        int s;
        if ((s = x.size()) > 0) {
            file.write((char *)&i, sizeof(i));
            file.write("f", 1);
            file.write((char *)&s, sizeof(s));
            file.write((char *)x.data(), s * sizeof(MLabel));
        }
        x = label_backward[i];
        if ((s = x.size()) > 0) {
            file.write((char *)&i, sizeof(i));
            file.write("b", 1);
            file.write((char *)&s, sizeof(s));
            file.write((char *)x.data(), s * sizeof(MLabel));
        }
    }
    file.close();
}
