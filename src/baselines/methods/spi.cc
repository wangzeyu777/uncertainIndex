#include "spi.h"
#include <chrono>

ReachabilityLabel::ReachabilityLabel() {}
ReachabilityLabel::ReachabilityLabel(int root, int rank, double distance, double reachability):
    Label(root, rank, distance), r(reachability) {}

ShortestPathIndexing::ShortestPathIndexing(Graph &graph, int sample_n, IndexOrder order, bool struct_opt, double t_r, int t_k): order(order),
graph(graph), struct_opt(struct_opt), label_backward(NULL), label_forward(NULL), sample_n(sample_n), threshold_hop(t_k), threshold_r(t_r) {}

ShortestPathIndexing::~ShortestPathIndexing() {
    delete_index();
}

void ShortestPathIndexing::delete_index() {
    delete []label_forward;
    delete []label_backward;
}

void ShortestPathIndexing::index() {
    delete_index();
    const int v_n = graph.vertex_number();
    label_forward = new RLabels[v_n];
    label_backward = new RLabels[v_n];
    int *counter = new int[v_n]();
    bool *valid = new bool[v_n]();
    int *vertices = new int[v_n+1]();
    int num = 0;

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
    } else if (order == PATH_COVER) {
        int k = log(v_n) / log(2);
        if (k < 1) k = 1;
        if (k >= v_n) k = v_n;
        ShortestPathTree::maximum_size = 10 * k * v_n;
        logger_info << "Number of nodes in trees: ~" << ShortestPathTree::maximum_size << std::endl;
        ShortestPathTree::add_trees(graph, dijk_search, true);
    }

    long long total_size = 0;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int rank = 0; rank < v_n; ++rank) {
        int root;
        if (order == BY_DEGREE) {
            root = v[rank].v_pos;
            logger_info << "round " << rank << " select root " << graph.id({root}) << std::endl;
        } else if (order == BY_INPUT_ORDER) {
            root = {rank};
        } else if (order == PATH_COVER) {
            root = ShortestPathTree::max_path_cover_node();
            if (root != -1) {
                logger_info << "round " << rank << " select root " << root << " path_n = " << ShortestPathTree::path_n[root]<< std::endl;
                ShortestPathTree::cut_all(root);
                ShortestPathTree::add_trees(graph, dijk_search);
            } else {
                // TODO: error
                logger_error << "SPI: No next root found while indexing." << std::endl;
            }
        }
        Distance d;
        for (int direction = 0; direction < 2; ++direction) {
            Graph *sub_graph = new Graph(0, 0, false, true);
            // logger_info << "round " << rank << " index " << root << " " << (direction ? "backward" : "forward") << std::endl;
            dijk_search.reset({{root}}, direction);
            dijk_search.next();
            RLabels *labels = direction ? label_forward : label_backward;
            while (graph.valid((d = dijk_search.front()).node)) {
                double queried = query_distance({root}, d.node, direction);
                // logger_info << " " << graph.id(d.node) << " (" << graph.id({root}) << " " << d.distance << ")" << " queried " << queried << std::endl;
                if (d.distance < queried + EPS || d.distance <= threshold_r - EPS || d.hop > threshold_hop ) {
                    dijk_search.next(true);
                } else {
                    labels[d.node.v_pos].emplace_back(root, rank, d.distance, 0);
                    valid[d.node.v_pos] = true;
                    vertices[num++] = d.node.v_pos;
                    // logger_info << "push node " << graph.id(d.node) << " dis " << d.distance << std::endl;
                    dijk_search.next();
                    total_size += 1;
                    if (total_size > 2ll * tree_bound * v_n) {
                        break;
                    }
                }
            }
            vertices[num] = -1;
            num = 0;
            labels[root].emplace_back(root, rank, 1, 1);
            // collect edge for sub graph
            auto first_edge = direction ? &Graph::first_in : &Graph::first_out;
            auto next = direction ? &Graph::next_in : &Graph::next_out;
            auto get_vertex = direction ? &Graph::get_source : &Graph::get_target;
            for (int *v = vertices; *v != -1; ++v) {
                for (auto e=(graph.*first_edge)({*v}); graph.valid(e); e=(graph.*next)(e)) {
                    int y = (graph.*get_vertex)(e).v_pos;
                    if (valid[y]) {
                        EdgeLoad el = graph.get_edge_load(e);
                        if (direction)
                            sub_graph->add_edge(y, *v, el);
                        else {
                            sub_graph->add_edge(*v, y, el);
                        }
                    }
                }
            }
            // estimate reachability
            naive_monte_carlo({root}, sample_n, graph, direction, valid, counter);
            // put reliability in labels
            for (int *v = vertices; *v != -1; ++v) {
                labels[*v].back().r = counter[*v] / (double) sample_n;
                // logger_info << graph.id({*v}) << " r=" <<counter[*v] / (double) sample_n << std::endl;
                // clear flags and counters
                valid[*v] = 0, counter[*v] = 0;
            }
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            logger_info << "Time Elapse: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " s" << std::endl;
            delete sub_graph;
        }

        logger_info << "total size " << total_size <<  " avg size " << total_size / (double) (rank+1) << std::endl;
    }
    logger_info << "total size " << total_size <<  " avg size " << total_size / (double) graph.vertex_number() << std::endl;

    delete []v;
}

double ShortestPathIndexing::query(const DistanceQuery &q) {
    return struct_opt ? ShortestPathIndexing::query_reachability_struct_opt(q.s, q.t) : ShortestPathIndexing::query_reachability(q.s, q.t);
}

double ShortestPathIndexing::query_distance(Vertex s, Vertex t, int direction) {
    if (direction) std::swap(s, t);
    RLabels & forward_s = label_forward[s.v_pos];
    RLabels & backward_t = label_backward[t.v_pos];
    auto i = forward_s.begin(), j = backward_t.begin();
    double res = 0;
    for (; i != forward_s.end() && j != backward_t.end(); ++i) {
        while (j != backward_t.end() && j->rank < i->rank) {
            ++j;
        }
        if (j != backward_t.end() && j->node == i->node) {
            res = std::max(res, j->distance * i->distance);
        }
    }
    return res;
}

void print(const ReachabilityLabel &label) {
    logger_info << "root " << label.node << " dis " << label.distance << " r " << label.r << std::endl;
}

double ShortestPathIndexing::query_reachability(Vertex s, Vertex t) {
    // logger_info << "here" << std::endl;
    RLabels & forward_s = label_forward[s.v_pos];
    RLabels & backward_t = label_backward[t.v_pos];
    auto i = forward_s.begin(), j = backward_t.begin();
    double res = 0;
    for (; i != forward_s.end() && j != backward_t.end(); ++i) {
        while (j != backward_t.end() && j->rank < i->rank) {
            ++j;
        }
        if (j != backward_t.end() && j->node == i->node) {
            // logger_info << graph.id(s) << " forward "; 
            // print(*i);
            // logger_info << graph.id(t) << " backward "; 
            // print(*j);
            
            res = std::max(res, j->r * i->r);
        }
    }
    return res;
}

typedef std::vector<int> EdgeSet;
bool intersect(const EdgeSet* a, const EdgeSet *b) {
    auto j = b->begin();
    for (auto i: *a) {
        while (*j < i) {
            j++;
        }
        if (*j == i) {
            return true;
        }
    }
    return false;
}

double maximum_clique_dfs(std::vector<double> & values, std::vector<EdgeSet *> & nodes, std::vector<int> & selected_nodes, int i) {
    if (i == values.size()) {
        double res = 1;
        for (auto id : selected_nodes) {
            res *= 1-values[id];
        }
        return 1-res;
    }
    double res = maximum_clique_dfs(values, nodes, selected_nodes, i+1);
    for (auto j: selected_nodes) {
        if (intersect(nodes[i], nodes[j])) {
            return res;
        }
    }
    selected_nodes.push_back(i);
    res = std::max(res, maximum_clique_dfs(values, nodes, selected_nodes, i+1));
    selected_nodes.pop_back();
    return res;
}
double maximum_clique(std::vector<double> & values, std::vector<EdgeSet *> & nodes) {
    std::vector<int> selected;
    return maximum_clique_dfs(values, nodes, selected, 0);
}

double ShortestPathIndexing::query_reachability_struct_opt(Vertex s, Vertex t) {
    RLabels & forward_s = label_forward[s.v_pos];
    RLabels & backward_t = label_backward[t.v_pos];
    auto i = forward_s.begin(), j = backward_t.begin();

    std::vector<double> reachabilities;
    std::vector<EdgeSet *> edge_sets;
    // bfs
    static Vertex *q = new Vertex[graph.vertex_number()];
    static int *visited = new int[graph.vertex_number()]();
    static int visit_flag = 0;

    for (; i != forward_s.end() && j != backward_t.end(); ++i) {
        while (j != backward_t.end() && j->rank < i->rank) {
            ++j;
        }
        if (j != backward_t.end() && j->node == i->node) {
            // res = std::max(res, j->r * i->r);
            reachabilities.push_back(j->r * i->r);
            EdgeSet *es = new EdgeSet(); // TODO: delete
            // use bfs to calculate edge set
            for (int direction = 0; direction < 2; ++direction) {
                // prepare bfs
                if (visit_flag >= BreadthFirstSearch::MAX_FLAG) {
                    memset(visited, 0, sizeof(int) * graph.vertex_number());
                    visit_flag = 0;
                }
                visit_flag++;
                int op = 0, cl = 0;
                visited[(q[cl++] = direction ? t : s).v_pos] = visit_flag;

                auto first = direction ? &Graph::first_in : &Graph::first_out;
                auto next = direction ? &Graph::next_in : &Graph::next_out;
                auto get_node = direction ? &Graph:: get_source : &Graph::get_target;
                auto labels = direction ? label_backward : label_forward;
                while (op < cl) {
                    Vertex v = q[op++], u;
                    // printf("node %d\n", graph.id(v));
                    for (auto e = (graph.*first)(v); graph.valid(e); e = (graph.*next)(e)) {
                        u = (graph.*get_node)(e);
                        // printf("edge %d->%d\n", graph.id(v), graph.id(u));
                        // printf("u %d\n", u.v_pos);
                        if (exist_rank(labels[u.v_pos], i->rank)) {
                            // printf("u %d\n", u.v_pos);
                            es->push_back(e.e_pos);
                            // printf("u %d\n", u.v_pos);
                            if (visited[u.v_pos] != visit_flag) {
                                visited[u.v_pos] = visit_flag;
                                // printf("cl %d v_n %d\n", cl, graph.vertex_number());
                                q[cl++] = u;
                            }
                        }
                    }
                }
            }
            std::sort(es->begin(), es->end());
            edge_sets.push_back(es);
        }
    }

    return maximum_clique(reachabilities, edge_sets);
}

void ShortestPathIndexing::load_index(const char *filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        logger_error << "Cannot Open File " << filename << std::endl;
    }
    delete_index();
    label_forward = new RLabels[graph.vertex_number()];
    label_backward = new RLabels[graph.vertex_number()];

    int i, s;
    char c;
    ReachabilityLabel label;
    RLabels *labels;
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
        RLabels &l = labels[i];
        // logger_info << i << " " << c << " size " << s << std::endl;
        l.reserve(s);
        while (s--) {
            file.read((char *)&label, sizeof(label));
            // logger_info << label.node << " " << label.rank << " " << label.distance << " "<< label.r << std::endl;
            l.emplace_back(label);
        }
    }
    file.close();
}


void ShortestPathIndexing::save_index(const char* filename) {
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
            file.write((char *)x.data(), s * sizeof(ReachabilityLabel));
        }
        x = label_backward[i];
        if ((s = x.size()) > 0) {
            file.write((char *)&i, sizeof(i));
            file.write("b", 1);
            file.write((char *)&s, sizeof(s));
            file.write((char *)x.data(), s * sizeof(ReachabilityLabel));
        }
    }
    file.close();
}

bool operator < (const ReachabilityLabel &a, const ReachabilityLabel &b) {
    return a.rank < b.rank;
}

bool exist_rank(const RLabels &labels, int rank) {
    ReachabilityLabel l(0, rank, 0, 0);
    auto it = std::lower_bound(labels.begin(), labels.end(), l);
    return it != labels.end() && it->rank == rank;
}

void ShortestPathIndexing::set_tree_bound(int bound) {
    tree_bound = bound;
}