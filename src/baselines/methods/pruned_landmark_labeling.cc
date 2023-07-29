#include "pruned_landmark_labeling.h"

// bool degree_cmp(const Vertex *a, const Vertex *b) {
//     return a && (!b || get_degree() > b->get_degree());
// }

PrunedLandmarkLabeling::PrunedLandmarkLabeling(Graph &graph, IndexOrder order): order(order), graph(graph), label_backward(NULL), label_forward(NULL) {}

PrunedLandmarkLabeling::~PrunedLandmarkLabeling() {
    delete_index();
}

void PrunedLandmarkLabeling::delete_index() {
    delete []label_forward;
    delete []label_backward;
}

void PrunedLandmarkLabeling::index() {
    delete_index();
    const int v_n = graph.vertex_number();
    label_forward = new Labels[v_n];
    label_backward = new Labels[v_n];

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
                logger_error << "PrunedLandmarkLabaling: No next root found while indexing." << std::endl;
            }
        }
        Distance d;
        for (int direction = 0; direction < 2; ++direction) {
            // logger_info << "round " << rank << " index " << root << " " << (direction ? "backward" : "forward") << std::endl;
            dijk_search.reset({{root}}, direction);
            dijk_search.next();
            Labels *labels = direction ? label_forward : label_backward;
            while (graph.valid((d = dijk_search.front()).node)) {
                double queried = query({root}, d.node, direction);
                // logger_info << " " << graph.id(d.node) << " (" << graph.id({root}) << " " << d.distance << ")" << " queried " << queried << std::endl;
                if (d.distance < queried + EPS) {
                    dijk_search.next(true);
                } else {
                    labels[d.node.v_pos].emplace_back(root, rank, d.distance);
                    // logger_info << "push node " << graph.id(d.node) << " dis " << d.distance << std::endl;
                    dijk_search.next();
                    total_size += 1;
                }
            }
            // logger_info << "index " << (direction ? "forward" : "backward") << std::endl;
        }
        label_backward[root].emplace_back(root, rank, 1);
        label_forward[root].emplace_back(root, rank, 1);
        // logger_info << "total size " << total_size <<  " avg size " << total_size / (double) (rank+1) << std::endl;
    }
    logger_info << "total size " << total_size <<  " avg size " << total_size / (double) graph.vertex_number() << std::endl;

    delete []v;
}

double PrunedLandmarkLabeling::query(const DistanceQuery &q) {
    return PrunedLandmarkLabeling::query(q.s, q.t);
}

double PrunedLandmarkLabeling::query(Vertex s, Vertex t, int direction) {
    if (direction) std::swap(s, t);
    Labels & forward_s = label_forward[s.v_pos];
    Labels & backward_t = label_backward[t.v_pos];
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

void PrunedLandmarkLabeling::load_index(const char *filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        logger_error << "Cannot Open File " << filename << std::endl;
    }
    delete_index();
    label_forward = new Labels[graph.vertex_number()];
    label_backward = new Labels[graph.vertex_number()];

    int i, s;
    char c;
    Label label;
    Labels *labels;
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
        Labels &l = labels[i];
        // logger_info << i << " " << c << " size " << s << std::endl;
        l.reserve(s);
        while (s--) {
            file.read((char *)&label, sizeof(label));
            // logger_info << label.node << " " << label.rank << " " << label.distance << std::endl;
            l.emplace_back(label);
        }
    }
    file.close();
}


void PrunedLandmarkLabeling::save_index(const char* filename) {
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
            file.write((char *)x.data(), s * sizeof(Label));
        }
        x = label_backward[i];
        if ((s = x.size()) > 0) {
            file.write((char *)&i, sizeof(i));
            file.write("b", 1);
            file.write((char *)&s, sizeof(s));
            file.write((char *)x.data(), s * sizeof(Label));
        }
    }
    file.close();
}

long long ShortestPathTree::maximum_size;
long long ShortestPathTree::total_size;
Fib<std::pair<long long, int>> ShortestPathTree::heap;
long long* ShortestPathTree::path_n = NULL;
bool* ShortestPathTree::cut = NULL;
std::vector<ShortestPathTree *> ShortestPathTree::trees;

int ShortestPathTree::new_tree(Graph &graph, DijkstraBreadthFirstSearch &dijk_search, int root) {
    logger_info << "new tree rooted at "<< root << std::endl;
    auto tree = new ShortestPathTree(graph, dijk_search); // TODO: delete this
    trees.push_back(tree);
    return tree->generate(root);
}

void ShortestPathTree::cut_all(int root) {
    int temp = 0;
    for (auto iter = trees.begin(); iter != trees.end();) {
        ShortestPathTree *tree = *iter;
        tree->cut_subtree(root);
        if (!tree->size()) {
            logger_info << "remove tree " << root << std::endl;
            delete tree;
            iter = trees.erase(iter);
        } else {
            ++iter;
        }
    }
    cut[root] = true;
    logger_info << trees.size() << " trees, " << total_size << " nodes in total." << std::endl;
}

int ShortestPathTree::max_path_cover_node() {
    while (!heap.empty() && (cut[heap.root_value().second] || heap.root_value().first + path_n[heap.root_value().second])) {
        heap.delete_min();
    }
    return heap.empty() ? -1 : heap.root_value().second;
}

void ShortestPathTree::add_trees(Graph &graph, DijkstraBreadthFirstSearch &dijk_search, bool initialize_order) {
    // TODO: delete this
    if (!cut) cut = new bool[graph.vertex_number()]();
    if (!path_n) path_n = new long long[graph.vertex_number()]();

    static int *v, *cur;
    if (initialize_order) {
        // TODO: delete
        v = new int[graph.vertex_number()+1];
        for (int i = 0; i < graph.vertex_number(); ++i) {
            v[i] = i;
        }
        v[graph.vertex_number()] = -1;
        std::random_shuffle(v, v+graph.vertex_number());
        cur = v;
    }
    while (total_size < maximum_size) {
        int root;
        for (root = *cur; root != -1 && !graph.valid((Vertex){root}); root = *++cur);
        // randomly select a root
        int inc = 0;
        if (root != -1) {
            // create tree
            if (!cut[root]) {
                inc = ShortestPathTree::new_tree(graph, dijk_search, root);
            }
            ++cur;
        }
        if (!inc) {
            break;
        } else {
            total_size += inc;
        }
    }
}

ShortestPathTree::ShortestPathTree(Graph &graph, DijkstraBreadthFirstSearch &dijk_search): graph(graph), dijk_search(dijk_search) {
    tree = new Graph();
}

ShortestPathTree::~ShortestPathTree() {
    delete tree;
}

int ShortestPathTree::generate(int root) {
    this->root = root;
    // int *queue = new int[graph.vertex_number()]();
    std::vector<Vertex> queue; // nodes on tree

    dijk_search.reset({{root}});
    dijk_search.next();
    Distance d;
    // construct tree
    Vertex r_v = tree->new_vertex(root);
    while (graph.valid((d = dijk_search.front()).node)) {
        if (cut[d.node.v_pos]) {
            dijk_search.next(1);
        } else {
            dijk_search.next();
            // logger_info << "add tree edge " << d.parent.v_pos << " " << d.node.v_pos << std::endl;
            tree->add_edge(d.parent.v_pos, d.node.v_pos, 0);
            queue.push_back(tree->get_vertex_by_id(d.node.v_pos));
        }
    }
    subtree_size.resize(tree->vertex_number(), 0);

    // count sub-tree size
    for (int v = queue.size()-1; ~v; --v) {
        Vertex x = queue[v];
        int x_pos_in_graph = tree->get_vertex_load(x).id;
        // logger_info << "check node " << x_pos_in_graph << " tree pos " << x.v_pos << std::endl;
        subtree_size[x.v_pos] += 1;
        path_n[x_pos_in_graph] += subtree_size[x.v_pos];
        heap.insert({-path_n[x_pos_in_graph], x_pos_in_graph});
        for (auto e=tree->first_in(x); tree->valid(e); e=tree->next_in(e)) {
            subtree_size[tree->get_source(e).v_pos] += subtree_size[x.v_pos];   
        }
    }
    subtree_size[r_v.v_pos] += 1;
    path_n[root] += subtree_size[r_v.v_pos];
    heap.insert({-path_n[root], root});
    return subtree_size[r_v.v_pos];
}

void ShortestPathTree::cut_subtree(int cut_id) {
    Vertex cutting = tree->get_vertex_by_id(cut_id);
    if (!tree->valid(cutting) || !subtree_size[cutting.v_pos]) return;
    // logger_info << "cutting " << cut_id << " on tree rooted at " << root << std::endl;
    // logger_info << "id pos in tree: " << cutting.v_pos << "/" << subtree_size.size() << std::endl;
    // update number of nodes in all trees
    int cutted_size = subtree_size[cutting.v_pos];
    total_size -= cutted_size;
    // update subtree size of ancestors
    BreadthFirstSearch bfs(*tree);
    bfs.reset({cutting}, 1);
    bfs.next();
    Vertex v;
    while (tree->valid(v = bfs.next())) {
        int v_id = tree->id(v);
        path_n[v_id] -= cutted_size;
        subtree_size[v.v_pos] -= cutted_size;
        heap.insert({-path_n[v_id], v_id});
    }
    // update subtree size of descendants
    bfs.reset({cutting});
    while (tree->valid(v = bfs.next())) {
        int v_id = tree->id(v);
        path_n[v_id] -= subtree_size[v.v_pos];
        subtree_size[v.v_pos] = 0;
        heap.insert({-path_n[v_id], v_id});
    }
    // cut the tree
    if (cut_id != root) {
        auto e=tree->first_in(cutting);
        // logger_info << "cutting " << cutting.v_pos << " edge " << tree->get_source(e).v_pos << "->" << tree->get_target(e).v_pos << std::endl;
        tree->delete_edge(e);
    }
}

int ShortestPathTree::size() {
    return subtree_size[0];
}

Label::Label() {}
Label::Label(int node, int rank, double distance): node(node), rank(rank), distance(distance) {}
// Labels::Labels(): n(0), labels(NULL) {}

// Labels::Labels(Distance *_labels, int n): n(n) {
//     labels = (Distance *)malloc(sizeof(Distance) * (n+1));
//     memcpy(labels, _labels, sizeof(Distance) * n);
//     labels[n] = {};
// }

// Labels::~Labels() {
//     free(labels);
// }
const char *to_string(IndexOrder order) {
    switch (order) {
        case BY_INPUT_ORDER: return "io";
        case RANDOM: return "random";
        case BY_DEGREE: return "d";
        case PATH_COVER: return "pc";
    }
    return "";
}
