#ifndef __PLL_H__
#define __PLL_H__

#include <string.h>
#include <algorithm>
#include <vector>
#include "math.h"
#include "struct.h"
#include "query.h"
#include "breadth_first_search.h"
#include "method.h"
#include "config.h"

// bool degree_cmp(const Vertex *a, const Vertex *b);

struct Label {
    int node, rank;
    double distance;
    Label();
    Label(int node, int rank, double distance);
};
// struct Labels {
//     Distance *labels;
//     int n;
//     Labels();
//     Labels(Distance *_labels, int n);
//     ~Labels();
// };
typedef std::vector<Label> Labels;

const char *to_string(IndexOrder order);

struct ShortestPathTree {
    static long long maximum_size;
    static long long total_size;
    static long long* path_n;
    static bool* cut;
    static Fib<std::pair<long long, int>> heap; // <path_n, v_id>, TODO: a better heap
    static std::vector<ShortestPathTree *> trees;
    static int new_tree(Graph &graph, DijkstraBreadthFirstSearch &dijk_search, int root);
    static void cut_all(int root);
    static int max_path_cover_node();
    static void add_trees(Graph &, DijkstraBreadthFirstSearch &, bool f=false);

    Graph &graph;
    DijkstraBreadthFirstSearch &dijk_search;

    Graph *tree;
    int root; // id of tree root at original graph
    std::vector<int> subtree_size;

    ShortestPathTree(Graph &graph, DijkstraBreadthFirstSearch &dijk_search);
    ~ShortestPathTree();

    int generate(int root);
    void cut_subtree(int root);
    int size();
};

struct PrunedLandmarkLabeling: Method {
    // TODO: tunable EPS
    const double EPS = 1e-9;
    Graph &graph;
    IndexOrder order;

    Labels *label_forward;
    Labels *label_backward;
    PrunedLandmarkLabeling(Graph &graph, IndexOrder order=BY_DEGREE);
    ~PrunedLandmarkLabeling();
    void index();
    void save_index(const char *filename);
    void delete_index();
    double query(const DistanceQuery &q);

    double query(Vertex s, Vertex t, int direction=0);
    // long long size();
    void load_index(const char *filename);
};

#endif