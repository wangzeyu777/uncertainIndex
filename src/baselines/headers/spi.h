#ifndef __SPI_H__
#define __SPI_H__

#include <string.h>
#include <algorithm>
#include <vector>
#include <set>
#include "math.h"
#include "struct.h"
#include "query.h"
#include "breadth_first_search.h"
#include "method.h"
#include "config.h"
#include "pruned_landmark_labeling.h"
#include "naive_monte_carlo.h"

struct ReachabilityLabel: Label {
    double r;
    ReachabilityLabel();
    ReachabilityLabel(int root, int rank, double distance, double reachability);    
};

bool operator < (const ReachabilityLabel &a, const ReachabilityLabel &b);

typedef std::vector<ReachabilityLabel> RLabels; // TODO: use pointer
bool exist_rank(const RLabels &labels, int rank);

struct ShortestPathIndexing: Method {
    // TODO: tunable EPS
    const double EPS = 1e-9;
    Graph &graph;
    IndexOrder order;
    int sample_n;
    bool struct_opt;
    double threshold_r;
    int threshold_hop;

    RLabels *label_forward;
    RLabels *label_backward;
    ShortestPathIndexing(Graph &graph, int sample_n, IndexOrder order=BY_DEGREE, bool struct_opt=false, double threshold_r=0, int threshold_hop=2e9);
    ~ShortestPathIndexing();

    void index();
    void save_index(const char *filename);
    void delete_index();
    double query(const DistanceQuery &q);
    double query_distance(Vertex s, Vertex t, int direction=0);
    double query_reachability_struct_opt(Vertex s, Vertex t);

    double query_reachability(Vertex s, Vertex t);
    // long long size();
    void load_index(const char *filename);

    int tree_bound = 2e9;
    void set_tree_bound(int bound);
};

#endif