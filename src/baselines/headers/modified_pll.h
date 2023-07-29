#ifndef __MPLL_H__
#define __MPLL_H__

#include <string.h>
#include <algorithm>
#include <vector>
#include "math.h"
#include "struct.h"
#include "query.h"
#include "breadth_first_search.h"
#include "method.h"
#include "config.h"
#include "pruned_landmark_labeling.h"

struct MLabel {
    int node, rank;
    double distance, prob;
    int parent, parent_pos;
    MLabel();
    MLabel(int node, int parent, int parent_pos, int rank, double distance);
};

typedef std::vector<MLabel> MLabels;

struct MPLL: Method {
    // TODO: tunable EPS
    const double EPS = 1e-9;
    Graph &graph;
    IndexOrder order;

    MLabels *label_forward;
    MLabels *label_backward;
    MPLL(Graph &graph, IndexOrder order=BY_DEGREE);
    ~MPLL();
    void index();
    void save_index(const char *filename);
    void delete_index();
    double query(const DistanceQuery &q);

    double query(Vertex s, Vertex t, int direction=0);
    // long long size();
    void load_index(const char *filename);

    double path_length(MLabel *, Vertex, MLabel *, Vertex, int direction=0); // 0 forward, 1 backward
};

#endif