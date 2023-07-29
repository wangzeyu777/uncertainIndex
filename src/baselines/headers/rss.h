#ifndef __RSS_H__
#define __RSS_H__

#include "struct.h"
#include "recursive_sampling.h"

double recursive_stratified_sampling_II(Graph &graph,
                                        const DistanceQuery &query,
                                        ReachedVertexStack &edge_stack,
                                        bool *invalid_vertex,
                                        int sample_n,
                                        int depth);

double recursive_stratified_sampling_II(Graph &graph, const DistanceQuery &query, int sample_n);

struct RecursiveStratifiedSampling: Method {
    Graph &graph;
    int sample_n;
    static int max_edge_numger;
    RecursiveStratifiedSampling(Graph &graph, int sample_n);
    void set_rss_r(int);
    double query(const DistanceQuery &query);
};

#endif