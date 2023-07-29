#ifndef __RSS_I_H__
#define __RSS_I_H__

#include "struct.h"
#include "recursive_sampling.h"

double recursive_stratified_sampling_I(Graph &graph,
                                        const DistanceQuery &query,
                                        ReachedVertexStack &edge_stack,
                                        bool *invalid_vertex,
                                        int sample_n,
                                        int depth);

double recursive_stratified_sampling_I(Graph &graph, const DistanceQuery &query, int sample_n);

struct RecursiveStratifiedSamplingI: Method {
    Graph &graph;
    int sample_n;
    RecursiveStratifiedSamplingI(Graph &graph, int sample_n);
    double query(const DistanceQuery &query);
};

#endif