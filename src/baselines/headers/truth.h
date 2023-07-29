#ifndef BF_H
#define BF_H
#include "method.h"

struct BruteForce: Method {
    Graph &graph;
    int sample_n;
    BruteForce(Graph &graph);
    double query(const DistanceQuery &query);
};

#endif