#include <random>
#include "query.h"
#include "struct.h"
#include "method.h"
#include "naive_monte_carlo.h"

struct SampledEdge {
    int round;
    Edge e;
    // bool operator < (const SampledEdge &that);
};

struct LazyPropagationSampling: Method {
    static std::default_random_engine generator;
    Graph &graph;
    int sample_n;
    LazyPropagationSampling(Graph &graph, int sample_n);
    int geometric_dist(double p);
    double query(const DistanceQuery &q);
};