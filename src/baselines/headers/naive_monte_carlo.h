#ifndef NMC
#define NMC
#include <queue>
#include "method.h"

double naive_monte_carlo(const DistanceQuery &query, int sample_n, Graph &graph);
double naive_monte_carlo(const DistanceQuery &query, int sample_n, Graph &graph, bool *valid);
void naive_monte_carlo(Vertex s, int sample_n, Graph &graph, int direction, const bool *valid, int *counter);
double naive_monte_carlo_logging(const DistanceQuery &query, int iter_num, Graph &graph, bool *valid, bool silence=false);
int naive_monte_carlo_logging_int(const DistanceQuery &query, int iter_num, Graph &graph, bool *valid, bool silence);

struct NaiveMonteCarlo: Method {
    static const double EPS;
    Graph &graph;
    int sample_n;
    NaiveMonteCarlo(Graph &graph, int sample_n);
    double query(const DistanceQuery &query);
};

bool sample_prob(double p);

#endif