#ifndef SP_DIJK
#define SP_DIJK
#include "method.h"
#include "breadth_first_search.h"

struct ShortestPathDijkstra: Method {
    Graph &graph;
    ShortestPathDijkstra(Graph &graph);
    double query(const DistanceQuery &query);
};
#endif