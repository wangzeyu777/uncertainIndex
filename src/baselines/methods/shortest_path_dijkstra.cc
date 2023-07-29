#include "shortest_path_dijkstra.h"

ShortestPathDijkstra::ShortestPathDijkstra(Graph &graph): graph(graph) {}

double ShortestPathDijkstra::query(const DistanceQuery &query) {
    DijkstraBreadthFirstSearch dijk_search(graph);
    // logger_info << "query " << query.s.v_pos << " " << query.s_id << std::endl;
    dijk_search.reset({query.s});
    Distance d;
    while (graph.valid((d = dijk_search.next()).node)) {
        if (d.node.v_pos == query.t.v_pos) {
            return d.distance;
        }
    }
    return 0;
}