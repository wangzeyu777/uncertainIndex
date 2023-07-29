#include "truth.h"
inline double kahanSum(std::vector<double> &fa) {
    double sum = 0.0;
    double c = 0.0;
    for(double f : fa) {
        double y = f - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}
BruteForce::BruteForce(Graph &graph): graph(graph) {}
double BruteForce::query(const DistanceQuery &query) {
    double result = -1;
    const int max_edge_number = 24;
    int edge_num = graph.edge_number();
    if (edge_num < max_edge_number) {
        result = 0;
        std::vector<double> p_list;
        int visited[graph.vertex_number()];
        memset(visited, 0, graph.vertex_number() * sizeof(int));
        for (int i = 0; i < (1 << edge_num); ++i) {
            visited[query.s.v_pos] = i+1;
            Vertex q[graph.vertex_number()];
            q[0] = query.s;
            int tail = 1, head = 0;
            while (head < tail) {
                Vertex x = q[head];
                for (auto edge=graph.first_out(x); graph.valid(edge); edge=graph.next_out(edge)) {
                    if ((1 << edge.e_pos) & i) {
                        Vertex y = graph.get_target(edge);
                        if (i+1 != visited[y.v_pos]) {
                            visited[y.v_pos] = i+1;
                            q[tail++] = {y};
                            if (y.v_pos == query.t.v_pos) {
                                break;
                            }
                        }
                    }
                }
                if (visited[query.t.v_pos] == i+1) {
                    break;
                }
                head++;
            }
            double p = 1;
            for (int e = 0; e < graph.edge_number(); ++e) {
                double pe = graph.get_edge_load({e});
                p *= (i & (1 << e)) ? pe : 1-pe;
            }
            if (visited[query.t.v_pos] == i+1) {
                p_list.push_back(p);
            }
        }
        result = kahanSum(p_list);
    }
    return result;
}
