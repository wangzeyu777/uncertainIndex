#include "rss1.h"

double recursive_stratified_sampling_I(Graph &graph,
                                        const DistanceQuery &query,
                                        ReachedVertexStack &edge_stack,
                                        bool *invalid_vertex,
                                        int sample_n,
                                        int depth) {
    static const int r = 5; // TODO: make tunable
    static const int sample_threshold = 5;
    static const int max_depth = 1024;
    if (sample_n <= 0) {
        return 0;
    }
    if (sample_n < sample_threshold || depth >= max_depth) { // TODO: should depth be constrained?
        if (depth >= max_depth) {
            Logger::Get(logINFO) << "Recursive Stratified Sampling: Exceed Max Depth " << max_depth << std::endl;
        }
        return naive_monte_carlo(query, sample_n, graph, invalid_vertex);
    } else {
        // TODO: fastly check 1
        if (edge_stack.find(query.t)) {
            return 1;
        }
        // TODO: fastly check 0
        // TODO: use invalid_vertex to prune vertices
        // TODO: how to choose good edges?
        edge_stack.checkpoint();
        std::vector<Edge> edges;
        for (int i = 0; i < r; ++i) {
            // TODO: skip edges pointing to reached vertices
            Edge edge = edge_stack.next_edge();
            // printf("next_edge %d\n", edge.e_pos);
            if (graph.valid(edge)) {
                edges.push_back(edge);
            } else {
                break;
            }
        }
        int e_n = edges.size();
        if (!e_n) {
            // different from single edge choice because in that case all edges are sampled
            // printf("nmc sample_n=%d\n", sample_n);
            edge_stack.restore();
            // graph.print();
            double res = naive_monte_carlo(query, sample_n, graph, invalid_vertex);
            // printf("nmc res=%f\n", res);
            return res;
        }
        double edge_probs[e_n];
        // edge_stack.print();
        // printf("edges ");
        for (int i = 0; i < e_n; ++i) {
            edge_probs[i] = graph.get_edge_load(edges[i]);
            // printf(" %d->%d", graph.id(graph.get_source(edges[i])), graph.id(graph.get_target(edges[i])));
        }
        // printf("\n");
        double result = 0;
        for (int mask = 0; mask < (1 << e_n); ++mask) {
            double world_prob = 1;
            edge_stack.checkpoint();
            for (int i = 0; i < e_n; ++i) {
                if (mask & (1<<i)) {
                    world_prob *=  edge_probs[i];
                    graph.update_edge(edges[i], 1);
                    Vertex v = graph.get_target(edges[i]);
                    if (!edge_stack.find(v)) {
                        edge_stack.push_vertex(v);
                    }
                } else {
                    world_prob *= 1-edge_probs[i];
                    graph.update_edge(edges[i], 0);
                }
            }

            // printf("edges ");
            // for (int i = 0; i < e_n; ++i) {
            //     printf(" %d->%d", graph.id(graph.get_source(edges[i])), graph.id(graph.get_target(edges[i])));
            // }
            // printf("\n");

            // printf("mask %d prob %f into recursion\n", mask, world_prob);
            // for (int i = 0; i < edge_stack.vertices.size(); ++i) {
            //     printf("%d ", edge_stack.vertices[i].vertex.v_pos);
            // }
            // printf("\n");
            double res = recursive_stratified_sampling_I(graph, query, edge_stack, invalid_vertex, sample_n * world_prob, depth+1);
            // printf("mask %d prob %f res %f\n", mask, world_prob, res);
            result += world_prob * res;
            edge_stack.restore();
        }
        
        // backtracking
        edge_stack.restore();
        for (int i = 0; i < e_n; ++i) {            
            graph.update_edge(edges[i], edge_probs[i]);
        }
        return result;
    }
    return 0;
}

double recursive_stratified_sampling_I(Graph &graph, const DistanceQuery &query, int sample_n) {
    ReachedVertexStack edge_stack(graph);
    edge_stack.push_vertex(query.s);
    bool *invalid_vertex = new bool[graph.vertex_number()]();
    double result = recursive_stratified_sampling_I(graph, query, edge_stack, invalid_vertex, sample_n, 0);
    delete []invalid_vertex;
    return result;
}

RecursiveStratifiedSamplingI::RecursiveStratifiedSamplingI(Graph &graph, int sample_n): graph(graph), sample_n(sample_n) {}
double RecursiveStratifiedSamplingI::query(const DistanceQuery &query) {
    return recursive_stratified_sampling_I(graph, query, sample_n);
}