#include <math.h>
// #include <algorithm>
#include "rss.h"

double recursive_stratified_sampling_II(Graph &graph,
                                        const DistanceQuery &query,
                                        ReachedVertexStack &edge_stack,
                                        bool *invalid_vertex,
                                        int sample_n,
                                        int depth) {
    static const int sample_threshold = 10;
    static const int max_depth = 1024;
    if (sample_n <= 0) {
        return 0;
    }
    if (sample_n <= sample_threshold || depth >= max_depth) { // TODO: should depth be constrained?
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
        for (int i = 0; i < RecursiveStratifiedSampling::max_edge_numger; ++i) {
            // TODO: skip edges pointing to reached vertices
            Edge edge = edge_stack.next_edge();
            if (graph.valid(edge)) {
                edges.push_back(edge);
            }
        }

        if (!edges.size()) {
            // different from single edge choice because in that case all edges are sampled
            // printf("nmc sample_n=%d\n", sample_n);
            edge_stack.restore();
            // graph.print();
            double res = naive_monte_carlo(query, sample_n, graph, invalid_vertex);
            // printf("nmc res=%f\n", res);
            return res;
        }
        double edge_probs[edges.size()];
        double prob_prefix[edges.size()];
        double prefix = 1;
        
        // printf("edges: ");
        for (int i = 0; i < edges.size(); ++i) {
            edge_probs[i] = graph.get_edge_load(edges[i]);
            graph.update_edge(edges[i], 0);
            prob_prefix[i] = prefix;
            prefix *= 1-edge_probs[i];
            // printf(" %d->%d", graph.id(graph.get_source(edges[i])), graph.id(graph.get_target(edges[i])));
        }
        // printf("\n");

        // all edges are not sampled
        double res = recursive_stratified_sampling_II(graph, query, edge_stack, invalid_vertex,
            std::max(1, (int)round(sample_n * prefix)), depth+1);
        double result = prefix * res;
        // printf("depth %d, world_prob=%f, res=%f \n", depth, prefix, res);

        for (int i = edges.size()-1; ~i; --i) {
            Edge edge = edges[i];
            double world_prob = prob_prefix[i] * edge_probs[i];
            // recurse
            graph.update_edge(edge, 1);
            bool push_vertex = false;
            if (!edge_stack.find(graph.get_target(edge))) {
                // logger_info << "push vertex " << graph.id(graph.get_target(edge)) << std::endl;
                edge_stack.push_vertex(graph.get_target(edge));
                push_vertex = true;
            }
            // printf("sample edge %d->%d\n", graph.id(graph.get_source(edge)), graph.id(graph.get_target(edge)));
            double res = recursive_stratified_sampling_II(graph, query, edge_stack, invalid_vertex,
                std::max(1, (int) round(sample_n * world_prob)), depth+1); // TODO: sample number unconsistent(in total)
            // printf("depth %d, world_prob=%f, res=%f \n", depth, world_prob, res);
            result += world_prob * res;
            if (push_vertex) {
                // logger_info << "pop vertex " << graph.id(graph.get_target(edge)) << std::endl;
                edge_stack.pop_vertex();
            }
            // printf("restore edge %d->%d\n", graph.id(graph.get_source(edge)), graph.id(graph.get_target(edge)));
            // set for next edge
            graph.update_edge(edge, edge_probs[i]);
        }
        // backtracking
        edge_stack.restore();
        // for (int i = 0; i < edges.size(); ++i) {            
        //     graph.update_edge(edges[i], edge_probs[i]);
        // }
        return result;
    }
    return 0;
}

double recursive_stratified_sampling_II(Graph &graph, const DistanceQuery &query, int sample_n) {
    ReachedVertexStack edge_stack(graph);
    edge_stack.push_vertex(query.s);
    bool *invalid_vertex = new bool[graph.vertex_number()]();
    double result = recursive_stratified_sampling_II(graph, query, edge_stack, invalid_vertex, sample_n, 0);
    delete []invalid_vertex;
    return result;
}

RecursiveStratifiedSampling::RecursiveStratifiedSampling(Graph &graph, int sample_n): graph(graph), sample_n(sample_n) {}
double RecursiveStratifiedSampling::query(const DistanceQuery &query) {
    return recursive_stratified_sampling_II(graph, query, sample_n);
}

int RecursiveStratifiedSampling::max_edge_numger = 50;
void RecursiveStratifiedSampling::set_rss_r(int r) { max_edge_numger = r;}