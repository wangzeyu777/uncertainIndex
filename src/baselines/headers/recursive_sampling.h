#ifndef RECRUSIVE
#define RECRUSIVE

#include "struct.h"
#include "query.h"
#include "logger.h" 
#include "naive_monte_carlo.h"

struct ReachedVertexStack {
    struct ReachedVertex {
        Vertex vertex;
        Edge edge;
    };
    struct Checkpoint {
        int vertex_index;
        Edge edge;
    };
    Graph &graph;
    std::vector<ReachedVertex> vertices;
    static std::vector<ReachedVertexStack> checkpoint_stack;
    int vertex_index;

    Edge next_edge();
    void restore_last_edge();
    bool find(Vertex vertex);
    void push_vertex(Vertex v);
    void pop_vertex();
    void checkpoint();
    void restore();
    void print();
    ReachedVertexStack(Graph &graph);
    ~ReachedVertexStack();
};

struct RecursiveSampling: Method {
    Graph &graph;
    int sample_n;
    RecursiveSampling(Graph &graph, int sample_n);
    double query(const DistanceQuery &query);
};

double recursive_sampling(Graph &graph, const DistanceQuery &query, int sample_n);
double recursive_sampling(Graph &graph, const DistanceQuery &query, ReachedVertexStack &edge_stack,
                          int sample_n, int depth);

#endif