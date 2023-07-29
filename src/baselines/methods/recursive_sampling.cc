#include "recursive_sampling.h"
#include "math.h"

std::vector<ReachedVertexStack> ReachedVertexStack::checkpoint_stack;

Edge ReachedVertexStack::next_edge() {
    // find next valid edge
    // while (vertex_index < vertices.size() && !graph.valid(vertices[vertex_index].edge)) {
    //     vertex_index++;
    // }
    // if (vertex_index >= vertices.size()) return {-1}; // TODO: static const invalid edge

    // Edge result = vertices[vertex_index].edge;
    // // if not last edge of current vertex
    // if (graph.last_out(vertices[vertex_index].vertex).e_pos != result.e_pos) {
    //     // move to next edge of current vertex
    //     vertices[vertex_index].edge = graph.next_out(vertices[vertex_index].edge);
    // } else {
    //     // move to next vertex
    //     vertex_index++;
    //     while (vertex_index < vertices.size() && !graph.valid(graph.last_out(vertices[vertex_index].vertex))) {
    //         vertex_index++;
    //     }
    // }
    // return result;
    while (vertex_index < vertices.size() && !graph.valid(vertices[vertex_index].edge)) {
        vertex_index++;
    }
    if (vertex_index >= vertices.size()) return {-1}; // TODO: static const invalid edge
    Edge result = vertices[vertex_index].edge;
    vertices[vertex_index].edge = graph.next_out(vertices[vertex_index].edge);
    return result;
}

void ReachedVertexStack::print() {
    printf("vertex_index %d\n", vertex_index);
    for (auto v: vertices) {
        printf("v %d e %d->%d\n", graph.id(v.vertex), graph.id(graph.get_source(v.edge)), graph.id(graph.get_target(v.edge)));
    }
}

void ReachedVertexStack::restore_last_edge() {
    if (vertex_index >= vertices.size()) {
        vertex_index = vertices.size() - 1;
        while (vertex_index >= 0 && !graph.valid(graph.last_out(vertices[vertex_index].vertex))) {
            vertex_index--;
        }
    } else {
        if (vertices[vertex_index].edge.e_pos != graph.first_out(vertices[vertex_index].vertex).e_pos) {// TODO: Edge class equal
            vertices[vertex_index].edge = graph.last_out(vertices[vertex_index].edge);
        } else {
            vertex_index--;
            while (vertex_index >= 0 && !graph.valid(graph.last_out(vertices[vertex_index].vertex))) {
                vertex_index--;
            }
        }
    }
    // printf("restore edge %d %d\n", vertex_index, vertices[vertex_index].edge.e_pos);
}

void ReachedVertexStack::checkpoint() {
    // Edge e = vertex_index < vertices.size() ? vertices[vertex_index].edge : (Edge){-1};
    // checkpoint_stack.push_back({vertex_index, e});
    // printf("check point ");
    // for (int i = 0; i < vertices.size(); ++i) {
    //     printf("%d ",vertices[i].vertex.v_pos);
    // }
    // printf(" index %d\n", vertex_index);
    checkpoint_stack.push_back(*this);
}

void ReachedVertexStack::restore() {
    vertex_index = checkpoint_stack.back().vertex_index;
    vertices = checkpoint_stack.back().vertices;
    // printf("check point restore ");
    // for (int i = 0; i < vertices.size(); ++i) {
    //     printf("%d ",vertices[i].vertex.v_pos);
    // }
    // printf(" index %d\n", vertex_index);
    checkpoint_stack.pop_back();
    // vertex_index = checkpoint_stack.back().vertex_index;
    // if (vertex_index < vertices.size()) vertices[vertex_index].edge = checkpoint_stack.back().edge;
    // checkpoint_stack.pop_back();
}

bool ReachedVertexStack::find(Vertex vertex) {
    for (auto v: vertices) {
        if (v.vertex.v_pos == vertex.v_pos) { // TODO: equal for Vertex class
            return true;
        }
    }
    return false;
}

void ReachedVertexStack::push_vertex(Vertex v) {
    // Logger::Get(logINFO) << "push vertex " << graph.id(v) << std::endl;
    vertices.push_back({v, graph.first_out(v)});
}

void ReachedVertexStack::pop_vertex() {
    // TODO: pop error if vertex has unrestored edges
    // Logger::Get(logINFO) << "pop vertex " << vertex_index << " " << vertices.size() << std::endl;
    vertices.pop_back();
    // Logger::Get(logINFO) << "after pop vertex " << vertex_index << " " << vertices.size() << std::endl;
}

ReachedVertexStack::ReachedVertexStack(Graph &graph): vertex_index(0), graph(graph) {}
ReachedVertexStack::~ReachedVertexStack() {}


double recursive_sampling(Graph &graph,
                          const DistanceQuery &query,
                          ReachedVertexStack &edge_stack,
                          int sample_n,
                          int depth) {
    static const int sample_threshold = 10; // TODO: make tunable
    static const int max_depth = 1024;
    // logger_info << "recursion" << std::endl;
    if (sample_n <= 0) {
        return 0;
    }
    if (sample_n <= sample_threshold || depth >= max_depth) { // TODO: should depth be constrained?
        if (depth >= max_depth) {
            logger_info << "Recursive Sampling: Exceed Max Depth " << max_depth << std::endl;
        }
        return naive_monte_carlo(query, sample_n, graph);
    } else {
        // TODO: fastly check 1
        if (edge_stack.find(query.t)) {
            return 1;
        }
        // TODO: fastly check 0
        // TODO: how to choose a good edge?
        Edge edge;
        edge_stack.checkpoint();
        // edge_stack.print();
        do {
            edge = edge_stack.next_edge();
        } while (graph.valid(edge) && edge_stack.find(graph.get_target(edge)));
        // if (graph.valid(edge)) logger_info << "find edge " << graph.id(graph.get_source(edge)) << " " << graph.id(graph.get_target(edge)) << std::endl;
        // else logger_info << "find edge -1"<< std::endl;
        if (!graph.valid(edge)) {
            edge_stack.restore();
            return 0;
        }
        // sample an edge
        double p = graph.get_edge_load(edge);
        double result = 0;
        bool have_pushed_vertex = false;
        // edge_stack.checkpoint();
        if (!edge_stack.find(graph.get_target(edge))) {
            have_pushed_vertex = true;
            edge_stack.push_vertex(graph.get_target(edge));
        } else {
            // TODO: skip edges pointing to reached vertex. They do not influence the result.
        }
        // recursion
        graph.update_edge(edge, 1);
        // logger_info << "into recursion sample edge " << graph.id(graph.get_source(edge)) << "->" << graph.id(graph.get_target(edge)) << std::endl;
        double r = p * recursive_sampling(graph, query, edge_stack, round(sample_n * p), depth+1);
        // logger_info << "backtracking sample edge " << graph.id(graph.get_source(edge)) << "->" << graph.id(graph.get_target(edge)) << std::endl;
        // logger_info << "result " << r << std::endl;
        result += r;
        // backtracking: remove reached vertex
        // edge_stack.restore();
        if (have_pushed_vertex) {
            edge_stack.pop_vertex();
        }
        // recursion
        // logger_info << "into recursion remove edge " << graph.id(graph.get_source(edge)) << "->" << graph.id(graph.get_target(edge)) << std::endl;
        graph.update_edge(edge, 0);
        r = (1-p) * recursive_sampling(graph, query, edge_stack, round(sample_n * (1-p)), depth+1);
        // logger_info << "backtracking remove edge " << graph.id(graph.get_source(edge)) << "->" << graph.id(graph.get_target(edge)) << std::endl;
        // logger_info << "result " << r << std::endl;
        result += r;
        // backtracking
        graph.update_edge(edge, p);
        edge_stack.restore();
        return result;
    }
}

double recursive_sampling(Graph &graph, const DistanceQuery &query, int sample_n) {
    ReachedVertexStack edge_stack(graph);
    edge_stack.push_vertex(query.s);
    return recursive_sampling(graph, query, edge_stack, sample_n, 0);
}

RecursiveSampling::RecursiveSampling(Graph &graph, int sample_n): graph(graph), sample_n(sample_n) {}

double RecursiveSampling::query(const DistanceQuery &query) {
    return recursive_sampling(graph, query, sample_n);
}