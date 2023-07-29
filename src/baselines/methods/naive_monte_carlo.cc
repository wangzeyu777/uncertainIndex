#include "naive_monte_carlo.h"

const double NaiveMonteCarlo::EPS = 1e-8;

bool sample_prob(double p) {
    if (p < NaiveMonteCarlo::EPS) return 0;
    if (p > 1-NaiveMonteCarlo::EPS) return 1;
    return rand() <= (int)(p * RAND_MAX);
}

NaiveMonteCarlo::NaiveMonteCarlo(Graph &graph, int sample_n): graph(graph), sample_n(sample_n) {}
double NaiveMonteCarlo::query(const DistanceQuery &query) {
    return naive_monte_carlo_logging(query, sample_n, graph, NULL);
}

// sample once
double naive_monte_carlo_once(const DistanceQuery &query, Graph &graph, bool *invalid) {
    if (query.s.v_pos == query.t.v_pos) {return 1;}
    static unsigned int *visited = new unsigned int[graph.vertex_number()](); // TODO: should this be deleted? or make global struct
    static unsigned int flag = 0;
    static const unsigned int reset_flag_round = 4e9;

    if ((++flag) == reset_flag_round) {
        flag = 1;
        memset(visited, 0, sizeof(unsigned int) * graph.vertex_number());
    }

    std::queue<Vertex> q;
    q.push(query.s);
    visited[query.s.v_pos] = flag;
    // printf("sample from %d\n", graph.id(query.s));
    while (!q.empty() && flag != visited[query.t.v_pos]) {
        Vertex x = q.front();
        // for (auto e: x->edges) {
        for (auto e=graph.first_out(x); graph.valid(e); e=graph.next_out(e)) { // iterate neighbors
            Vertex y = graph.get_target(e);
            if (flag != visited[y.v_pos] && (!invalid || !invalid[y.v_pos])) {
                if (sample_prob(graph.get_edge_load(e))) {
                    if (y.v_pos == query.t.v_pos) return 1;
                    q.push(y);
                    // printf("get %d ", graph.id(y));
                    visited[y.v_pos] = flag;
                }
            }
        }
        q.pop();
    }

    return flag == visited[query.t.v_pos];
}

// sample multiple times
double naive_monte_carlo(const DistanceQuery &query, int iter_num, Graph &graph, bool *valid) {
    int result = 0;
    for (int iter = 1; iter <= iter_num; ++iter) {
        result += naive_monte_carlo_once(query, graph, valid);
    }
    return result / (double) iter_num;
}

int naive_monte_carlo_logging_int(const DistanceQuery &query, int iter_num, Graph &graph, bool *valid, bool silence) {
    int result = 0;
    for (int iter = 1; iter <= iter_num; ++iter) {
        result += naive_monte_carlo_once(query, graph, valid);
        if (iter % 50 == 0 && !silence) {
            logger_info << "sample_n: " << iter << " result: " << result << std::endl;
        }
    }
    return result;
}

double naive_monte_carlo_logging(const DistanceQuery &query, int iter_num, Graph &graph, bool *valid, bool silence) {
    return naive_monte_carlo_logging_int(query, iter_num, graph, valid, silence) / (double) iter_num;
}

double naive_monte_carlo(const DistanceQuery &query, int iter_num, Graph &graph) {
    return naive_monte_carlo(query, iter_num, graph, NULL);
}

void naive_monte_carlo(Vertex s, int sample_n, Graph &graph, int direction, const bool *valid, int *counter) {
    // TODO: put following thing to namespace in this file
    static unsigned int *visited = new unsigned int[graph.vertex_number()](); // TODO: should this be deleted? or make global struct
    static unsigned int flag = 0;
    static const unsigned int reset_flag_round = 4e9;

    counter[s.v_pos] = sample_n;
    auto first_edge = direction ? &Graph::first_in : &Graph::first_out;
    auto next = direction ? &Graph::next_in : &Graph::next_out;
    auto get_vertex = direction ? &Graph::get_source : &Graph::get_target;

    while (sample_n--) {
        if ((++flag) == reset_flag_round) {
            flag = 1;
            memset(visited, 0, sizeof(unsigned int) * graph.vertex_number());
        }

        std::queue<Vertex> q;
        q.push(s);
        visited[s.v_pos] = flag;
        while (!q.empty()) {
            Vertex x = q.front();
            // for (auto e: x->edges) { 
            for (auto e=(graph.*first_edge)(x); graph.valid(e); e=(graph.*next)(e)) { // iterate neighbors
                Vertex y = (graph.*get_vertex)(e);
                if (valid && !valid[y.v_pos]) {
                    continue;
                }
                int bound = graph.get_edge_load(e) * RAND_MAX;
                int sample = rand();
                if (sample < bound) {
                    if (flag != visited[y.v_pos]) {
                        q.push(y);
                        visited[y.v_pos] = flag;
                        ++counter[y.v_pos];
                    }
                }
            }
            q.pop();
        }
    }
}
