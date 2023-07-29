#include "query.h"
#include "struct.h"
#include "lazy_propagation_sampling.h"
#include "fib.h"

std::default_random_engine LazyPropagationSampling::generator;

bool operator < (const SampledEdge &lhs, const SampledEdge &rhs) {
    return lhs.round < rhs.round;
}

template<class T>
struct Heap {
    Storage<T> values;
    void push(const T&x) {
        values.resize(values.size + 1);
        values.size++;
        int i = values.size - 1, j = (i-1) / 2;
        while (i != 0 && x < *values.get(j)) {
            values.set(i, *values.get(j));
            i = j;
            j = (j-1) / 2;
        }
        values.set(i, x);
    }
    void increase_root(const T&x) {
        int i = 0, j = i*2+1;
        while (j < values.size) {
            if (j+1 < values.size && *values.get(j+1) < *values.get(j)) {
                j = j+1;
            }
            if (*values.get(j) < x) {
                values.set(i, *values.get(j));
                i = j;
                j = j*2+1;
            } else {
                break;
            }
        }
        values.set(i, x);
    }
    const T&root() {
        return *values.get(0);
    }
    bool empty() {
        return values.size == 0;
    }
};

LazyPropagationSampling::LazyPropagationSampling(Graph &graph, int iter_num) : graph(graph), sample_n(iter_num) {}

int LazyPropagationSampling::geometric_dist(double p) {
    if (p > 1-NaiveMonteCarlo::EPS) {
        return 0;
    }
    if (p < NaiveMonteCarlo::EPS) {
        return 2e9;
    }
    std::geometric_distribution<int> dist(p);
    return dist(generator);
}

double LazyPropagationSampling::query(const DistanceQuery &q) {
    // TODO: new in constructor
    int *visited = new int[graph.vertex_number()]();
    int *queue = new int[graph.vertex_number()];
    int *visit_count = new int[graph.vertex_number()]();
    // TODO: use heap with smaller memory
    Heap<SampledEdge> *heap = new Heap<SampledEdge>[graph.vertex_number()]();
    int result = 0;

    for (int round = 1; round <= sample_n; ++round) {
        int op = 0, cl = 1;
        queue[0] = q.s.v_pos;
        visited[q.s.v_pos] = round;
        int reached = 0;
        while (op < cl) {
            int v = queue[op++];
            Vertex x = {v};
            if (!visit_count[v]) {
                // init
                visit_count[v] = 1;
                for (auto e = graph.first_out(x); graph.valid(e); e = graph.next_out(e)) {
                    heap[v].push({geometric_dist(graph.get_edge_load(e)) + visit_count[v], e});
                }
            }
            int v_round;
            while (!heap[v].empty() && visit_count[v] >= (v_round = heap[v].root().round)) {
                Edge e = heap[v].root().e;
                Vertex u = graph.get_target(e);
                heap[v].increase_root({geometric_dist(graph.get_edge_load(e)) + v_round + 1, e});
                if (visited[u.v_pos] != round && v_round == visit_count[v]) {
                    if (u.v_pos == q.t.v_pos) {
                        reached = 1;
                        break;
                    }
                    visited[u.v_pos] = round;
                    queue[cl++] = u.v_pos;
                }
            }
            visit_count[v]++;
            if (reached) {
                result += 1;
                break;
            }
        }
        // logger_info << "round " << round << std::endl;
        // for (int i = 0; i < cl; ++i) {
        //     logger_info << queue[i] << " ";
        // }
        // if (reached) logger_info << q.t->id;
        // logger_info << std::endl;
        // char c = getchar();
        if (round % 1000 == 0) {
            logger_info << "sample_n: " << round << " result: " << result << std::endl;
        }
    }

    delete []heap;
    delete []queue;
    delete []visited;
    delete []visit_count;

    return result / (double) sample_n;
}
