// Independent complete-matrix synthesis extracted from the authors' anytime linear engine.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <cmath>
#include <atomic>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
using Row = std::uint64_t;
using Matrix = std::vector<Row>;
struct FreshWorkLimit {};
static std::uint64_t fresh_work_limit = 100000000, fresh_work_used = 0;
static void consume_fresh_work() {
    if (fresh_work_used >= fresh_work_limit) throw FreshWorkLimit{};
    ++fresh_work_used;
}

using RowOp = std::pair<int, int>;

enum class OpSide {
    Row,
    Column
};

static int popcount64(Row x) {
#if defined(__GNUG__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    int c = 0;
    while (x != 0) {
        x &= x - 1;
        ++c;
    }
    return c;
#endif
}

static Row mask_for(int n) {
    if (n <= 0 || n > 64) {
        throw std::runtime_error("matrix dimension must be in [1, 64]");
    }
    if (n == 64) return ~Row(0);
    return (Row(1) << n) - 1;
}

static Matrix identity_matrix(int n) {
    Matrix out(n);
    for (int i = 0; i < n; ++i) out[i] = Row(1) << i;
    return out;
}

static int check_square_rows(const Matrix& rows) {
    const int n = static_cast<int>(rows.size());
    const Row mask = mask_for(n);
    for (Row r : rows) {
        if ((r & ~mask) != 0) throw std::runtime_error("row bitset exceeds matrix dimension");
    }
    return n;
}

static void row_add(Matrix& rows, int control, int target) {
    if (control == target) throw std::runtime_error("control and target must differ");
    rows[target] ^= rows[control];
}

static void column_add(Matrix& rows, int control, int target) {
    if (control == target) throw std::runtime_error("control and target must differ");
    const Row control_bit = Row(1) << control;
    const Row target_bit = Row(1) << target;
    for (Row& row : rows) {
        if ((row & target_bit) != 0) row ^= control_bit;
    }
}

static Row column_bits(const Matrix& rows, int column) {
    Row bits = 0;
    const Row column_bit = Row(1) << column;
    for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
        if ((rows[r] & column_bit) != 0) bits |= Row(1) << r;
    }
    return bits;
}

static Matrix invert_matrix(const Matrix& rows) {
    const int n = check_square_rows(rows);
    Matrix left = rows;
    Matrix right = identity_matrix(n);
    for (int col = 0; col < n; ++col) {
        int pivot = -1;
        for (int r = col; r < n; ++r) {
            if (((left[r] >> col) & 1U) != 0) {
                pivot = r;
                break;
            }
        }
        if (pivot < 0) throw std::runtime_error("matrix is singular over GF(2)");
        if (pivot != col) {
            std::swap(left[pivot], left[col]);
            std::swap(right[pivot], right[col]);
        }
        for (int r = 0; r < n; ++r) {
            if (r != col && ((left[r] >> col) & 1U) != 0) {
                left[r] ^= left[col];
                right[r] ^= right[col];
            }
        }
    }
    return right;
}

static int gf2_rank(std::vector<Row> rows) {
    if (rows.empty()) return 0;
    const int n = check_square_rows(rows);
    int rank = 0;
    for (int col = 0; col < n; ++col) {
        int pivot = -1;
        for (int r = rank; r < n; ++r) {
            if (((rows[static_cast<std::size_t>(r)] >> col) & 1U) != 0) {
                pivot = r;
                break;
            }
        }
        if (pivot < 0) continue;
        std::swap(rows[static_cast<std::size_t>(rank)], rows[static_cast<std::size_t>(pivot)]);
        for (int r = 0; r < n; ++r) {
            if (r != rank && ((rows[static_cast<std::size_t>(r)] >> col) & 1U) != 0) {
                rows[static_cast<std::size_t>(r)] ^= rows[static_cast<std::size_t>(rank)];
            }
        }
        ++rank;
    }
    return rank;
}

static std::vector<RowOp> gaussian_reduction_ops(const Matrix& rows) {
    const int n = check_square_rows(rows);
    Matrix work = rows;
    std::vector<RowOp> ops;
    for (int col = 0; col < n; ++col) {
        if (((work[col] >> col) & 1U) == 0) {
            int pivot = -1;
            for (int r = col + 1; r < n; ++r) {
                if (((work[r] >> col) & 1U) != 0) {
                    pivot = r;
                    break;
                }
            }
            if (pivot < 0) throw std::runtime_error("matrix is singular over GF(2)");
            row_add(work, pivot, col);
            ops.emplace_back(pivot, col);
        }
        for (int r = 0; r < n; ++r) {
            if (r != col && ((work[r] >> col) & 1U) != 0) {
                row_add(work, col, r);
                ops.emplace_back(col, r);
            }
        }
    }
    if (work != identity_matrix(n)) throw std::runtime_error("Gaussian fallback failed");
    return ops;
}

static Matrix permute_rows(const Matrix& rows, const std::vector<int>& perm) {
    if (rows.size() != perm.size()) throw std::runtime_error("permutation length mismatch");
    Matrix out;
    out.reserve(rows.size());
    for (int idx : perm) out.push_back(rows.at(static_cast<std::size_t>(idx)));
    return out;
}

static bool is_permutation_matrix(const Matrix& rows) {
    const int n = check_square_rows(rows);
    Row seen = 0;
    for (Row row : rows) {
        if (popcount64(row) != 1) return false;
        if ((seen & row) != 0) return false;
        seen |= row;
    }
    return seen == mask_for(n);
}

static std::vector<int> permutation_columns(const Matrix& rows) {
    if (!is_permutation_matrix(rows)) throw std::runtime_error("matrix is not a permutation");
    std::vector<int> out;
    out.reserve(rows.size());
    for (Row row : rows) out.push_back(popcount64(row - 1));
    return out;
}

static std::vector<int> inverse_permutation(const std::vector<int>& perm) {
    std::vector<int> inv(perm.size(), -1);
    for (int i = 0; i < static_cast<int>(perm.size()); ++i) {
        const int v = perm[static_cast<std::size_t>(i)];
        if (v < 0 || v >= static_cast<int>(perm.size()) || inv[static_cast<std::size_t>(v)] != -1) {
            throw std::runtime_error("invalid permutation");
        }
        inv[static_cast<std::size_t>(v)] = i;
    }
    return inv;
}

static bool infer_output_permutation(
    const Matrix& target_rows,
    const Matrix& final_rows,
    std::vector<int>& output_permutation
) {
    const int n = check_square_rows(target_rows);
    if (static_cast<int>(final_rows.size()) != n) return false;
    check_square_rows(final_rows);
    output_permutation.assign(static_cast<std::size_t>(n), -1);
    std::vector<char> used(static_cast<std::size_t>(n), 0);
    for (int r = 0; r < n; ++r) {
        int match = -1;
        for (int t = 0; t < n; ++t) {
            if (!used[static_cast<std::size_t>(t)] && final_rows[r] == target_rows[t]) {
                match = t;
                break;
            }
        }
        if (match < 0) return false;
        output_permutation[static_cast<std::size_t>(r)] = match;
        used[static_cast<std::size_t>(match)] = 1;
    }
    return true;
}

static bool logical_depth_one_ops(const Matrix& rows, std::vector<RowOp>& ops) {
    const int n = check_square_rows(rows);
    Row occupied_weight2_columns = 0;
    ops.clear();
    for (int i = 0; i < n; ++i) {
        const int w = popcount64(rows[i]);
        if (w > 2 || w == 0) return false;
        if (w == 2) {
            if ((occupied_weight2_columns & rows[i]) != 0) return false;
            occupied_weight2_columns |= rows[i];
        }
    }

    for (int i = 0; i < n; ++i) {
        if (popcount64(rows[i]) != 2) continue;
        int control = -1;
        for (int t = 0; t < n; ++t) {
            if (i == t || popcount64(rows[t]) != 1) continue;
            if ((rows[i] & rows[t]) != 0) {
                control = t;
                break;
            }
        }
        if (control < 0) return false;
        ops.emplace_back(control, i);
    }
    return !ops.empty();
}

struct Layout {
    int n = 0;
    int data_rows = 0;
    int data_cols = 0;
    int grid_rows = 0;
    int grid_cols = 0;
    std::vector<int> data_pos;
    std::vector<char> occupied_data;

    static std::pair<int, int> balanced_shape(int n) {
        int best_rows = 1;
        int best_cols = n;
        std::tuple<int, int, int> best{n, n, n};
        for (int rows = 1; rows * rows <= n + rows; ++rows) {
            const int cols = (n + rows - 1) / rows;
            const int excess = rows * cols - n;
            const auto score = std::make_tuple(excess, std::abs(cols - rows), rows * cols);
            if (score < best) {
                best = score;
                best_rows = rows;
                best_cols = cols;
            }
        }
        return {best_rows, best_cols};
    }

    static Layout for_data_shape(int n, int dr, int dc) {
        if (n <= 0) throw std::runtime_error("layout needs a positive qubit count");
        if (dr <= 0 && dc <= 0) throw std::runtime_error("data_rows or data_cols must be positive");
        if (dr <= 0) dr = (n + dc - 1) / dc;
        if (dc <= 0) dc = (n + dr - 1) / dr;
        if (dr <= 0 || dc <= 0 || dr * dc < n) {
            throw std::runtime_error("forced data layout does not have enough data slots");
        }
        Layout layout;
        layout.n = n;
        layout.data_rows = dr;
        layout.data_cols = dc;
        layout.grid_rows = 2 * dr + 1;
        layout.grid_cols = 2 * dc + 1;
        const int vertices = layout.grid_rows * layout.grid_cols;
        layout.occupied_data.assign(vertices, 0);
        for (int r = 0; r < dr; ++r) {
            for (int c = 0; c < dc; ++c) {
                if (static_cast<int>(layout.data_pos.size()) >= n) break;
                const int id = layout.vertex_id(2 * r + 1, 2 * c + 1);
                layout.data_pos.push_back(id);
                layout.occupied_data[id] = 1;
            }
        }
        return layout;
    }

    static Layout for_qubits(int n) {
        auto [dr, dc] = balanced_shape(n);
        return for_data_shape(n, dr, dc);
    }

    int vertex_id(int r, int c) const { return r * grid_cols + c; }
    int row(int id) const { return id / grid_cols; }
    int col(int id) const { return id % grid_cols; }

    std::vector<int> neighbors(int id) const {
        const int r = row(id);
        const int c = col(id);
        std::vector<int> out;
        if (r > 0) out.push_back(vertex_id(r - 1, c));
        if (r + 1 < grid_rows) out.push_back(vertex_id(r + 1, c));
        if (c > 0) out.push_back(vertex_id(r, c - 1));
        if (c + 1 < grid_cols) out.push_back(vertex_id(r, c + 1));
        return out;
    }

    int manhattan_between_qubits(int a, int b) const {
        const int va = data_pos[a];
        const int vb = data_pos[b];
        return std::abs(row(va) - row(vb)) + std::abs(col(va) - col(vb));
    }

    static std::uint64_t edge_key(int a, int b) {
        if (a > b) std::swap(a, b);
        return (std::uint64_t(a) << 32U) | std::uint32_t(b);
    }

    std::vector<std::uint64_t> path_edges(const std::vector<int>& path) const {
        std::vector<std::uint64_t> edges;
        if (path.size() < 2) return edges;
        edges.reserve(path.size() - 1);
        for (std::size_t i = 0; i + 1 < path.size(); ++i) {
            edges.push_back(edge_key(path[i], path[i + 1]));
        }
        return edges;
    }

    bool are_vertex_disjoint(const std::vector<std::vector<int>>& paths) const {
        std::vector<char> seen(grid_rows * grid_cols, 0);
        for (const auto& path : paths) {
            for (int v : path) {
                if (seen[v]) return false;
                seen[v] = 1;
            }
        }
        return true;
    }

    std::vector<int> find_operator_path(
        int control,
        int target,
        const std::unordered_set<std::uint64_t>& banned_edges,
        const std::vector<char>* banned_vertices
    ) const {
        if (control == target) throw std::runtime_error("control and target must differ");
        const int start = data_pos[control];
        const int goal = data_pos[target];
        if (banned_vertices != nullptr && ((*banned_vertices)[start] || (*banned_vertices)[goal])) {
            return {};
        }

        std::vector<int> parent(grid_rows * grid_cols, -2);
        std::deque<int> q;
        parent[start] = -1;
        q.push_back(start);
        while (!q.empty()) {
            const int cur = q.front();
            q.pop_front();
            if (cur == goal) break;
            for (int nb : neighbors(cur)) {
                if (banned_edges.find(edge_key(cur, nb)) != banned_edges.end()) continue;
                if (parent[nb] != -2) continue;
                if (occupied_data[nb] && nb != start && nb != goal) continue;
                if (banned_vertices != nullptr && (*banned_vertices)[nb]) continue;
                if (cur == start && col(nb) != col(start)) continue;
                if (nb == goal && row(cur) != row(goal)) continue;
                parent[nb] = cur;
                q.push_back(nb);
            }
        }
        if (parent[goal] == -2) return {};
        std::vector<int> path;
        for (int cur = goal; cur != -1; cur = parent[cur]) path.push_back(cur);
        std::reverse(path.begin(), path.end());
        return path;
    }

    std::vector<std::vector<int>> find_operator_paths(
        int control,
        int target,
        const std::unordered_set<std::uint64_t>& banned_edges,
        const std::vector<char>* banned_vertices,
        int max_paths,
        int max_extra_edges
    ) const {
        std::vector<std::vector<int>> paths;
        if (max_paths <= 0) return paths;
        std::vector<int> first = find_operator_path(control, target, banned_edges, banned_vertices);
        if (first.empty()) return paths;
        paths.push_back(first);
        if (max_paths == 1) return paths;

        const int start = data_pos[control];
        const int goal = data_pos[target];
        const int min_edges = static_cast<int>(first.size()) - 1;
        const int max_edges = min_edges + std::max(0, max_extra_edges);
        std::vector<char> visited(static_cast<std::size_t>(grid_rows * grid_cols), 0);
        std::vector<int> path;
        path.reserve(static_cast<std::size_t>(max_edges + 1));
        path.push_back(start);
        visited[static_cast<std::size_t>(start)] = 1;

        auto heuristic = [&](int v) {
            return std::abs(row(v) - row(goal)) + std::abs(col(v) - col(goal));
        };
        auto have_path = [&](const std::vector<int>& candidate) {
            return std::any_of(paths.begin(), paths.end(), [&](const std::vector<int>& existing) {
                return existing == candidate;
            });
        };

        int expanded = 0;
        const int expansion_limit = std::max(4000, max_paths * 8000);
        int edge_limit = min_edges;
        std::function<void(int, int)> dfs = [&](int cur, int depth_edges) {
            if (static_cast<int>(paths.size()) >= max_paths) return;
            if (++expanded > expansion_limit) return;
            if (cur == goal) {
                if (depth_edges == edge_limit && !have_path(path)) paths.push_back(path);
                return;
            }
            if (depth_edges >= edge_limit) return;
            if (depth_edges + heuristic(cur) > edge_limit) return;

            std::vector<int> nbs = neighbors(cur);
            std::sort(nbs.begin(), nbs.end(), [&](int a, int b) {
                const int ha = heuristic(a);
                const int hb = heuristic(b);
                if (ha != hb) return ha < hb;
                return std::make_tuple(row(a), col(a)) < std::make_tuple(row(b), col(b));
            });

            for (int nb : nbs) {
                if (banned_edges.find(edge_key(cur, nb)) != banned_edges.end()) continue;
                if (visited[static_cast<std::size_t>(nb)]) continue;
                if (occupied_data[static_cast<std::size_t>(nb)] && nb != start && nb != goal) continue;
                if (banned_vertices != nullptr && (*banned_vertices)[static_cast<std::size_t>(nb)]) continue;
                if (cur == start && col(nb) != col(start)) continue;
                if (nb == goal && row(cur) != row(goal)) continue;
                if (depth_edges + 1 + heuristic(nb) > edge_limit) continue;

                visited[static_cast<std::size_t>(nb)] = 1;
                path.push_back(nb);
                dfs(nb, depth_edges + 1);
                path.pop_back();
                visited[static_cast<std::size_t>(nb)] = 0;
                if (static_cast<int>(paths.size()) >= max_paths || expanded > expansion_limit) return;
            }
        };
        // Grid paths have fixed parity. Iterative deepening by exact path
        // length keeps candidate sets nested when k or path-extra increases.
        for (edge_limit = min_edges; edge_limit <= max_edges; edge_limit += 2) {
            dfs(start, 0);
            if (static_cast<int>(paths.size()) >= max_paths || expanded > expansion_limit) break;
        }
        std::sort(paths.begin(), paths.end(), [](const auto& a, const auto& b) {
            if (a.size() != b.size()) return a.size() < b.size();
            return a < b;
        });
        if (static_cast<int>(paths.size()) > max_paths) paths.resize(static_cast<std::size_t>(max_paths));
        return paths;
    }
};

struct SearchWeights {
    double row_identity = 10.0;
    double column_identity = 7.0;
    double row_shape = 0.7;
    double column_shape = 0.5;
    double route_penalty = 0.03;
    double surface_depth_penalty = 0.80;
    bool allow_columns = true;
};

struct AlgebraicChoice {
    OpSide side = OpSide::Row;
    int control = 0;
    int target = 0;
    double delta = 0.0;
    double priority = 0.0;
};

struct GreedyCandidate {
    int control = 0;
    int target = 0;
    int gain = 0;
    double priority = 0.0;
};

enum class ShiFengCostKind {
    Square,
    Product
};

struct ShiFengSurfaceConfig {
    ShiFengCostKind cost_kind = ShiFengCostKind::Square;
    double route_penalty = 0.03;
    double surface_depth_penalty = 0.45;
    bool allow_columns = true;
    double logical_layer_penalty = 8.0;
    bool rank_global_cost = true;
    bool endpoint_balanced_rank = false;
    bool residual_demand_rank = false;
    bool layer_matrix_rank = false;
};

struct LogicalSearchParams {
    std::size_t max_search_targets = 4;
    std::uint64_t seeds_per_target = 64;
    std::size_t choice_window = 1;
    int depth_goal = 10;
    int cnot_goal = 131;
    bool use_layer_matrix_rank = false;
    bool use_layer_matrix_only = false;
    int layer_matrix_extra_orders = 0;
    int layer_matrix_order_variants = 1;
    int layer_matrix_candidate_multiplier = 1;
    int layer_matrix_beam_multiplier = 1;
    int layer_matrix_tail_cap = 0;
    int layer_matrix_layer_cap_variants = 1;
    double layer_matrix_nonimproving_slack = 0.0;
    bool layer_matrix_goal_only = false;
    bool layer_matrix_endpoint_first = false;
    // Reuse one candidate pool per beam state for tail checks and expansion.
    // This trades repeated randomized pools for substantially less work.
    bool layer_matrix_reuse_state_candidates = false;
    bool layer_matrix_vdp_hard = false;
    int logical_reschedule_orders = 0;
    std::uint64_t seed_offset = 0;
};

struct SurfaceFinalPackParams {
    int k_paths = 1;
    int path_extra = 0;
    int repack_limit = 0;
    int state_limit = 0;
    int conditional_state_limit = 0;

    bool enabled() const {
        return k_paths > 1 && repack_limit > 1;
    }
};

static LogicalSearchParams g_logical_search_params;

static SurfaceFinalPackParams g_surface_final_pack_params;

static const Layout* g_layer_matrix_surface_layout = nullptr;

static std::string g_layer_matrix_surface_mode;

static bool logical_depth_pressure_goal() {
    return g_logical_search_params.depth_goal > 0
        && g_logical_search_params.depth_goal <= 10;
}

static std::uint64_t effective_seeds_per_target(std::size_t config_count) {
    const std::uint64_t requested = std::max<std::uint64_t>(1, g_logical_search_params.seeds_per_target);
    (void)config_count;
    return requested;
}

struct RoutedCNOT {
    int control = 0;
    int target = 0;
    std::vector<int> path;
    int gain = 0;
};

struct ShiFengChoice {
    OpSide side = OpSide::Row;
    int control = 0;
    int target = 0;
    double cost_after = 0.0;
    double rank_cost_after = 0.0;
    double priority = 0.0;
    int layer_depth_delta = 0;
    int endpoint_goal_excess_after = 0;
    int endpoint_max_after = 0;
    long long endpoint_energy_after = 0;
    int residual_endpoint_goal_excess_after = 0;
    int residual_endpoint_max_after = 0;
    long long residual_endpoint_energy_after = 0;
    RoutedCNOT routed;
};

struct Layer {
    std::string mode;
    std::vector<RoutedCNOT> ops;
    int logical_depth = 0;
};

struct Result {
    int n = 0;
    std::string mode;
    Layout layout;
    Matrix target_rows;
    std::vector<int> output_permutation;
    std::vector<Layer> layers;
    Matrix final_rows;

    int cnot_count() const {
        int total = 0;
        for (const auto& layer : layers) total += static_cast<int>(layer.ops.size());
        return total;
    }

    int layer_count() const { return static_cast<int>(layers.size()); }

    int surface_depth() const {
        int total = 0;
        for (const auto& layer : layers) total += layer.logical_depth;
        return total;
    }

    bool verify() const {
        return final_rows == permute_rows(target_rows, output_permutation);
    }
};

static void validate_permutation_vector(const std::vector<int>& perm, int n) {
    if (static_cast<int>(perm.size()) != n) throw std::runtime_error("output permutation size mismatch");
    std::vector<char> seen(static_cast<std::size_t>(n), 0);
    for (int v : perm) {
        if (v < 0 || v >= n) throw std::runtime_error("output permutation entry out of range");
        if (seen[static_cast<std::size_t>(v)]) throw std::runtime_error("output permutation has duplicate entries");
        seen[static_cast<std::size_t>(v)] = 1;
    }
}

static bool logical_result_meets_goals(const Result& result) {
    const bool depth_ok = g_logical_search_params.depth_goal <= 0
        || result.layer_count() <= g_logical_search_params.depth_goal;
    const bool cnot_ok = g_logical_search_params.cnot_goal <= 0
        || result.cnot_count() <= g_logical_search_params.cnot_goal;
    return depth_ok && cnot_ok;
}

static bool try_route_into_layer(
    const Layout& layout,
    const std::string& mode,
    const Layer& layer,
    int control,
    int target,
    RoutedCNOT& routed
);

static int layer_depth(const std::string& mode, const Layout& layout, const std::vector<RoutedCNOT>& ops);

struct SurfaceScheduleState {
    const Layout* layout = nullptr;
    const std::string* mode = nullptr;
    std::vector<Layer> layers;
    std::vector<int> last_touch;
    int logical_depth = 0;

    SurfaceScheduleState(const Layout& layout_ref, const std::string& mode_ref)
        : layout(&layout_ref), mode(&mode_ref), last_touch(layout_ref.n, -1) {}

    int incremental_depth(int control, int target) const;
    void append(int control, int target);
};

static std::string matrix_key(const Matrix& m) {
    std::string out;
    out.reserve(m.size() * sizeof(Row));
    for (Row row : m) {
        for (int i = 0; i < 8; ++i) {
            out.push_back(static_cast<char>((row >> (8 * i)) & 0xffU));
        }
    }
    return out;
}

static double line_shape_cost(int weight) {
    if (weight <= 0) return 64.0;
    return static_cast<double>(weight * weight) + 4.0 * std::log2(static_cast<double>(weight));
}

static double matrix_cost(const Matrix& rows, const SearchWeights& weights) {
    const int n = check_square_rows(rows);
    double cost = 0.0;
    for (int i = 0; i < n; ++i) {
        const int d = popcount64(rows[i] ^ (Row(1) << i));
        const int w = popcount64(rows[i]);
        cost += weights.row_identity * static_cast<double>(d * d + d)
              + weights.row_shape * line_shape_cost(w);
    }
    for (int c = 0; c < n; ++c) {
        const Row col = column_bits(rows, c);
        const int d = popcount64(col ^ (Row(1) << c));
        const int w = popcount64(col);
        cost += weights.column_identity * static_cast<double>(d * d + d)
              + weights.column_shape * line_shape_cost(w);
    }
    return cost;
}

static Matrix transpose_matrix(const Matrix& rows) {
    const int n = check_square_rows(rows);
    Matrix out(n, 0);
    for (int r = 0; r < n; ++r) {
        Row bits = rows[r];
        while (bits != 0) {
            const Row lsb = bits & (~bits + 1);
            const int c = popcount64(lsb - 1);
            out[c] |= Row(1) << r;
            bits ^= lsb;
        }
    }
    return out;
}

static double h_square_cost(const Matrix& rows) {
    check_square_rows(rows);
    double cost = 0.0;
    for (Row row : rows) {
        const int w = popcount64(row);
        cost += static_cast<double>(w * w);
    }
    return cost;
}

static double h_product_cost(const Matrix& rows) {
    check_square_rows(rows);
    double cost = 0.0;
    for (Row row : rows) {
        const int w = popcount64(row);
        if (w <= 0) return std::numeric_limits<double>::infinity();
        cost += std::log2(static_cast<double>(w));
    }
    return cost;
}

static double shi_base_cost(const Matrix& rows, ShiFengCostKind kind) {
    return kind == ShiFengCostKind::Square ? h_square_cost(rows) : h_product_cost(rows);
}

static double shi_row_cost(const Matrix& rows, ShiFengCostKind kind) {
    const Matrix inv = invert_matrix(rows);
    return shi_base_cost(rows, kind) + shi_base_cost(transpose_matrix(inv), kind);
}

static double shi_column_cost(const Matrix& rows, ShiFengCostKind kind) {
    const Matrix inv = invert_matrix(rows);
    return shi_base_cost(transpose_matrix(rows), kind) + shi_base_cost(inv, kind);
}

static double shi_row_cost_with_inverse(const Matrix& rows, const Matrix& inverse, ShiFengCostKind kind) {
    return shi_base_cost(rows, kind) + shi_base_cost(transpose_matrix(inverse), kind);
}

static double shi_column_cost_with_inverse(const Matrix& rows, const Matrix& inverse, ShiFengCostKind kind) {
    return shi_base_cost(transpose_matrix(rows), kind) + shi_base_cost(inverse, kind);
}

static double shi_current_cost(const Matrix& rows, ShiFengCostKind kind) {
    return std::max(shi_row_cost(rows, kind), shi_column_cost(rows, kind));
}

static double shi_current_cost_with_inverse(const Matrix& rows, const Matrix& inverse, ShiFengCostKind kind) {
    return std::max(
        shi_row_cost_with_inverse(rows, inverse, kind),
        shi_column_cost_with_inverse(rows, inverse, kind)
    );
}

static double shi_side_cost_with_inverse(
    const Matrix& rows,
    const Matrix& inverse,
    OpSide side,
    ShiFengCostKind kind
) {
    return side == OpSide::Row
        ? shi_row_cost_with_inverse(rows, inverse, kind)
        : shi_column_cost_with_inverse(rows, inverse, kind);
}

static double shi_after_cost(const Matrix& rows, const AlgebraicChoice& choice, ShiFengCostKind kind) {
    Matrix after = rows;
    if (choice.side == OpSide::Row) {
        row_add(after, choice.control, choice.target);
    } else {
        column_add(after, choice.control, choice.target);
    }
    return choice.side == OpSide::Row ? shi_row_cost(after, kind) : shi_column_cost(after, kind);
}

struct ShiFengAfterCosts {
    double side_after = 0.0;
    double global_after = 0.0;
};

static double shi_after_side_cost_with_inverse(
    const Matrix& rows,
    const Matrix& inverse,
    const AlgebraicChoice& choice,
    ShiFengCostKind kind
) {
    Matrix after = rows;
    Matrix after_inverse = inverse;
    if (choice.side == OpSide::Row) {
        row_add(after, choice.control, choice.target);
        column_add(after_inverse, choice.control, choice.target);
    } else {
        column_add(after, choice.control, choice.target);
        row_add(after_inverse, choice.control, choice.target);
    }
    return shi_side_cost_with_inverse(after, after_inverse, choice.side, kind);
}

static ShiFengAfterCosts shi_after_costs_with_inverse(
    const Matrix& rows,
    const Matrix& inverse,
    const AlgebraicChoice& choice,
    ShiFengCostKind kind
) {
    Matrix after = rows;
    Matrix after_inverse = inverse;
    if (choice.side == OpSide::Row) {
        row_add(after, choice.control, choice.target);
        column_add(after_inverse, choice.control, choice.target);
    } else {
        column_add(after, choice.control, choice.target);
        row_add(after_inverse, choice.control, choice.target);
    }
    return {
        shi_side_cost_with_inverse(after, after_inverse, choice.side, kind),
        shi_current_cost_with_inverse(after, after_inverse, kind)
    };
}

static void apply_algebraic_choice(Matrix& rows, const AlgebraicChoice& choice) {
    if (choice.side == OpSide::Row) {
        row_add(rows, choice.control, choice.target);
    } else {
        column_add(rows, choice.control, choice.target);
    }
}

static void apply_algebraic_choice_with_inverse(
    Matrix& rows,
    Matrix& inverse,
    const AlgebraicChoice& choice
) {
    if (choice.side == OpSide::Row) {
        row_add(rows, choice.control, choice.target);
        column_add(inverse, choice.control, choice.target);
    } else {
        column_add(rows, choice.control, choice.target);
        row_add(inverse, choice.control, choice.target);
    }
}

static Matrix after_algebraic_choice(const Matrix& rows, const AlgebraicChoice& choice) {
    Matrix out = rows;
    apply_algebraic_choice(out, choice);
    return out;
}

static bool cnot_commute(const RowOp& a, const RowOp& b) {
    const int ac = a.first;
    const int at = a.second;
    const int bc = b.first;
    const int bt = b.second;
    return at != bc && bt != ac;
}

static std::vector<RowOp> optimize_cnot_sequence(const std::vector<RowOp>& ops) {
    std::vector<RowOp> out;
    out.reserve(ops.size());
    constexpr int window = 32;
    for (const RowOp& op : ops) {
        bool cancelled = false;
        const int start = std::max(0, static_cast<int>(out.size()) - window);
        for (int i = static_cast<int>(out.size()) - 1; i >= start; --i) {
            if (!cnot_commute(op, out[i])) break;
            if (out[i] == op) {
                out.erase(out.begin() + i);
                cancelled = true;
                break;
            }
        }
        if (!cancelled) out.push_back(op);
    }
    return out;
}

struct LocalReductionDatabase {
    int k = 0;
    int max_len = 0;
    std::unordered_map<std::string, std::vector<std::vector<RowOp>>> by_matrix;
    std::vector<RowOp> generators;

    LocalReductionDatabase(int local_qubits, int max_sequence_len)
        : k(local_qubits), max_len(max_sequence_len) {
        if (k < 2 || k > 6) throw std::runtime_error("local reduction arity out of range");
        generators.reserve(static_cast<std::size_t>(k * (k - 1)));
        for (int control = 0; control < k; ++control) {
            for (int target = 0; target < k; ++target) {
                if (control != target) generators.emplace_back(control, target);
            }
        }
        Matrix current = identity_matrix(k);
        std::vector<RowOp> seq;
        enumerate(current, seq, 0);
        for (auto& item : by_matrix) {
            auto& seqs = item.second;
            std::sort(seqs.begin(), seqs.end(), [](const std::vector<RowOp>& a, const std::vector<RowOp>& b) {
                if (a.size() != b.size()) return a.size() < b.size();
                return a < b;
            });
            seqs.erase(std::unique(seqs.begin(), seqs.end()), seqs.end());
        }
    }

    void enumerate(const Matrix& current, std::vector<RowOp>& seq, int depth) {
        by_matrix[matrix_key(current)].push_back(seq);
        if (depth >= max_len) return;
        for (const RowOp& op : generators) {
            if (!seq.empty() && seq.back() == op) continue;
            Matrix next = current;
            row_add(next, op.first, op.second);
            seq.push_back(op);
            enumerate(next, seq, depth + 1);
            seq.pop_back();
        }
    }

    const std::vector<std::vector<RowOp>>& equivalents(const Matrix& matrix) const {
        static const std::vector<std::vector<RowOp>> empty;
        const auto it = by_matrix.find(matrix_key(matrix));
        return it == by_matrix.end() ? empty : it->second;
    }
};

static const LocalReductionDatabase& local_reduction_db(int k, int max_len) {
    static std::mutex cache_mutex;
    static std::unordered_map<int, std::shared_ptr<const LocalReductionDatabase>> cache;
    static thread_local std::unordered_map<int, const LocalReductionDatabase*> thread_cache;
    const int key = 100 * k + max_len;
    const auto local = thread_cache.find(key);
    if (local != thread_cache.end()) return *local->second;

    const LocalReductionDatabase* database = nullptr;
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        auto it = cache.find(key);
        if (it == cache.end()) {
            it = cache.emplace(
                key,
                std::make_shared<const LocalReductionDatabase>(k, max_len)
            ).first;
        }
        database = it->second.get();
    }
    thread_cache.emplace(key, database);
    return *database;
}

static std::vector<AlgebraicChoice> algebraic_candidates(
    const Matrix& residual,
    const Layout& layout,
    const SearchWeights& weights,
    bool allow_lookahead,
    const SurfaceScheduleState* row_schedule,
    const SurfaceScheduleState* column_schedule
) {
    consume_fresh_work();
    const int n = static_cast<int>(residual.size());
    const double before = matrix_cost(residual, weights);
    std::vector<AlgebraicChoice> direct;
    direct.reserve(static_cast<std::size_t>(2 * n * (n - 1)));

    auto surface_penalty = [&](OpSide side, int control, int target) {
        (void)side;
        (void)column_schedule;
        if (row_schedule == nullptr) return 0.0;
        const int depth_delta = row_schedule->incremental_depth(control, target);
        if (depth_delta <= 0) return 0.0;
        return weights.surface_depth_penalty * static_cast<double>(depth_delta);
    };

    auto consider = [&](OpSide side, int control, int target, double boost) {
        AlgebraicChoice choice{side, control, target, 0.0, 0.0};
        Matrix after = after_algebraic_choice(residual, choice);
        choice.delta = before - matrix_cost(after, weights);
        if (choice.delta <= 1.0e-9) return;
        const double route = weights.route_penalty * layout.manhattan_between_qubits(control, target);
        const double parallel = surface_penalty(side, control, target);
        const double side_bias = (side == OpSide::Column) ? 0.001 : 0.0;
        choice.priority = choice.delta + boost + side_bias - route - parallel;
        direct.push_back(choice);
    };

    for (int control = 0; control < n; ++control) {
        for (int target = 0; target < n; ++target) {
            if (control == target) continue;
            consider(OpSide::Row, control, target, 0.0);
            if (weights.allow_columns) consider(OpSide::Column, control, target, 0.0);
        }
    }

    if (!direct.empty() || !allow_lookahead) {
        std::sort(direct.begin(), direct.end(), [](const AlgebraicChoice& a, const AlgebraicChoice& b) {
            if (a.priority != b.priority) return a.priority > b.priority;
            if (a.delta != b.delta) return a.delta > b.delta;
            if (a.control != b.control) return a.control < b.control;
            return a.target < b.target;
        });
        return direct;
    }

    if (n > 24) return {};

    struct FirstStep {
        AlgebraicChoice choice;
        double cost_after = 0.0;
        double rank = 0.0;
    };
    std::vector<FirstStep> pool;
    pool.reserve(static_cast<std::size_t>(2 * n * (n - 1)));
    const int side_count = weights.allow_columns ? 2 : 1;
    for (int first_side = 0; first_side < side_count; ++first_side) {
        for (int control = 0; control < n; ++control) {
            for (int target = 0; target < n; ++target) {
                if (control == target) continue;
                AlgebraicChoice first{
                    first_side == 0 ? OpSide::Row : OpSide::Column,
                    control,
                    target,
                    0.0,
                    0.0
                };
                Matrix one = after_algebraic_choice(residual, first);
                const double cost_after = matrix_cost(one, weights);
                const double route = weights.route_penalty * layout.manhattan_between_qubits(control, target);
                const double parallel = surface_penalty(first.side, control, target);
                pool.push_back({first, cost_after, cost_after + route + parallel});
            }
        }
    }

    std::sort(pool.begin(), pool.end(), [](const FirstStep& a, const FirstStep& b) {
        if (a.rank != b.rank) return a.rank < b.rank;
        if (a.choice.control != b.choice.control) return a.choice.control < b.choice.control;
        return a.choice.target < b.choice.target;
    });
    const int lookahead_limit = std::min(static_cast<int>(pool.size()), std::max(48, 6 * n));

    std::vector<AlgebraicChoice> lookahead;
    lookahead.reserve(static_cast<std::size_t>(lookahead_limit));
    for (int i = 0; i < lookahead_limit; ++i) {
        const AlgebraicChoice first = pool[i].choice;
        Matrix one = after_algebraic_choice(residual, first);
        double best_second = pool[i].cost_after;
                for (int second_side = 0; second_side < side_count; ++second_side) {
                    for (int c2 = 0; c2 < n; ++c2) {
                        for (int t2 = 0; t2 < n; ++t2) {
                            if (c2 == t2) continue;
                            AlgebraicChoice second{
                                second_side == 0 ? OpSide::Row : OpSide::Column,
                                c2,
                                t2,
                                0.0,
                                0.0
                            };
                            Matrix two = after_algebraic_choice(one, second);
                            best_second = std::min(best_second, matrix_cost(two, weights));
                        }
                    }
                }
        const double combined_delta = before - best_second;
        if (combined_delta <= 1.0e-9) continue;
        AlgebraicChoice choice = first;
        const double route = weights.route_penalty * layout.manhattan_between_qubits(choice.control, choice.target);
        const double parallel = surface_penalty(choice.side, choice.control, choice.target);
        choice.delta = combined_delta;
        choice.priority = combined_delta - route - parallel - 0.25;
        lookahead.push_back(choice);
    }
    std::sort(lookahead.begin(), lookahead.end(), [](const AlgebraicChoice& a, const AlgebraicChoice& b) {
        if (a.priority != b.priority) return a.priority > b.priority;
        if (a.delta != b.delta) return a.delta > b.delta;
        if (a.control != b.control) return a.control < b.control;
        return a.target < b.target;
    });
    return lookahead;
}

static std::vector<GreedyCandidate> greedy_candidate_list(const Matrix& residual, const Layout& layout) {
    consume_fresh_work();
    const int n = static_cast<int>(residual.size());
    std::vector<int> row_cost(n);
    for (int i = 0; i < n; ++i) row_cost[i] = popcount64(residual[i] ^ (Row(1) << i));

    std::vector<GreedyCandidate> candidates;
    candidates.reserve(static_cast<std::size_t>(n * (n - 1)));
    for (int control = 0; control < n; ++control) {
        const Row c_row = residual[control];
        for (int target = 0; target < n; ++target) {
            if (control == target) continue;
            const Row after = residual[target] ^ c_row;
            const int after_cost = popcount64(after ^ (Row(1) << target));
            const int gain = row_cost[target] - after_cost;
            if (gain <= 0) continue;
            const double distance_penalty = 0.05 * layout.manhattan_between_qubits(control, target);
            const double overlap_bonus = 0.01 * popcount64(residual[target] & c_row);
            candidates.push_back({control, target, gain, 100.0 * gain + overlap_bonus - distance_penalty});
        }
    }
    std::sort(candidates.begin(), candidates.end(), [](const GreedyCandidate& a, const GreedyCandidate& b) {
        return a.priority > b.priority;
    });
    return candidates;
}

static bool layer_paths_vertex_disjoint(const Layout& layout, const std::vector<RoutedCNOT>& ops) {
    std::vector<std::vector<int>> paths;
    paths.reserve(ops.size());
    for (const auto& op : ops) paths.push_back(op.path);
    return layout.are_vertex_disjoint(paths);
}

static int layer_depth(const std::string& mode, const Layout& layout, const std::vector<RoutedCNOT>& ops) {
    if (mode == "logical") {
        (void)layout;
        (void)ops;
        return 1;
    }
    if (mode == "vdp") return 2;
    return layer_paths_vertex_disjoint(layout, ops) ? 2 : 4;
}

static std::vector<RoutedCNOT> greedy_pack_layer(
    const Layout& layout,
    const std::string& mode,
    const std::vector<GreedyCandidate>& candidates,
    int max_ops = std::numeric_limits<int>::max()
) {
    std::vector<char> used_qubits(layout.n, 0);
    std::vector<char> used_vertices(layout.grid_rows * layout.grid_cols, 0);
    std::unordered_set<std::uint64_t> used_edges;
    std::vector<RoutedCNOT> routed;

    for (const GreedyCandidate& cand : candidates) {
        if (used_qubits[cand.control] || used_qubits[cand.target]) continue;
        std::vector<int> path;
        if (mode != "logical") {
            const std::vector<char>* banned_vertices =
                (mode == "vdp" || mode == "edp") ? &used_vertices : nullptr;
            path = layout.find_operator_path(cand.control, cand.target, used_edges, banned_vertices);
            if (path.empty() && mode == "edp") {
                path = layout.find_operator_path(cand.control, cand.target, used_edges, nullptr);
            }
            if (path.empty()) continue;
        }

        routed.push_back({cand.control, cand.target, std::move(path), cand.gain});
        used_qubits[cand.control] = 1;
        used_qubits[cand.target] = 1;
        for (auto e : layout.path_edges(routed.back().path)) used_edges.insert(e);
        if (mode == "vdp" || mode == "edp") {
            for (int v : routed.back().path) used_vertices[v] = 1;
        }
        if (static_cast<int>(routed.size()) >= max_ops) break;
    }
    return routed;
}

static bool try_route_into_layer(
    const Layout& layout,
    const std::string& mode,
    const Layer& layer,
    int control,
    int target,
    RoutedCNOT& routed
) {
    consume_fresh_work();
    std::vector<char> used_qubits(layout.n, 0);
    std::vector<char> used_vertices(layout.grid_rows * layout.grid_cols, 0);
    std::unordered_set<std::uint64_t> used_edges;
    for (const auto& op : layer.ops) {
        used_qubits[op.control] = 1;
        used_qubits[op.target] = 1;
        for (auto e : layout.path_edges(op.path)) used_edges.insert(e);
        if (mode == "vdp" || mode == "edp") {
            for (int v : op.path) used_vertices[v] = 1;
        }
    }
    if (used_qubits[control] || used_qubits[target]) return false;
    if (mode == "logical") {
        routed = {control, target, {}, 0};
        return true;
    }
    const std::vector<char>* banned_vertices = (mode == "vdp" || mode == "edp") ? &used_vertices : nullptr;
    std::vector<int> path = layout.find_operator_path(control, target, used_edges, banned_vertices);
    if (path.empty() && mode == "edp") {
        path = layout.find_operator_path(control, target, used_edges, nullptr);
    }
    if (path.empty()) return false;
    routed = {control, target, std::move(path), 0};
    return true;
}

static bool paths_conflict(
    const Layout& layout,
    const std::string& mode,
    const std::vector<int>& path,
    const std::unordered_set<std::uint64_t>& used_edges,
    const std::vector<char>& used_vertices
) {
    for (auto edge : layout.path_edges(path)) {
        if (used_edges.find(edge) != used_edges.end()) return true;
    }
    if (mode == "vdp") {
        for (int v : path) {
            if (used_vertices[static_cast<std::size_t>(v)]) return true;
        }
    }
    return false;
}

static bool pack_layer_ops_with_k_paths(
    const Layout& layout,
    const std::string& mode,
    const std::vector<RowOp>& ops,
    std::vector<RoutedCNOT>& routed_ops
) {
    consume_fresh_work();
    if (mode == "logical") {
        routed_ops.clear();
        routed_ops.reserve(ops.size());
        for (const auto& op : ops) routed_ops.push_back({op.first, op.second, {}, 0});
        return true;
    }
    if (static_cast<int>(ops.size()) > g_surface_final_pack_params.repack_limit) return false;
    if (mode == "edp") {
        // A vertex-disjoint EDP layer costs two cycles instead of four.  Search
        // that strict subset first, then fall back to general edge-disjoint paths.
        std::vector<RoutedCNOT> vertex_disjoint;
        if (pack_layer_ops_with_k_paths(layout, "vdp", ops, vertex_disjoint)) {
            routed_ops = std::move(vertex_disjoint);
            return true;
        }
    }
    std::vector<char> used_qubits(static_cast<std::size_t>(layout.n), 0);
    for (const auto& op : ops) {
        if (op.first == op.second) return false;
        if (used_qubits[static_cast<std::size_t>(op.first)]
            || used_qubits[static_cast<std::size_t>(op.second)]) {
            return false;
        }
        used_qubits[static_cast<std::size_t>(op.first)] = 1;
        used_qubits[static_cast<std::size_t>(op.second)] = 1;
    }

    const std::unordered_set<std::uint64_t> no_edges;
    std::vector<std::vector<std::vector<int>>> candidates(ops.size());
    static std::mutex path_cache_mutex;
    static std::unordered_map<std::string, std::vector<std::vector<int>>> path_cache;
    for (std::size_t i = 0; i < ops.size(); ++i) {
        std::ostringstream key_stream;
        key_stream << layout.n << '|' << layout.data_rows << 'x' << layout.data_cols
                   << '|' << ops[i].first << '>' << ops[i].second
                   << "|k=" << g_surface_final_pack_params.k_paths
                   << "|extra=" << g_surface_final_pack_params.path_extra;
        const std::string key = key_stream.str();
        {
            std::lock_guard<std::mutex> lock(path_cache_mutex);
            const auto it = path_cache.find(key);
            if (it != path_cache.end()) candidates[i] = it->second;
        }
        if (candidates[i].empty()) {
            candidates[i] = layout.find_operator_paths(
                ops[i].first,
                ops[i].second,
                no_edges,
                nullptr,
                g_surface_final_pack_params.k_paths,
                g_surface_final_pack_params.path_extra
            );
            std::lock_guard<std::mutex> lock(path_cache_mutex);
            path_cache.emplace(key, candidates[i]);
        }
        if (candidates[i].empty()) return false;
    }

    std::vector<std::vector<int>> chosen(ops.size());
    std::unordered_set<std::uint64_t> used_edges;
    std::vector<char> used_vertices(static_cast<std::size_t>(layout.grid_rows * layout.grid_cols), 0);
    bool found = false;
    int visited_states = 0;
    const int state_limit = g_surface_final_pack_params.state_limit > 0
        ? g_surface_final_pack_params.state_limit
        : std::max(20000, g_surface_final_pack_params.k_paths * 50000);
    Row remaining = ops.size() == 64
        ? ~Row(0)
        : (Row(1) << static_cast<int>(ops.size())) - 1;

    std::function<void(Row)> dfs = [&](Row unresolved) {
        if (found) return;
        if (++visited_states > state_limit) return;
        if (unresolved == 0) {
            found = true;
            return;
        }

        int op_idx = -1;
        std::vector<const std::vector<int>*> compatible;
        for (std::size_t candidate_index = 0; candidate_index < ops.size(); ++candidate_index) {
            if ((unresolved & (Row(1) << candidate_index)) == 0) continue;
            std::vector<const std::vector<int>*> available;
            for (const auto& path : candidates[candidate_index]) {
                if (!paths_conflict(layout, mode, path, used_edges, used_vertices)) {
                    available.push_back(&path);
                }
            }
            if (available.empty()) return;
            if (op_idx < 0
                || available.size() < compatible.size()
                || (available.size() == compatible.size()
                    && candidates[candidate_index].front().size()
                       > candidates[static_cast<std::size_t>(op_idx)].front().size())) {
                op_idx = static_cast<int>(candidate_index);
                compatible = std::move(available);
            }
        }
        std::sort(compatible.begin(), compatible.end(), [](const auto* left, const auto* right) {
            if (left->size() != right->size()) return left->size() < right->size();
            return *left < *right;
        });
        for (const auto* path_ptr : compatible) {
            const auto& path = *path_ptr;
            std::vector<std::uint64_t> added_edges = layout.path_edges(path);
            for (auto edge : added_edges) used_edges.insert(edge);
            if (mode == "vdp") {
                for (int v : path) used_vertices[static_cast<std::size_t>(v)] = 1;
            }
            chosen[static_cast<std::size_t>(op_idx)] = path;
            dfs(unresolved & ~(Row(1) << op_idx));
            if (found) return;
            chosen[static_cast<std::size_t>(op_idx)].clear();
            if (mode == "vdp") {
                for (int v : path) used_vertices[static_cast<std::size_t>(v)] = 0;
            }
            for (auto edge : added_edges) used_edges.erase(edge);
        }
    };
    dfs(remaining);

    if (!found && mode == "vdp" && ops.size() <= 8) {
        // A dynamically routed shortest path can still find a detour absent
        // from the bounded static path lists. Limit this fallback because this
        // routine is also called inside the logical beam-search hot loop.
        int conditional_states = 0;
        const int conditional_state_limit =
            g_surface_final_pack_params.conditional_state_limit > 0
            ? g_surface_final_pack_params.conditional_state_limit
            : 50000;
        std::vector<RoutedCNOT> conditional(ops.size());
        std::function<void(Row, const Layer&)> conditional_dfs =
            [&](Row unresolved, const Layer& layer) {
                if (found || ++conditional_states > conditional_state_limit) return;
                if (unresolved == 0) {
                    routed_ops = conditional;
                    found = true;
                    return;
                }
                std::vector<std::pair<int, RoutedCNOT>> choices;
                for (std::size_t op_index = 0; op_index < ops.size(); ++op_index) {
                    if ((unresolved & (Row(1) << op_index)) == 0) continue;
                    RoutedCNOT routed;
                    if (try_route_into_layer(
                            layout,
                            mode,
                            layer,
                            ops[op_index].first,
                            ops[op_index].second,
                            routed
                        )) {
                        choices.emplace_back(static_cast<int>(op_index), std::move(routed));
                    }
                }
                std::sort(choices.begin(), choices.end(), [](const auto& left, const auto& right) {
                    if (left.second.path.size() != right.second.path.size()) {
                        return left.second.path.size() > right.second.path.size();
                    }
                    return left.first < right.first;
                });
                for (auto& choice : choices) {
                    Layer next = layer;
                    next.ops.push_back(choice.second);
                    conditional[static_cast<std::size_t>(choice.first)] = choice.second;
                    conditional_dfs(
                        unresolved & ~(Row(1) << choice.first),
                        next
                    );
                    if (found) return;
                }
            };
        Layer empty_layer;
        empty_layer.mode = mode;
        conditional_dfs(remaining, empty_layer);
    }
    if (!found) return false;
    if (!routed_ops.empty()) return true;

    routed_ops.clear();
    routed_ops.reserve(ops.size());
    for (std::size_t i = 0; i < ops.size(); ++i) {
        routed_ops.push_back({ops[i].first, ops[i].second, std::move(chosen[i]), 0});
    }
    return true;
}

static bool repack_layer_with_op(
    const Layout& layout,
    const std::string& mode,
    const Layer& layer,
    int control,
    int target,
    Layer& replacement
) {
    if (!g_surface_final_pack_params.enabled()) return false;
    std::vector<RowOp> ops;
    ops.reserve(layer.ops.size() + 1);
    for (const auto& op : layer.ops) ops.emplace_back(op.control, op.target);
    ops.emplace_back(control, target);

    std::vector<RoutedCNOT> routed_ops;
    if (!pack_layer_ops_with_k_paths(layout, mode, ops, routed_ops)) return false;
    replacement.mode = mode;
    replacement.ops = std::move(routed_ops);
    replacement.logical_depth = layer_depth(mode, layout, replacement.ops);
    return true;
}

static bool try_route_into_search_node(
    const Layout& layout,
    const std::string& mode,
    const Layer& current_layer,
    int control,
    int target,
    RoutedCNOT& routed,
    int& depth_delta
) {
    if (!try_route_into_layer(layout, mode, current_layer, control, target, routed)) return false;
    std::vector<RoutedCNOT> ops = current_layer.ops;
    ops.push_back(routed);
    const int new_depth = layer_depth(mode, layout, ops);
    depth_delta = std::max(0, new_depth - current_layer.logical_depth);
    return true;
}

static void append_to_search_node(
    Layer& current_layer,
    const Layout& layout,
    const std::string& mode,
    RoutedCNOT routed
) {
    current_layer.mode = mode;
    current_layer.ops.push_back(std::move(routed));
    current_layer.logical_depth = layer_depth(mode, layout, current_layer.ops);
}

int SurfaceScheduleState::incremental_depth(int control, int target) const {
    const int earliest = std::max(last_touch[control], last_touch[target]) + 1;
    int best_delta = std::numeric_limits<int>::max() / 4;
    for (int li = earliest; li < static_cast<int>(layers.size()); ++li) {
        RoutedCNOT routed;
        if (!try_route_into_layer(*layout, *mode, layers[li], control, target, routed)) continue;
        std::vector<RoutedCNOT> ops = layers[li].ops;
        ops.push_back(std::move(routed));
        const int new_layer_depth = layer_depth(*mode, *layout, ops);
        const int delta = std::max(0, new_layer_depth - layers[li].logical_depth);
        best_delta = std::min(best_delta, delta);
        if (best_delta == 0) return 0;
    }

    Layer layer;
    layer.mode = *mode;
    RoutedCNOT routed;
    if (!try_route_into_layer(*layout, *mode, layer, control, target, routed)) {
        return std::numeric_limits<int>::max() / 4;
    }
    layer.ops.push_back(std::move(routed));
    return std::min(best_delta, layer_depth(*mode, *layout, layer.ops));
}

void SurfaceScheduleState::append(int control, int target) {
    const int earliest = std::max(last_touch[control], last_touch[target]) + 1;
    int best_li = -1;
    int best_delta = std::numeric_limits<int>::max();
    int best_new_depth = std::numeric_limits<int>::max();
    int best_path_len = std::numeric_limits<int>::max();
    RoutedCNOT best_routed;
    for (int li = earliest; li < static_cast<int>(layers.size()); ++li) {
        RoutedCNOT routed;
        if (!try_route_into_layer(*layout, *mode, layers[li], control, target, routed)) continue;
        std::vector<RoutedCNOT> ops = layers[li].ops;
        ops.push_back(routed);
        const int new_depth = layer_depth(*mode, *layout, ops);
        const int delta = std::max(0, new_depth - layers[li].logical_depth);
        const int path_len = static_cast<int>(routed.path.size());
        const auto key = std::make_tuple(delta, li, new_depth, path_len);
        const auto best_key = std::make_tuple(best_delta, best_li, best_new_depth, best_path_len);
        if (best_li < 0 || key < best_key) {
            best_li = li;
            best_delta = delta;
            best_new_depth = new_depth;
            best_path_len = path_len;
            best_routed = std::move(routed);
        }
    }

    if (best_li >= 0) {
        const int old_depth = layers[best_li].logical_depth;
        layers[best_li].ops.push_back(std::move(best_routed));
        layers[best_li].logical_depth = layer_depth(*mode, *layout, layers[best_li].ops);
        logical_depth += layers[best_li].logical_depth - old_depth;
        last_touch[control] = best_li;
        last_touch[target] = best_li;
        return;
    }

    Layer layer;
    layer.mode = *mode;
    RoutedCNOT routed;
    if (!try_route_into_layer(*layout, *mode, layer, control, target, routed)) {
        throw std::runtime_error("could not route surface schedule hint CNOT");
    }
    layer.ops.push_back(std::move(routed));
    layer.logical_depth = layer_depth(*mode, *layout, layer.ops);
    logical_depth += layer.logical_depth;
    layers.push_back(std::move(layer));
    const int li = static_cast<int>(layers.size()) - 1;
    last_touch[control] = li;
    last_touch[target] = li;
}

static bool cnots_commute(const RowOp& a, const RowOp& b) {
    return a.first != b.second && b.first != a.second;
}

struct LogicalSequenceBounds {
    int capacity_bound = 0;
    int endpoint_bound = 0;
    int critical_bound = 0;
    int lower_bound = 0;
};

static void build_logical_precedence(
    const std::vector<RowOp>& ops,
    std::vector<std::vector<int>>& succ,
    std::vector<int>& pred_count
) {
    const int m = static_cast<int>(ops.size());
    succ.assign(static_cast<std::size_t>(m), {});
    pred_count.assign(static_cast<std::size_t>(m), 0);
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            if (cnots_commute(ops[static_cast<std::size_t>(i)], ops[static_cast<std::size_t>(j)])) continue;
            succ[static_cast<std::size_t>(i)].push_back(j);
            ++pred_count[static_cast<std::size_t>(j)];
        }
    }
}

static std::vector<int> logical_critical_lengths(const std::vector<std::vector<int>>& succ) {
    const int m = static_cast<int>(succ.size());
    std::vector<int> critical(static_cast<std::size_t>(m), 1);
    for (int i = m - 1; i >= 0; --i) {
        for (int j : succ[static_cast<std::size_t>(i)]) {
            critical[static_cast<std::size_t>(i)] = std::max(
                critical[static_cast<std::size_t>(i)],
                1 + critical[static_cast<std::size_t>(j)]
            );
        }
    }
    return critical;
}

static LogicalSequenceBounds logical_sequence_bounds(const std::vector<RowOp>& ops, int n) {
    if (n <= 1) throw std::runtime_error("logical bounds need at least two qubits");
    std::vector<std::vector<int>> succ;
    std::vector<int> pred_count;
    build_logical_precedence(ops, succ, pred_count);
    const std::vector<int> critical = logical_critical_lengths(succ);

    std::vector<int> endpoint_count(static_cast<std::size_t>(n), 0);
    for (const auto& op : ops) {
        ++endpoint_count[static_cast<std::size_t>(op.first)];
        ++endpoint_count[static_cast<std::size_t>(op.second)];
    }

    LogicalSequenceBounds bounds;
    bounds.capacity_bound = static_cast<int>((ops.size() + static_cast<std::size_t>(n / 2) - 1)
        / static_cast<std::size_t>(n / 2));
    bounds.endpoint_bound = endpoint_count.empty()
        ? 0
        : *std::max_element(endpoint_count.begin(), endpoint_count.end());
    bounds.critical_bound = critical.empty()
        ? 0
        : *std::max_element(critical.begin(), critical.end());
    bounds.lower_bound = std::max(bounds.capacity_bound, std::max(bounds.endpoint_bound, bounds.critical_bound));
    return bounds;
}

static std::vector<Layer> schedule_logical_commuting_ops_once(
    const std::vector<RowOp>& ops,
    const Layout& layout,
    int order_kind
) {
    const int m = static_cast<int>(ops.size());
    std::vector<std::vector<int>> succ;
    std::vector<int> pred_count;
    build_logical_precedence(ops, succ, pred_count);
    const std::vector<int> critical = logical_critical_lengths(succ);

    std::vector<char> scheduled(static_cast<std::size_t>(m), 0);
    int remaining = m;
    std::vector<Layer> layers;
    while (remaining > 0) {
        std::vector<int> endpoint_remaining(static_cast<std::size_t>(layout.n), 0);
        for (int i = 0; i < m; ++i) {
            if (scheduled[static_cast<std::size_t>(i)]) continue;
            const auto& op = ops[static_cast<std::size_t>(i)];
            ++endpoint_remaining[static_cast<std::size_t>(op.first)];
            ++endpoint_remaining[static_cast<std::size_t>(op.second)];
        }

        std::vector<int> ready;
        ready.reserve(static_cast<std::size_t>(remaining));
        for (int i = 0; i < m; ++i) {
            if (!scheduled[static_cast<std::size_t>(i)] && pred_count[static_cast<std::size_t>(i)] == 0) {
                ready.push_back(i);
            }
        }
        if (ready.empty()) throw std::runtime_error("logical CNOT precedence graph has a cycle");

        auto endpoint_sum = [&](int idx) {
            const auto& op = ops[static_cast<std::size_t>(idx)];
            return endpoint_remaining[static_cast<std::size_t>(op.first)]
                 + endpoint_remaining[static_cast<std::size_t>(op.second)];
        };
        auto endpoint_max = [&](int idx) {
            const auto& op = ops[static_cast<std::size_t>(idx)];
            return std::max(
                endpoint_remaining[static_cast<std::size_t>(op.first)],
                endpoint_remaining[static_cast<std::size_t>(op.second)]
            );
        };
        auto pseudo_random_key = [&](int idx) {
            std::uint64_t x = 0x9e3779b97f4a7c15ULL
                ^ (static_cast<std::uint64_t>(order_kind) * 0xbf58476d1ce4e5b9ULL)
                ^ (static_cast<std::uint64_t>(idx) * 0x94d049bb133111ebULL);
            x ^= x >> 30;
            x *= 0xbf58476d1ce4e5b9ULL;
            x ^= x >> 27;
            x *= 0x94d049bb133111ebULL;
            x ^= x >> 31;
            return x;
        };

        std::sort(ready.begin(), ready.end(), [&](int a, int b) {
            if (order_kind == 1) {
                if (endpoint_sum(a) != endpoint_sum(b)) return endpoint_sum(a) > endpoint_sum(b);
            } else if (order_kind == 2) {
                if (endpoint_max(a) != endpoint_max(b)) return endpoint_max(a) > endpoint_max(b);
            } else if (order_kind == 3) {
                if (succ[static_cast<std::size_t>(a)].size() != succ[static_cast<std::size_t>(b)].size()) {
                    return succ[static_cast<std::size_t>(a)].size() > succ[static_cast<std::size_t>(b)].size();
                }
            } else if (order_kind == 4) {
                if (endpoint_sum(a) != endpoint_sum(b)) return endpoint_sum(a) < endpoint_sum(b);
            } else if (order_kind == 5) {
                return a > b;
            } else if (order_kind == 6) {
                const auto& oa = ops[static_cast<std::size_t>(a)];
                const auto& ob = ops[static_cast<std::size_t>(b)];
                const auto ak = std::make_tuple(oa.first, oa.second, -critical[static_cast<std::size_t>(a)]);
                const auto bk = std::make_tuple(ob.first, ob.second, -critical[static_cast<std::size_t>(b)]);
                if (ak != bk) return ak < bk;
            } else if (order_kind == 7) {
                const auto& oa = ops[static_cast<std::size_t>(a)];
                const auto& ob = ops[static_cast<std::size_t>(b)];
                const auto ak = std::make_tuple(oa.second, oa.first, -critical[static_cast<std::size_t>(a)]);
                const auto bk = std::make_tuple(ob.second, ob.first, -critical[static_cast<std::size_t>(b)]);
                if (ak != bk) return ak < bk;
            } else if (order_kind >= 8) {
                const auto ak = pseudo_random_key(a);
                const auto bk = pseudo_random_key(b);
                if (ak != bk) return ak < bk;
            }
            if (critical[static_cast<std::size_t>(a)] != critical[static_cast<std::size_t>(b)]) {
                return critical[static_cast<std::size_t>(a)] > critical[static_cast<std::size_t>(b)];
            }
            if (endpoint_sum(a) != endpoint_sum(b)) return endpoint_sum(a) > endpoint_sum(b);
            if (succ[static_cast<std::size_t>(a)].size() != succ[static_cast<std::size_t>(b)].size()) {
                return succ[static_cast<std::size_t>(a)].size() > succ[static_cast<std::size_t>(b)].size();
            }
            return a < b;
        });

        Layer layer;
        layer.mode = "logical";
        std::vector<int> selected;
        for (int idx : ready) {
            const auto& op = ops[static_cast<std::size_t>(idx)];
            RoutedCNOT routed;
            if (!try_route_into_layer(layout, "logical", layer, op.first, op.second, routed)) continue;
            layer.ops.push_back(std::move(routed));
            scheduled[static_cast<std::size_t>(idx)] = 1;
            selected.push_back(idx);
        }
        if (layer.ops.empty()) throw std::runtime_error("could not schedule ready logical CNOT");

        for (int idx : selected) {
            for (int j : succ[static_cast<std::size_t>(idx)]) {
                --pred_count[static_cast<std::size_t>(j)];
            }
            --remaining;
        }

        layer.logical_depth = 1;
        layers.push_back(std::move(layer));
    }
    return layers;
}

static bool better_logical_schedule(
    const std::vector<Layer>& candidate,
    const std::vector<Layer>& incumbent
) {
    if (candidate.size() != incumbent.size()) return candidate.size() < incumbent.size();
    if (candidate.empty()) return false;
    std::vector<int> candidate_sizes;
    std::vector<int> incumbent_sizes;
    candidate_sizes.reserve(candidate.size());
    incumbent_sizes.reserve(incumbent.size());
    for (const auto& layer : candidate) candidate_sizes.push_back(static_cast<int>(layer.ops.size()));
    for (const auto& layer : incumbent) incumbent_sizes.push_back(static_cast<int>(layer.ops.size()));
    std::sort(candidate_sizes.begin(), candidate_sizes.end());
    std::sort(incumbent_sizes.begin(), incumbent_sizes.end());
    return candidate_sizes > incumbent_sizes;
}

static std::vector<Layer> schedule_logical_commuting_ops(
    const std::vector<RowOp>& ops,
    const Layout& layout
) {
    if (ops.empty()) return {};
    std::vector<Layer> best = schedule_logical_commuting_ops_once(ops, layout, 0);
    const int order_count = (logical_depth_pressure_goal() ? 28 : 10)
        + std::max(0, g_logical_search_params.logical_reschedule_orders);
    const LogicalSequenceBounds bounds = logical_sequence_bounds(ops, layout.n);
    for (int order_kind = 1; order_kind < order_count; ++order_kind) {
        std::vector<Layer> candidate = schedule_logical_commuting_ops_once(ops, layout, order_kind);
        if (better_logical_schedule(candidate, best)) best = std::move(candidate);
        if (static_cast<int>(best.size()) == bounds.lower_bound) break;
    }
    return best;
}

static std::vector<Layer> schedule_ops(
    const std::vector<RowOp>& ops,
    const Layout& layout,
    const std::string& mode
) {
    if (mode == "logical") return schedule_logical_commuting_ops(ops, layout);

    std::vector<Layer> layers;
    std::vector<int> last_touch(layout.n, -1);
    for (const auto& [control, target] : ops) {
        const int earliest = std::max(last_touch[control], last_touch[target]) + 1;
        int best_li = -1;
        int best_delta = std::numeric_limits<int>::max();
        int best_new_depth = std::numeric_limits<int>::max();
        int best_path_len = std::numeric_limits<int>::max();
        RoutedCNOT best_routed;
        Layer best_replacement;
        bool best_uses_replacement = false;
        for (int li = earliest; li < static_cast<int>(layers.size()); ++li) {
            RoutedCNOT routed;
            std::vector<RoutedCNOT> trial_ops;
            Layer replacement;
            bool uses_replacement = false;
            if (try_route_into_layer(layout, mode, layers[li], control, target, routed)) {
                trial_ops = layers[li].ops;
                trial_ops.push_back(routed);
            } else if (repack_layer_with_op(layout, mode, layers[li], control, target, replacement)) {
                trial_ops = replacement.ops;
                uses_replacement = true;
            } else {
                continue;
            }
            const int new_depth = layer_depth(mode, layout, trial_ops);
            const int delta = std::max(0, new_depth - layers[li].logical_depth);
            const int path_len = uses_replacement
                ? static_cast<int>(trial_ops.back().path.size())
                : static_cast<int>(routed.path.size());
            const auto key = std::make_tuple(delta, li, new_depth, path_len);
            const auto best_key = std::make_tuple(best_delta, best_li, best_new_depth, best_path_len);
            if (best_li < 0 || key < best_key) {
                best_li = li;
                best_delta = delta;
                best_new_depth = new_depth;
                best_path_len = path_len;
                best_uses_replacement = uses_replacement;
                if (uses_replacement) {
                    best_replacement = std::move(replacement);
                } else {
                    best_routed = std::move(routed);
                }
            }
        }
        if (best_li >= 0) {
            if (best_uses_replacement) {
                layers[best_li] = std::move(best_replacement);
            } else {
                layers[best_li].ops.push_back(std::move(best_routed));
                layers[best_li].logical_depth = layer_depth(mode, layout, layers[best_li].ops);
            }
            last_touch[control] = best_li;
            last_touch[target] = best_li;
        } else {
            Layer layer;
            layer.mode = mode;
            RoutedCNOT routed;
            if (!try_route_into_layer(layout, mode, layer, control, target, routed)) {
                throw std::runtime_error("could not route scheduled CNOT");
            }
            layer.ops.push_back(std::move(routed));
            layer.logical_depth = layer_depth(mode, layout, layer.ops);
            layers.push_back(std::move(layer));
            const int li = static_cast<int>(layers.size()) - 1;
            last_touch[control] = li;
            last_touch[target] = li;
        }
    }
    for (auto& layer : layers) {
        layer.logical_depth = layer_depth(mode, layout, layer.ops);
    }
    return layers;
}

static std::vector<RowOp> topological_reorder_ops(
    const std::vector<RowOp>& ops,
    const Layout& layout,
    int order_kind
) {
    const int m = static_cast<int>(ops.size());
    std::vector<std::vector<int>> succ;
    std::vector<int> pred_count;
    build_logical_precedence(ops, succ, pred_count);
    const std::vector<int> critical = logical_critical_lengths(succ);

    std::vector<char> emitted(static_cast<std::size_t>(m), 0);
    std::vector<RowOp> out;
    out.reserve(ops.size());
    while (static_cast<int>(out.size()) < m) {
        std::vector<int> ready;
        for (int i = 0; i < m; ++i) {
            if (!emitted[static_cast<std::size_t>(i)] && pred_count[static_cast<std::size_t>(i)] == 0) {
                ready.push_back(i);
            }
        }
        if (ready.empty()) throw std::runtime_error("CNOT precedence graph has a cycle");

        auto distance = [&](int idx) {
            const auto& op = ops[static_cast<std::size_t>(idx)];
            return layout.manhattan_between_qubits(op.first, op.second);
        };
        auto qubit_row = [&](int q) { return layout.row(layout.data_pos[static_cast<std::size_t>(q)]); };
        auto qubit_col = [&](int q) { return layout.col(layout.data_pos[static_cast<std::size_t>(q)]); };
        auto pseudo_random_key = [&](int idx) {
            std::uint64_t x = 0x9e3779b97f4a7c15ULL
                ^ (static_cast<std::uint64_t>(order_kind) * 0xbf58476d1ce4e5b9ULL)
                ^ (static_cast<std::uint64_t>(idx) * 0x94d049bb133111ebULL);
            x ^= x >> 30;
            x *= 0xbf58476d1ce4e5b9ULL;
            x ^= x >> 27;
            x *= 0x94d049bb133111ebULL;
            x ^= x >> 31;
            return x;
        };
        std::sort(ready.begin(), ready.end(), [&](int a, int b) {
            if (order_kind == 1) {
                const int da = distance(a);
                const int db = distance(b);
                if (da != db) return da < db;
            } else if (order_kind == 2) {
                const int da = distance(a);
                const int db = distance(b);
                if (da != db) return da > db;
            } else if (order_kind == 3) {
                const auto& oa = ops[static_cast<std::size_t>(a)];
                const auto& ob = ops[static_cast<std::size_t>(b)];
                const auto ak = std::make_tuple(qubit_row(oa.first), qubit_col(oa.first), qubit_row(oa.second), qubit_col(oa.second));
                const auto bk = std::make_tuple(qubit_row(ob.first), qubit_col(ob.first), qubit_row(ob.second), qubit_col(ob.second));
                if (ak != bk) return ak < bk;
            } else if (order_kind == 4) {
                const auto& oa = ops[static_cast<std::size_t>(a)];
                const auto& ob = ops[static_cast<std::size_t>(b)];
                const auto ak = std::make_tuple(qubit_row(oa.second), qubit_col(oa.second), qubit_row(oa.first), qubit_col(oa.first));
                const auto bk = std::make_tuple(qubit_row(ob.second), qubit_col(ob.second), qubit_row(ob.first), qubit_col(ob.first));
                if (ak != bk) return ak < bk;
            } else if (order_kind == 5) {
                const auto& oa = ops[static_cast<std::size_t>(a)];
                const auto& ob = ops[static_cast<std::size_t>(b)];
                const int span_a = std::abs(qubit_col(oa.first) - qubit_col(oa.second));
                const int span_b = std::abs(qubit_col(ob.first) - qubit_col(ob.second));
                if (span_a != span_b) return span_a > span_b;
            } else if (order_kind >= 6) {
                const auto ak = pseudo_random_key(a);
                const auto bk = pseudo_random_key(b);
                if (ak != bk) return ak < bk;
            }
            if (critical[static_cast<std::size_t>(a)] != critical[static_cast<std::size_t>(b)]) {
                return critical[static_cast<std::size_t>(a)] > critical[static_cast<std::size_t>(b)];
            }
            if (succ[static_cast<std::size_t>(a)].size() != succ[static_cast<std::size_t>(b)].size()) {
                return succ[static_cast<std::size_t>(a)].size() > succ[static_cast<std::size_t>(b)].size();
            }
            return a < b;
        });

        const int chosen = ready.front();
        emitted[static_cast<std::size_t>(chosen)] = 1;
        out.push_back(ops[static_cast<std::size_t>(chosen)]);
        for (int j : succ[static_cast<std::size_t>(chosen)]) {
            --pred_count[static_cast<std::size_t>(j)];
        }
    }
    return out;
}

static std::vector<Layer> schedule_surface_commuting_ops(
    const std::vector<RowOp>& ops,
    const Layout& layout,
    const std::string& mode,
    int order_kind
) {
    if (mode == "logical") return schedule_logical_commuting_ops(ops, layout);

    const int m = static_cast<int>(ops.size());
    std::vector<std::vector<int>> succ;
    std::vector<int> pred_count;
    build_logical_precedence(ops, succ, pred_count);
    const std::vector<int> critical = logical_critical_lengths(succ);

    std::vector<char> scheduled(static_cast<std::size_t>(m), 0);
    int remaining = m;
    std::vector<Layer> layers;
    while (remaining > 0) {
        std::vector<int> ready;
        ready.reserve(static_cast<std::size_t>(remaining));
        for (int i = 0; i < m; ++i) {
            if (!scheduled[static_cast<std::size_t>(i)] && pred_count[static_cast<std::size_t>(i)] == 0) {
                ready.push_back(i);
            }
        }
        if (ready.empty()) throw std::runtime_error("surface CNOT precedence graph has a cycle");

        auto distance = [&](int idx) {
            const auto& op = ops[static_cast<std::size_t>(idx)];
            return layout.manhattan_between_qubits(op.first, op.second);
        };
        auto qubit_row = [&](int q) { return layout.row(layout.data_pos[static_cast<std::size_t>(q)]); };
        auto qubit_col = [&](int q) { return layout.col(layout.data_pos[static_cast<std::size_t>(q)]); };
        std::sort(ready.begin(), ready.end(), [&](int a, int b) {
            if (order_kind == 1) {
                const int da = distance(a);
                const int db = distance(b);
                if (da != db) return da < db;
            } else if (order_kind == 2) {
                const int da = distance(a);
                const int db = distance(b);
                if (da != db) return da > db;
            } else if (order_kind == 3) {
                const auto& oa = ops[static_cast<std::size_t>(a)];
                const auto& ob = ops[static_cast<std::size_t>(b)];
                const auto ak = std::make_tuple(
                    std::min(qubit_row(oa.first), qubit_row(oa.second)),
                    std::min(qubit_col(oa.first), qubit_col(oa.second)),
                    distance(a)
                );
                const auto bk = std::make_tuple(
                    std::min(qubit_row(ob.first), qubit_row(ob.second)),
                    std::min(qubit_col(ob.first), qubit_col(ob.second)),
                    distance(b)
                );
                if (ak != bk) return ak < bk;
            } else if (order_kind == 4) {
                const auto& oa = ops[static_cast<std::size_t>(a)];
                const auto& ob = ops[static_cast<std::size_t>(b)];
                const int span_a = std::abs(qubit_row(oa.first) - qubit_row(oa.second));
                const int span_b = std::abs(qubit_row(ob.first) - qubit_row(ob.second));
                if (span_a != span_b) return span_a > span_b;
            }
            if (critical[static_cast<std::size_t>(a)] != critical[static_cast<std::size_t>(b)]) {
                return critical[static_cast<std::size_t>(a)] > critical[static_cast<std::size_t>(b)];
            }
            if (succ[static_cast<std::size_t>(a)].size() != succ[static_cast<std::size_t>(b)].size()) {
                return succ[static_cast<std::size_t>(a)].size() > succ[static_cast<std::size_t>(b)].size();
            }
            return a < b;
        });

        Layer layer;
        layer.mode = mode;
        std::vector<int> selected;
        for (int idx : ready) {
            const auto& op = ops[static_cast<std::size_t>(idx)];
            RoutedCNOT routed;
            if (try_route_into_layer(layout, mode, layer, op.first, op.second, routed)) {
                layer.ops.push_back(std::move(routed));
            } else {
                Layer replacement;
                if (!repack_layer_with_op(layout, mode, layer, op.first, op.second, replacement)) continue;
                layer = std::move(replacement);
            }
            scheduled[static_cast<std::size_t>(idx)] = 1;
            selected.push_back(idx);
        }
        if (layer.ops.empty()) throw std::runtime_error("could not schedule ready surface CNOT");

        for (int idx : selected) {
            for (int j : succ[static_cast<std::size_t>(idx)]) {
                --pred_count[static_cast<std::size_t>(j)];
            }
            --remaining;
        }

        layer.logical_depth = layer_depth(mode, layout, layer.ops);
        layers.push_back(std::move(layer));
    }
    return layers;
}

struct SequenceCost {
    int surface_depth = 0;
    int layers = 0;
    int cnots = 0;
};

static SequenceCost cost_of_layers(const std::vector<Layer>& layers) {
    SequenceCost cost;
    cost.layers = static_cast<int>(layers.size());
    for (const auto& layer : layers) {
        cost.surface_depth += layer.logical_depth;
        cost.cnots += static_cast<int>(layer.ops.size());
    }
    return cost;
}

static bool better_sequence_cost(const SequenceCost& a, const SequenceCost& b) {
    return std::make_tuple(a.surface_depth, a.layers, a.cnots)
         < std::make_tuple(b.surface_depth, b.layers, b.cnots);
}

static SequenceCost cost_of_sequence(
    const std::vector<RowOp>& ops,
    const Layout& layout,
    const std::string& mode
) {
    return cost_of_layers(schedule_ops(ops, layout, mode));
}

static std::vector<RowOp> flatten_layers(const std::vector<Layer>& layers) {
    std::vector<RowOp> out;
    for (const auto& layer : layers) {
        for (const auto& op : layer.ops) out.emplace_back(op.control, op.target);
    }
    return out;
}

static std::vector<RowOp> flatten_rowop_layers(
    const std::vector<std::vector<RowOp>>& layers,
    bool reverse_layers,
    bool reverse_ops
) {
    std::vector<RowOp> out;
    if (reverse_layers) {
        for (auto layer_it = layers.rbegin(); layer_it != layers.rend(); ++layer_it) {
            if (reverse_ops) {
                out.insert(out.end(), layer_it->rbegin(), layer_it->rend());
            } else {
                out.insert(out.end(), layer_it->begin(), layer_it->end());
            }
        }
    } else {
        for (const auto& layer : layers) {
            if (reverse_ops) {
                out.insert(out.end(), layer.rbegin(), layer.rend());
            } else {
                out.insert(out.end(), layer.begin(), layer.end());
            }
        }
    }
    return out;
}

static std::vector<RowOp> transform_ops_by_permutation(
    const std::vector<RowOp>& ops,
    const std::vector<int>& perm
) {
    std::vector<RowOp> out;
    out.reserve(ops.size());
    for (const auto& op : ops) out.emplace_back(perm[op.first], perm[op.second]);
    return out;
}

static std::vector<RowOp> swap_controls_targets(const std::vector<RowOp>& ops) {
    std::vector<RowOp> out;
    out.reserve(ops.size());
    for (const auto& op : ops) out.emplace_back(op.second, op.first);
    return out;
}

static Matrix apply_ops_to_identity(int n, const std::vector<RowOp>& ops) {
    Matrix current = identity_matrix(n);
    for (const auto& op : ops) row_add(current, op.first, op.second);
    return current;
}

static std::vector<int> active_qubits_for_segment(
    const std::vector<RowOp>& ops,
    int begin,
    int len
) {
    std::vector<int> active;
    active.reserve(static_cast<std::size_t>(2 * len));
    for (int i = 0; i < len; ++i) {
        active.push_back(ops[begin + i].first);
        active.push_back(ops[begin + i].second);
    }
    std::sort(active.begin(), active.end());
    active.erase(std::unique(active.begin(), active.end()), active.end());
    return active;
}

static Matrix segment_matrix(
    const std::vector<RowOp>& ops,
    int begin,
    int len,
    const std::vector<int>& active,
    int n
) {
    std::vector<int> local_of(static_cast<std::size_t>(n), -1);
    for (int i = 0; i < static_cast<int>(active.size()); ++i) local_of[active[i]] = i;
    Matrix matrix = identity_matrix(static_cast<int>(active.size()));
    for (int i = 0; i < len; ++i) {
        const auto& op = ops[begin + i];
        row_add(matrix, local_of[op.first], local_of[op.second]);
    }
    return matrix;
}

static std::vector<RowOp> map_local_sequence(
    const std::vector<RowOp>& local_ops,
    const std::vector<int>& active
) {
    std::vector<RowOp> mapped;
    mapped.reserve(local_ops.size());
    for (const auto& op : local_ops) mapped.emplace_back(active[op.first], active[op.second]);
    return mapped;
}

static bool segment_equals(
    const std::vector<RowOp>& ops,
    int begin,
    int len,
    const std::vector<RowOp>& replacement
) {
    if (static_cast<int>(replacement.size()) != len) return false;
    for (int i = 0; i < len; ++i) {
        if (ops[begin + i] != replacement[static_cast<std::size_t>(i)]) return false;
    }
    return true;
}

static std::vector<RowOp> replace_segment(
    const std::vector<RowOp>& ops,
    int begin,
    int len,
    const std::vector<RowOp>& replacement
) {
    std::vector<RowOp> out;
    out.reserve(ops.size() - static_cast<std::size_t>(len) + replacement.size());
    out.insert(out.end(), ops.begin(), ops.begin() + begin);
    out.insert(out.end(), replacement.begin(), replacement.end());
    out.insert(out.end(), ops.begin() + begin + len, ops.end());
    return out;
}

static std::vector<RowOp> recombine_cnot_sequence(
    std::vector<RowOp> ops,
    const Layout& layout,
    const std::string& mode,
    int max_segment_len = 5,
    int max_active_qubits = 4,
    int max_passes = 6
) {
    if (ops.size() < 2) return ops;
    ops = optimize_cnot_sequence(ops);
    SequenceCost best_cost = cost_of_sequence(ops, layout, mode);

    for (int pass = 0; pass < max_passes; ++pass) {
        bool improved = false;
        const int largest = std::min(max_segment_len, static_cast<int>(ops.size()));
        for (int len = largest; len >= 2 && !improved; --len) {
            for (int begin = 0; begin + len <= static_cast<int>(ops.size()) && !improved; ++begin) {
                const std::vector<int> active = active_qubits_for_segment(ops, begin, len);
                const int k = static_cast<int>(active.size());
                if (k < 2 || k > max_active_qubits) continue;

                const Matrix matrix = segment_matrix(ops, begin, len, active, layout.n);
                const auto& candidates = local_reduction_db(k, max_segment_len).equivalents(matrix);
                for (const auto& local_candidate : candidates) {
                    const int replacement_len = static_cast<int>(local_candidate.size());
                    if (replacement_len > len) continue;
                    const std::vector<RowOp> replacement = map_local_sequence(local_candidate, active);
                    if (segment_equals(ops, begin, len, replacement)) continue;

                    std::vector<RowOp> trial = replace_segment(ops, begin, len, replacement);
                    trial = optimize_cnot_sequence(trial);
                    if (trial == ops) continue;

                    SequenceCost trial_cost;
                    try {
                        trial_cost = cost_of_sequence(trial, layout, mode);
                    } catch (const std::exception&) {
                        continue;
                    }
                    if (better_sequence_cost(trial_cost, best_cost)) {
                        ops = std::move(trial);
                        best_cost = trial_cost;
                        improved = true;
                        break;
                    }
                }
            }
        }
        if (!improved) break;
    }
    return ops;
}

static void apply_layer(Matrix& rows, const std::vector<RoutedCNOT>& ops) {
    for (const auto& op : ops) row_add(rows, op.control, op.target);
}

static void validate_routed_path(
    const Result& result,
    const Layer& layer,
    const RoutedCNOT& op
) {
    if (op.control < 0 || op.control >= result.n || op.target < 0 || op.target >= result.n
        || op.control == op.target) {
        throw std::runtime_error("CNOT endpoint out of range");
    }

    if (layer.mode == "logical") {
        if (!op.path.empty()) throw std::runtime_error("logical CNOT unexpectedly has a routed path");
        return;
    }

    if (op.path.size() < 2) throw std::runtime_error("surface CNOT has an empty routed path");
    const int start = result.layout.data_pos[static_cast<std::size_t>(op.control)];
    const int goal = result.layout.data_pos[static_cast<std::size_t>(op.target)];
    if (op.path.front() != start || op.path.back() != goal) {
        throw std::runtime_error("routed path endpoints do not match CNOT endpoints");
    }

    std::unordered_set<int> seen_vertices;
    for (std::size_t i = 0; i < op.path.size(); ++i) {
        const int v = op.path[i];
        if (v < 0 || v >= result.layout.grid_rows * result.layout.grid_cols) {
            throw std::runtime_error("routed path vertex out of layout range");
        }
        if (!seen_vertices.insert(v).second) throw std::runtime_error("routed path is not simple");
        if (i > 0) {
            const int prev = op.path[i - 1];
            const int dr = std::abs(result.layout.row(prev) - result.layout.row(v));
            const int dc = std::abs(result.layout.col(prev) - result.layout.col(v));
            if (dr + dc != 1) throw std::runtime_error("routed path has non-adjacent vertices");
        }
        if (i > 0 && i + 1 < op.path.size() && result.layout.occupied_data[static_cast<std::size_t>(v)]) {
            throw std::runtime_error("routed path passes through an occupied data qubit");
        }
    }

    if (result.layout.col(op.path[1]) != result.layout.col(start)) {
        throw std::runtime_error("routed path does not leave the control vertically");
    }
    if (result.layout.row(op.path[op.path.size() - 2]) != result.layout.row(goal)) {
        throw std::runtime_error("routed path does not enter the target horizontally");
    }
}

static int checked_layer_depth(const Result& result, const Layer& layer) {
    std::vector<char> used_qubits(static_cast<std::size_t>(result.n), 0);
    std::unordered_set<int> used_vertices;
    std::unordered_set<std::uint64_t> used_edges;
    bool vertex_disjoint = true;

    for (const auto& op : layer.ops) {
        validate_routed_path(result, layer, op);
        if (used_qubits[static_cast<std::size_t>(op.control)]
            || used_qubits[static_cast<std::size_t>(op.target)]) {
            throw std::runtime_error("layer reuses a logical endpoint");
        }
        used_qubits[static_cast<std::size_t>(op.control)] = 1;
        used_qubits[static_cast<std::size_t>(op.target)] = 1;

        if (layer.mode == "logical") continue;
        for (int v : op.path) {
            if (!used_vertices.insert(v).second) vertex_disjoint = false;
        }
        for (auto e : result.layout.path_edges(op.path)) {
            if (!used_edges.insert(e).second) throw std::runtime_error("layer reuses a routing edge");
        }
    }

    if (layer.mode == "logical") return 1;
    if (layer.mode == "vdp") {
        if (!vertex_disjoint) throw std::runtime_error("VDP layer has intersecting routed paths");
        return 2;
    }
    if (layer.mode == "edp") return vertex_disjoint ? 2 : 4;
    throw std::runtime_error("unknown layer mode in result");
}

static bool verify_geometry_and_stats(const Result& result) {
    try {
        if (result.n <= 0 || result.layout.n != result.n) return false;
        check_square_rows(result.target_rows);
        check_square_rows(result.final_rows);
        validate_permutation_vector(result.output_permutation, result.n);

        Matrix rows = identity_matrix(result.n);
        int cnots = 0;
        int surface_depth = 0;
        for (const auto& layer : result.layers) {
            const int actual_depth = checked_layer_depth(result, layer);
            if (layer.logical_depth != actual_depth) return false;
            apply_layer(rows, layer.ops);
            cnots += static_cast<int>(layer.ops.size());
            surface_depth += actual_depth;
        }
        if (cnots != result.cnot_count()) return false;
        if (surface_depth != result.surface_depth()) return false;
        return rows == result.final_rows && result.verify();
    } catch (const std::exception&) {
        return false;
    }
}
static std::vector<Result> completed;
static Matrix fresh_target;
static std::vector<int> fresh_output;
static Layout fresh_layout;
static bool export_verified_anytime(const Result& candidate, const char*) {
    if (candidate.mode != "vdp" || !verify_geometry_and_stats(candidate)) return false;
    if (candidate.target_rows != fresh_target || candidate.output_permutation != fresh_output
        || candidate.layout.data_pos != fresh_layout.data_pos
        || candidate.layout.occupied_data != fresh_layout.occupied_data) return false;
    if (!completed.empty() && std::make_pair(candidate.surface_depth(), candidate.cnot_count()) >=
        std::make_pair(completed.back().surface_depth(), completed.back().cnot_count())) return false;
    completed.push_back(candidate);
    if (completed.size() > 8) completed.erase(completed.begin());
    return true;
}

static std::vector<ShiFengSurfaceConfig> shi_feng_surface_configs() {
    return {
        {ShiFengCostKind::Square, 0.03, 0.45, true, 0.0, false, false, false},
        {ShiFengCostKind::Product, 0.04, 0.40, true, 0.0, false, false, false},
        {ShiFengCostKind::Square, 0.05, 0.25, false, 0.0, false, false, false},
        {ShiFengCostKind::Product, 0.05, 0.20, false, 0.0, false, false, false},
        {ShiFengCostKind::Square, 0.03, 0.45, true, 8.0, true, false, false},
        {ShiFengCostKind::Product, 0.04, 0.40, true, 0.75, true, false, false},
        {ShiFengCostKind::Square, 0.05, 0.25, false, 8.0, true, false, false},
        {ShiFengCostKind::Product, 0.05, 0.20, false, 0.75, true, false, false}
    };
}

static std::vector<ShiFengSurfaceConfig> shi_feng_logical_configs() {
    std::vector<ShiFengSurfaceConfig> configs = {
        {ShiFengCostKind::Square, 0.03, 0.45, true, 0.0, false, false, false},
        {ShiFengCostKind::Product, 0.04, 0.40, true, 0.0, false, false, false},
        {ShiFengCostKind::Square, 0.05, 0.25, false, 0.0, false, false, false},
        {ShiFengCostKind::Product, 0.05, 0.20, false, 0.0, false, false, false}
    };

    if (g_logical_search_params.use_layer_matrix_rank) {
        if (g_logical_search_params.use_layer_matrix_only) {
            return {
                {ShiFengCostKind::Square, 0.03, 0.45, true, 0.0, true, false, false, true},
                {ShiFengCostKind::Product, 0.04, 0.40, true, 0.0, true, false, false, true},
                {ShiFengCostKind::Square, 0.03, 0.45, true, 0.0, true, true, false, true},
                {ShiFengCostKind::Product, 0.04, 0.40, true, 0.0, true, true, false, true}
            };
        }
        configs.push_back({ShiFengCostKind::Square, 0.03, 0.45, true, 0.0, true, false, false, true});
        configs.push_back({ShiFengCostKind::Product, 0.04, 0.40, true, 0.0, true, false, false, true});
        configs.push_back({ShiFengCostKind::Square, 0.03, 0.45, true, 0.0, true, true, false, true});
        configs.push_back({ShiFengCostKind::Product, 0.04, 0.40, true, 0.0, true, true, false, true});
    }
    return configs;
}

static Result synthesize_shi_feng_surface_once(
    const Matrix& target_rows,
    const std::string& mode,
    const Layout& layout,
    const std::vector<int>& output_permutation,
    const ShiFengSurfaceConfig& config,
    int max_layers = 10000
) {
    const int n = check_square_rows(target_rows);
    if (layout.n != n) throw std::runtime_error("layout size mismatch");
    if (static_cast<int>(output_permutation.size()) != n) throw std::runtime_error("permutation size mismatch");

    const Matrix desired = permute_rows(target_rows, output_permutation);
    Matrix residual = invert_matrix(desired);
    const Matrix identity = identity_matrix(n);
    std::vector<RowOp> row_ops;
    std::vector<RowOp> column_ops;
    Layer current_row_layer{mode, {}, 0};
    Layer current_column_layer{mode, {}, 0};
    int closed_search_layers = 0;
    int closed_surface_depth = 0;
    const int step_budget = std::max(128, 16 * n * n);

    auto flush_search_layers = [&]() {
        bool flushed = false;
        if (!current_row_layer.ops.empty()) {
            ++closed_search_layers;
            closed_surface_depth += current_row_layer.logical_depth;
            current_row_layer.ops.clear();
            current_row_layer.logical_depth = 0;
            flushed = true;
        }
        if (!current_column_layer.ops.empty()) {
            ++closed_search_layers;
            closed_surface_depth += current_column_layer.logical_depth;
            current_column_layer.ops.clear();
            current_column_layer.logical_depth = 0;
            flushed = true;
        }
        return flushed;
    };

    for (int iter = 0; iter < step_budget; ++iter) {
        if (residual == identity) break;
        if (closed_search_layers >= max_layers || closed_surface_depth >= 2 * max_layers) break;

        const double current_cost = shi_current_cost(residual, config.cost_kind);
        std::vector<ShiFengChoice> choices;
        choices.reserve(static_cast<std::size_t>(2 * n * (n - 1)));

        auto consider = [&](OpSide side, int control, int target) {
            Layer& node = (side == OpSide::Row) ? current_row_layer : current_column_layer;
            RoutedCNOT routed;
            int depth_delta = 0;
            if (!try_route_into_search_node(layout, mode, node, control, target, routed, depth_delta)) {
                return;
            }

            AlgebraicChoice choice{side, control, target, 0.0, 0.0};
            const double after_cost = shi_after_cost(residual, choice, config.cost_kind);
            if (!(after_cost + 1.0e-9 < current_cost)) return;

            const double route = config.route_penalty * layout.manhattan_between_qubits(control, target);
            const double depth = config.surface_depth_penalty * static_cast<double>(depth_delta);
            const double delta = current_cost - after_cost;
            const double side_bias = (side == OpSide::Column) ? 0.0005 : 0.0;
            const double rank_cost_after = config.rank_global_cost
                ? after_cost + route + depth
                : after_cost;
            choices.push_back({
                side,
                control,
                target,
                after_cost,
                rank_cost_after,
                delta + side_bias - route - depth,
                depth_delta,
                0,
                0,
                0,
                0,
                0,
                0,
                std::move(routed)
            });
        };

        for (int control = 0; control < n; ++control) {
            for (int target = 0; target < n; ++target) {
                if (control == target) continue;
                consider(OpSide::Row, control, target);
                if (config.allow_columns) consider(OpSide::Column, control, target);
            }
        }

        if (choices.empty()) {
            if (flush_search_layers()) continue;
            break;
        }

        std::sort(choices.begin(), choices.end(), [](const ShiFengChoice& a, const ShiFengChoice& b) {
            if (a.rank_cost_after != b.rank_cost_after) return a.rank_cost_after < b.rank_cost_after;
            if (a.cost_after != b.cost_after) return a.cost_after < b.cost_after;
            if (a.layer_depth_delta != b.layer_depth_delta) return a.layer_depth_delta < b.layer_depth_delta;
            if (a.priority != b.priority) return a.priority > b.priority;
            if (a.control != b.control) return a.control < b.control;
            if (a.target != b.target) return a.target < b.target;
            return static_cast<int>(a.side) < static_cast<int>(b.side);
        });

        ShiFengChoice choice = std::move(choices.front());
        AlgebraicChoice algebraic{choice.side, choice.control, choice.target, 0.0, choice.priority};
        apply_algebraic_choice(residual, algebraic);
        if (choice.side == OpSide::Row) {
            row_ops.emplace_back(choice.control, choice.target);
            append_to_search_node(current_row_layer, layout, mode, std::move(choice.routed));
        } else {
            column_ops.emplace_back(choice.control, choice.target);
            append_to_search_node(current_column_layer, layout, mode, std::move(choice.routed));
        }
    }

    flush_search_layers();

    if (residual != identity) {
        const auto fallback = gaussian_reduction_ops(residual);
        for (const RowOp& op : fallback) {
            row_add(residual, op.first, op.second);
            row_ops.push_back(op);
        }
    }
    if (residual != identity) throw std::runtime_error("Shi-Feng surface synthesis failed");

    std::vector<RowOp> circuit_ops = row_ops;
    for (auto it = column_ops.rbegin(); it != column_ops.rend(); ++it) {
        circuit_ops.push_back(*it);
    }
    circuit_ops = optimize_cnot_sequence(circuit_ops);
    std::vector<Layer> layers = schedule_ops(circuit_ops, layout, mode);

    Matrix current = identity_matrix(n);
    for (const auto& layer : layers) apply_layer(current, layer.ops);

    Result result;
    result.n = n;
    result.mode = mode;
    result.layout = layout;
    result.target_rows = target_rows;
    result.output_permutation = output_permutation;
    result.layers = std::move(layers);
    result.final_rows = current;
    if (!verify_geometry_and_stats(result)) throw std::runtime_error("Shi-Feng surface verification failed");
    return result;
}

static Result synthesize_shi_feng_logical(
    const Matrix& target_rows,
    const Layout& layout,
    int max_layers,
    const std::vector<int>* required_output_permutation = nullptr
);

static Result synthesize_cached_layer_matrix_logical(
    const Matrix& target_rows,
    const Layout& layout,
    int max_layers,
    const std::vector<int>* required_output_permutation = nullptr
) {
    static std::mutex cache_mutex;
    static std::unordered_map<std::string, Result> cache;
    std::ostringstream key_stream;
    key_stream << matrix_key(target_rows)
               << "|n=" << layout.n
               << "|targets=" << g_logical_search_params.max_search_targets
               << "|seeds=" << g_logical_search_params.seeds_per_target
               << "|window=" << g_logical_search_params.choice_window
               << "|depth=" << g_logical_search_params.depth_goal
               << "|cnot=" << g_logical_search_params.cnot_goal
               << "|layer=" << g_logical_search_params.use_layer_matrix_rank
               << "|layer_only=" << g_logical_search_params.use_layer_matrix_only
               << "|lm_extra_orders=" << g_logical_search_params.layer_matrix_extra_orders
               << "|lm_order_variants=" << g_logical_search_params.layer_matrix_order_variants
               << "|lm_candidate_multiplier=" << g_logical_search_params.layer_matrix_candidate_multiplier
               << "|lm_beam_multiplier=" << g_logical_search_params.layer_matrix_beam_multiplier
               << "|lm_tail_cap=" << g_logical_search_params.layer_matrix_tail_cap
               << "|lm_layer_cap_variants=" << g_logical_search_params.layer_matrix_layer_cap_variants
               << "|lm_nonimproving_slack=" << g_logical_search_params.layer_matrix_nonimproving_slack
               << "|lm_goal_only=" << g_logical_search_params.layer_matrix_goal_only
               << "|lm_endpoint_first=" << g_logical_search_params.layer_matrix_endpoint_first
               << "|lm_vdp_hard=" << g_logical_search_params.layer_matrix_vdp_hard
               << "|logical_reschedule_orders=" << g_logical_search_params.logical_reschedule_orders
               << "|max=" << max_layers
               << "|required_output=";
    if (required_output_permutation == nullptr) {
        key_stream << "any";
    } else {
        for (int value : *required_output_permutation) key_stream << value << ',';
    }
    const std::string key = key_stream.str();

    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        const auto it = cache.find(key);
        if (it != cache.end()) return it->second;
    }

    Result result = synthesize_shi_feng_logical(
        target_rows,
        layout,
        max_layers,
        required_output_permutation
    );
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        cache[key] = result;
    }
    return result;
}

static Result synthesize_shi_feng_surface(
    const Matrix& target_rows,
    const std::string& mode,
    const Layout& layout,
    const std::vector<int>& output_permutation,
    int max_layers = 10000,
    const std::vector<int>* required_output_permutation = nullptr
) {
    if (mode == "logical") {
        (void)output_permutation;
        return synthesize_shi_feng_logical(
            target_rows,
            layout,
            max_layers,
            required_output_permutation
        );
    }

    bool have_best = false;
    Result best;
    for (const auto& config : shi_feng_surface_configs()) {
        try {
            Result result = synthesize_shi_feng_surface_once(
                target_rows,
                mode,
                layout,
                output_permutation,
                config,
                max_layers
            );
            export_verified_anytime(result, "shi_feng_complete_surface_configuration");
            const auto key = std::make_tuple(result.surface_depth(), result.layer_count(), result.cnot_count());
            const auto best_key = std::make_tuple(best.surface_depth(), best.layer_count(), best.cnot_count());
            if (!have_best || key < best_key) {
                best = std::move(result);
                have_best = true;
            }
        } catch (const std::exception&) {
        }
    }
    if (!have_best) throw std::runtime_error("no Shi-Feng surface result");
    return best;
}

static Result make_result_from_logical_ops(
    const Matrix& target_rows,
    const Layout& layout,
    std::vector<RowOp> ops
) {
    const int n = check_square_rows(target_rows);
    ops = optimize_cnot_sequence(ops);
    const Matrix final_rows = apply_ops_to_identity(n, ops);
    std::vector<int> output_permutation;
    if (!infer_output_permutation(target_rows, final_rows, output_permutation)) {
        throw std::runtime_error("logical circuit is not an output permutation of target");
    }
    std::vector<Layer> layers = schedule_ops(ops, layout, "logical");

    Result result;
    result.n = n;
    result.mode = "logical";
    result.layout = layout;
    result.target_rows = target_rows;
    result.output_permutation = std::move(output_permutation);
    result.layers = std::move(layers);
    result.final_rows = final_rows;
    if (!verify_geometry_and_stats(result)) throw std::runtime_error("logical verification failed");
    return result;
}

static Result make_result_from_logical_layers(
    const Matrix& target_rows,
    const Layout& layout,
    const std::vector<std::vector<RowOp>>& rowop_layers
) {
    const int n = check_square_rows(target_rows);
    Matrix final_rows = identity_matrix(n);
    std::vector<Layer> layers;
    layers.reserve(rowop_layers.size());
    for (const auto& rowop_layer : rowop_layers) {
        if (rowop_layer.empty()) continue;
        Layer layer;
        layer.mode = "logical";
        for (const auto& op : rowop_layer) {
            RoutedCNOT routed;
            if (!try_route_into_layer(layout, "logical", layer, op.first, op.second, routed)) {
                throw std::runtime_error("logical search layer is not parallel");
            }
            layer.ops.push_back(std::move(routed));
            row_add(final_rows, op.first, op.second);
        }
        layer.logical_depth = 1;
        layers.push_back(std::move(layer));
    }

    std::vector<int> output_permutation;
    if (!infer_output_permutation(target_rows, final_rows, output_permutation)) {
        throw std::runtime_error("logical layers are not an output permutation of target");
    }

    if (g_logical_search_params.logical_reschedule_orders > 0) {
        std::vector<RowOp> flat;
        for (const auto& layer_ops : rowop_layers) {
            flat.insert(flat.end(), layer_ops.begin(), layer_ops.end());
        }
        try {
            std::vector<Layer> rescheduled = schedule_logical_commuting_ops(flat, layout);
            if (better_logical_schedule(rescheduled, layers)) {
                Matrix replay = identity_matrix(n);
                for (const auto& layer : rescheduled) apply_layer(replay, layer.ops);
                if (replay == final_rows) layers = std::move(rescheduled);
            }
        } catch (const std::exception&) {
        }
    }

    Result result;
    result.n = n;
    result.mode = "logical";
    result.layout = layout;
    result.target_rows = target_rows;
    result.output_permutation = std::move(output_permutation);
    result.layers = std::move(layers);
    result.final_rows = final_rows;
    if (!verify_geometry_and_stats(result)) throw std::runtime_error("logical layer verification failed");
    return result;
}

static Result schedule_logical_result_on_surface(
    const Result& logical_result,
    const std::string& mode,
    const Layout& layout,
    const std::vector<int>& output_permutation
) {
    if (mode == "logical") return logical_result;
    if (logical_result.output_permutation != output_permutation) {
        throw std::runtime_error("logical candidate has incompatible output permutation");
    }

    std::vector<std::vector<RowOp>> grouped;
    grouped.reserve(logical_result.layers.size());
    for (const auto& layer : logical_result.layers) {
        std::vector<RowOp> ops;
        ops.reserve(layer.ops.size());
        for (const auto& op : layer.ops) ops.emplace_back(op.control, op.target);
        grouped.push_back(std::move(ops));
    }

    std::vector<std::vector<std::vector<RowOp>>> variants;
    variants.push_back(grouped);
    auto sorted_variant = [&](auto less) {
        std::vector<std::vector<RowOp>> variant = grouped;
        for (auto& layer : variant) std::sort(layer.begin(), layer.end(), less);
        variants.push_back(std::move(variant));
    };
    sorted_variant([&](const RowOp& a, const RowOp& b) {
        const int da = layout.manhattan_between_qubits(a.first, a.second);
        const int db = layout.manhattan_between_qubits(b.first, b.second);
        if (da != db) return da < db;
        return a < b;
    });
    sorted_variant([&](const RowOp& a, const RowOp& b) {
        const int da = layout.manhattan_between_qubits(a.first, a.second);
        const int db = layout.manhattan_between_qubits(b.first, b.second);
        if (da != db) return da > db;
        return a < b;
    });
    sorted_variant([](const RowOp& a, const RowOp& b) {
        return a < b;
    });
    sorted_variant([](const RowOp& a, const RowOp& b) {
        return std::make_tuple(a.second, a.first) < std::make_tuple(b.second, b.first);
    });
    for (int seed = 0; seed < 4; ++seed) {
        std::vector<std::vector<RowOp>> variant = grouped;
        std::uint64_t state = 0x243f6a8885a308d3ULL
            ^ (static_cast<std::uint64_t>(seed) * 0x9e3779b97f4a7c15ULL);
        for (auto& layer : variant) {
            state ^= static_cast<std::uint64_t>(layer.size()) + 0x9e3779b97f4a7c15ULL + (state << 6) + (state >> 2);
            std::mt19937_64 rng(state);
            std::shuffle(layer.begin(), layer.end(), rng);
        }
        variants.push_back(std::move(variant));
    }

    bool have_best = false;
    Result best;
    for (const auto& variant : variants) {
        std::vector<RowOp> flat;
        for (const auto& layer : variant) flat.insert(flat.end(), layer.begin(), layer.end());
        std::vector<std::vector<Layer>> scheduled_variants;
        scheduled_variants.push_back(schedule_ops(flat, layout, mode));
        for (int order_kind = 0; order_kind < 14; ++order_kind) {
            scheduled_variants.push_back(schedule_ops(
                topological_reorder_ops(flat, layout, order_kind),
                layout,
                mode
            ));
        }
        for (int order_kind = 0; order_kind < 5; ++order_kind) {
            scheduled_variants.push_back(schedule_surface_commuting_ops(flat, layout, mode, order_kind));
        }

        for (auto& layers : scheduled_variants) {
            Matrix current = identity_matrix(logical_result.n);
            for (const auto& layer : layers) apply_layer(current, layer.ops);

            Result result;
            result.n = logical_result.n;
            result.mode = mode;
            result.layout = layout;
            result.target_rows = logical_result.target_rows;
            result.output_permutation = output_permutation;
            result.layers = std::move(layers);
            result.final_rows = current;
            if (!verify_geometry_and_stats(result)) continue;

            const auto key = std::make_tuple(result.surface_depth(), result.layer_count(), result.cnot_count());
            const auto best_key = std::make_tuple(best.surface_depth(), best.layer_count(), best.cnot_count());
            if (!have_best || key < best_key) {
                best = std::move(result);
                have_best = true;
            }
        }
    }
    if (!have_best) throw std::runtime_error("surface schedule of logical candidate failed verification");
    return best;
}

static Result retarget_logical_result(Result result, const Matrix& target_rows) {
    std::vector<int> output_permutation;
    if (!infer_output_permutation(target_rows, result.final_rows, output_permutation)) {
        throw std::runtime_error("logical result is not a permutation of the requested target");
    }
    result.target_rows = target_rows;
    result.output_permutation = std::move(output_permutation);
    if (!verify_geometry_and_stats(result)) throw std::runtime_error("retargeted logical result failed verification");
    return result;
}

static std::tuple<int, int, int, int, int, int> logical_result_priority_key(const Result& item) {
    if (g_logical_search_params.use_layer_matrix_rank
        && logical_depth_pressure_goal()) {
        const LogicalSequenceBounds bounds = logical_sequence_bounds(flatten_layers(item.layers), item.n);
        const int misses_goal = item.layer_count() <= g_logical_search_params.depth_goal ? 0 : 1;
        if (misses_goal == 0) {
            return std::make_tuple(
                0,
                item.layer_count(),
                item.cnot_count(),
                bounds.critical_bound,
                bounds.endpoint_bound,
                0
            );
        }
        const int depth_or_bound = bounds.lower_bound;
        return std::make_tuple(
            misses_goal,
            depth_or_bound,
            bounds.critical_bound,
            bounds.endpoint_bound,
            item.layer_count(),
            item.cnot_count()
        );
    }
    return std::make_tuple(
        0,
        item.layer_count(),
        item.layer_count(),
        item.layer_count(),
        item.layer_count(),
        item.cnot_count()
    );
}

static Result reconstruct_logical_search_result(
    const Matrix& target_rows,
    const Layout& layout,
    const Matrix& residual,
    const std::vector<std::vector<RowOp>>& row_layers,
    const std::vector<std::vector<RowOp>>& column_layers
) {
    if (!is_permutation_matrix(residual)) {
        throw std::runtime_error("logical search did not reach a permutation");
    }

    const std::vector<int> p = permutation_columns(residual);
    const std::vector<int> inv_p = inverse_permutation(p);
    std::vector<std::vector<RowOp>> candidates;
    std::vector<std::vector<std::vector<RowOp>>> layer_candidates;

    for (bool col_reverse_layers : {false, true}) {
        for (bool col_swap : {false, true}) {
            for (int transform = 0; transform < 3; ++transform) {
                for (bool rows_first : {false, true}) {
                    std::vector<std::vector<RowOp>> layers;
                    auto append_columns = [&]() {
                        if (col_reverse_layers) {
                            for (auto it = column_layers.rbegin(); it != column_layers.rend(); ++it) {
                                layers.push_back(col_swap ? swap_controls_targets(*it) : *it);
                            }
                        } else {
                            for (const auto& layer : column_layers) {
                                layers.push_back(col_swap ? swap_controls_targets(layer) : layer);
                            }
                        }
                    };
                    auto append_rows = [&]() {
                        for (auto it = row_layers.rbegin(); it != row_layers.rend(); ++it) {
                            std::vector<RowOp> layer = *it;
                            if (transform == 1) layer = transform_ops_by_permutation(layer, p);
                            if (transform == 2) layer = transform_ops_by_permutation(layer, inv_p);
                            layers.push_back(std::move(layer));
                        }
                    };
                    if (rows_first) {
                        append_rows();
                        append_columns();
                    } else {
                        append_columns();
                        append_rows();
                    }
                    layer_candidates.push_back(std::move(layers));
                }
            }
        }
    }

    for (bool col_reverse_layers : {false, true}) {
        for (bool col_reverse_ops : {false, true}) {
            std::vector<RowOp> col = flatten_rowop_layers(column_layers, col_reverse_layers, col_reverse_ops);
            for (bool col_swap : {false, true}) {
                std::vector<RowOp> col_variant = col_swap ? swap_controls_targets(col) : col;
                for (bool row_reverse_ops : {false, true}) {
                    std::vector<RowOp> row = flatten_rowop_layers(row_layers, true, row_reverse_ops);
                    for (int transform = 0; transform < 3; ++transform) {
                        std::vector<RowOp> row_variant = row;
                        if (transform == 1) row_variant = transform_ops_by_permutation(row, p);
                        if (transform == 2) row_variant = transform_ops_by_permutation(row, inv_p);

                        std::vector<RowOp> seq;
                        seq.reserve(col_variant.size() + row_variant.size());
                        seq.insert(seq.end(), col_variant.begin(), col_variant.end());
                        seq.insert(seq.end(), row_variant.begin(), row_variant.end());
                        candidates.push_back(seq);

                        seq.clear();
                        seq.reserve(col_variant.size() + row_variant.size());
                        seq.insert(seq.end(), row_variant.begin(), row_variant.end());
                        seq.insert(seq.end(), col_variant.begin(), col_variant.end());
                        candidates.push_back(seq);
                    }
                }
            }
        }
    }

    bool have_best = false;
    Result best;
    auto layers_are_surface_packable = [&](const Result& item) {
        if (!g_logical_search_params.layer_matrix_vdp_hard) return true;
        if (g_layer_matrix_surface_layout == nullptr || g_layer_matrix_surface_mode != "vdp") {
            return false;
        }
        for (const Layer& layer : item.layers) {
            std::vector<RowOp> ops;
            ops.reserve(layer.ops.size());
            for (const RoutedCNOT& op : layer.ops) ops.emplace_back(op.control, op.target);
            std::vector<RoutedCNOT> routed;
            if (!pack_layer_ops_with_k_paths(
                    *g_layer_matrix_surface_layout,
                    g_layer_matrix_surface_mode,
                    ops,
                    routed)) {
                return false;
            }
        }
        return true;
    };
    auto reconstruction_key = [&](const Result& item) {
        return logical_result_priority_key(item);
    };
    for (const auto& layers : layer_candidates) {
        try {
            Result result = make_result_from_logical_layers(target_rows, layout, layers);
            if (!layers_are_surface_packable(result)) continue;
            if (!have_best || reconstruction_key(result) < reconstruction_key(best)) {
                best = std::move(result);
                have_best = true;
            }
        } catch (const std::exception&) {
        }
    }
    for (auto& ops : candidates) {
        try {
            Result result = make_result_from_logical_ops(target_rows, layout, std::move(ops));
            if (!layers_are_surface_packable(result)) continue;
            if (!have_best || reconstruction_key(result) < reconstruction_key(best)) {
                best = std::move(result);
                have_best = true;
            }
        } catch (const std::exception&) {
        }
    }
    if (!have_best) throw std::runtime_error("could not reconstruct logical circuit");
    return best;
}

static int permutation_excess(const Matrix& rows) {
    check_square_rows(rows);
    int excess = 0;
    for (Row row : rows) excess += std::abs(popcount64(row) - 1);
    const Matrix transposed = transpose_matrix(rows);
    for (Row col : transposed) excess += std::abs(popcount64(col) - 1);
    return excess;
}

static int identity_delta_rank(const Matrix& rows) {
    const int n = check_square_rows(rows);
    Matrix delta = rows;
    for (int i = 0; i < n; ++i) delta[static_cast<std::size_t>(i)] ^= Row(1) << i;
    return gf2_rank(std::move(delta));
}

static bool logical_depth_one_ops_on_side(
    const Matrix& residual,
    OpSide side,
    std::vector<RowOp>& ops
) {
    if (side == OpSide::Row) return logical_depth_one_ops(residual, ops);
    Matrix transposed = transpose_matrix(residual);
    return logical_depth_one_ops(transposed, ops);
}

static void apply_logical_layer_with_inverse(
    Matrix& rows,
    Matrix& inverse,
    OpSide side,
    const std::vector<RowOp>& ops
) {
    for (const RowOp& op : ops) {
        AlgebraicChoice choice{side, op.first, op.second, 0.0, 0.0};
        apply_algebraic_choice_with_inverse(rows, inverse, choice);
    }
}

struct LayerMatrixEdge {
    int control = 0;
    int target = 0;
    double side_delta = 0.0;
    double global_delta = 0.0;
    double side_after = 0.0;
    double global_after = 0.0;
    int permutation_excess_after = 0;
    int identity_rank_after = 0;
};

struct LayerMatrixCandidate {
    OpSide side = OpSide::Row;
    std::vector<RowOp> ops;
    Matrix residual_after;
    Matrix inverse_after;
    double global_after = 0.0;
    double side_after = 0.0;
    int permutation_excess_after = 0;
    int identity_rank_after = 0;
    int surface_pack_misses = 0;
    int surface_path_score = 0;
    int endpoint_goal_excess_after = 0;
    int endpoint_max_after = 0;
    long long endpoint_energy_after = 0;
};

static std::string layer_matrix_candidate_key(const LayerMatrixCandidate& candidate) {
    std::vector<RowOp> ops = candidate.ops;
    std::sort(ops.begin(), ops.end());
    std::string key;
    key.reserve(1 + 2 * ops.size());
    key.push_back(candidate.side == OpSide::Row ? 'r' : 'c');
    for (const RowOp& op : ops) {
        key.push_back(static_cast<char>(op.first));
        key.push_back(static_cast<char>(op.second));
    }
    return key;
}

static std::vector<LayerMatrixEdge> layer_matrix_edges(
    const Matrix& residual,
    const Matrix& residual_inverse,
    OpSide side,
    ShiFengCostKind kind,
    bool include_non_improving
) {
    const int n = check_square_rows(residual);
    const double current_global = shi_current_cost_with_inverse(residual, residual_inverse, kind);
    const double current_side = shi_side_cost_with_inverse(residual, residual_inverse, side, kind);
    const int current_permutation_excess = include_non_improving ? permutation_excess(residual) : 0;
    const int current_identity_rank = include_non_improving ? identity_delta_rank(residual) : 0;
    std::vector<LayerMatrixEdge> edges;
    edges.reserve(static_cast<std::size_t>(n * (n - 1)));
    for (int control = 0; control < n; ++control) {
        for (int target = 0; target < n; ++target) {
            if (control == target) continue;
            AlgebraicChoice choice{side, control, target, 0.0, 0.0};
            Matrix after = residual;
            Matrix after_inverse = residual_inverse;
            apply_algebraic_choice_with_inverse(after, after_inverse, choice);
            const double side_after = shi_side_cost_with_inverse(after, after_inverse, side, kind);
            const double global_after = shi_current_cost_with_inverse(after, after_inverse, kind);
            const double side_delta = current_side - side_after;
            const double global_delta = current_global - global_after;
            if (!include_non_improving && side_delta <= 1.0e-9 && global_delta <= 1.0e-9) continue;
            const int perm_after = permutation_excess(after);
            const int rank_after = identity_delta_rank(after);
            if (side_delta <= 1.0e-9 && global_delta <= 1.0e-9
                && g_logical_search_params.layer_matrix_nonimproving_slack <= 0.0) {
                if (perm_after >= current_permutation_excess && rank_after >= current_identity_rank) continue;
            }
            edges.push_back({
                control,
                target,
                side_delta,
                global_delta,
                side_after,
                global_after,
                perm_after,
                rank_after
            });
        }
    }
    return edges;
}

static std::pair<int, int> layer_matrix_surface_pack_score(const std::vector<RowOp>& ops) {
    if (g_layer_matrix_surface_layout == nullptr || g_layer_matrix_surface_mode.empty()) return {0, 0};
    const Layout& layout = *g_layer_matrix_surface_layout;
    if (g_logical_search_params.layer_matrix_vdp_hard) {
        if (g_layer_matrix_surface_mode != "vdp") {
            return {static_cast<int>(ops.size()), std::numeric_limits<int>::max() / 4};
        }

        std::vector<RowOp> canonical = ops;
        std::sort(canonical.begin(), canonical.end());
        std::ostringstream key_stream;
        key_stream << layout.data_rows << 'x' << layout.data_cols
                   << "|k=" << g_surface_final_pack_params.k_paths
                   << "|extra=" << g_surface_final_pack_params.path_extra
                   << "|limit=" << g_surface_final_pack_params.repack_limit << '|';
        for (const RowOp& op : canonical) key_stream << op.first << '>' << op.second << ',';
        const std::string key = key_stream.str();

        static std::mutex cache_mutex;
        static std::unordered_map<std::string, std::pair<int, int>> cache;
        {
            std::lock_guard<std::mutex> lock(cache_mutex);
            const auto it = cache.find(key);
            if (it != cache.end()) return it->second;
        }

        std::vector<RoutedCNOT> routed;
        std::pair<int, int> score{
            static_cast<int>(ops.size()),
            std::numeric_limits<int>::max() / 4
        };
        if (pack_layer_ops_with_k_paths(layout, "vdp", canonical, routed)) {
            int path_score = 0;
            for (const RoutedCNOT& op : routed) path_score += static_cast<int>(op.path.size());
            score = {0, path_score};
        }
        {
            std::lock_guard<std::mutex> lock(cache_mutex);
            cache.emplace(key, score);
        }
        return score;
    }

    auto score_order = [&](std::vector<RowOp> ordered) {
        Layer layer;
        layer.mode = g_layer_matrix_surface_mode;
        int routed_count = 0;
        int path_score = 0;
        for (const RowOp& op : ordered) {
            RoutedCNOT routed;
            if (!try_route_into_layer(layout, g_layer_matrix_surface_mode, layer, op.first, op.second, routed)) {
                continue;
            }
            path_score += static_cast<int>(routed.path.size());
            layer.ops.push_back(std::move(routed));
            ++routed_count;
        }
        return std::make_pair(static_cast<int>(ops.size()) - routed_count, path_score);
    };

    std::vector<RowOp> by_near = ops;
    std::sort(by_near.begin(), by_near.end(), [&](const RowOp& a, const RowOp& b) {
        const int da = layout.manhattan_between_qubits(a.first, a.second);
        const int db = layout.manhattan_between_qubits(b.first, b.second);
        if (da != db) return da < db;
        return a < b;
    });
    std::vector<RowOp> by_far = by_near;
    std::reverse(by_far.begin(), by_far.end());
    const auto near_score = score_order(std::move(by_near));
    const auto far_score = score_order(std::move(by_far));
    return std::min(near_score, far_score);
}

static std::vector<std::size_t> ordered_layer_edges(
    const std::vector<LayerMatrixEdge>& edges,
    int order_kind,
    std::mt19937_64& rng,
    const std::vector<int>& endpoint_load,
    bool endpoint_rank
) {
    std::vector<std::size_t> order(edges.size());
    for (std::size_t i = 0; i < edges.size(); ++i) order[i] = i;
    const bool endpoint_active = endpoint_rank
        || g_logical_search_params.layer_matrix_endpoint_first;
    auto endpoint_key = [&](const LayerMatrixEdge& edge) {
        const int control_load = endpoint_load.empty()
            ? 1
            : endpoint_load[static_cast<std::size_t>(edge.control)] + 1;
        const int target_load = endpoint_load.empty()
            ? 1
            : endpoint_load[static_cast<std::size_t>(edge.target)] + 1;
        const int projected_max = std::max(control_load, target_load);
        const int goal = g_logical_search_params.depth_goal;
        return std::make_tuple(
            goal > 0 ? std::max(0, projected_max - goal) : 0,
            projected_max,
            control_load + target_load
        );
    };
    auto deterministic_less = [&](std::size_t ai, std::size_t bi) {
        const auto& a = edges[ai];
        const auto& b = edges[bi];
        if (endpoint_active && order_kind == 3) {
            const auto ak = endpoint_key(a);
            const auto bk = endpoint_key(b);
            if (ak != bk) return ak < bk;
        }
        if (order_kind == 1 && std::abs(a.side_delta - b.side_delta) > 1.0e-9) {
            return a.side_delta > b.side_delta;
        }
        if (order_kind == 2) {
            const auto ak = std::make_tuple(a.permutation_excess_after, a.identity_rank_after);
            const auto bk = std::make_tuple(b.permutation_excess_after, b.identity_rank_after);
            if (ak != bk) return ak < bk;
        }
        if (std::abs(a.global_delta - b.global_delta) > 1.0e-9) return a.global_delta > b.global_delta;
        if (std::abs(a.side_delta - b.side_delta) > 1.0e-9) return a.side_delta > b.side_delta;
        if (std::abs(a.global_after - b.global_after) > 1.0e-9) return a.global_after < b.global_after;
        if (a.control != b.control) return a.control < b.control;
        return a.target < b.target;
    };
    std::sort(order.begin(), order.end(), deterministic_less);
    if (order_kind >= 3 && order.size() > 1) {
        const std::size_t hot = std::min<std::size_t>(
            order.size(),
            logical_depth_pressure_goal() ? 192 : 96
        );
        std::vector<std::pair<double, std::size_t>> scored;
        scored.reserve(hot);
        for (std::size_t i = 0; i < hot; ++i) {
            const auto& edge = edges[order[i]];
            const double jitter = static_cast<double>(rng() & 0xffffU) / 65536.0;
            double endpoint_penalty = 0.0;
            if (endpoint_active) {
                const auto key = endpoint_key(edge);
                endpoint_penalty =
                    4.0 * static_cast<double>(std::get<0>(key))
                    + 0.20 * static_cast<double>(std::get<1>(key))
                    + 0.02 * static_cast<double>(std::get<2>(key));
            }
            scored.emplace_back(edge.global_delta + 0.25 * edge.side_delta + jitter - endpoint_penalty, order[i]);
        }
        std::sort(scored.begin(), scored.end(), [&](const auto& a, const auto& b) {
            if (std::abs(a.first - b.first) > 1.0e-9) return a.first > b.first;
            return deterministic_less(a.second, b.second);
        });
        for (std::size_t i = 0; i < hot; ++i) order[i] = scored[i].second;
    }
    return order;
}

static bool build_layer_matrix_candidate(
    const Matrix& residual,
    const Matrix& residual_inverse,
    OpSide side,
    ShiFengCostKind kind,
    const std::vector<LayerMatrixEdge>& edges,
    const std::vector<std::size_t>& order,
    const std::vector<int>& endpoint_load,
    bool endpoint_rank,
    bool fill_layer,
    int max_layer_ops,
    LayerMatrixCandidate& candidate
) {
    const int n = check_square_rows(residual);
    const double current_global = shi_current_cost_with_inverse(residual, residual_inverse, kind);
    const double current_side = shi_side_cost_with_inverse(residual, residual_inverse, side, kind);
    std::vector<char> used(static_cast<std::size_t>(n), 0);
    std::vector<RowOp> ops;
    ops.reserve(static_cast<std::size_t>(n / 2));

    for (std::size_t edge_idx : order) {
        const auto& edge = edges[edge_idx];
        if (used[static_cast<std::size_t>(edge.control)]
            || used[static_cast<std::size_t>(edge.target)]) {
            continue;
        }
        if ((endpoint_rank || g_logical_search_params.layer_matrix_endpoint_first)
            && g_logical_search_params.depth_goal > 0) {
            if (endpoint_load[static_cast<std::size_t>(edge.control)] >= g_logical_search_params.depth_goal
                || endpoint_load[static_cast<std::size_t>(edge.target)] >= g_logical_search_params.depth_goal) {
                continue;
            }
        }
        if (!fill_layer && edge.side_delta <= 1.0e-9 && edge.global_delta <= 1.0e-9) continue;
        used[static_cast<std::size_t>(edge.control)] = 1;
        used[static_cast<std::size_t>(edge.target)] = 1;
        ops.emplace_back(edge.control, edge.target);
        if (max_layer_ops > 0 && static_cast<int>(ops.size()) >= max_layer_ops) break;
        if (static_cast<int>(ops.size()) == n / 2) break;
    }
    if (ops.empty()) return false;

    Matrix after = residual;
    Matrix after_inverse = residual_inverse;
    apply_logical_layer_with_inverse(after, after_inverse, side, ops);
    const double global_after = shi_current_cost_with_inverse(after, after_inverse, kind);
    const double side_after = shi_side_cost_with_inverse(after, after_inverse, side, kind);
    if (!(global_after + 1.0e-9 < current_global || side_after + 1.0e-9 < current_side)) {
        const double slack = g_logical_search_params.layer_matrix_nonimproving_slack;
        if (slack <= 0.0
            || (global_after > current_global + slack
                && side_after > current_side + slack)) {
            return false;
        }
    }

    candidate.side = side;
    candidate.ops = std::move(ops);
    candidate.residual_after = std::move(after);
    candidate.inverse_after = std::move(after_inverse);
    candidate.global_after = global_after;
    candidate.side_after = side_after;
    candidate.permutation_excess_after = permutation_excess(candidate.residual_after);
    candidate.identity_rank_after = identity_delta_rank(candidate.residual_after);
    const auto surface_score = layer_matrix_surface_pack_score(candidate.ops);
    candidate.surface_pack_misses = surface_score.first;
    candidate.surface_path_score = surface_score.second;
    if (g_logical_search_params.layer_matrix_vdp_hard
        && candidate.surface_pack_misses != 0) {
        return false;
    }
    if (!endpoint_load.empty()) {
        std::vector<int> projected_load = endpoint_load;
        long long endpoint_energy = 0;
        int endpoint_max = *std::max_element(projected_load.begin(), projected_load.end());
        for (const RowOp& op : candidate.ops) {
            const int control_load = projected_load[static_cast<std::size_t>(op.first)];
            const int target_load = projected_load[static_cast<std::size_t>(op.second)];
            endpoint_energy += static_cast<long long>(2 * control_load + 1);
            endpoint_energy += static_cast<long long>(2 * target_load + 1);
            ++projected_load[static_cast<std::size_t>(op.first)];
            ++projected_load[static_cast<std::size_t>(op.second)];
            endpoint_max = std::max(
                endpoint_max,
                std::max(
                    projected_load[static_cast<std::size_t>(op.first)],
                    projected_load[static_cast<std::size_t>(op.second)]
                )
            );
        }
        candidate.endpoint_max_after = endpoint_max;
        candidate.endpoint_energy_after = endpoint_energy;
        if (g_logical_search_params.depth_goal > 0) {
            candidate.endpoint_goal_excess_after =
                std::max(0, endpoint_max - g_logical_search_params.depth_goal);
        }
    }
    return true;
}

static std::vector<std::vector<std::size_t>> layer_matching_order_variants(
    const std::vector<LayerMatrixEdge>& edges,
    const std::vector<std::size_t>& base_order,
    std::mt19937_64& rng
) {
    const int variants = std::max(1, g_logical_search_params.layer_matrix_order_variants);
    std::vector<std::vector<std::size_t>> out;
    out.reserve(static_cast<std::size_t>(variants));
    out.push_back(base_order);
    if (variants <= 1 || base_order.size() <= 1) return out;

    const std::size_t hot = std::min<std::size_t>(
        base_order.size(),
        logical_depth_pressure_goal() ? 256 : 128
    );
    for (int variant = 1; variant < variants; ++variant) {
        std::vector<std::size_t> order = base_order;
        if (variant == 1) {
            std::rotate(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(hot / 3), order.begin() + static_cast<std::ptrdiff_t>(hot));
        } else if (variant == 2) {
            std::reverse(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(hot));
        } else {
            std::vector<std::size_t> prefix(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(hot));
            std::shuffle(prefix.begin(), prefix.end(), rng);
            std::copy(prefix.begin(), prefix.end(), order.begin());
        }

        const int endpoint_bias = variant % 4;
        if (endpoint_bias != 0) {
            const std::size_t prefix_size = std::min<std::size_t>(hot, 96);
            std::stable_sort(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(prefix_size), [&](std::size_t ai, std::size_t bi) {
                const auto& a = edges[ai];
                const auto& b = edges[bi];
                if (endpoint_bias == 1) {
                    const auto ak = std::make_tuple(std::min(a.control, a.target), std::max(a.control, a.target));
                    const auto bk = std::make_tuple(std::min(b.control, b.target), std::max(b.control, b.target));
                    if (ak != bk) return ak < bk;
                } else if (endpoint_bias == 2) {
                    const auto ak = std::make_tuple(std::max(a.control, a.target), std::min(a.control, a.target));
                    const auto bk = std::make_tuple(std::max(b.control, b.target), std::min(b.control, b.target));
                    if (ak != bk) return ak > bk;
                } else {
                    const auto ak = std::make_tuple((a.control + a.target) & 7, a.control, a.target);
                    const auto bk = std::make_tuple((b.control + b.target) & 7, b.control, b.target);
                    if (ak != bk) return ak < bk;
                }
                return ai < bi;
            });
        }
        out.push_back(std::move(order));
    }
    return out;
}

static std::vector<LayerMatrixCandidate> layer_matrix_candidates(
    const Matrix& residual,
    const Matrix& residual_inverse,
    const ShiFengSurfaceConfig& config,
    const std::vector<int>& endpoint_load,
    std::mt19937_64& rng
) {
    consume_fresh_work();
    const int n = check_square_rows(residual);
    std::vector<LayerMatrixCandidate> candidates;
    std::unordered_set<std::string> seen;
    const int side_count = config.allow_columns ? 2 : 1;
    const int random_orders = logical_depth_pressure_goal()
        ? std::max<int>(6, static_cast<int>(g_logical_search_params.choice_window) * 4)
            + std::max(0, g_logical_search_params.layer_matrix_extra_orders)
        : std::max<int>(2, static_cast<int>(g_logical_search_params.choice_window) * 2)
            + std::max(0, g_logical_search_params.layer_matrix_extra_orders / 2);
    for (int side_value = 0; side_value < side_count; ++side_value) {
        const OpSide side = side_value == 0 ? OpSide::Row : OpSide::Column;
        const bool include_non_improving = config.layer_matrix_rank
            && logical_depth_pressure_goal();
        const auto edges = layer_matrix_edges(
            residual,
            residual_inverse,
            side,
            config.cost_kind,
            include_non_improving
        );
        if (edges.empty()) continue;
        for (int order_kind = 0; order_kind < 3 + random_orders; ++order_kind) {
            const auto order = ordered_layer_edges(
                edges,
                order_kind,
                rng,
                endpoint_load,
                config.endpoint_balanced_rank
            );
            std::vector<int> layer_caps;
            const int cap_variants = std::min(
                n / 2,
                std::max(1, g_logical_search_params.layer_matrix_layer_cap_variants)
            );
            layer_caps.reserve(static_cast<std::size_t>(cap_variants));
            for (int offset = 0; offset < cap_variants; ++offset) {
                layer_caps.push_back(n / 2 - offset);
            }
            for (bool fill_layer : {false, true}) {
                for (int layer_cap : layer_caps) {
                    for (const auto& variant_order : layer_matching_order_variants(edges, order, rng)) {
                        LayerMatrixCandidate candidate;
                        if (!build_layer_matrix_candidate(
                                residual,
                                residual_inverse,
                                side,
                                config.cost_kind,
                                edges,
                                variant_order,
                                endpoint_load,
                                config.endpoint_balanced_rank,
                                fill_layer,
                                layer_cap,
                                candidate)) {
                            continue;
                        }
                        const std::string key = layer_matrix_candidate_key(candidate);
                        if (seen.insert(key).second) candidates.push_back(std::move(candidate));
                    }
                }
            }
        }
    }
    std::sort(candidates.begin(), candidates.end(), [](const LayerMatrixCandidate& a, const LayerMatrixCandidate& b) {
        if (logical_depth_pressure_goal() && g_logical_search_params.layer_matrix_endpoint_first) {
            const auto ak = std::make_tuple(
                a.endpoint_goal_excess_after,
                a.endpoint_max_after,
                std::min<long long>(a.endpoint_energy_after, std::numeric_limits<int>::max())
            );
            const auto bk = std::make_tuple(
                b.endpoint_goal_excess_after,
                b.endpoint_max_after,
                std::min<long long>(b.endpoint_energy_after, std::numeric_limits<int>::max())
            );
            if (ak != bk) return ak < bk;
        }
        const auto ak = logical_depth_pressure_goal()
            ? std::make_tuple(
                static_cast<double>(a.permutation_excess_after),
                static_cast<double>(a.identity_rank_after),
                static_cast<double>(a.surface_pack_misses),
                static_cast<double>(a.surface_path_score),
                a.global_after,
                -static_cast<double>(a.ops.size()),
                a.side_after,
                static_cast<int>(a.side)
            )
            : std::make_tuple(
                a.global_after,
                static_cast<double>(a.permutation_excess_after),
                static_cast<double>(a.identity_rank_after),
                static_cast<double>(a.surface_pack_misses),
                static_cast<double>(a.surface_path_score),
                -static_cast<double>(a.ops.size()),
                a.side_after,
                static_cast<int>(a.side)
            );
        const auto bk = logical_depth_pressure_goal()
            ? std::make_tuple(
                static_cast<double>(b.permutation_excess_after),
                static_cast<double>(b.identity_rank_after),
                static_cast<double>(b.surface_pack_misses),
                static_cast<double>(b.surface_path_score),
                b.global_after,
                -static_cast<double>(b.ops.size()),
                b.side_after,
                static_cast<int>(b.side)
            )
            : std::make_tuple(
                b.global_after,
                static_cast<double>(b.permutation_excess_after),
                static_cast<double>(b.identity_rank_after),
                static_cast<double>(b.surface_pack_misses),
                static_cast<double>(b.surface_path_score),
                -static_cast<double>(b.ops.size()),
                b.side_after,
                static_cast<int>(b.side)
            );
        if (ak != bk) return ak < bk;
        return a.ops < b.ops;
    });
    const std::size_t cap = std::min<std::size_t>(
        candidates.size(),
        std::max<int>(1, g_logical_search_params.layer_matrix_candidate_multiplier)
        * (logical_depth_pressure_goal()
            ? std::max<std::size_t>(16, 16 * g_logical_search_params.choice_window)
            : std::max<std::size_t>(8, 8 * g_logical_search_params.choice_window))
    );
    if (candidates.size() > cap) candidates.resize(cap);
    return candidates;
}

struct LayerMatrixState {
    Matrix residual;
    Matrix residual_inverse;
    std::vector<std::vector<RowOp>> row_layers;
    std::vector<std::vector<RowOp>> column_layers;
    std::vector<int> endpoint_load;
    double cost = 0.0;
    int permutation_excess_value = 0;
    int identity_rank_value = 0;
    int endpoint_max_load = 0;
    long long endpoint_energy = 0;
    int cnot_count = 0;
};

static std::string layer_matrix_state_seen_key(const LayerMatrixState& state) {
    std::string key = matrix_key(state.residual);
    if (!logical_depth_pressure_goal()) return key;

    key.push_back('|');
    const auto append_small = [&](int value) {
        key.push_back(static_cast<char>(std::max(0, std::min(255, value))));
    };
    append_small(static_cast<int>(state.row_layers.size()));
    append_small(static_cast<int>(state.column_layers.size()));
    append_small(state.endpoint_max_load);
    append_small(state.cnot_count);
    for (int load : state.endpoint_load) append_small(load);
    return key;
}

static std::tuple<int, int, int, long long, double, int, int, int> layer_matrix_state_key(
    const LayerMatrixState& state,
    bool endpoint_rank
) {
    const bool endpoint_active = endpoint_rank
        || g_logical_search_params.layer_matrix_endpoint_first;
    const int endpoint_goal_excess = endpoint_active && g_logical_search_params.depth_goal > 0
        ? std::max(0, state.endpoint_max_load - g_logical_search_params.depth_goal)
        : 0;
    if (logical_depth_pressure_goal()) {
        if (g_logical_search_params.layer_matrix_endpoint_first) {
            return std::make_tuple(
                is_permutation_matrix(state.residual) ? 0 : 1,
                endpoint_goal_excess,
                endpoint_active ? state.endpoint_max_load : 0,
                endpoint_active ? static_cast<int>(std::min<long long>(state.endpoint_energy, std::numeric_limits<int>::max())) : 0,
                state.permutation_excess_value,
                static_cast<long long>(state.identity_rank_value),
                state.cost,
                -state.cnot_count
            );
        }
        return std::make_tuple(
            is_permutation_matrix(state.residual) ? 0 : 1,
            endpoint_goal_excess,
            state.permutation_excess_value,
            static_cast<long long>(state.identity_rank_value),
            state.cost,
            endpoint_active ? state.endpoint_max_load : 0,
            endpoint_active ? static_cast<int>(std::min<long long>(state.endpoint_energy, std::numeric_limits<int>::max())) : 0,
            -state.cnot_count
        );
    }
    return std::make_tuple(
        is_permutation_matrix(state.residual) ? 0 : 1,
        endpoint_goal_excess,
        endpoint_active ? state.endpoint_max_load : 0,
        endpoint_active ? state.endpoint_energy : 0,
        state.cost,
        state.permutation_excess_value,
        state.identity_rank_value,
        state.cnot_count
    );
}

static std::size_t scaled_layer_matrix_width(std::size_t base) {
    const int multiplier = std::max(1, g_logical_search_params.layer_matrix_beam_multiplier);
    const std::size_t scaled = base * static_cast<std::size_t>(multiplier);
    const std::size_t hard_cap = logical_depth_pressure_goal() ? 4096 : 1024;
    return std::min<std::size_t>(hard_cap, scaled);
}

static Result synthesize_layer_matrix_logical_once(
    const Matrix& target_rows,
    const Layout& layout,
    const ShiFengSurfaceConfig& config,
    std::uint64_t tie_seed,
    int max_layers
) {
    check_square_rows(target_rows);
    LayerMatrixState initial;
    initial.residual = target_rows;
    initial.residual_inverse = invert_matrix(initial.residual);
    initial.cost = shi_current_cost_with_inverse(initial.residual, initial.residual_inverse, config.cost_kind);
    initial.permutation_excess_value = permutation_excess(initial.residual);
    initial.identity_rank_value = identity_delta_rank(initial.residual);
    initial.endpoint_load.assign(static_cast<std::size_t>(layout.n), 0);
    initial.endpoint_max_load = 0;
    initial.endpoint_energy = 0;
    initial.cnot_count = 0;

    std::mt19937_64 rng(tie_seed == std::numeric_limits<std::uint64_t>::max()
        ? 0x4d595df4d0f33173ULL
        : tie_seed);
    std::vector<LayerMatrixState> beam{std::move(initial)};
    const std::size_t beam_width = std::max<std::size_t>(
        8,
        scaled_layer_matrix_width(logical_depth_pressure_goal()
            ? std::min<std::size_t>(256, 32 * std::max<std::size_t>(1, g_logical_search_params.choice_window))
            : std::min<std::size_t>(96, 16 * std::max<std::size_t>(1, g_logical_search_params.choice_window)))
    );
    const int configured_depth_limit = g_logical_search_params.layer_matrix_goal_only
        ? std::max(1, g_logical_search_params.depth_goal)
        : std::max(1, std::max(g_logical_search_params.depth_goal, 12));
    const int depth_limit = std::min(max_layers, configured_depth_limit);
    bool have_best = false;
    Result best;

    auto keep_if_complete = [&](const LayerMatrixState& state) {
        if (!is_permutation_matrix(state.residual)) return;
        Result result = reconstruct_logical_search_result(
            target_rows,
            layout,
            state.residual,
            state.row_layers,
            state.column_layers
        );
        if (!have_best || logical_result_priority_key(result) < logical_result_priority_key(best)) {
            best = std::move(result);
            have_best = true;
        }
    };
    auto append_layer_to_state = [&](LayerMatrixState& state, OpSide side, const std::vector<RowOp>& ops) {
        if (side == OpSide::Row) {
            state.row_layers.push_back(ops);
        } else {
            state.column_layers.push_back(ops);
        }
        if (state.endpoint_load.empty()) return;
        for (const RowOp& op : ops) {
            const int control_load = state.endpoint_load[static_cast<std::size_t>(op.first)];
            const int target_load = state.endpoint_load[static_cast<std::size_t>(op.second)];
            state.endpoint_energy += static_cast<long long>(2 * control_load + 1);
            state.endpoint_energy += static_cast<long long>(2 * target_load + 1);
            ++state.endpoint_load[static_cast<std::size_t>(op.first)];
            ++state.endpoint_load[static_cast<std::size_t>(op.second)];
            state.endpoint_max_load = std::max(
                state.endpoint_max_load,
                std::max(
                    state.endpoint_load[static_cast<std::size_t>(op.first)],
                    state.endpoint_load[static_cast<std::size_t>(op.second)]
                )
            );
        }
        state.cnot_count += static_cast<int>(ops.size());
    };
    auto try_depth_one_tail = [&](const LayerMatrixState& state) {
        const int tail_side_count = config.allow_columns ? 2 : 1;
        for (int side_value = 0; side_value < tail_side_count; ++side_value) {
            const OpSide side = side_value == 0 ? OpSide::Row : OpSide::Column;
            std::vector<RowOp> tail;
            if (!logical_depth_one_ops_on_side(state.residual, side, tail)) continue;
            LayerMatrixState child = state;
            apply_logical_layer_with_inverse(child.residual, child.residual_inverse, side, tail);
            if (!is_permutation_matrix(child.residual)) continue;
            append_layer_to_state(child, side, tail);
            keep_if_complete(child);
        }
    };
    auto try_depth_two_tail = [&](
        const LayerMatrixState& state,
        const std::vector<LayerMatrixCandidate>* precomputed_candidates
    ) {
        std::vector<LayerMatrixCandidate> generated_candidates;
        const std::vector<LayerMatrixCandidate>* candidates = precomputed_candidates;
        if (candidates == nullptr) {
            generated_candidates = layer_matrix_candidates(
                state.residual,
                state.residual_inverse,
                config,
                state.endpoint_load,
                rng
            );
            candidates = &generated_candidates;
        }
        for (const auto& candidate : *candidates) {
            LayerMatrixState child = state;
            child.residual = candidate.residual_after;
            child.residual_inverse = candidate.inverse_after;
            append_layer_to_state(child, candidate.side, candidate.ops);
            const int tail_side_count = config.allow_columns ? 2 : 1;
            for (int side_value = 0; side_value < tail_side_count; ++side_value) {
                const OpSide side = side_value == 0 ? OpSide::Row : OpSide::Column;
                std::vector<RowOp> tail;
                if (!logical_depth_one_ops_on_side(candidate.residual_after, side, tail)) continue;
                LayerMatrixState tail_child = child;
                apply_logical_layer_with_inverse(tail_child.residual, tail_child.residual_inverse, side, tail);
                if (!is_permutation_matrix(tail_child.residual)) continue;
                append_layer_to_state(tail_child, side, tail);
                keep_if_complete(tail_child);
            }
        }
    };
    auto try_depth_three_tail = [&](
        const LayerMatrixState& state,
        const std::vector<LayerMatrixCandidate>* precomputed_candidates
    ) {
        std::vector<LayerMatrixCandidate> generated_candidates;
        const std::vector<LayerMatrixCandidate>* first_candidates = precomputed_candidates;
        if (first_candidates == nullptr) {
            generated_candidates = layer_matrix_candidates(
                state.residual,
                state.residual_inverse,
                config,
                state.endpoint_load,
                rng
            );
            first_candidates = &generated_candidates;
        }
        const std::size_t first_cap = std::min<std::size_t>(
            first_candidates->size(),
            static_cast<std::size_t>(std::max(0, g_logical_search_params.layer_matrix_tail_cap))
        );
        for (std::size_t i = 0; i < first_cap; ++i) {
            const auto& first = (*first_candidates)[i];
            LayerMatrixState child = state;
            child.residual = first.residual_after;
            child.residual_inverse = first.inverse_after;
            append_layer_to_state(child, first.side, first.ops);
            try_depth_two_tail(child, nullptr);
        }
    };

    for (int depth = 0; depth <= depth_limit; ++depth) {
        std::vector<std::vector<LayerMatrixCandidate>> cached_state_candidates;
        if (g_logical_search_params.layer_matrix_reuse_state_candidates
            && depth < depth_limit) {
            cached_state_candidates.reserve(beam.size());
            for (const auto& state : beam) {
                cached_state_candidates.push_back(layer_matrix_candidates(
                    state.residual,
                    state.residual_inverse,
                    config,
                    state.endpoint_load,
                    rng
                ));
            }
        }
        for (std::size_t state_index = 0; state_index < beam.size(); ++state_index) {
            const auto& state = beam[state_index];
            const auto* cached_candidates = cached_state_candidates.empty()
                ? nullptr
                : &cached_state_candidates[state_index];
            keep_if_complete(state);
            if (depth < depth_limit) try_depth_one_tail(state);
            if (depth + 1 < depth_limit) try_depth_two_tail(state, cached_candidates);
            if (g_logical_search_params.layer_matrix_tail_cap > 0
                && depth + 2 < depth_limit) {
                try_depth_three_tail(state, cached_candidates);
            }
        }
        if (have_best && logical_result_meets_goals(best)) return best;
        if (depth == depth_limit) break;

        std::vector<LayerMatrixState> next;
        std::unordered_set<std::string> seen;
        for (std::size_t state_index = 0; state_index < beam.size(); ++state_index) {
            const auto& state = beam[state_index];
            std::vector<LayerMatrixCandidate> generated_candidates;
            const std::vector<LayerMatrixCandidate>* candidates = cached_state_candidates.empty()
                ? nullptr
                : &cached_state_candidates[state_index];
            if (candidates == nullptr) {
                generated_candidates = layer_matrix_candidates(
                    state.residual,
                    state.residual_inverse,
                    config,
                    state.endpoint_load,
                    rng
                );
                candidates = &generated_candidates;
            }
            for (const auto& candidate : *candidates) {
                LayerMatrixState child;
                child.residual = candidate.residual_after;
                child.residual_inverse = candidate.inverse_after;
                child.row_layers = state.row_layers;
                child.column_layers = state.column_layers;
                child.endpoint_load = state.endpoint_load;
                child.endpoint_max_load = state.endpoint_max_load;
                child.endpoint_energy = state.endpoint_energy;
                if (candidate.side == OpSide::Row) {
                    child.row_layers.push_back(candidate.ops);
                } else {
                    child.column_layers.push_back(candidate.ops);
                }
                for (const RowOp& op : candidate.ops) {
                    const int control_load = child.endpoint_load[static_cast<std::size_t>(op.first)];
                    const int target_load = child.endpoint_load[static_cast<std::size_t>(op.second)];
                    child.endpoint_energy += static_cast<long long>(2 * control_load + 1);
                    child.endpoint_energy += static_cast<long long>(2 * target_load + 1);
                    ++child.endpoint_load[static_cast<std::size_t>(op.first)];
                    ++child.endpoint_load[static_cast<std::size_t>(op.second)];
                    child.endpoint_max_load = std::max(
                        child.endpoint_max_load,
                        std::max(
                            child.endpoint_load[static_cast<std::size_t>(op.first)],
                            child.endpoint_load[static_cast<std::size_t>(op.second)]
                        )
                    );
                }
                child.cost = candidate.global_after;
                child.permutation_excess_value = candidate.permutation_excess_after;
                child.identity_rank_value = candidate.identity_rank_after;
                child.cnot_count = state.cnot_count + static_cast<int>(candidate.ops.size());
                const std::string key = layer_matrix_state_seen_key(child);
                if (seen.insert(key).second) next.push_back(std::move(child));
            }
        }
        if (next.empty()) break;
        std::sort(next.begin(), next.end(), [&](const LayerMatrixState& a, const LayerMatrixState& b) {
            const auto ak = layer_matrix_state_key(a, config.endpoint_balanced_rank);
            const auto bk = layer_matrix_state_key(b, config.endpoint_balanced_rank);
            if (ak != bk) return ak < bk;
            return matrix_key(a.residual) < matrix_key(b.residual);
        });
        if (next.size() > beam_width) next.resize(beam_width);
        beam = std::move(next);
    }
    if (have_best) return best;
    throw std::runtime_error("layer-matrix logical search did not reach a permutation");
}

static Result synthesize_inverse_row_layer_matrix_logical_once(
    const Matrix& target_rows,
    const Layout& layout,
    const ShiFengSurfaceConfig& base_config,
    std::uint64_t tie_seed,
    int max_layers
) {
    check_square_rows(target_rows);
    ShiFengSurfaceConfig config = base_config;
    config.allow_columns = false;

    LayerMatrixState initial;
    initial.residual = invert_matrix(target_rows);
    initial.residual_inverse = target_rows;
    initial.cost = shi_current_cost_with_inverse(initial.residual, initial.residual_inverse, config.cost_kind);
    initial.permutation_excess_value = permutation_excess(initial.residual);
    initial.identity_rank_value = identity_delta_rank(initial.residual);
    initial.endpoint_load.assign(static_cast<std::size_t>(layout.n), 0);
    initial.endpoint_max_load = 0;
    initial.endpoint_energy = 0;
    initial.cnot_count = 0;

    std::mt19937_64 rng(tie_seed == std::numeric_limits<std::uint64_t>::max()
        ? 0xa0761d6478bd642fULL
        : (tie_seed ^ 0xa0761d6478bd642fULL));
    std::vector<LayerMatrixState> beam{std::move(initial)};
    const std::size_t beam_width = std::max<std::size_t>(
        8,
        scaled_layer_matrix_width(logical_depth_pressure_goal()
            ? std::min<std::size_t>(384, 32 * std::max<std::size_t>(1, g_logical_search_params.choice_window))
            : std::min<std::size_t>(128, 16 * std::max<std::size_t>(1, g_logical_search_params.choice_window)))
    );
    const int depth_limit = std::min(
        max_layers,
        std::max(1, g_logical_search_params.depth_goal)
    );
    bool have_best = false;
    Result best;

    auto append_row_layer_to_state = [&](LayerMatrixState& state, const std::vector<RowOp>& ops) {
        state.row_layers.push_back(ops);
        for (const RowOp& op : ops) {
            const int control_load = state.endpoint_load[static_cast<std::size_t>(op.first)];
            const int target_load = state.endpoint_load[static_cast<std::size_t>(op.second)];
            state.endpoint_energy += static_cast<long long>(2 * control_load + 1);
            state.endpoint_energy += static_cast<long long>(2 * target_load + 1);
            ++state.endpoint_load[static_cast<std::size_t>(op.first)];
            ++state.endpoint_load[static_cast<std::size_t>(op.second)];
            state.endpoint_max_load = std::max(
                state.endpoint_max_load,
                std::max(
                    state.endpoint_load[static_cast<std::size_t>(op.first)],
                    state.endpoint_load[static_cast<std::size_t>(op.second)]
                )
            );
        }
        state.cnot_count += static_cast<int>(ops.size());
    };
    auto keep_if_complete = [&](const LayerMatrixState& state) {
        if (!is_permutation_matrix(state.residual)) return;
        Result result = make_result_from_logical_layers(target_rows, layout, state.row_layers);
        if (!have_best || logical_result_priority_key(result) < logical_result_priority_key(best)) {
            best = std::move(result);
            have_best = true;
        }
    };
    auto try_depth_one_tail = [&](const LayerMatrixState& state) {
        std::vector<RowOp> tail;
        if (!logical_depth_one_ops(state.residual, tail)) return;
        LayerMatrixState child = state;
        apply_logical_layer_with_inverse(child.residual, child.residual_inverse, OpSide::Row, tail);
        if (!is_permutation_matrix(child.residual)) return;
        append_row_layer_to_state(child, tail);
        keep_if_complete(child);
    };
    auto try_depth_two_tail = [&](
        const LayerMatrixState& state,
        const std::vector<LayerMatrixCandidate>* precomputed_candidates
    ) {
        std::vector<LayerMatrixCandidate> generated_candidates;
        const std::vector<LayerMatrixCandidate>* candidates = precomputed_candidates;
        if (candidates == nullptr) {
            generated_candidates = layer_matrix_candidates(
                state.residual,
                state.residual_inverse,
                config,
                state.endpoint_load,
                rng
            );
            candidates = &generated_candidates;
        }
        for (const auto& candidate : *candidates) {
            LayerMatrixState child = state;
            child.residual = candidate.residual_after;
            child.residual_inverse = candidate.inverse_after;
            append_row_layer_to_state(child, candidate.ops);
            std::vector<RowOp> tail;
            if (!logical_depth_one_ops(child.residual, tail)) continue;
            apply_logical_layer_with_inverse(child.residual, child.residual_inverse, OpSide::Row, tail);
            if (!is_permutation_matrix(child.residual)) continue;
            append_row_layer_to_state(child, tail);
            keep_if_complete(child);
        }
    };
    auto try_depth_three_tail = [&](
        const LayerMatrixState& state,
        const std::vector<LayerMatrixCandidate>* precomputed_candidates
    ) {
        std::vector<LayerMatrixCandidate> generated_candidates;
        const std::vector<LayerMatrixCandidate>* first_candidates = precomputed_candidates;
        if (first_candidates == nullptr) {
            generated_candidates = layer_matrix_candidates(
                state.residual,
                state.residual_inverse,
                config,
                state.endpoint_load,
                rng
            );
            first_candidates = &generated_candidates;
        }
        const std::size_t first_cap = std::min<std::size_t>(
            first_candidates->size(),
            static_cast<std::size_t>(std::max(0, g_logical_search_params.layer_matrix_tail_cap))
        );
        for (std::size_t i = 0; i < first_cap; ++i) {
            const auto& first = (*first_candidates)[i];
            LayerMatrixState child = state;
            child.residual = first.residual_after;
            child.residual_inverse = first.inverse_after;
            append_row_layer_to_state(child, first.ops);
            try_depth_two_tail(child, nullptr);
        }
    };

    for (int depth = 0; depth <= depth_limit; ++depth) {
        std::vector<std::vector<LayerMatrixCandidate>> cached_state_candidates;
        if (g_logical_search_params.layer_matrix_reuse_state_candidates
            && depth < depth_limit) {
            cached_state_candidates.reserve(beam.size());
            for (const auto& state : beam) {
                cached_state_candidates.push_back(layer_matrix_candidates(
                    state.residual,
                    state.residual_inverse,
                    config,
                    state.endpoint_load,
                    rng
                ));
            }
        }
        for (std::size_t state_index = 0; state_index < beam.size(); ++state_index) {
            const auto& state = beam[state_index];
            const auto* cached_candidates = cached_state_candidates.empty()
                ? nullptr
                : &cached_state_candidates[state_index];
            keep_if_complete(state);
            if (depth < depth_limit) try_depth_one_tail(state);
            if (depth + 1 < depth_limit) try_depth_two_tail(state, cached_candidates);
            if (g_logical_search_params.layer_matrix_tail_cap > 0
                && depth + 2 < depth_limit) {
                try_depth_three_tail(state, cached_candidates);
            }
        }
        if (have_best && logical_result_meets_goals(best)) return best;
        if (depth == depth_limit) break;

        std::vector<LayerMatrixState> next;
        std::unordered_set<std::string> seen;
        for (std::size_t state_index = 0; state_index < beam.size(); ++state_index) {
            const auto& state = beam[state_index];
            std::vector<LayerMatrixCandidate> generated_candidates;
            const std::vector<LayerMatrixCandidate>* candidates = cached_state_candidates.empty()
                ? nullptr
                : &cached_state_candidates[state_index];
            if (candidates == nullptr) {
                generated_candidates = layer_matrix_candidates(
                    state.residual,
                    state.residual_inverse,
                    config,
                    state.endpoint_load,
                    rng
                );
                candidates = &generated_candidates;
            }
            for (const auto& candidate : *candidates) {
                LayerMatrixState child = state;
                child.residual = candidate.residual_after;
                child.residual_inverse = candidate.inverse_after;
                append_row_layer_to_state(child, candidate.ops);
                child.cost = candidate.global_after;
                child.permutation_excess_value = candidate.permutation_excess_after;
                child.identity_rank_value = candidate.identity_rank_after;
                const std::string key = layer_matrix_state_seen_key(child);
                if (seen.insert(key).second) next.push_back(std::move(child));
            }
        }
        if (next.empty()) break;
        std::sort(next.begin(), next.end(), [&](const LayerMatrixState& a, const LayerMatrixState& b) {
            const auto ak = layer_matrix_state_key(a, config.endpoint_balanced_rank);
            const auto bk = layer_matrix_state_key(b, config.endpoint_balanced_rank);
            if (ak != bk) return ak < bk;
            return matrix_key(a.residual) < matrix_key(b.residual);
        });
        if (next.size() > beam_width) next.resize(beam_width);
        beam = std::move(next);
    }
    if (have_best) return best;
    throw std::runtime_error("inverse row layer-matrix search did not reach a permutation");
}

static Result synthesize_shi_feng_logical_once(
    const Matrix& target_rows,
    const Layout& layout,
    const ShiFengSurfaceConfig& config,
    std::uint64_t tie_seed,
    int max_layers = 10000
) {
    if (config.layer_matrix_rank) {
        bool have_best = false;
        Result best;
        if (logical_depth_pressure_goal()) {
            try {
                best = synthesize_inverse_row_layer_matrix_logical_once(
                    target_rows,
                    layout,
                    config,
                    tie_seed,
                    max_layers
                );
                have_best = true;
                if (logical_result_meets_goals(best)) return best;
            } catch (const std::exception&) {
            }
        }
        try {
            Result result = synthesize_layer_matrix_logical_once(
                target_rows,
                layout,
                config,
                tie_seed,
                max_layers
            );
            if (!have_best || logical_result_priority_key(result) < logical_result_priority_key(best)) {
                best = std::move(result);
                have_best = true;
            }
        } catch (const std::exception&) {
        }
        if (have_best) return best;
        throw std::runtime_error("layer-matrix logical search did not reach a permutation");
    }
    const int n = check_square_rows(target_rows);
    Matrix residual = target_rows;
    Matrix residual_inverse = invert_matrix(residual);
    std::vector<std::vector<RowOp>> row_layers;
    std::vector<std::vector<RowOp>> column_layers;
    Layer current_row_layer{"logical", {}, 0};
    Layer current_column_layer{"logical", {}, 0};
    const int step_budget = std::max(128, 16 * n * n);
    int closed_layers = 0;
    std::vector<int> endpoint_load(static_cast<std::size_t>(n), 0);
    int endpoint_max_load = 0;
    long long endpoint_energy = 0;
    std::mt19937_64 tie_rng(tie_seed);

    auto flush_search_layers = [&]() {
        bool flushed = false;
        if (!current_row_layer.ops.empty()) {
            std::vector<RowOp> layer_ops;
            layer_ops.reserve(current_row_layer.ops.size());
            for (const auto& op : current_row_layer.ops) layer_ops.emplace_back(op.control, op.target);
            row_layers.push_back(std::move(layer_ops));
            current_row_layer.ops.clear();
            current_row_layer.logical_depth = 0;
            ++closed_layers;
            flushed = true;
        }
        if (!current_column_layer.ops.empty()) {
            std::vector<RowOp> layer_ops;
            layer_ops.reserve(current_column_layer.ops.size());
            for (const auto& op : current_column_layer.ops) layer_ops.emplace_back(op.control, op.target);
            column_layers.push_back(std::move(layer_ops));
            current_column_layer.ops.clear();
            current_column_layer.logical_depth = 0;
            ++closed_layers;
            flushed = true;
        }
        return flushed;
    };

    for (int iter = 0; iter < step_budget; ++iter) {
        if (is_permutation_matrix(residual)) break;
        if (closed_layers >= max_layers) break;

        std::vector<RowOp> final_depth_one;
        if (logical_depth_one_ops(residual, final_depth_one)) {
            bool fits_current_row = true;
            std::vector<RoutedCNOT> routed_ops;
            routed_ops.reserve(final_depth_one.size());
            Layer trial_layer = current_row_layer;
            for (const RowOp& op : final_depth_one) {
                RoutedCNOT routed;
                int depth_delta = 0;
                if (!try_route_into_search_node(layout, "logical", trial_layer, op.first, op.second, routed, depth_delta)) {
                    fits_current_row = false;
                    break;
                }
                append_to_search_node(trial_layer, layout, "logical", routed);
                routed_ops.push_back(std::move(routed));
            }
            if (!fits_current_row) {
                flush_search_layers();
            }
            for (const RowOp& op : final_depth_one) {
                RoutedCNOT routed;
                int depth_delta = 0;
                if (!try_route_into_search_node(layout, "logical", current_row_layer, op.first, op.second, routed, depth_delta)) {
                    throw std::runtime_error("logical depth-one layer is not parallel");
                }
                row_add(residual, op.first, op.second);
                append_to_search_node(current_row_layer, layout, "logical", std::move(routed));
            }
            break;
        }

        const double current_cost = shi_current_cost_with_inverse(residual, residual_inverse, config.cost_kind);
        std::vector<ShiFengChoice> choices;
        choices.reserve(static_cast<std::size_t>(2 * n * (n - 1)));
        Matrix residual_transposed;
        Matrix inverse_transposed;
        if (config.residual_demand_rank) {
            residual_transposed = transpose_matrix(residual);
            inverse_transposed = transpose_matrix(residual_inverse);
        }

        auto line_demand = [](Row line) {
            return std::max(0, popcount64(line) - 1);
        };

        auto residual_endpoint_projection = [&](OpSide side, int control, int target) {
            constexpr int demand_divisor = 4;
            const Row control_bit = Row(1) << control;
            const Row target_bit = Row(1) << target;
            int endpoint_max = 0;
            long long demand_energy = 0;
            for (int q = 0; q < n; ++q) {
                const Row q_bit = Row(1) << q;
                Row residual_row = residual[static_cast<std::size_t>(q)];
                Row residual_col = residual_transposed[static_cast<std::size_t>(q)];
                Row inverse_row = residual_inverse[static_cast<std::size_t>(q)];
                Row inverse_col = inverse_transposed[static_cast<std::size_t>(q)];

                if (side == OpSide::Row) {
                    if (q == target) residual_row ^= residual[static_cast<std::size_t>(control)];
                    if ((residual[static_cast<std::size_t>(control)] & q_bit) != 0) {
                        residual_col ^= target_bit;
                    }
                    if ((residual_inverse[static_cast<std::size_t>(q)] & target_bit) != 0) {
                        inverse_row ^= control_bit;
                    }
                    if (q == control) inverse_col ^= inverse_transposed[static_cast<std::size_t>(target)];
                } else {
                    if ((residual[static_cast<std::size_t>(q)] & target_bit) != 0) {
                        residual_row ^= control_bit;
                    }
                    if (q == control) residual_col ^= residual_transposed[static_cast<std::size_t>(target)];
                    if (q == target) inverse_row ^= residual_inverse[static_cast<std::size_t>(control)];
                    if ((residual_inverse[static_cast<std::size_t>(control)] & q_bit) != 0) {
                        inverse_col ^= target_bit;
                    }
                }

                const int demand =
                    line_demand(residual_row)
                    + line_demand(residual_col)
                    + line_demand(inverse_row)
                    + line_demand(inverse_col);
                const int projected = endpoint_load[static_cast<std::size_t>(q)]
                    + (demand + demand_divisor - 1) / demand_divisor;
                endpoint_max = std::max(endpoint_max, projected);
                demand_energy += static_cast<long long>(projected) * projected;
            }
            return std::make_tuple(
                std::max(0, endpoint_max - g_logical_search_params.depth_goal),
                endpoint_max,
                demand_energy
            );
        };

        auto consider = [&](OpSide side, int control, int target) {
            Layer& node = (side == OpSide::Row) ? current_row_layer : current_column_layer;
            RoutedCNOT routed;
            int depth_delta = 0;
            if (!try_route_into_search_node(layout, "logical", node, control, target, routed, depth_delta)) {
                return;
            }

            AlgebraicChoice choice{side, control, target, 0.0, 0.0};
            double side_after = 0.0;
            double global_after = 0.0;
            if (config.rank_global_cost) {
                const ShiFengAfterCosts after_costs = shi_after_costs_with_inverse(
                    residual,
                    residual_inverse,
                    choice,
                    config.cost_kind
                );
                side_after = after_costs.side_after;
                global_after = after_costs.global_after;
            } else {
                side_after = shi_after_side_cost_with_inverse(
                    residual,
                    residual_inverse,
                    choice,
                    config.cost_kind
                );
                global_after = side_after;
            }
            if (!(side_after + 1.0e-9 < current_cost)) return;

            const double delta = current_cost - global_after;
            const double side_bias = (side == OpSide::Column) ? 0.0005 : 0.0;
            const double base_rank_cost = config.rank_global_cost ? global_after : side_after;
            const double rank_cost_after =
                base_rank_cost
                + config.logical_layer_penalty * static_cast<double>(depth_delta);
            const int next_control_load = endpoint_load[static_cast<std::size_t>(control)] + 1;
            const int next_target_load = endpoint_load[static_cast<std::size_t>(target)] + 1;
            const int endpoint_max_after = std::max(
                endpoint_max_load,
                std::max(next_control_load, next_target_load)
            );
            const long long endpoint_energy_after =
                endpoint_energy
                + static_cast<long long>(2 * endpoint_load[static_cast<std::size_t>(control)] + 1)
                + static_cast<long long>(2 * endpoint_load[static_cast<std::size_t>(target)] + 1);
            const int endpoint_goal_excess_after = std::max(
                0,
                endpoint_max_after - g_logical_search_params.depth_goal
            );
            auto residual_endpoint_key = std::make_tuple(0, 0, 0LL);
            if (config.residual_demand_rank) {
                residual_endpoint_key = residual_endpoint_projection(side, control, target);
            }
            choices.push_back({
                side,
                control,
                target,
                side_after,
                rank_cost_after,
                delta + side_bias,
                depth_delta,
                endpoint_goal_excess_after,
                endpoint_max_after,
                endpoint_energy_after,
                std::get<0>(residual_endpoint_key),
                std::get<1>(residual_endpoint_key),
                std::get<2>(residual_endpoint_key),
                std::move(routed)
            });
        };

        for (int control = 0; control < n; ++control) {
            for (int target = 0; target < n; ++target) {
                if (control == target) continue;
                consider(OpSide::Row, control, target);
                if (config.allow_columns) consider(OpSide::Column, control, target);
            }
        }

        if (choices.empty()) {
            if (flush_search_layers()) continue;
            break;
        }

        auto choice_rank_less = [&](const ShiFengChoice& a, const ShiFengChoice& b) {
            if (config.residual_demand_rank && g_logical_search_params.depth_goal > 0) {
                const auto a_demand_key = std::make_tuple(
                    a.residual_endpoint_goal_excess_after,
                    a.residual_endpoint_max_after,
                    a.residual_endpoint_energy_after
                );
                const auto b_demand_key = std::make_tuple(
                    b.residual_endpoint_goal_excess_after,
                    b.residual_endpoint_max_after,
                    b.residual_endpoint_energy_after
                );
                if (a_demand_key != b_demand_key) return a_demand_key < b_demand_key;
            }
            if (config.endpoint_balanced_rank && g_logical_search_params.depth_goal > 0) {
                const auto a_endpoint_key = std::make_tuple(
                    a.endpoint_goal_excess_after,
                    a.endpoint_max_after,
                    a.endpoint_energy_after
                );
                const auto b_endpoint_key = std::make_tuple(
                    b.endpoint_goal_excess_after,
                    b.endpoint_max_after,
                    b.endpoint_energy_after
                );
                if (a_endpoint_key != b_endpoint_key) return a_endpoint_key < b_endpoint_key;
            }
            if (std::abs(a.rank_cost_after - b.rank_cost_after) > 1.0e-9) {
                return a.rank_cost_after < b.rank_cost_after;
            }
            if (std::abs(a.cost_after - b.cost_after) > 1.0e-9) return a.cost_after < b.cost_after;
            if (a.layer_depth_delta != b.layer_depth_delta) return a.layer_depth_delta < b.layer_depth_delta;
            if (std::abs(a.priority - b.priority) > 1.0e-9) return a.priority > b.priority;
            if (a.control != b.control) return a.control < b.control;
            if (a.target != b.target) return a.target < b.target;
            return static_cast<int>(a.side) < static_cast<int>(b.side);
        };

        std::size_t selected = 0;
        if (g_logical_search_params.choice_window <= 1) {
            const auto best_it = std::min_element(
                choices.begin(),
                choices.end(),
                choice_rank_less
            );
            std::vector<std::size_t> best_indices;
            for (std::size_t idx = 0; idx < choices.size(); ++idx) {
                if (!choice_rank_less(choices[idx], *best_it)
                    && !choice_rank_less(*best_it, choices[idx])) {
                    best_indices.push_back(idx);
                }
            }

            selected = best_indices.front();
            if (tie_seed != std::numeric_limits<std::uint64_t>::max() && best_indices.size() > 1) {
                selected = best_indices[static_cast<std::size_t>(tie_rng() % best_indices.size())];
            } else {
                for (std::size_t idx : best_indices) {
                    const auto& a = choices[idx];
                    const auto& b = choices[selected];
                    if (std::make_tuple(a.layer_depth_delta, -a.priority, a.control, a.target, static_cast<int>(a.side))
                        < std::make_tuple(b.layer_depth_delta, -b.priority, b.control, b.target, static_cast<int>(b.side))) {
                        selected = idx;
                    }
                }
            }
        } else {
            std::vector<std::size_t> ranked;
            ranked.reserve(choices.size());
            for (std::size_t idx = 0; idx < choices.size(); ++idx) ranked.push_back(idx);
            std::sort(ranked.begin(), ranked.end(), [&](std::size_t ai, std::size_t bi) {
                return choice_rank_less(choices[ai], choices[bi]);
            });
            const std::size_t cap = std::min(g_logical_search_params.choice_window, ranked.size());
            std::vector<std::size_t> pool(ranked.begin(), ranked.begin() + static_cast<std::ptrdiff_t>(cap));
            selected = pool.front();
            if (tie_seed != std::numeric_limits<std::uint64_t>::max() && pool.size() > 1) {
                selected = pool[static_cast<std::size_t>(tie_rng() % pool.size())];
            }
        }

        ShiFengChoice choice = std::move(choices[selected]);
        AlgebraicChoice algebraic{choice.side, choice.control, choice.target, 0.0, choice.priority};
        apply_algebraic_choice_with_inverse(residual, residual_inverse, algebraic);
        if (choice.side == OpSide::Row) {
            append_to_search_node(current_row_layer, layout, "logical", std::move(choice.routed));
        } else {
            append_to_search_node(current_column_layer, layout, "logical", std::move(choice.routed));
        }
        endpoint_energy =
            choice.endpoint_energy_after;
        ++endpoint_load[static_cast<std::size_t>(choice.control)];
        ++endpoint_load[static_cast<std::size_t>(choice.target)];
        endpoint_max_load = choice.endpoint_max_after;
    }

    flush_search_layers();
    return reconstruct_logical_search_result(target_rows, layout, residual, row_layers, column_layers);
}

static Result synthesize_shi_feng_logical(
    const Matrix& target_rows,
    const Layout& layout,
    int max_layers,
    const std::vector<int>* required_output_permutation
) {
    bool have_best = false;
    Result best;
    std::vector<Matrix> search_targets;
    std::unordered_set<std::string> queued_targets;
    auto queue_search_target = [&](const Matrix& rows) {
        const std::string key = matrix_key(rows);
        if (queued_targets.insert(key).second) search_targets.push_back(rows);
    };
    queue_search_target(target_rows);
    if (required_output_permutation != nullptr) {
        queue_search_target(permute_rows(target_rows, *required_output_permutation));
    }

    auto keep_best = [&](Result result) {
        result = retarget_logical_result(std::move(result), target_rows);
        queue_search_target(permute_rows(target_rows, result.output_permutation));
        if (required_output_permutation != nullptr
            && result.output_permutation != *required_output_permutation) {
            return;
        }
        auto logical_best_key = [&](const Result& item) {
            return logical_result_priority_key(item);
        };
        if (!have_best || logical_best_key(result) < logical_best_key(best)) {
            best = std::move(result);
            have_best = true;
        }
    };

    const std::size_t max_search_targets = required_output_permutation == nullptr
        ? std::max<std::size_t>(1, g_logical_search_params.max_search_targets)
        : std::max<std::size_t>(2, g_logical_search_params.max_search_targets);
    const auto configs = shi_feng_logical_configs();
    const std::uint64_t seeds_per_target = effective_seeds_per_target(configs.size());

    for (std::size_t target_idx = 0; target_idx < search_targets.size() && target_idx < max_search_targets; ++target_idx) {
        const Matrix search_target = search_targets[target_idx];
        struct LogicalTask {
            std::size_t order = 0;
            std::size_t config_idx = 0;
            std::uint64_t seed = 0;
        };
        struct LogicalTaskResult {
            std::size_t order = 0;
            Result result;
        };

        std::vector<LogicalTask> tasks;
        tasks.reserve(configs.size() * static_cast<std::size_t>(seeds_per_target));
        for (std::size_t config_idx = 0; config_idx < configs.size(); ++config_idx) {
            for (std::uint64_t seed = 0; seed < seeds_per_target; ++seed) {
                tasks.push_back({tasks.size(), config_idx, seed});
            }
        }

        std::vector<LogicalTaskResult> results;
        std::mutex results_mutex;
        std::atomic<std::size_t> next_task{0};
        auto run_worker = [&]() {
            while (true) {
                const std::size_t task_idx = next_task.fetch_add(1);
                if (task_idx >= tasks.size()) break;
                const LogicalTask task = tasks[task_idx];
                const bool deterministic_first =
                    g_logical_search_params.seed_offset == 0
                    && task.seed == 0
                    && target_idx == 0;
                const std::uint64_t tie_seed = deterministic_first
                    ? std::numeric_limits<std::uint64_t>::max()
                    : (0x9e3779b97f4a7c15ULL
                       ^ (g_logical_search_params.seed_offset * 0xd6e8feb86659fd93ULL)
                       ^ (task.seed * 0xbf58476d1ce4e5b9ULL)
                       ^ (target_idx * 0x94d049bb133111ebULL));
                try {
                    Result result = synthesize_shi_feng_logical_once(
                        search_target,
                        layout,
                        configs[task.config_idx],
                        tie_seed,
                        max_layers
                    );
                    std::lock_guard<std::mutex> lock(results_mutex);
                    results.push_back({task.order, std::move(result)});
                } catch (const std::exception&) {
                }
            }
        };

        run_worker();

        std::sort(results.begin(), results.end(), [](const LogicalTaskResult& a, const LogicalTaskResult& b) {
            return a.order < b.order;
        });
        for (auto& item : results) {
            keep_best(std::move(item.result));
            if (have_best && logical_result_meets_goals(best)) return best;
        }
    }
    if (!have_best) throw std::runtime_error("no logical Shi-Feng result");
    return best;
}

static Result synthesize_layer_matrix_surface_pool(
    const Matrix& target_rows,
    const std::string& mode,
    const Layout& layout,
    int max_layers,
    const std::vector<int>* required_output_permutation = nullptr
) {
    if (mode == "logical") throw std::runtime_error("surface pool needs a surface mode");

    struct SurfaceRankScope {
        const Layout* previous_layout = nullptr;
        std::string previous_mode;
        SurfaceRankScope(const Layout& layout_ref, const std::string& mode_ref)
            : previous_layout(g_layer_matrix_surface_layout),
              previous_mode(g_layer_matrix_surface_mode) {
            g_layer_matrix_surface_layout = &layout_ref;
            g_layer_matrix_surface_mode = mode_ref;
        }
        ~SurfaceRankScope() {
            g_layer_matrix_surface_layout = previous_layout;
            g_layer_matrix_surface_mode = previous_mode;
        }
    } surface_rank_scope(layout, mode);

    std::vector<ShiFengSurfaceConfig> configs;
    for (const auto& config : shi_feng_logical_configs()) {
        if (config.layer_matrix_rank) configs.push_back(config);
    }
    if (configs.empty()) throw std::runtime_error("no layer-matrix configs available");

    bool have_best = false;
    Result best;
    std::vector<Matrix> search_targets;
    std::unordered_set<std::string> queued_targets;
    auto queue_search_target = [&](const Matrix& rows) {
        const std::string key = matrix_key(rows);
        if (queued_targets.insert(key).second) search_targets.push_back(rows);
    };
    queue_search_target(target_rows);
    if (required_output_permutation != nullptr) {
        queue_search_target(permute_rows(target_rows, *required_output_permutation));
    }

    const std::size_t max_search_targets = required_output_permutation == nullptr
        ? std::max<std::size_t>(1, g_logical_search_params.max_search_targets)
        : std::max<std::size_t>(2, g_logical_search_params.max_search_targets);
    const std::uint64_t seeds_per_target = effective_seeds_per_target(configs.size());

    for (std::size_t target_idx = 0; target_idx < search_targets.size() && target_idx < max_search_targets; ++target_idx) {
        const Matrix search_target = search_targets[target_idx];
        struct SurfacePoolTask {
            std::size_t order = 0;
            std::size_t config_idx = 0;
            std::uint64_t seed = 0;
        };
        struct SurfacePoolResult {
            std::size_t order = 0;
            Result result;
        };

        std::vector<SurfacePoolTask> tasks;
        tasks.reserve(configs.size() * static_cast<std::size_t>(seeds_per_target));
        for (std::size_t config_idx = 0; config_idx < configs.size(); ++config_idx) {
            for (std::uint64_t seed = 0; seed < seeds_per_target; ++seed) {
                tasks.push_back({tasks.size(), config_idx, seed});
            }
        }

        std::vector<SurfacePoolResult> results;
        std::mutex results_mutex;
        std::atomic<std::size_t> next_task{0};
        auto run_worker = [&]() {
            while (true) {
                const std::size_t task_idx = next_task.fetch_add(1);
                if (task_idx >= tasks.size()) break;
                const SurfacePoolTask task = tasks[task_idx];
                const bool deterministic_first =
                    g_logical_search_params.seed_offset == 0
                    && task.seed == 0
                    && target_idx == 0;
                const std::uint64_t tie_seed = deterministic_first
                    ? std::numeric_limits<std::uint64_t>::max()
                    : (0x9e3779b97f4a7c15ULL
                       ^ (g_logical_search_params.seed_offset * 0xd6e8feb86659fd93ULL)
                       ^ (task.seed * 0xbf58476d1ce4e5b9ULL)
                       ^ (target_idx * 0x94d049bb133111ebULL));
                try {
                    Result logical = synthesize_shi_feng_logical_once(
                        search_target,
                        layout,
                        configs[task.config_idx],
                        tie_seed,
                        max_layers
                    );
                    logical = retarget_logical_result(std::move(logical), target_rows);
                    Result surface = schedule_logical_result_on_surface(
                            logical,
                            mode,
                            layout,
                            logical.output_permutation
                        );
                    export_verified_anytime(surface, "layer_matrix_complete_surface_configuration");
                    std::lock_guard<std::mutex> lock(results_mutex);
                    results.push_back({task.order, std::move(surface)});
                } catch (const std::exception&) {
                }
            }
        };

        run_worker();

        std::sort(results.begin(), results.end(), [](const SurfacePoolResult& a, const SurfacePoolResult& b) {
            return a.order < b.order;
        });
        for (auto& item : results) {
            queue_search_target(permute_rows(target_rows, item.result.output_permutation));
            if (required_output_permutation != nullptr
                && item.result.output_permutation != *required_output_permutation) {
                continue;
            }
            const auto key = std::make_tuple(item.result.surface_depth(), item.result.layer_count(), item.result.cnot_count());
            const auto best_key = std::make_tuple(best.surface_depth(), best.layer_count(), best.cnot_count());
            if (!have_best || key < best_key) {
                best = std::move(item.result);
                have_best = true;
            }
        }
    }

    if (!have_best) throw std::runtime_error("no layer-matrix surface-pool result");
    return best;
}

static Result synthesize_routed_greedy(
    const Matrix& target_rows,
    const std::string& mode,
    const Layout& layout,
    const std::vector<int>& output_permutation,
    int max_layers = 10000
) {
    const int n = check_square_rows(target_rows);
    const Matrix desired = permute_rows(target_rows, output_permutation);
    Matrix residual = invert_matrix(desired);
    Matrix current = identity_matrix(n);
    std::vector<Layer> layers;
    std::unordered_set<std::string> seen;
    const Matrix identity = identity_matrix(n);
    const int heuristic_layer_budget = std::min(max_layers, 4 * n * n);

    for (int iter = 0; iter < heuristic_layer_budget; ++iter) {
        if (residual == identity) break;
        const std::string key = matrix_key(residual);
        if (seen.find(key) != seen.end()) break;
        seen.insert(key);

        std::vector<RoutedCNOT> routed = greedy_pack_layer(layout, mode, greedy_candidate_list(residual, layout));
        if (routed.empty()) break;
        apply_layer(residual, routed);
        apply_layer(current, routed);
        layers.push_back({mode, std::move(routed), 0});
        layers.back().logical_depth = layer_depth(mode, layout, layers.back().ops);
    }

    if (residual != identity) {
        const auto fallback_ops = gaussian_reduction_ops(residual);
        std::vector<Layer> fallback_layers = schedule_ops(fallback_ops, layout, mode);
        for (const auto& layer : fallback_layers) {
            apply_layer(residual, layer.ops);
            apply_layer(current, layer.ops);
            layers.push_back(layer);
        }
    }
    if (residual != identity) throw std::runtime_error("greedy synthesis failed");

    std::vector<RowOp> circuit_ops = recombine_cnot_sequence(flatten_layers(layers), layout, mode);
    layers = schedule_ops(circuit_ops, layout, mode);
    current = identity_matrix(n);
    for (const auto& layer : layers) apply_layer(current, layer.ops);

    Result result;
    result.n = n;
    result.mode = mode;
    result.layout = layout;
    result.target_rows = target_rows;
    result.output_permutation = output_permutation;
    result.layers = std::move(layers);
    result.final_rows = current;
    if (!verify_geometry_and_stats(result)) throw std::runtime_error("greedy synthesis verification failed");
    return result;
}

static Result synthesize_once(
    const Matrix& target_rows,
    const std::string& mode,
    const Layout& layout,
    const std::vector<int>& output_permutation,
    const SearchWeights& weights,
    int max_layers = 10000
) {
    const int n = check_square_rows(target_rows);
    if (layout.n != n) throw std::runtime_error("layout size mismatch");
    if (static_cast<int>(output_permutation.size()) != n) throw std::runtime_error("permutation size mismatch");

    const Matrix desired = permute_rows(target_rows, output_permutation);
    const Matrix target_inv = invert_matrix(desired);
    Matrix residual = target_inv;
    std::unordered_set<std::string> seen;
    std::vector<RowOp> row_ops;
    std::vector<RowOp> column_ops;
    SurfaceScheduleState row_schedule(layout, mode);
    const int heuristic_step_budget = std::min(max_layers, std::max(64, 12 * n * n));

    for (int iter = 0; iter < heuristic_step_budget; ++iter) {
        if (residual == identity_matrix(n)) break;
        const std::string key = matrix_key(residual);
        if (seen.find(key) != seen.end()) break;
        seen.insert(key);

        auto choices = algebraic_candidates(
            residual,
            layout,
            weights,
            true,
            &row_schedule,
            nullptr
        );
        if (choices.empty()) break;
        const AlgebraicChoice choice = choices.front();
        apply_algebraic_choice(residual, choice);
        if (choice.side == OpSide::Row) {
            row_ops.emplace_back(choice.control, choice.target);
            row_schedule.append(choice.control, choice.target);
        } else {
            column_ops.emplace_back(choice.control, choice.target);
        }
    }

    if (residual != identity_matrix(n)) {
        const auto fallback = gaussian_reduction_ops(residual);
        for (const RowOp& op : fallback) {
            row_add(residual, op.first, op.second);
            row_ops.push_back(op);
            row_schedule.append(op.first, op.second);
        }
    }
    if (residual != identity_matrix(n)) throw std::runtime_error("matrix reduction failed");

    std::vector<RowOp> circuit_ops = row_ops;
    for (auto it = column_ops.rbegin(); it != column_ops.rend(); ++it) {
        circuit_ops.push_back(*it);
    }
    circuit_ops = optimize_cnot_sequence(circuit_ops);
    circuit_ops = recombine_cnot_sequence(circuit_ops, layout, mode);
    std::vector<Layer> layers = schedule_ops(circuit_ops, layout, mode);

    Matrix current = identity_matrix(n);
    for (const auto& layer : layers) apply_layer(current, layer.ops);

    Result result;
    result.n = n;
    result.mode = mode;
    result.layout = layout;
    result.target_rows = target_rows;
    result.output_permutation = output_permutation;
    result.layers = std::move(layers);
    result.final_rows = current;
    if (!verify_geometry_and_stats(result)) throw std::runtime_error("synthesis verification failed");
    return result;
}

static std::vector<SearchWeights> search_presets() {
    return {
        {12.0, 7.0, 0.70, 0.50, 0.03, 0.00, true},
        {8.0, 12.0, 0.50, 0.70, 0.03, 0.00, true},
        {10.0, 10.0, 1.20, 1.20, 0.04, 0.00, true},
        {16.0, 4.0, 0.35, 0.20, 0.05, 0.00, false},
        {12.0, 7.0, 0.70, 0.50, 0.03, 0.35, true},
        {8.0, 12.0, 0.50, 0.70, 0.03, 0.40, true},
        {10.0, 10.0, 1.20, 1.20, 0.04, 0.30, true},
        {16.0, 4.0, 0.35, 0.20, 0.05, 0.25, false}
    };
}
