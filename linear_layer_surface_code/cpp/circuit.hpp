#pragma once
#include "algebra.hpp"

namespace linear_surface {
struct Layout {
    int n = 0;
    int data_rows = 0;
    int data_cols = 0;
    int grid_rows = 0;
    int grid_cols = 0;
    std::vector<int> data_pos;
    std::vector<char> occupied_data;

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

    int vertex_id(int r, int c) const { return r * grid_cols + c; }
    int row(int id) const { return id / grid_cols; }
    int col(int id) const { return id % grid_cols; }

    // Return the same deterministic neighbor order as the historical vector
    // helper, but write into caller-owned storage.  Path routing invokes this
    // routine millions of times; avoiding a temporary vector and allocator
    // traffic materially improves throughput without changing BFS/DFS order.
    int neighbor_ids(int id, std::array<int, 4>& out) const noexcept {
        const int r = row(id);
        const int c = col(id);
        int count = 0;
        if (r > 0) out[static_cast<std::size_t>(count++)] = vertex_id(r - 1, c);
        if (r + 1 < grid_rows) out[static_cast<std::size_t>(count++)] = vertex_id(r + 1, c);
        if (c > 0) out[static_cast<std::size_t>(count++)] = vertex_id(r, c - 1);
        if (c + 1 < grid_cols) out[static_cast<std::size_t>(count++)] = vertex_id(r, c + 1);
        return count;
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
            std::array<int, 4> neighbors{};
            const int neighbor_count = neighbor_ids(cur, neighbors);
            for (int neighbor_index = 0; neighbor_index < neighbor_count; ++neighbor_index) {
                const int nb = neighbors[static_cast<std::size_t>(neighbor_index)];
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

            std::array<int, 4> nbs{};
            const int neighbor_count = neighbor_ids(cur, nbs);
            std::sort(nbs.begin(), nbs.begin() + neighbor_count, [&](int a, int b) {
                const int ha = heuristic(a);
                const int hb = heuristic(b);
                if (ha != hb) return ha < hb;
                return std::make_tuple(row(a), col(a)) < std::make_tuple(row(b), col(b));
            });

            for (int neighbor_index = 0; neighbor_index < neighbor_count; ++neighbor_index) {
                const int nb = nbs[static_cast<std::size_t>(neighbor_index)];
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

struct RoutedCNOT {
    int control = 0;
    int target = 0;
    std::vector<int> path;
};

struct Layer {
    std::string mode;
    std::vector<RoutedCNOT> ops;
    int cycles = 0;
};

struct Circuit {
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

    int cycles() const {
        int total = 0;
        for (const auto& layer : layers) total += layer.cycles;
        return total;
    }

    bool verify() const {
        return final_rows == permute_rows(target_rows, output_permutation);
    }
};

inline void apply_layer(Matrix& rows, const std::vector<RoutedCNOT>& ops) {
    for (const auto& op : ops) row_add(rows, op.control, op.target);
}

inline void validate_routed_path(
    const Circuit& result,
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

inline int checked_layer_depth(const Circuit& result, const Layer& layer) {
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
    throw std::runtime_error("unknown layer mode in result");
}

inline bool verify_geometry_and_stats(const Circuit& result) {
    try {
        if (result.n <= 0 || result.layout.n != result.n) return false;
        check_square_rows(result.target_rows);
        check_square_rows(result.final_rows);
        validate_permutation_vector(result.output_permutation, result.n);

        Matrix rows = identity_matrix(result.n);
        int cnots = 0;
        int cycles = 0;
        for (const auto& layer : result.layers) {
            const int actual_depth = checked_layer_depth(result, layer);
            if (layer.cycles != actual_depth) return false;
            apply_layer(rows, layer.ops);
            cnots += static_cast<int>(layer.ops.size());
            cycles += actual_depth;
        }
        if (cnots != result.cnot_count()) return false;
        if (cycles != result.cycles()) return false;
        return rows == result.final_rows && result.verify();
    } catch (const std::exception&) {
        return false;
    }
}

inline std::vector<RowOp> flatten_layers(const std::vector<Layer>& layers) {
    std::vector<RowOp> out;
    for (const auto& layer : layers) {
        for (const auto& op : layer.ops) out.emplace_back(op.control, op.target);
    }
    return out;
}

} // namespace linear_surface
