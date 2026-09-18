#pragma once
#include "circuit.hpp"
#include "policy.hpp"

namespace linear_surface {

inline RoutingPolicy routing_policy;

inline int layer_depth(const std::string& mode, const Layout&, const std::vector<RoutedCNOT>&) {
    if (mode == "logical") return 1;
    if (mode == "vdp") return 2;
    throw std::runtime_error("unsupported layer mode");
}

inline bool try_route_into_layer(
    const Layout& layout,
    const std::string& mode,
    const Layer& layer,
    int control,
    int target,
    RoutedCNOT& routed
) {
    std::vector<char> used_qubits(layout.n, 0);
    std::vector<char> used_vertices(layout.grid_rows * layout.grid_cols, 0);
    std::unordered_set<std::uint64_t> used_edges;
    for (const auto& op : layer.ops) {
        used_qubits[op.control] = 1;
        used_qubits[op.target] = 1;
        for (auto e : layout.path_edges(op.path)) used_edges.insert(e);
        if (mode == "vdp") {
            for (int v : op.path) used_vertices[v] = 1;
        }
    }
    if (used_qubits[control] || used_qubits[target]) return false;
    if (mode == "logical") {
        routed = {control, target, {}};
        return true;
    }
    const std::vector<char>* banned_vertices = (mode == "vdp") ? &used_vertices : nullptr;
    std::vector<int> path = layout.find_operator_path(control, target, used_edges, banned_vertices);

    if (path.empty()) return false;
    routed = {control, target, std::move(path)};
    return true;
}

inline bool paths_conflict(
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

inline bool pack_layer_ops_with_k_paths(
    const Layout& layout,
    const std::string& mode,
    const std::vector<RowOp>& ops,
    std::vector<RoutedCNOT>& routed_ops
) {
    if (mode == "logical") {
        routed_ops.clear();
        routed_ops.reserve(ops.size());
        for (const auto& op : ops) routed_ops.push_back({op.first, op.second, {}});
        return true;
    }
    if (static_cast<int>(ops.size()) > routing_policy.repack_limit) return false;

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

    static std::unordered_map<std::string, std::vector<std::vector<int>>> path_cache;
    for (std::size_t i = 0; i < ops.size(); ++i) {
        std::ostringstream key_stream;
        key_stream << layout.n << '|' << layout.data_rows << 'x' << layout.data_cols
                   << '|' << ops[i].first << '>' << ops[i].second
                   << "|k=" << routing_policy.k_paths
                   << "|extra=" << routing_policy.path_extra;
        const std::string key = key_stream.str();
        {

            const auto it = path_cache.find(key);
            if (it != path_cache.end()) candidates[i] = it->second;
        }
        if (candidates[i].empty()) {
            candidates[i] = layout.find_operator_paths(
                ops[i].first,
                ops[i].second,
                no_edges,
                nullptr,
                routing_policy.k_paths,
                routing_policy.path_extra
            );

            path_cache.emplace(key, candidates[i]);
        }
        if (candidates[i].empty()) return false;
    }

    std::vector<std::vector<int>> chosen(ops.size());
    std::unordered_set<std::uint64_t> used_edges;
    std::vector<char> used_vertices(static_cast<std::size_t>(layout.grid_rows * layout.grid_cols), 0);
    bool found = false;
    int visited_states = 0;
    const int state_limit = std::max(20000, routing_policy.k_paths * 50000);
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
        constexpr int conditional_state_limit = 50000;
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
        routed_ops.push_back({ops[i].first, ops[i].second, std::move(chosen[i])});
    }
    return true;
}

inline bool repack_layer_with_op(
    const Layout& layout,
    const std::string& mode,
    const Layer& layer,
    int control,
    int target,
    Layer& replacement
) {
    if (!routing_policy.enabled()) return false;
    std::vector<RowOp> ops;
    ops.reserve(layer.ops.size() + 1);
    for (const auto& op : layer.ops) ops.emplace_back(op.control, op.target);
    ops.emplace_back(control, target);

    std::vector<RoutedCNOT> routed_ops;
    if (!pack_layer_ops_with_k_paths(layout, mode, ops, routed_ops)) return false;
    replacement.mode = mode;
    replacement.ops = std::move(routed_ops);
    replacement.cycles = layer_depth(mode, layout, replacement.ops);
    return true;
}

struct LogicalSequenceBounds {
    int capacity_bound = 0;
    int endpoint_bound = 0;
    int critical_bound = 0;
    int lower_bound = 0;
};

inline void build_logical_precedence(
    const std::vector<RowOp>& ops,
    std::vector<std::vector<int>>& succ,
    std::vector<int>& pred_count
) {
    const int m = static_cast<int>(ops.size());
    succ.assign(static_cast<std::size_t>(m), {});
    pred_count.assign(static_cast<std::size_t>(m), 0);
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            if (cnot_commute(ops[static_cast<std::size_t>(i)], ops[static_cast<std::size_t>(j)])) continue;
            succ[static_cast<std::size_t>(i)].push_back(j);
            ++pred_count[static_cast<std::size_t>(j)];
        }
    }
}

inline std::vector<int> logical_critical_lengths(const std::vector<std::vector<int>>& succ) {
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

inline LogicalSequenceBounds logical_sequence_bounds(const std::vector<RowOp>& ops, int n) {
    const auto fast = mixcolumn::logical_sequence_bounds_linear(ops, n);
    return {fast.capacity_bound, fast.endpoint_bound, fast.critical_bound, fast.lower_bound};
}

inline std::vector<Layer> schedule_logical_commuting_ops_once(
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

        layer.cycles = 1;
        layers.push_back(std::move(layer));
    }
    return layers;
}

inline bool better_logical_schedule(
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

inline std::vector<Layer> schedule_logical_commuting_ops(
    const std::vector<RowOp>& ops,
    const Layout& layout
) {
    if (ops.empty()) return {};
    std::vector<Layer> best = schedule_logical_commuting_ops_once(ops, layout, 0);
    constexpr int order_count = 28;
    const LogicalSequenceBounds bounds = logical_sequence_bounds(ops, layout.n);
    for (int order_kind = 1; order_kind < order_count; ++order_kind) {
        std::vector<Layer> candidate = schedule_logical_commuting_ops_once(ops, layout, order_kind);
        if (better_logical_schedule(candidate, best)) best = std::move(candidate);
        if (static_cast<int>(best.size()) == bounds.lower_bound) break;
    }
    return best;
}

inline std::vector<Layer> schedule_ops(
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
            const int delta = std::max(0, new_depth - layers[li].cycles);
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
                layers[best_li].cycles = layer_depth(mode, layout, layers[best_li].ops);
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
            layer.cycles = layer_depth(mode, layout, layer.ops);
            layers.push_back(std::move(layer));
            const int li = static_cast<int>(layers.size()) - 1;
            last_touch[control] = li;
            last_touch[target] = li;
        }
    }
    for (auto& layer : layers) {
        layer.cycles = layer_depth(mode, layout, layer.ops);
    }
    return layers;
}

inline std::vector<RowOp> topological_reorder_ops(
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

inline std::vector<Layer> schedule_surface_commuting_ops(
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

        layer.cycles = layer_depth(mode, layout, layer.ops);
        layers.push_back(std::move(layer));
    }
    return layers;
}

} // namespace linear_surface
