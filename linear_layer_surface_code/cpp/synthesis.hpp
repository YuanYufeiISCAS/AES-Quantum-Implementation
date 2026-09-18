#pragma once
#include "routing.hpp"

namespace linear_surface {
struct AlgebraicChoice {
    OpSide side = OpSide::Row;
    int control = 0;
    int target = 0;
    double delta = 0.0;
    double priority = 0.0;
};

enum class ShiFengCostKind { Square };

struct MatrixSearchConfig {
    static constexpr ShiFengCostKind cost_kind = ShiFengCostKind::Square;
    static constexpr bool endpoint_balanced_rank = false;
    bool allow_columns = true;
};

inline const Layout* synthesis_layout = nullptr;

inline double h_square_cost(const Matrix& rows) {
    check_square_rows(rows);
    double cost = 0.0;
    for (Row row : rows) {
        const int w = popcount64(row);
        cost += static_cast<double>(w * w);
    }
    return cost;
}

inline double shi_base_cost(const Matrix& rows, ShiFengCostKind) { return h_square_cost(rows); }

inline double shi_row_cost_with_inverse(const Matrix& rows, const Matrix& inverse, ShiFengCostKind kind) {
    return shi_base_cost(rows, kind) + shi_base_cost(transpose_matrix(inverse), kind);
}

inline double shi_column_cost_with_inverse(const Matrix& rows, const Matrix& inverse, ShiFengCostKind kind) {
    return shi_base_cost(transpose_matrix(rows), kind) + shi_base_cost(inverse, kind);
}

inline double shi_current_cost_with_inverse(const Matrix& rows, const Matrix& inverse, ShiFengCostKind kind) {
    return std::max(
        shi_row_cost_with_inverse(rows, inverse, kind),
        shi_column_cost_with_inverse(rows, inverse, kind)
    );
}

struct ShiFengAfterCosts {
    double side_after = 0.0;
    double global_after = 0.0;
};

inline ShiFengAfterCosts shi_side_and_global_costs_with_inverse(
    const Matrix& rows, const Matrix& inverse, OpSide side, ShiFengCostKind kind
) {
    const double row = shi_row_cost_with_inverse(rows, inverse, kind);
    const double column = shi_column_cost_with_inverse(rows, inverse, kind);
    return {side == OpSide::Row ? row : column, std::max(row, column)};
}

inline void apply_algebraic_choice_with_inverse(
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

inline Circuit make_result_from_logical_ops(
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

    Circuit result;
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

inline Circuit make_result_from_logical_layers(
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
        layer.cycles = 1;
        layers.push_back(std::move(layer));
    }

    std::vector<int> output_permutation;
    if (!infer_output_permutation(target_rows, final_rows, output_permutation)) {
        throw std::runtime_error("logical layers are not an output permutation of target");
    }

    Circuit result;
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

inline Circuit schedule_logical_result_on_surface_basic(
    const Circuit& logical_result,
    const std::string& mode,
    const Layout& layout,
    const std::string& output_mode = ""
);

inline Circuit schedule_logical_result_on_surface_basic(
    const Circuit& logical_result,
    const std::string& mode,
    const Layout& layout,
    const std::string& output_mode
) {
    if (mode != "vdp") {
        throw std::runtime_error("surface routing requires VDP");
    }

    std::vector<Layer> surface_layers;
    Matrix current = identity_matrix(logical_result.n);
    for (const auto& logical_layer : logical_result.layers) {
        std::vector<RowOp> ops;
        ops.reserve(logical_layer.ops.size());
        for (const auto& op : logical_layer.ops) ops.emplace_back(op.control, op.target);

        std::vector<Layer> scheduled = schedule_ops(ops, layout, mode);
        for (auto& layer : scheduled) {
            apply_layer(current, layer.ops);
            surface_layers.push_back(std::move(layer));
        }
    }

    Circuit result;
    result.n = logical_result.n;
    result.mode = output_mode.empty() ? "basic-" + mode : output_mode;
    result.layout = layout;
    result.target_rows = logical_result.target_rows;
    result.output_permutation = logical_result.output_permutation;
    result.layers = std::move(surface_layers);
    result.final_rows = current;
    if (!verify_geometry_and_stats(result)) throw std::runtime_error("basic logical-to-surface verification failed");
    return result;
}

inline std::tuple<int, int, int, int, int, int> logical_result_priority_key(const Circuit& item) {
    const auto bounds = logical_sequence_bounds(flatten_layers(item.layers), item.n);
    if (item.layer_count() <= 10) {
        return {0, item.layer_count(), item.cnot_count(), bounds.critical_bound, bounds.endpoint_bound, 0};
    }
    return {1, bounds.lower_bound, bounds.critical_bound, bounds.endpoint_bound,
            item.layer_count(), item.cnot_count()};
}

// Bounded candidate archive, alternating four priorities: candidate depth,
// endpoint pressure, geometric distance, and CNOT count.
inline std::string surface_logical_candidate_key(const Circuit& result) {
    std::string key;
    key.reserve(2U + result.output_permutation.size()
        + static_cast<std::size_t>(2 * result.cnot_count()));
    key.push_back(static_cast<char>(result.n));
    for (int value : result.output_permutation) key.push_back(static_cast<char>(value));
    key.push_back(static_cast<char>(255));
    for (const Layer& layer : result.layers) {
        for (const RoutedCNOT& op : layer.ops) {
            key.push_back(static_cast<char>(op.control));
            key.push_back(static_cast<char>(op.target));
        }
    }
    return key;
}

class SurfaceLogicalCandidateCollector {
    using CheapKey = std::array<long long, 6>;
    struct Entry {
        Circuit result;
        std::string word_key;
        std::array<CheapKey, 4> ranks;
    };

    std::size_t capacity_;
    std::vector<Entry> entries_;

    static Entry make_entry(const Circuit& result, std::string word_key) {
        const auto ops = flatten_layers(result.layers);
        const LogicalSequenceBounds bounds = result.n > 1
            ? logical_sequence_bounds(ops, result.n) : LogicalSequenceBounds{};
        const auto logical_key = result.n > 1 ? logical_result_priority_key(result)
            : std::make_tuple(0, result.layer_count(), result.layer_count(),
                              result.layer_count(), result.layer_count(), result.cnot_count());
        long long distance = 0;
        for (const RowOp& op : ops) {
            distance += result.layout.manhattan_between_qubits(op.first, op.second);
        }
        const long long cnots = static_cast<long long>(ops.size());
        const long long layers = result.layer_count();
        Entry entry{result, std::move(word_key), {}};
        entry.ranks[0] = {
            std::get<0>(logical_key), std::get<1>(logical_key),
            std::get<2>(logical_key), std::get<3>(logical_key),
            std::get<4>(logical_key), std::get<5>(logical_key)
        };
        entry.ranks[1] = {
            bounds.endpoint_bound, bounds.critical_bound, distance, cnots, layers, 0
        };
        entry.ranks[2] = {
            distance, bounds.endpoint_bound, bounds.critical_bound, cnots, layers, 0
        };
        entry.ranks[3] = {
            cnots, bounds.critical_bound, bounds.endpoint_bound, layers, distance, 0
        };
        return entry;
    }

    // Alternate logical depth, endpoint pressure, geometric distance and gate
    // count rather than reducing all candidates to one logical priority key.
    // The archive has at most capacity+1 entries even while admitting a new
    // candidate.  Evicted words need not stay in an unbounded deduplication set.
    std::vector<std::size_t> diverse_indices(std::size_t limit) const {
        std::vector<std::size_t> selected;
        selected.reserve(std::min(limit, entries_.size()));
        std::vector<char> used(entries_.size(), 0);
        for (std::size_t round = 0; selected.size() < limit
             && selected.size() < entries_.size(); ++round) {
            const std::size_t profile = round % 4;
            std::size_t best = entries_.size();
            for (std::size_t i = 0; i < entries_.size(); ++i) {
                if (used[i]) continue;
                if (best == entries_.size()
                    || entries_[i].ranks[profile] < entries_[best].ranks[profile]
                    || (entries_[i].ranks[profile] == entries_[best].ranks[profile]
                        && entries_[i].word_key < entries_[best].word_key)) {
                    best = i;
                }
            }
            used[best] = 1;
            selected.push_back(best);
        }
        return selected;
    }

public:
    explicit SurfaceLogicalCandidateCollector(std::size_t capacity)
        : capacity_(std::max<std::size_t>(1, capacity)) {}

    std::size_t size() const { return entries_.size(); }

    bool consider(const Circuit& result) {
        if (result.mode != "logical"
            || std::any_of(result.layers.begin(), result.layers.end(), [](const Layer& layer) {
                return layer.mode != "logical";
            })
            || !verify_geometry_and_stats(result)) {
            return false;
        }
        std::string key = surface_logical_candidate_key(result);
        for (Entry& entry : entries_) {
            if (entry.word_key != key) continue;
            // Different grouping of one gate word does not occupy a second
            // slot.  The exact primary grouping is restored when it is routed.
            if (result.layer_count() < entry.result.layer_count()) {
                entry = make_entry(result, std::move(key));
            }
            return false;
        }
        entries_.push_back(make_entry(result, std::move(key)));
        if (entries_.size() > capacity_) {
            const auto selected = diverse_indices(capacity_);
            std::vector<Entry> retained;
            retained.reserve(capacity_);
            for (std::size_t i : selected) retained.push_back(std::move(entries_[i]));
            entries_ = std::move(retained);
        }
        return true;
    }

    std::vector<Circuit> take_with_best(Circuit primary) {
        const std::string best_word_key = surface_logical_candidate_key(primary);
        std::vector<Circuit> selected;
        selected.reserve(capacity_);
        selected.push_back(std::move(primary));
        const auto order = diverse_indices(entries_.size());
        for (std::size_t i : order) {
            if (selected.size() == capacity_) break;
            if (entries_[i].word_key == best_word_key) continue;
            selected.push_back(std::move(entries_[i].result));
        }
        entries_.clear();
        return selected;
    }
};

inline SurfaceLogicalCandidateCollector* g_surface_logical_candidate_collector = nullptr;

class SurfaceLogicalCandidateScope {
    SurfaceLogicalCandidateCollector* previous_;
public:
    explicit SurfaceLogicalCandidateScope(SurfaceLogicalCandidateCollector* collector)
        : previous_(g_surface_logical_candidate_collector) {
        g_surface_logical_candidate_collector = collector;
    }
    ~SurfaceLogicalCandidateScope() { g_surface_logical_candidate_collector = previous_; }
    SurfaceLogicalCandidateScope(const SurfaceLogicalCandidateScope&) = delete;
    SurfaceLogicalCandidateScope& operator=(const SurfaceLogicalCandidateScope&) = delete;
};

inline void collect_surface_logical_candidate(const Circuit& result) {
    if (g_surface_logical_candidate_collector != nullptr) {
        g_surface_logical_candidate_collector->consider(result);
    }
}

inline Circuit reconstruct_logical_search_result(
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
    Circuit best;

    auto reconstruction_key = [&](const Circuit& item) {
        return logical_result_priority_key(item);
    };
    for (const auto& layers : layer_candidates) {
        try {
            Circuit result = make_result_from_logical_layers(target_rows, layout, layers);
            collect_surface_logical_candidate(result);
            if (!have_best || reconstruction_key(result) < reconstruction_key(best)) {
                best = std::move(result);
                have_best = true;
            }
        } catch (const std::exception&) {
        }
    }
    for (auto& ops : candidates) {
        try {
            Circuit result = make_result_from_logical_ops(target_rows, layout, std::move(ops));
            collect_surface_logical_candidate(result);
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

inline int permutation_excess(const Matrix& rows) {
    (void)check_square_rows(rows);
    // For a non-negative integer weight w, |w-1| = w-1 when w>0 and 1
    // otherwise.  Summing this identity over rows and columns gives
    //
    //   E = 2 * ( total_popcount + zero_rows - nonzero_columns ).
    //
    // The OR of all rows exposes the nonzero columns directly.  This is
    // algebraically identical to the previous O(n^2) column scan, while
    // reducing every heuristic evaluation to one O(n) pass for n <= 64.
    int total_popcount = 0;
    int zero_rows = 0;
    Row nonzero_columns = 0;
    for (Row row : rows) {
        total_popcount += popcount64(row);
        if (row == 0) ++zero_rows;
        nonzero_columns |= row;
    }
    return 2 * (total_popcount + zero_rows - popcount64(nonzero_columns));
}

inline int identity_delta_rank(const Matrix& rows) {
    const int n = check_square_rows(rows);
    std::array<Row, 64> delta{};
    for (int i = 0; i < n; ++i) delta[static_cast<std::size_t>(i)] = rows[static_cast<std::size_t>(i)] ^ (Row(1) << i);
    return mixcolumn::gf2_rank_valid_rows(delta.data(), static_cast<std::size_t>(n));
}

inline bool logical_depth_one_ops_on_side(
    const Matrix& residual,
    OpSide side,
    std::vector<RowOp>& ops
) {
    if (side == OpSide::Row) return logical_depth_one_ops(residual, ops);
    Matrix transposed = transpose_matrix(residual);
    return logical_depth_one_ops(transposed, ops);
}

inline void apply_logical_layer_with_inverse(
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

inline std::string layer_matrix_candidate_key(const LayerMatrixCandidate& candidate) {
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

inline std::vector<LayerMatrixEdge> layer_matrix_edges(
    const Matrix& residual,
    const Matrix& residual_inverse,
    OpSide side,
    ShiFengCostKind kind,
    bool include_non_improving
) {
    const int n = check_square_rows(residual);
    const auto current_costs = shi_side_and_global_costs_with_inverse(residual, residual_inverse, side, kind);
    const double current_global = current_costs.global_after;
    const double current_side = current_costs.side_after;
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
            const auto after_costs = shi_side_and_global_costs_with_inverse(after, after_inverse, side, kind);
            const double side_after = after_costs.side_after;
            const double global_after = after_costs.global_after;
            const double side_delta = current_side - side_after;
            const double global_delta = current_global - global_after;
            if (!include_non_improving && side_delta <= 1.0e-9 && global_delta <= 1.0e-9) continue;
            const int perm_after = permutation_excess(after);
            const int rank_after = identity_delta_rank(after);
            if (side_delta <= 1.0e-9 && global_delta <= 1.0e-9
                && 0.0 <= 0.0) {
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

// Geometry-aware matching score: retain the better of near-first and
// far-first greedy insertion orders. This score only ranks beam candidates.
inline std::pair<int, int> layer_matrix_surface_pack_score_uncached(const std::vector<RowOp>& ops) {
    if (synthesis_layout == nullptr) return {0, 0};
    const Layout& layout = *synthesis_layout;

    auto score_order = [&](std::vector<RowOp> ordered) {
        Layer layer;
        int routed_count = 0;
        int path_score = 0;
        for (const RowOp& op : ordered) {
            RoutedCNOT routed;
            if (!try_route_into_layer(layout, "vdp", layer, op.first, op.second, routed)) {
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

inline std::unordered_map<std::string, std::pair<int, int>> surface_score_cache;
inline std::deque<std::string> surface_score_order;

inline std::pair<int, int> layer_matrix_surface_pack_score(const std::vector<RowOp>& ops) {
    auto canonical = ops;
    std::sort(canonical.begin(), canonical.end());
    std::string key;
    for (const auto& op : canonical) {
        key.push_back(static_cast<char>(op.first));
        key.push_back(static_cast<char>(op.second));
    }
    const auto found = surface_score_cache.find(key);
    if (found != surface_score_cache.end()) return found->second;
    const auto value = layer_matrix_surface_pack_score_uncached(ops);
    if (surface_score_cache.size() == 4096) {
        surface_score_cache.erase(surface_score_order.front());
        surface_score_order.pop_front();
    }
    surface_score_order.push_back(key);
    surface_score_cache.emplace(std::move(key), value);
    return value;
}

inline std::vector<std::size_t> ordered_layer_edges(
    const std::vector<LayerMatrixEdge>& edges,
    int order_kind,
    std::mt19937_64& rng
) {
    std::vector<std::size_t> order(edges.size());
    for (std::size_t i = 0; i < edges.size(); ++i) order[i] = i;

    auto deterministic_less = [&](std::size_t ai, std::size_t bi) {
        const auto& a = edges[ai];
        const auto& b = edges[bi];

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
            192
        );
        std::vector<std::pair<double, std::size_t>> scored;
        scored.reserve(hot);
        for (std::size_t i = 0; i < hot; ++i) {
            const auto& edge = edges[order[i]];
            const double jitter = static_cast<double>(rng() & 0xffffU) / 65536.0;

            scored.emplace_back(edge.global_delta + 0.25 * edge.side_delta + jitter, order[i]);
        }
        std::sort(scored.begin(), scored.end(), [&](const auto& a, const auto& b) {
            if (std::abs(a.first - b.first) > 1.0e-9) return a.first > b.first;
            return deterministic_less(a.second, b.second);
        });
        for (std::size_t i = 0; i < hot; ++i) order[i] = scored[i].second;
    }
    return order;
}

inline bool build_layer_matrix_candidate(
    const Matrix& residual,
    const Matrix& residual_inverse,
    OpSide side,
    ShiFengCostKind kind,
    const std::vector<LayerMatrixEdge>& edges,
    const std::vector<std::size_t>& order,
    bool fill_layer,
    int max_layer_ops,
    LayerMatrixCandidate& candidate
) {
    const int n = check_square_rows(residual);
    const auto current_costs = shi_side_and_global_costs_with_inverse(residual, residual_inverse, side, kind);
    const double current_global = current_costs.global_after;
    const double current_side = current_costs.side_after;
    std::vector<char> used(static_cast<std::size_t>(n), 0);
    std::vector<RowOp> ops;
    ops.reserve(static_cast<std::size_t>(n / 2));

    for (std::size_t edge_idx : order) {
        const auto& edge = edges[edge_idx];
        if (used[static_cast<std::size_t>(edge.control)]
            || used[static_cast<std::size_t>(edge.target)]) {
            continue;
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
    const auto after_costs = shi_side_and_global_costs_with_inverse(after, after_inverse, side, kind);
    const double global_after = after_costs.global_after;
    const double side_after = after_costs.side_after;
    if (!(global_after + 1.0e-9 < current_global || side_after + 1.0e-9 < current_side)) return false;

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

    return true;
}

inline std::vector<LayerMatrixCandidate> layer_matrix_candidates(
    const Matrix& residual, const Matrix& residual_inverse,
    const MatrixSearchConfig& config, const std::vector<int>&,
    std::mt19937_64& rng
) {
    const int n = check_square_rows(residual);
    std::vector<LayerMatrixCandidate> candidates;
    std::unordered_set<std::string> seen;
    for (int side_value = 0; side_value < (config.allow_columns ? 2 : 1); ++side_value) {
        const OpSide side = side_value == 0 ? OpSide::Row : OpSide::Column;
        const auto edges = layer_matrix_edges(residual, residual_inverse, side, config.cost_kind, true);
        if (edges.empty()) continue;
        // Three deterministic orders, then eight seeded jittered orders.
        for (int order_kind = 0; order_kind < 11; ++order_kind) {
            const auto order = ordered_layer_edges(edges, order_kind, rng);
            for (bool fill_layer : {false, true}) {
                LayerMatrixCandidate candidate;
                if (!build_layer_matrix_candidate(residual, residual_inverse, side, config.cost_kind,
                                                  edges, order, fill_layer, n / 2, candidate)) continue;
                if (seen.insert(layer_matrix_candidate_key(candidate)).second)
                    candidates.push_back(std::move(candidate));
            }
        }
    }
    const auto key = [](const LayerMatrixCandidate& a) {
        return std::make_tuple(double(a.permutation_excess_after), double(a.identity_rank_after),
            double(a.surface_pack_misses), double(a.surface_path_score), a.global_after,
            -double(a.ops.size()), a.side_after, static_cast<int>(a.side));
    };
    std::sort(candidates.begin(), candidates.end(), [&](const auto& a, const auto& b) {
        if (key(a) != key(b)) return key(a) < key(b);
        return a.ops < b.ops;
    });
    if (candidates.size() > SynthesisPolicy::candidate_cap) candidates.resize(SynthesisPolicy::candidate_cap);
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

inline std::string layer_matrix_state_seen_key(const LayerMatrixState& state) {
    std::string key = matrix_key(state.residual);

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

inline bool layer_matrix_state_tie_less(const LayerMatrixState& a, const LayerMatrixState& b) {
    return matrix_key(a.residual) < matrix_key(b.residual);
}

inline auto layer_matrix_state_key(const LayerMatrixState& state, bool) {
    return std::make_tuple(is_permutation_matrix(state.residual) ? 0 : 1,
        state.permutation_excess_value, state.identity_rank_value, state.cost, -state.cnot_count);
}

inline Circuit synthesize_layer_matrix_logical_once(
    const Matrix& target_rows,
    const Layout& layout,
    const MatrixSearchConfig& config,
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
    const std::size_t beam_width = SynthesisPolicy::beam_width;
    const int depth_limit = std::min(max_layers, SynthesisPolicy::depth_limit);
    bool have_best = false;
    Circuit best;

    auto keep_if_complete = [&](const LayerMatrixState& state) {
        if (!is_permutation_matrix(state.residual)) return;
        Circuit result = reconstruct_logical_search_result(
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

    for (int depth = 0; depth <= depth_limit; ++depth) {

        for (std::size_t state_index = 0; state_index < beam.size(); ++state_index) {
            const auto& state = beam[state_index];
            keep_if_complete(state);
            if (depth < depth_limit) try_depth_one_tail(state);
            if (depth + 1 < depth_limit) try_depth_two_tail(state, nullptr);

        }
        if (depth == depth_limit) break;

        std::vector<LayerMatrixState> next;
        std::unordered_set<std::string> seen;
        for (std::size_t state_index = 0; state_index < beam.size(); ++state_index) {
            const auto& state = beam[state_index];
            std::vector<LayerMatrixCandidate> generated_candidates;
            const std::vector<LayerMatrixCandidate>* candidates = nullptr;
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
            return layer_matrix_state_tie_less(a, b);
        });
        if (next.size() > beam_width) next.resize(beam_width);
        beam = std::move(next);
    }
    if (have_best) return best;
    throw std::runtime_error("layer-matrix logical search did not reach a permutation");
}

inline Circuit synthesize_inverse_row_layer_matrix_logical_once(
    const Matrix& target_rows,
    const Layout& layout,
    const MatrixSearchConfig& base_config,
    std::uint64_t tie_seed,
    int max_layers
) {
    check_square_rows(target_rows);
    MatrixSearchConfig config = base_config;
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
    const std::size_t beam_width = SynthesisPolicy::beam_width;
    const int depth_limit = std::min(max_layers, SynthesisPolicy::depth_limit);
    bool have_best = false;
    Circuit best;

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
        Circuit result = make_result_from_logical_layers(target_rows, layout, state.row_layers);
        collect_surface_logical_candidate(result);
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

    for (int depth = 0; depth <= depth_limit; ++depth) {

        for (std::size_t state_index = 0; state_index < beam.size(); ++state_index) {
            const auto& state = beam[state_index];
            keep_if_complete(state);
            if (depth < depth_limit) try_depth_one_tail(state);
            if (depth + 1 < depth_limit) try_depth_two_tail(state, nullptr);

        }
        if (depth == depth_limit) break;

        std::vector<LayerMatrixState> next;
        std::unordered_set<std::string> seen;
        for (std::size_t state_index = 0; state_index < beam.size(); ++state_index) {
            const auto& state = beam[state_index];
            std::vector<LayerMatrixCandidate> generated_candidates;
            const std::vector<LayerMatrixCandidate>* candidates = nullptr;
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
            return layer_matrix_state_tie_less(a, b);
        });
        if (next.size() > beam_width) next.resize(beam_width);
        beam = std::move(next);
    }
    if (have_best) return best;
    throw std::runtime_error("inverse row layer-matrix search did not reach a permutation");
}

inline Circuit schedule_layer_matrix_surface_bounded(
    const Circuit& logical, const std::string& mode, const Layout& layout,
    int trials, std::uint64_t seed
) {
    if (trials < 1 || trials > 64) throw std::runtime_error("bounded surface trials must be in [1,64]");
    const auto flat = flatten_layers(logical.layers);
    auto expected_multiset = flat;
    std::sort(expected_multiset.begin(), expected_multiset.end());
    bool have_best = false;
    Circuit best;
    const int lower = mode == "vdp" ? 2 * logical_sequence_bounds(flat, logical.n).lower_bound : -1;
    for (int trial = 0; trial < trials; ++trial) {
        Circuit candidate;
        std::string route;
        try {
            if (trial == 0) {
                candidate = schedule_logical_result_on_surface_basic(logical, mode, layout, mode);
                route = "basic";
            } else {
                std::vector<Layer> layers;
                if (trial == 1) {
                    layers = schedule_ops(flat, layout, mode);
                    route = "flat";
                } else if (trial % 2 == 0) {
                    const int kind = static_cast<int>((seed + static_cast<std::uint64_t>(trial / 2)) % 5);
                    layers = schedule_surface_commuting_ops(flat, layout, mode, kind);
                    route = "commuting-" + std::to_string(kind);
                } else {
                    const int kind = static_cast<int>((seed + static_cast<std::uint64_t>(trial / 2)) % 14);
                    layers = schedule_ops(topological_reorder_ops(flat, layout, kind), layout, mode);
                    route = "topological-" + std::to_string(kind);
                }
                candidate = logical;
                candidate.mode = mode;
                candidate.layout = layout;
                candidate.layers = std::move(layers);
                candidate.final_rows = identity_matrix(logical.n);
                for (const auto& layer : candidate.layers) apply_layer(candidate.final_rows, layer.ops);
            }
            auto routed_multiset = flatten_layers(candidate.layers);
            std::sort(routed_multiset.begin(), routed_multiset.end());
            if (candidate.final_rows != logical.final_rows || routed_multiset != expected_multiset
                || candidate.output_permutation != logical.output_permutation
                || !verify_geometry_and_stats(candidate)) {
                throw std::runtime_error("bounded archive route failed full verification");
            }
        } catch (const std::exception&) {
            // One optional routing failure consumes one trial, never the incumbent.
            continue;
        }
        const auto key = std::make_tuple(candidate.cycles(), candidate.layer_count(), candidate.cnot_count());
        const auto best_key = std::make_tuple(best.cycles(), best.layer_count(), best.cnot_count());
        if (!have_best || key < best_key) {
            // Publish outside ordinary routing catches. A timeout after this
            // point can retain the independently checked complete circuit.
            best = std::move(candidate);
            have_best = true;
        }
        if (have_best && lower >= 0 && best.cycles() == lower) break;
    }
    if (!have_best) throw std::runtime_error("no complete bounded archive route");
    return best;
}

// Run exactly one configuration/seed on the original target.  In particular,
// this does not expand target permutations or invoke other synthesis families.
inline Circuit synthesize_matrix(const Matrix& target, const Layout& layout, std::uint64_t seed) {
    (void)invert_matrix(target);
    synthesis_layout = &layout;
    surface_score_cache.clear();
    surface_score_order.clear();
    const std::uint64_t tie_seed = seed == 0 ? std::numeric_limits<std::uint64_t>::max()
        : 0x9e3779b97f4a7c15ULL ^ (seed * 0xd6e8feb86659fd93ULL);
    SurfaceLogicalCandidateCollector collector(SynthesisPolicy::archive_capacity);
    bool have_logical = false;
    Circuit best_logical;
    const MatrixSearchConfig config;
    {
        SurfaceLogicalCandidateScope scope(&collector);
        for (int branch = 0; branch < 2; ++branch) {
            try {
                Circuit candidate = branch == 0
                    ? synthesize_inverse_row_layer_matrix_logical_once(target, layout, config, tie_seed, SynthesisPolicy::depth_limit)
                    : synthesize_layer_matrix_logical_once(target, layout, config, tie_seed, SynthesisPolicy::depth_limit);
                if (!have_logical || logical_result_priority_key(candidate) < logical_result_priority_key(best_logical)) {
                    best_logical = std::move(candidate);
                    have_logical = true;
                }
            } catch (const std::runtime_error&) {
                // A beam may exhaust its depth without completing a circuit.
            }
        }
    }
    if (!have_logical) throw std::runtime_error("matrix search found no complete candidate");
    auto candidates = collector.take_with_best(std::move(best_logical));
    Circuit best;
    bool have_best = false;
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        auto candidate = schedule_layer_matrix_surface_bounded(candidates[index], "vdp", layout, SynthesisPolicy::route_trials, tie_seed + index);
        const auto key = [](const Circuit& value) {
            return std::make_tuple(value.cycles(), value.layer_count(), value.cnot_count());
        };
        if (!have_best || key(candidate) < key(best)) {
            best = std::move(candidate);
            have_best = true;
        }
    }
    return best;
}

} // namespace linear_surface
