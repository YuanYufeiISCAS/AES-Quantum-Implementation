#pragma once
#include "routing.hpp"
#include "protocol.hpp"
#include "gl4.hpp"
#include "free_output.hpp"

namespace linear_surface {
struct LocalReductionDatabase {
    int k = 0;
    int max_len = 0;
    std::unordered_map<std::string, std::vector<std::vector<RowOp>>> by_matrix;
    std::vector<RowOp> generators;

    LocalReductionDatabase(int local_qubits, int max_sequence_len)
        : k(local_qubits), max_len(max_sequence_len) {
        if (k < 2 || k > 4 || max_len != 5) throw std::runtime_error("local reduction arity out of range");
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

inline const LocalReductionDatabase& local_reduction_db(int k, int max_len) {
    static std::unordered_map<int, std::unique_ptr<const LocalReductionDatabase>> tables;
    const int key = 100 * k + max_len;
    auto found = tables.find(key);
    if (found == tables.end())
        found = tables.emplace(key, std::make_unique<const LocalReductionDatabase>(k, max_len)).first;
    return *found->second;
}

struct RefinementInput {
    Problem problem;
    Circuit incumbent;
    bool free_output = false;
};
struct RefineStatistics {
    std::uint64_t proposals = 0;
    std::uint64_t route_calls = 0;
    std::uint64_t rounds = 0;
    std::uint64_t walk_accepts = 0;
};
inline std::string refine_word_key(const std::vector<RowOp>& ops) {
    std::string key;
    key.reserve(2 * ops.size());
    for (const auto& op : ops) {
        key.push_back(static_cast<char>(op.first));
        key.push_back(static_cast<char>(op.second));
    }
    return key;
}

inline std::vector<RowOp> refine_gate_multiset(std::vector<RowOp> ops) {
    std::sort(ops.begin(), ops.end());
    return ops;
}

inline std::vector<RowOp> refine_cancel(std::vector<RowOp> ops) {
    for (int pass = 0; pass < 4; ++pass) {
        auto next = optimize_cnot_sequence(ops);
        if (next == ops) break;
        ops = std::move(next);
    }
    return ops;
}

// Expose a small-support subsequence by swapping only proven commuting gates.
// This changes no algebraic operation; an ensuing local matrix rewrite must
// change the gate multiset before any expensive surface route is requested.
inline void refine_gather_support(std::vector<RowOp>& ops, std::mt19937_64& rng) {
    if (ops.size() < 3) return;
    const std::size_t begin = static_cast<std::size_t>(rng() % ops.size());
    const std::size_t end = std::min(ops.size(), begin + 48);
    std::vector<int> support{ops[begin].first, ops[begin].second};
    for (std::size_t index = begin + 1; index < end && support.size() < 4; ++index) {
        const auto& op = ops[index];
        const bool control = std::find(support.begin(), support.end(), op.first) != support.end();
        const bool target = std::find(support.begin(), support.end(), op.second) != support.end();
        if (control != target) support.push_back(control ? op.second : op.first);
    }
    std::size_t insertion = begin + 1;
    for (std::size_t index = begin + 1; index < end; ++index) {
        const RowOp op = ops[index];
        if (std::find(support.begin(), support.end(), op.first) == support.end()
            || std::find(support.begin(), support.end(), op.second) == support.end()) continue;
        bool can_move = true;
        for (std::size_t barrier = insertion; barrier < index; ++barrier) {
            if (!cnot_commute(op, ops[barrier])) { can_move = false; break; }
        }
        if (!can_move) continue;
        std::rotate(ops.begin() + static_cast<std::ptrdiff_t>(insertion),
                    ops.begin() + static_cast<std::ptrdiff_t>(index),
                    ops.begin() + static_cast<std::ptrdiff_t>(index + 1));
        ++insertion;
    }
}

class Refiner {

    const RefinementInput& input_;
    const RefinementBudget& options_;
    RefineStatistics stats_;
    Circuit best_;
    std::vector<RowOp> best_word_;
    std::vector<RowOp> best_multiset_;
    std::vector<int> best_permutation_;
    Matrix best_final_rows_;
    const std::vector<RowOp> input_multiset_;
    std::unordered_set<std::string> seen_words_;
    std::mt19937_64 rng_;
    Circuit walk_result_;
    std::vector<RowOp> walk_word_;
    std::vector<int> walk_permutation_;
    Matrix walk_final_rows_;
    std::unordered_set<std::string> exploration_keys_;

    bool exhausted() const {
        return stats_.proposals >= static_cast<std::uint64_t>(options_.max_proposals)
            || stats_.route_calls >= static_cast<std::uint64_t>(options_.max_routes);
    }
    bool routing_exhausted() const {
        return stats_.route_calls >= static_cast<std::uint64_t>(options_.max_routes);
    }
    static std::tuple<int, int, int> score(const Circuit& result) {
        return {result.cycles(), result.layer_count(), result.cnot_count()};
    }

    void emit(const Circuit& result, const std::string& source, bool) {
        verify_result(input_.problem, result);
        std::cout << "{\"kind\":\"candidate\",\"result\":";
        write_json_to_stream(result, std::cout);
        std::cout << "}\n" << std::flush;
        std::cerr << source << "  " << result.cnot_count() << " CNOT / "
                  << result.layer_count() << " layers\n";
    }

    // Export a small set of accepted exploratory states for population
    // restarts. Export never consumes randomness or changes the local walk.
    void emit_exploration(const Circuit& result) {
        if (exploration_keys_.size() == 8) return;
        const auto multiset = refine_gate_multiset(flatten_layers(result.layers));
        if ((multiset == input_multiset_ && result.output_permutation == input_.incumbent.output_permutation)
            || (multiset == best_multiset_ && result.output_permutation == best_permutation_)) return;
        std::string key = refine_word_key(multiset);
        key.push_back(static_cast<char>(255));
        for (int wire : result.output_permutation) key.push_back(static_cast<char>(wire));
        if (!exploration_keys_.insert(std::move(key)).second) return;
        verify_result(input_.problem, result);
        std::cout << "{\"kind\":\"exploration\",\"result\":";
        write_json_to_stream(result, std::cout);
        std::cout << "}\n" << std::flush;
    }

    void consider_surface(std::vector<Layer> layers, const std::vector<RowOp>& expected_word,
                          const std::vector<int>& expected_permutation,
                          const std::string& source) {
        Circuit candidate = best_;
        candidate.layers = std::move(layers);
        candidate.output_permutation = expected_permutation;
        candidate.final_rows = identity_matrix(candidate.n);
        for (const auto& layer : candidate.layers) apply_layer(candidate.final_rows, layer.ops);
        verify_result(input_.problem, candidate);
        if (refine_gate_multiset(flatten_layers(candidate.layers)) != refine_gate_multiset(expected_word)) {
            throw std::runtime_error("routing did not preserve the proposed algebraically changed gate multiset");
        }

        const bool improved = score(candidate) < score(best_);
        if (improved) {

            // Keep the full verified best independently of any exploratory state.
            best_ = candidate;
            best_word_ = flatten_layers(best_.layers);
            best_multiset_ = refine_gate_multiset(best_word_);
            best_permutation_ = best_.output_permutation;
            best_final_rows_ = best_.final_rows;
            emit(best_, source, true);
        }
        if (options_.walk_slack > 0
            && candidate.cycles() <= best_.cycles() + options_.walk_slack
            && candidate.cnot_count() <= best_.cnot_count() + options_.max_growth) {
            // Every accepted state is itself an exact, fully routed circuit.
            // Neutral/uphill states are search diagnostics, never emitted as winners.
            const bool downhill = score(candidate) < score(walk_result_);
            const int depth_delta = candidate.cycles() - walk_result_.cycles();
            const int gate_delta = candidate.cnot_count() - walk_result_.cnot_count();
            const double energy_delta = double(depth_delta) + .05 * double(gate_delta);
            const double temperature = .75 + .75 * double(stats_.rounds % 8) / 7.0;
            const double probability = std::min(.35, .20 * std::exp(-std::max(0.0, energy_delta) / temperature));
            const bool accept = improved || downhill
                || std::generate_canonical<double, 53>(rng_) < probability;
            if (accept) {
                walk_result_ = std::move(candidate);
                walk_word_ = flatten_layers(walk_result_.layers);
                walk_permutation_ = walk_result_.output_permutation;
                walk_final_rows_ = walk_result_.final_rows;
                ++stats_.walk_accepts;
                if (!improved) emit_exploration(walk_result_);
            }
        }
    }
    void route_changed_word(const std::vector<RowOp>& word,
                            const std::vector<int>& permutation,
                            const std::string& source) {
        if (routing_exhausted()) return;
        ++stats_.route_calls;

        auto layers = schedule_surface_commuting_ops(word, input_.problem.layout, "vdp", 0);
        consider_surface(std::move(layers), word, permutation, source + "/route=commuting-0");
    }

    void propose(std::vector<RowOp> word, std::vector<int> permutation, const std::string& source) {
        if (exhausted()) return;
        ++stats_.proposals;
        word = refine_cancel(std::move(word));
        if (word.size() > best_word_.size() + static_cast<std::size_t>(options_.max_growth)) return;
        const auto multiset = refine_gate_multiset(word);
        const bool unchanged_best = multiset == best_multiset_
            && (!input_.free_output || permutation == best_permutation_);
        const bool unchanged_input = multiset == input_multiset_
            && (!input_.free_output || permutation == input_.incumbent.output_permutation);
        if (unchanged_best || unchanged_input) {

            return;
        }
        std::string key = refine_word_key(word);
        if (input_.free_output) {
            key.push_back(static_cast<char>(255));
            for (int wire : permutation) key.push_back(static_cast<char>(wire));
        }
        // Finite attempt budget plus a bounded seen set; no unbounded history.
        if (seen_words_.find(key) != seen_words_.end()) { return; }
        if (seen_words_.size() == 8192) seen_words_.clear();
        seen_words_.insert(key);
        const Matrix expected_final = input_.free_output
            ? permute_rows(input_.problem.target, permutation) : input_.incumbent.final_rows;
        if (apply_ops_to_identity(input_.problem.layout.n, word) != expected_final) {

            return;
        }

        if (!word.empty()) {
            const auto bound = logical_sequence_bounds(word, input_.problem.layout.n);
            const int surface_lower_bound = 2 * bound.lower_bound;
            if (surface_lower_bound > best_.cycles() + options_.walk_slack
                || (options_.walk_slack == 0 && surface_lower_bound == best_.cycles()
                    && word.size() >= static_cast<std::size_t>(best_.cnot_count()))) {

                return;
            }
        }
        route_changed_word(word, permutation, source);
    }

public:
    Refiner(const RefinementInput& input, const RefinementBudget& options)
        : input_(input), options_(options), best_(input.incumbent),
          best_word_(flatten_layers(best_.layers)), best_multiset_(refine_gate_multiset(best_word_)),
          best_permutation_(best_.output_permutation), best_final_rows_(best_.final_rows),
          input_multiset_(best_multiset_),
          rng_(options.seed), walk_result_(best_), walk_word_(best_word_),
          walk_permutation_(best_.output_permutation), walk_final_rows_(best_.final_rows) {}

    void run() {
        emit(best_, "input-incumbent", false);
        for (int round = 0; round < options_.rounds && !exhausted(); ++round) {
            ++stats_.rounds;
            if (options_.walk_slack > 0 && round > 0 && round % options_.walk_restart_rounds == 0) {
                walk_result_ = best_;
                walk_word_ = best_word_;
                walk_permutation_ = best_permutation_;
                walk_final_rows_ = best_final_rows_;

            }
            const std::uint64_t generation = stats_.walk_accepts;
            std::vector<RowOp> base = options_.walk_slack > 0 ? walk_word_ : best_word_;
            // A greedy improvement may change best P while this round still
            // visits windows of the previous base. Bind every proposal to the
            // base's own P rather than whichever best/walk happens to be current.
            const std::vector<int> base_permutation = options_.walk_slack > 0
                ? walk_permutation_ : best_permutation_;
            if (input_.free_output) {
                const Matrix& base_final = options_.walk_slack > 0 ? walk_final_rows_ : best_final_rows_;
                if (base_final != permute_rows(input_.problem.target, base_permutation)
                    || apply_ops_to_identity(input_.problem.layout.n, base) != base_final) {
                    throw std::runtime_error("free-output search state lost its circuit/permutation invariant");
                }
            }
            if (round > 0) {
                const int kind = int(rng_() % 64);
                base = topological_reorder_ops(base, input_.problem.layout, kind);
                for (int gather = 0; gather < 4; ++gather) refine_gather_support(base, rng_);
            }
            propose(base, base_permutation, "commute-cancel/round=" + std::to_string(round));
            if (options_.walk_slack > 0 && stats_.walk_accepts != generation) continue;

            struct Window { int begin; int length; std::vector<int> active; };
            std::vector<Window> windows;
            for (int begin = 0; begin < static_cast<int>(base.size()); ++begin) {
                if (exhausted()) break;
                for (int length = 2; length <= options_.max_window
                     && begin + length <= static_cast<int>(base.size()); ++length) {
                    auto active = active_qubits_for_segment(base, begin, length);
                    if (active.size() > 4) break;
                    if (active.size() >= 2) windows.push_back({begin, length, std::move(active)});
                }
            }
            std::shuffle(windows.begin(), windows.end(), rng_);
            for (const auto& window : windows) {
                if (exhausted()) break;
                const Matrix local = segment_matrix(base, window.begin, window.length,
                                                     window.active, input_.problem.layout.n);
                if (window.active.size() == 4) {

                    const auto encoded = gl4::encode_rows(local);
                    if (input_.free_output) {

                        using Permutation = gl4::Permutation;
                        std::vector<Permutation> permutations;
                        Permutation q{0, 1, 2, 3};
                        do { permutations.push_back(q); }
                        while (std::next_permutation(q.begin(), q.end()));
                        // Logical length only prioritizes a bounded route
                        // budget. Acceptance always uses final surface cost.
                        std::stable_sort(permutations.begin(), permutations.end(), [&](const auto& a, const auto& b) {
                            return gl4::table().distance(gl4::left_row_permute(encoded, a))
                                 < gl4::table().distance(gl4::left_row_permute(encoded, b));
                        });
                        for (const auto& permutation : permutations) {
                            if (exhausted()) break;

                            const auto replacement = gl4::table().shortest(
                                gl4::left_row_permute(encoded, permutation));
                            if (replacement.size() > static_cast<std::size_t>(window.length + options_.max_growth)) continue;
                            auto trial = free_output::replace_window(
                                base, window.begin, window.length, window.active,
                                replacement, permutation, base_permutation);
                            if (trial.word == base && trial.output_permutation == base_permutation) continue;
                            std::string q_text;
                            for (int wire : permutation) q_text.push_back(static_cast<char>('0' + wire));
                            const std::string source = "free-output-gl4/round=" + std::to_string(round)
                                + "/begin=" + std::to_string(window.begin) + "/old=" + std::to_string(window.length)
                                + "/new=" + std::to_string(replacement.size()) + "/arity=4/q=" + q_text;

                            propose(std::move(trial.word), std::move(trial.output_permutation), source);
                            if ((options_.walk_slack > 0 && stats_.walk_accepts != generation)) break;
                        }

                        if (exhausted() || (options_.walk_slack > 0 && stats_.walk_accepts != generation)) break;
                    }
                    const auto representative_seed = rng_();
                    const auto full_representatives = gl4::shortest_variants(encoded, representative_seed, options_.gl4_variants);
                    for (std::size_t variant = 0; variant < full_representatives.size(); ++variant) {
                        if (exhausted()) break;
                        const auto replacement = map_local_sequence(full_representatives[variant], window.active);
                        if (replacement.size() > static_cast<std::size_t>(window.length + options_.max_growth)
                            || segment_equals(base, window.begin, window.length, replacement)) continue;
                        const auto trial = replace_segment(base, window.begin, window.length, replacement);
                        const std::string source = std::string("full-gl4/round=") + std::to_string(round)
                            + "/begin=" + std::to_string(window.begin) + "/old=" + std::to_string(window.length)
                            + "/new=" + std::to_string(replacement.size()) + "/arity=4/variant=" + std::to_string(variant);

                        propose(trial, base_permutation, source);
                        if ((options_.walk_slack > 0 && stats_.walk_accepts != generation)) break;
                    }

                    if (exhausted() || (options_.walk_slack > 0 && stats_.walk_accepts != generation)) break;
                }
                // Fixed 2..4 wires / <=5 gates: at most 193261 enumerated words
                // for the largest table, unlike unbounded full-circuit search.
                const auto& equivalents = local_reduction_db(static_cast<int>(window.active.size()), 5)
                                              .equivalents(local);
                std::vector<std::size_t> alternatives;
                for (std::size_t index = 0; index < equivalents.size(); ++index) {
                    if (equivalents[index].size()
                        <= static_cast<std::size_t>(window.length + options_.max_growth)) {
                        alternatives.push_back(index);
                    }
                }
                // Shortest representative plus sampled equal/longer words:
                // a CNOT-count tie is not a reason to discard a surface candidate.
                if (alternatives.size() > 1) {
                    std::shuffle(alternatives.begin() + 1, alternatives.end(), rng_);
                }
                if (alternatives.size() > 16) alternatives.resize(16);
                for (std::size_t index : alternatives) {
                    if (exhausted()) break;
                    const auto replacement = map_local_sequence(equivalents[index], window.active);
                    if (segment_equals(base, window.begin, window.length, replacement)) continue;
                    const auto trial = replace_segment(base, window.begin, window.length, replacement);
                    const std::string source = "local-matrix/round=" + std::to_string(round)
                        + "/begin=" + std::to_string(window.begin) + "/old=" + std::to_string(window.length)
                        + "/new=" + std::to_string(replacement.size()) + "/arity=" + std::to_string(window.active.size());
                    propose(trial, base_permutation, source);
                    if ((options_.walk_slack > 0 && stats_.walk_accepts != generation)) break;
                }

                if (options_.walk_slack > 0 && stats_.walk_accepts != generation) break;
            }
            if (base.empty()) break;

        }
        std::cout << "{\"kind\":\"complete\",\"proposals\":" << stats_.proposals
                  << ",\"routes\":" << stats_.route_calls
                  << ",\"rounds\":" << stats_.rounds << "}\n";
    }
};

} // namespace linear_surface
