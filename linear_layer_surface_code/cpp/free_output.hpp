#ifndef MIXCOLUMN_REFINE_FREE_OUTPUT_H
#define MIXCOLUMN_REFINE_FREE_OUTPUT_H

// Algebraic free-output GL4 window rewrite. This changes a CNOT word and its
// output permutation only: logical wire IDs still denote the same fixed
// physical data positions. No SWAP, mapping update, or extra qubit is added.
#include "gl4.hpp"

namespace free_output {

struct Candidate {
    gl4::Word word;
    std::vector<int> output_permutation;
};

inline Candidate replace_window(
    const gl4::Word& word,
    int begin,
    int length,
    const std::vector<int>& active,
    const gl4::Word& replacement,
    const gl4::Permutation& q,
    const std::vector<int>& output_permutation
) {
    const int n = static_cast<int>(output_permutation.size());
    if (n < 4 || n > 64 || begin < 0 || length < 0
        || static_cast<std::size_t>(begin) > word.size()
        || static_cast<std::size_t>(length) > word.size() - static_cast<std::size_t>(begin)
        || active.size() != 4) {
        throw std::runtime_error("invalid free-output window bounds");
    }
    std::vector<bool> used(static_cast<std::size_t>(n), false);
    for (int wire : output_permutation) {
        if (wire < 0 || wire >= n || used[static_cast<std::size_t>(wire)]) {
            throw std::runtime_error("invalid free-output incumbent permutation");
        }
        used[static_cast<std::size_t>(wire)] = true;
    }
    std::fill(used.begin(), used.end(), false);
    for (int wire : active) {
        if (wire < 0 || wire >= n || used[static_cast<std::size_t>(wire)]) {
            throw std::runtime_error("invalid free-output active support");
        }
        used[static_cast<std::size_t>(wire)] = true;
    }
    // Reuse the GL4 validator for q, including noninvolutive permutations.
    (void)gl4::left_row_permute(gl4::identity, q);
    std::vector<int> endpoint_map(static_cast<std::size_t>(n));
    for (int wire = 0; wire < n; ++wire) endpoint_map[static_cast<std::size_t>(wire)] = wire;
    Candidate candidate;
    candidate.output_permutation = output_permutation;
    for (int local = 0; local < 4; ++local) {
        const int global = active[static_cast<std::size_t>(local)];
        const int previous = active[static_cast<std::size_t>(q[static_cast<std::size_t>(local)])];
        endpoint_map[static_cast<std::size_t>(previous)] = global; // q^{-1}
        candidate.output_permutation[static_cast<std::size_t>(global)] =
            output_permutation[static_cast<std::size_t>(previous)];
    }
    candidate.word.reserve(word.size() - static_cast<std::size_t>(length) + replacement.size());
    candidate.word.insert(candidate.word.end(), word.begin(), word.begin() + begin);
    for (const auto& gate : replacement) {
        if (gate.first < 0 || gate.first >= 4 || gate.second < 0 || gate.second >= 4
            || gate.first == gate.second) {
            throw std::runtime_error("invalid free-output local replacement CNOT");
        }
        candidate.word.emplace_back(active[static_cast<std::size_t>(gate.first)],
                                    active[static_cast<std::size_t>(gate.second)]);
    }
    for (std::size_t index = static_cast<std::size_t>(begin + length); index < word.size(); ++index) {
        const auto& gate = word[index];
        if (gate.first < 0 || gate.first >= n || gate.second < 0 || gate.second >= n
            || gate.first == gate.second) {
            throw std::runtime_error("invalid free-output suffix CNOT");
        }
        candidate.word.emplace_back(endpoint_map[static_cast<std::size_t>(gate.first)],
                                    endpoint_map[static_cast<std::size_t>(gate.second)]);
    }
    // For F=S M T and replacement C=Q M, relabeled suffix endpoints give
    // S'=Q S Q^{-1}; hence S' C T = Q F = P_new A. The caller must replay the
    // complete candidate against P_new A before any surface route or output.
    return candidate;
}

} // namespace free_output

#endif
