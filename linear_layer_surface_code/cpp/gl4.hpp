#ifndef MIXCOLUMN_REFINE_GL4_H
#define MIXCOLUMN_REFINE_GL4_H

// Complete four-wire CNOT synthesis, used by the local rewriting stage.
// A 4x4 binary matrix is stored in four row-major nibbles. The identity is
// 0x8421, and appending CNOT c->t XORs nibble c into nibble t. All twelve
// directed CNOTs are self-inverse, so one unweighted BFS gives an exact
// shortest gate word for every member of GL(4,2), without a length cutoff.
#include <algorithm>
#include <array>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gl4 {

using Code = std::uint16_t;
using Gate = std::pair<int, int>;
using Word = std::vector<Gate>;
using Permutation = std::array<int, 4>;
constexpr Code identity = 0x8421;
constexpr std::uint8_t unreachable = 255;

inline Code append_cnot(Code matrix, int control, int target) {
    if (control < 0 || control >= 4 || target < 0 || target >= 4 || control == target) {
        throw std::runtime_error("GL4 CNOT has invalid endpoints");
    }
    const unsigned control_row = (unsigned(matrix) >> (4 * control)) & 15U;
    return static_cast<Code>(unsigned(matrix) ^ (control_row << (4 * target)));
}

inline Code replay(const Word& word) {
    Code matrix = identity;
    for (const auto& gate : word) matrix = append_cnot(matrix, gate.first, gate.second);
    return matrix;
}

template<class Rows>
Code encode_rows(const Rows& rows) {
    if (rows.size() != 4) throw std::runtime_error("GL4 encoding requires exactly four rows");
    unsigned code = 0;
    for (int row = 0; row < 4; ++row) {
        const std::uint64_t bits = static_cast<std::uint64_t>(rows[static_cast<std::size_t>(row)]);
        if (bits > 15) {
            throw std::runtime_error("GL4 matrix row exceeds four bits");
        }
        code |= unsigned(bits) << (4 * row);
    }
    return static_cast<Code>(code);
}

class ShortestTable {
    std::array<Code, 65536> predecessor_{};
    std::array<std::uint8_t, 65536> edge_{};
    std::array<std::uint8_t, 65536> distance_{};
    std::size_t states_ = 0;
    int diameter_ = 0;

public:
    ShortestTable() {
        distance_.fill(unreachable);
        edge_.fill(unreachable);
        // Queue size is bounded by all binary 4x4 matrices, not by words.
        std::array<Code, 65536> queue{};
        std::size_t head = 0, tail = 0;
        queue[tail++] = identity;
        distance_[identity] = 0;
        predecessor_[identity] = identity;
        while (head < tail) {
            const Code matrix = queue[head++];
            for (int control = 0; control < 4; ++control) {
                for (int target = 0; target < 4; ++target) {
                    if (control == target) continue;
                    const Code next = append_cnot(matrix, control, target);
                    if (distance_[next] != unreachable) continue;
                    predecessor_[next] = matrix;
                    edge_[next] = static_cast<std::uint8_t>(4 * control + target);
                    distance_[next] = static_cast<std::uint8_t>(distance_[matrix] + 1);
                    diameter_ = std::max(diameter_, int(distance_[next]));
                    queue[tail++] = next;
                }
            }
        }
        states_ = tail;
        if (states_ != 20160) throw std::runtime_error("GL4 BFS did not cover exactly 20160 invertible matrices");
    }

    std::size_t state_count() const { return states_; }
    int diameter() const { return diameter_; }
    bool contains(Code matrix) const { return distance_[matrix] != unreachable; }
    int distance(Code matrix) const {
        if (!contains(matrix)) throw std::runtime_error("GL4 target is singular");
        return distance_[matrix];
    }
    Word shortest(Code matrix) const {
        const Code requested = matrix;
        Word word;
        word.reserve(static_cast<std::size_t>(distance(matrix)));
        while (matrix != identity) {
            const int edge = edge_[matrix];
            word.emplace_back(edge / 4, edge % 4);
            matrix = predecessor_[matrix];
        }
        std::reverse(word.begin(), word.end());
        if (replay(word) != requested) throw std::runtime_error("GL4 shortest representative failed exact replay");
        return word;
    }

};

inline const ShortestTable& table() {
    // The immutable exact table is constructed once per process.
    static const ShortestTable value;
    return value;
}

inline Code conjugate_for_labels(Code matrix, const Permutation& labels) {
    unsigned used = 0;
    for (int wire : labels) {
        if (wire < 0 || wire >= 4 || (used & (1U << wire))) {
            throw std::runtime_error("GL4 relabeling is not a permutation");
        }
        used |= 1U << wire;
    }
    // B[i,j] = A[labels[i],labels[j]]. Relabeling a B word's endpoints by
    // labels therefore implements A, including noninvolutive permutations.
    unsigned result = 0;
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            const unsigned bit = (unsigned(matrix) >> (4 * labels[row] + labels[column])) & 1U;
            result |= bit << (4 * row + column);
        }
    }
    return static_cast<Code>(result);
}

inline Code left_row_permute(Code matrix, const Permutation& rows) {
    unsigned used = 0;
    for (int wire : rows) {
        if (wire < 0 || wire >= 4 || (used & (1U << wire))) {
            throw std::runtime_error("GL4 row permutation is not a permutation");
        }
        used |= 1U << wire;
    }
    const unsigned input = unsigned(matrix);
    unsigned output = 0;
    for (int row = 0; row < 4; ++row) {
        output |= ((input >> (4 * rows[static_cast<std::size_t>(row)])) & 15U) << (4 * row);
    }
    return static_cast<Code>(output);
}

inline std::vector<Word> shortest_variants(Code matrix, std::uint64_t seed, int limit = 8) {
    if (limit < 1 || limit > 24) throw std::runtime_error("GL4 variant limit must be in [1,24]");
    const auto& database = table();
    const int expected_length = database.distance(matrix);
    std::vector<Permutation> permutations;
    Permutation labels{0, 1, 2, 3};
    do { permutations.push_back(labels); }
    while (std::next_permutation(labels.begin(), labels.end()));
    // Identity labeling first always supplies the canonical shortest word.
    // The seed only orders the other 23 physically different labelings.
    std::mt19937_64 rng(seed);
    std::shuffle(permutations.begin() + 1, permutations.end(), rng);
    std::vector<Word> variants;
    variants.reserve(static_cast<std::size_t>(limit));
    for (const auto& permutation : permutations) {
        Word word = database.shortest(conjugate_for_labels(matrix, permutation));
        for (auto& gate : word) {
            gate.first = permutation[static_cast<std::size_t>(gate.first)];
            gate.second = permutation[static_cast<std::size_t>(gate.second)];
        }
        if (static_cast<int>(word.size()) != expected_length || replay(word) != matrix) {
            throw std::runtime_error("GL4 wire-conjugation representative failed exact replay/length");
        }
        if (std::find(variants.begin(), variants.end(), word) == variants.end()) {
            variants.push_back(std::move(word));
            if (static_cast<int>(variants.size()) == limit) break;
        }
    }
    return variants;
}

}  // namespace gl4

#endif
