#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace mixcolumn {
inline int least_set_bit_index_portable_nonzero(std::uint64_t word) noexcept {
    int index = 0;
    while ((word & std::uint64_t(1)) == 0) {
        word >>= 1;
        ++index;
    }
    return index;
}

inline int least_set_bit_index_nonzero(std::uint64_t word) noexcept {
#if defined(__GNUG__) || defined(__clang__)
    return __builtin_ctzll(word);
#else
    return least_set_bit_index_portable_nonzero(word);
#endif
}

inline int gf2_rank_valid_rows(const std::uint64_t* rows, std::size_t count) noexcept {
    std::array<std::uint64_t, 64> basis{};
    int rank = 0;
    for (std::size_t row_index = 0; row_index < count; ++row_index) {
        std::uint64_t word = rows[row_index];
        while (word != 0) {
            const int pivot = least_set_bit_index_nonzero(word);
            std::uint64_t& previous = basis[static_cast<std::size_t>(pivot)];
            if (previous == 0) {
                previous = word;
                ++rank;
                break;
            }
            word ^= previous;
        }
    }
    return rank;
}

using LogicalRowOp = std::pair<int, int>;

struct LinearLogicalSequenceBounds {
    int capacity_bound = 0;
    int endpoint_bound = 0;
    int critical_bound = 0;
    int lower_bound = 0;
};

// Time is O(m+n), auxiliary storage is O(n), and the input word is untouched.
// Like that implementation, this internal helper assumes in-range endpoints
// and counts representable by int.  It does not validate or cancel gates.
// The forward longest-path maximum equals the reverse longest-path maximum
// returned by the old full DAG even though their per-gate depths differ.
inline LinearLogicalSequenceBounds logical_sequence_bounds_linear(
    const std::vector<LogicalRowOp>& ops, int n
) {
    if (n <= 1) throw std::runtime_error("logical bounds need at least two qubits");

    const std::size_t wire_count = static_cast<std::size_t>(n);
    std::vector<int> endpoint_count(wire_count, 0);
    std::vector<int> prior_control_depth(wire_count, 0);
    std::vector<int> prior_target_depth(wire_count, 0);
    int critical_bound = 0;

    for (const LogicalRowOp& op : ops) {
        const int control = op.first;
        const int target = op.second;

        ++endpoint_count[static_cast<std::size_t>(control)];
        ++endpoint_count[static_cast<std::size_t>(target)];

        // A prior (pc, pt) is non-commuting with (control, target) exactly
        // when pc == target or pt == control.  The two arrays retain the
        // maximum critical depth of those two predecessor classes.
        const int critical_depth = 1 + std::max(
            prior_target_depth[static_cast<std::size_t>(control)],
            prior_control_depth[static_cast<std::size_t>(target)]
        );
        prior_control_depth[static_cast<std::size_t>(control)] = std::max(
            prior_control_depth[static_cast<std::size_t>(control)], critical_depth);
        prior_target_depth[static_cast<std::size_t>(target)] = std::max(
            prior_target_depth[static_cast<std::size_t>(target)], critical_depth);
        critical_bound = std::max(critical_bound, critical_depth);
    }

    const std::size_t half = static_cast<std::size_t>(n / 2);
    LinearLogicalSequenceBounds bounds;
    bounds.capacity_bound = static_cast<int>((ops.size() + half - 1) / half);
    bounds.endpoint_bound = endpoint_count.empty()
        ? 0
        : *std::max_element(endpoint_count.begin(), endpoint_count.end());
    bounds.critical_bound = critical_bound;
    bounds.lower_bound = std::max(bounds.capacity_bound,
                                  std::max(bounds.endpoint_bound, bounds.critical_bound));
    return bounds;
}

}

namespace linear_surface {
using Row = std::uint64_t;
using Matrix = std::vector<Row>;
using RowOp = std::pair<int,int>;

enum class OpSide {
    Row,
    Column
};

inline int popcount64(Row x) {
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

inline Row mask_for(int n) {
    if (n <= 0 || n > 64) {
        throw std::runtime_error("matrix dimension must be in [1, 64]");
    }
    if (n == 64) return ~Row(0);
    return (Row(1) << n) - 1;
}

inline Matrix identity_matrix(int n) {
    Matrix out(n);
    for (int i = 0; i < n; ++i) out[i] = Row(1) << i;
    return out;
}

inline int check_square_rows(const Matrix& rows) {
    const int n = static_cast<int>(rows.size());
    const Row mask = mask_for(n);
    for (Row r : rows) {
        if ((r & ~mask) != 0) throw std::runtime_error("row bitset exceeds matrix dimension");
    }
    return n;
}

inline void row_add(Matrix& rows, int control, int target) {
    if (control == target) throw std::runtime_error("control and target must differ");
    rows[target] ^= rows[control];
}

inline void column_add(Matrix& rows, int control, int target) {
    if (control == target) throw std::runtime_error("control and target must differ");
    const Row control_bit = Row(1) << control;
    const Row target_bit = Row(1) << target;
    for (Row& row : rows) {
        if ((row & target_bit) != 0) row ^= control_bit;
    }
}

inline Matrix matmul(const Matrix& a, const Matrix& b) {
    const int n = check_square_rows(a);
    if (static_cast<int>(b.size()) != n) throw std::runtime_error("dimension mismatch");
    check_square_rows(b);
    Matrix out;
    out.reserve(n);
    for (Row row : a) {
        Row acc = 0;
        Row bits = row;
        while (bits != 0) {
            const Row lsb = bits & (~bits + 1);
            const int idx = mixcolumn::least_set_bit_index_nonzero(bits);
            acc ^= b[idx];
            bits ^= lsb;
        }
        out.push_back(acc);
    }
    return out;
}

inline Matrix invert_matrix(const Matrix& rows) {
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

inline Matrix permute_rows(const Matrix& rows, const std::vector<int>& perm) {
    if (rows.size() != perm.size()) throw std::runtime_error("permutation length mismatch");
    Matrix out;
    out.reserve(rows.size());
    for (int idx : perm) out.push_back(rows.at(static_cast<std::size_t>(idx)));
    return out;
}

inline bool is_permutation_matrix(const Matrix& rows) {
    const int n = check_square_rows(rows);
    Row seen = 0;
    for (Row row : rows) {
        if (popcount64(row) != 1) return false;
        if ((seen & row) != 0) return false;
        seen |= row;
    }
    return seen == mask_for(n);
}

inline std::vector<int> permutation_columns(const Matrix& rows) {
    if (!is_permutation_matrix(rows)) throw std::runtime_error("matrix is not a permutation");
    std::vector<int> out;
    out.reserve(rows.size());
    for (Row row : rows) out.push_back(popcount64(row - 1));
    return out;
}

inline std::vector<int> inverse_permutation(const std::vector<int>& perm) {
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

inline bool infer_output_permutation(
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

inline bool logical_depth_one_ops(const Matrix& rows, std::vector<RowOp>& ops) {
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

inline void validate_permutation_vector(const std::vector<int>& perm, int n) {
    if (static_cast<int>(perm.size()) != n) throw std::runtime_error("output permutation size mismatch");
    std::vector<char> seen(static_cast<std::size_t>(n), 0);
    for (int v : perm) {
        if (v < 0 || v >= n) throw std::runtime_error("output permutation entry out of range");
        if (seen[static_cast<std::size_t>(v)]) throw std::runtime_error("output permutation has duplicate entries");
        seen[static_cast<std::size_t>(v)] = 1;
    }
}

inline std::string matrix_key(const Matrix& m) {
    std::string out;
    out.reserve(m.size() * sizeof(Row));
    for (Row row : m) {
        for (int i = 0; i < 8; ++i) {
            out.push_back(static_cast<char>((row >> (8 * i)) & 0xffU));
        }
    }
    return out;
}

inline Matrix transpose_matrix(const Matrix& rows) {
    const int n = check_square_rows(rows);
    Matrix out(n, 0);
    for (int r = 0; r < n; ++r) {
        Row bits = rows[r];
        while (bits != 0) {
            const Row lsb = bits & (~bits + 1);
            const int c = mixcolumn::least_set_bit_index_nonzero(bits);
            out[c] |= Row(1) << r;
            bits ^= lsb;
        }
    }
    return out;
}

inline bool cnot_commute(const RowOp& a, const RowOp& b) {
    const int ac = a.first;
    const int at = a.second;
    const int bc = b.first;
    const int bt = b.second;
    return at != bc && bt != ac;
}

inline std::vector<RowOp> optimize_cnot_sequence(const std::vector<RowOp>& ops) {
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

inline std::vector<RowOp> flatten_rowop_layers(
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

inline std::vector<RowOp> transform_ops_by_permutation(
    const std::vector<RowOp>& ops,
    const std::vector<int>& perm
) {
    std::vector<RowOp> out;
    out.reserve(ops.size());
    for (const auto& op : ops) out.emplace_back(perm[op.first], perm[op.second]);
    return out;
}

inline std::vector<RowOp> swap_controls_targets(const std::vector<RowOp>& ops) {
    std::vector<RowOp> out;
    out.reserve(ops.size());
    for (const auto& op : ops) out.emplace_back(op.second, op.first);
    return out;
}

inline Matrix apply_ops_to_identity(int n, const std::vector<RowOp>& ops) {
    Matrix current = identity_matrix(n);
    for (const auto& op : ops) row_add(current, op.first, op.second);
    return current;
}

inline std::vector<int> active_qubits_for_segment(
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

inline Matrix segment_matrix(
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

inline std::vector<RowOp> map_local_sequence(
    const std::vector<RowOp>& local_ops,
    const std::vector<int>& active
) {
    std::vector<RowOp> mapped;
    mapped.reserve(local_ops.size());
    for (const auto& op : local_ops) mapped.emplace_back(active[op.first], active[op.second]);
    return mapped;
}

inline bool segment_equals(
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

inline std::vector<RowOp> replace_segment(
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

} // namespace linear_surface
