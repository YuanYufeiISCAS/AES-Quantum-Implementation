#include "synthesis.hpp"
#include "rewrite.hpp"

using namespace linear_surface;

static void check(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

static void test_gl4() {
    const auto& table = gl4::table();
    check(table.state_count() == 20160, "GL(4,2) state count");
    for (unsigned code = 0; code < 65536; ++code) {
        const auto matrix = static_cast<gl4::Code>(code);
        if (!table.contains(matrix)) continue;
        const auto word = table.shortest(matrix);
        check(gl4::replay(word) == matrix, "GL4 exact replay");
        check(static_cast<int>(word.size()) == table.distance(matrix), "GL4 exact distance");
        if (code % 137 == 0) {
            for (const auto& variant : gl4::shortest_variants(matrix, code, 8)) {
                check(gl4::replay(variant) == matrix, "GL4 variant replay");
                check(variant.size() == word.size(), "GL4 variant length");
            }
        }
    }
}

static void test_free_output() {
    std::mt19937_64 rng(731);
    for (int trial = 0; trial < 128; ++trial) {
        const std::vector<int> active{1, 3, 4, 6};
        std::vector<RowOp> local;
        for (int i = 0; i < 5; ++i) {
            int c = int(rng() % 4), t = int(rng() % 3);
            if (t >= c) ++t;
            local.emplace_back(c, t);
        }
        std::vector<RowOp> word{{0, 2}, {2, 7}};
        const auto global = map_local_sequence(local, active);
        word.insert(word.end(), global.begin(), global.end());
        word.insert(word.end(), {{6, 0}, {1, 2}, {4, 7}, {0, 3}});
        const Matrix target = apply_ops_to_identity(8, word);
        std::vector<int> p{0, 1, 2, 3, 4, 5, 6, 7};
        gl4::Permutation q{0, 1, 2, 3};
        do {
            const auto replacement = gl4::table().shortest(gl4::left_row_permute(gl4::replay(local), q));
            const auto candidate = free_output::replace_window(word, 2, 5, active, replacement, q, p);
            check(apply_ops_to_identity(8, candidate.word) == permute_rows(target, candidate.output_permutation),
                  "free-output suffix relabeling");
        } while (std::next_permutation(q.begin(), q.end()));
    }
}

static void test_routes_and_reorders() {
    routing_policy = refinement_routing;
    const auto layout = Layout::for_data_shape(8, 2, 4);
    std::mt19937_64 rng(812);
    for (int trial = 0; trial < 20; ++trial) {
        std::vector<RowOp> word;
        for (int i = 0; i < 24; ++i) {
            int c = int(rng() % 8), t = int(rng() % 7);
            if (t >= c) ++t;
            word.emplace_back(c, t);
        }
        const Matrix target = apply_ops_to_identity(8, word);
        check(matmul(target, invert_matrix(target)) == identity_matrix(8), "matrix inverse");
        check(apply_ops_to_identity(8, optimize_cnot_sequence(word)) == target, "commuting cancellation");
        for (int order = 0; order < 64; ++order)
            check(apply_ops_to_identity(8, topological_reorder_ops(word, layout, order)) == target,
                  "topological reorder equivalence");
        Circuit circuit;
        circuit.n = 8;
        circuit.mode = "vdp";
        circuit.layout = layout;
        circuit.target_rows = target;
        circuit.output_permutation = {0, 1, 2, 3, 4, 5, 6, 7};
        circuit.layers = schedule_surface_commuting_ops(word, layout, "vdp", 0);
        circuit.final_rows = identity_matrix(8);
        for (const auto& layer : circuit.layers) apply_layer(circuit.final_rows, layer.ops);
        check(verify_geometry_and_stats(circuit), "routed word verification");
    }
}

int main() {
    try {
        test_gl4();
        test_free_output();
        test_routes_and_reorders();
        std::cout << "PASS: exhaustive GL4, free-output rewrites, algebra and VDP routing\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
