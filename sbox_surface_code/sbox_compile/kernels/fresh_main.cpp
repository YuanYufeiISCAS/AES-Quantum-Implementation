// Fixed native-interface front end for the independent synthesis portfolio.
#include "fresh.cpp"

static void emit_result(const Result& r) {
    if (!verify_geometry_and_stats(r)) throw std::runtime_error("fresh result replay failed");
    auto& out = std::cout;
    out << "{\"n\":" << r.n << ",\"verified\":true,\"layout\":{\"data_rows\":"
        << r.layout.data_rows << ",\"data_cols\":" << r.layout.data_cols
        << ",\"grid_rows\":" << r.layout.grid_rows << ",\"grid_cols\":" << r.layout.grid_cols
        << "},\"output_permutation\":[";
    for (int i = 0; i < r.n; ++i) { if (i) out << ','; out << r.output_permutation[i]; }
    out << "],\"stats\":{\"surface_depth\":" << r.surface_depth()
        << ",\"layers\":" << r.layer_count() << ",\"cnots\":" << r.cnot_count()
        << "},\"layers\":[";
    for (std::size_t i = 0; i < r.layers.size(); ++i) {
        if (i) out << ',';
        const auto& layer = r.layers[i];
        out << "{\"mode\":\"vdp\",\"logical_depth\":2,\"operations\":[";
        for (std::size_t j = 0; j < layer.ops.size(); ++j) {
            if (j) out << ',';
            const auto& op = layer.ops[j];
            out << "{\"control\":" << op.control << ",\"target\":" << op.target << ",\"path\":[";
            for (std::size_t k = 0; k < op.path.size(); ++k) {
                if (k) out << ',';
                out << '[' << r.layout.row(op.path[k]) << ',' << r.layout.col(op.path[k]) << ']';
            }
            out << "]}";
        }
        out << "]}";
    }
    out << "]}";
}

int main() {
    try {
        int n, rows, cols, depth_goal, profile;
        std::uint64_t seed, work;
        if (!(std::cin >> n >> rows >> cols >> seed >> depth_goal >> profile >> work)
            || n < 1 || n > 32 || rows < 1 || rows > 5 || cols < 1 || n > rows * cols
            || depth_goal < 0 || depth_goal > 2048 || profile < 0 || profile > 2
            || work < 1 || work > 100000000)
            throw std::runtime_error("invalid fresh request");
        Matrix target(n);
        for (Row& row : target) {
            if (!(std::cin >> row) || row >= (Row(1) << n)) throw std::runtime_error("invalid matrix row");
        }
        std::string extra;
        if (std::cin >> extra) throw std::runtime_error("trailing fresh request data");
        (void)invert_matrix(target);
        const Layout layout = Layout::for_data_shape(n, rows, cols);
        std::vector<int> identity(n);
        for (int i = 0; i < n; ++i) identity[i] = i;
        fresh_target = target;
        fresh_output = identity;
        fresh_layout = layout;
        auto& p = g_logical_search_params;
        p.max_search_targets = 2;
        p.seeds_per_target = 1;
        p.choice_window = profile == 0 ? 1 : profile == 1 ? 2 : 3;
        p.depth_goal = depth_goal;
        p.cnot_goal = 0;
        p.seed_offset = seed;
        p.use_layer_matrix_rank = profile != 0;
        p.use_layer_matrix_only = profile != 0;
        p.layer_matrix_vdp_hard = profile != 0;
        p.layer_matrix_reuse_state_candidates = profile != 0;
        p.layer_matrix_extra_orders = profile == 0 ? 0 : 2;
        p.layer_matrix_order_variants = profile == 0 ? 1 : 2;
        p.layer_matrix_candidate_multiplier = profile == 0 ? 1 : 2;
        p.layer_matrix_beam_multiplier = profile == 0 ? 1 : 2;
        p.layer_matrix_tail_cap = profile == 0 ? 0 : 16;
        p.layer_matrix_layer_cap_variants = profile == 0 ? 1 : 2;
        p.layer_matrix_nonimproving_slack = profile == 2 ? 1 : 0;
        p.layer_matrix_goal_only = profile == 1;
        p.layer_matrix_endpoint_first = profile == 2;
        p.logical_reschedule_orders = profile == 0 ? 0 : 2;
        auto& routing = g_surface_final_pack_params;
        routing.k_paths = 4;
        routing.path_extra = 2;
        routing.repack_limit = profile == 0 ? 32 : 64;
        routing.state_limit = profile == 0 ? 128 : 256;
        fresh_work_limit = work;
        auto consider = [&](Result r) { export_verified_anytime(r, "complete_method"); };
        try {
            if (profile == 0) {
                try { consider(synthesize_shi_feng_surface(target, "vdp", layout, identity, 2048, &identity)); }
                catch (const std::exception&) { }
            } else {
                try { consider(synthesize_layer_matrix_surface_pool(target, "vdp", layout, 2048, &identity)); }
                catch (const std::exception&) { }
                try {
                    auto logical = synthesize_cached_layer_matrix_logical(target, layout, 2048, &identity);
                    consider(schedule_logical_result_on_surface(logical, "vdp", layout, logical.output_permutation));
                } catch (const std::exception&) { }
            }
            try { consider(synthesize_routed_greedy(target, "vdp", layout, identity, 2048)); }
            catch (const std::exception&) { }
            for (const auto& weights : search_presets()) {
                try { consider(synthesize_once(target, "vdp", layout, identity, weights, 2048)); }
                catch (const std::exception&) { }
            }
        } catch (const FreshWorkLimit&) { }
        std::cout << "{\"candidates\":[";
        for (std::size_t i = 0; i < completed.size(); ++i) {
            if (i) std::cout << ',';
            emit_result(completed[i]);
        }
        std::cout << "],\"work\":" << fresh_work_used << ",\"exhausted\":"
                  << (fresh_work_used >= work ? "true" : "false") << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "fresh synthesis: " << error.what() << '\n';
        return 2;
    }
}
