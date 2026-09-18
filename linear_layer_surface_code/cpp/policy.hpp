#pragma once
#include <cstdint>

namespace linear_surface {

// Fixed compact-layout algorithm policy. Only random seeds and finite work
// budgets are exposed by the Python interface and reproduction recipe.
struct RoutingPolicy {
    int k_paths = 4;
    int path_extra = 4;
    int repack_limit = 16;
    bool enabled() const { return k_paths > 1 && repack_limit > 1; }
};
inline constexpr RoutingPolicy synthesis_routing{2, 1, 4};
inline constexpr RoutingPolicy refinement_routing{4, 4, 16};

struct SynthesisPolicy {
    static constexpr int beam_width = 32;
    static constexpr int depth_limit = 28;
    static constexpr int candidate_cap = 24;
    static constexpr int archive_capacity = 4;
    static constexpr int route_trials = 12;
};

struct RefinementBudget {
    std::uint64_t seed = 307;
    int max_proposals = 65536;
    int max_routes = 512;
    int rounds = 64;
    int max_growth = 0;
    int walk_slack = 0;
    static constexpr int max_window = 12;
    static constexpr int walk_restart_rounds = 64;
    static constexpr int gl4_variants = 8;
};

} // namespace linear_surface
