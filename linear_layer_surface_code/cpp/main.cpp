#include "synthesis.hpp"
#include "rewrite.hpp"

using namespace linear_surface;

// Private, whitespace-delimited protocol. The Python entry point validates
// the public JSON schema before starting this process.
inline RefinementInput read_seed(std::istream& input, bool free_output) {
    std::ostringstream matrix_text;
    for (int line = 0; line < 4; ++line) {
        std::string text;
        if (!std::getline(input, text)) throw std::runtime_error("missing matrix line");
        matrix_text << text << '\n';
    }
    std::istringstream matrix_input(matrix_text.str());
    RefinementInput request;
    request.problem = read_problem(matrix_input);
    request.free_output = free_output;
    if (request.problem.fixed_output == free_output)
        throw std::runtime_error("seed output policy disagrees with stage policy");
    auto& seed = request.incumbent;
    seed.n = request.problem.layout.n;
    seed.mode = "vdp";
    seed.layout = request.problem.layout;
    seed.target_rows = request.problem.target;
    seed.output_permutation = request.problem.initial_permutation;
    seed.final_rows = identity_matrix(seed.n);
    int count = -1;
    if (!(input >> count) || count < 0 || count > 4096)
        throw std::runtime_error("invalid layer count");
    for (int index = 0; index < count; ++index) {
        int gates = -1;
        if (!(input >> gates) || gates < 1 || gates > seed.n / 2)
            throw std::runtime_error("invalid layer size");
        Layer layer{"vdp", {}, 2};
        for (int i = 0; i < gates; ++i) {
            RoutedCNOT gate;
            int length = -1;
            if (!(input >> gate.control >> gate.target >> length)
                || gate.control < 0 || gate.control >= seed.n
                || gate.target < 0 || gate.target >= seed.n
                || gate.control == gate.target || length < 2
                || length > seed.layout.grid_rows * seed.layout.grid_cols)
                throw std::runtime_error("invalid gate or path length");
            for (int j = 0; j < length; ++j) {
                int vertex = -1;
                if (!(input >> vertex) || vertex < 0
                    || vertex >= seed.layout.grid_rows * seed.layout.grid_cols)
                    throw std::runtime_error("invalid path vertex");
                gate.path.push_back(vertex);
            }
            row_add(seed.final_rows, gate.control, gate.target);
            layer.ops.push_back(std::move(gate));
        }
        seed.layers.push_back(std::move(layer));
    }
    std::string extra;
    if (input >> extra) throw std::runtime_error("extra input after circuit");
    if (free_output && !infer_output_permutation(seed.target_rows, seed.final_rows, seed.output_permutation))
        throw std::runtime_error("seed does not implement target");
    verify_result(request.problem, seed);
    return request;
}

int main(int argc, char** argv) {
    try {
        if (argc < 2) throw std::runtime_error("use python -m linear_surface --help");
        const std::string command = argv[1];
        if (command == "synth" && argc == 3) {
            routing_policy = synthesis_routing;
            const Problem problem = read_problem(std::cin);
            const auto result = synthesize_matrix(problem.target, problem.layout,
                                                   parse_unsigned(argv[2], "seed"));
            verify_result(problem, result);
            std::cout << "{\"kind\":\"candidate\",\"result\":";
            write_json_to_stream(result, std::cout);
            std::cout << "}\n";
        } else if (command == "refine" && argc == 9) {
            RefinementBudget options;
            options.seed = parse_unsigned(argv[2], "seed");
            options.max_proposals = parse_integer(argv[3], "proposals", 1);
            options.max_routes = parse_integer(argv[4], "routes", 1);
            options.rounds = parse_integer(argv[5], "rounds", 1);
            options.max_growth = parse_integer(argv[6], "growth");
            options.walk_slack = parse_integer(argv[7], "slack");
            const int policy = parse_integer(argv[8], "free");
            if (options.max_proposals > 262144 || options.max_routes > 4096
                || options.rounds > 4096 || options.max_growth > 8
                || options.walk_slack > 8 || policy > 1)
                throw std::runtime_error("refinement budget exceeds supported bounds");
            const bool free_output = policy != 0;
            routing_policy = refinement_routing;
            const auto input = read_seed(std::cin, free_output);
            Refiner(input, options).run();
        } else {
            throw std::runtime_error("invalid backend command");
        }
        if (!std::cout) throw std::runtime_error("output write failed");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "linear-surface: " << error.what() << '\n';
        return 1;
    }
}
