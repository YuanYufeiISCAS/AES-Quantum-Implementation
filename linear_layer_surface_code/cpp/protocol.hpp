#pragma once
#include "routing.hpp"

namespace linear_surface {
struct Problem {
    Matrix target;
    Layout layout;
    std::string mode;
    bool fixed_output = false;
    std::vector<int> initial_permutation;
};

inline std::vector<std::string> read_tokens(std::istream& input, const std::string& name) {
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("missing input line: " + name);
    std::istringstream stream(line);
    return {std::istream_iterator<std::string>(stream), std::istream_iterator<std::string>()};
}

inline std::uint64_t parse_unsigned(const std::string& token, const std::string& name) {
    if (token.empty() || token.find_first_not_of("0123456789") != std::string::npos) {
        throw std::runtime_error(name + " must be an unsigned decimal integer");
    }
    std::size_t consumed = 0;
    const std::uint64_t value = std::stoull(token, &consumed, 10);
    if (consumed != token.size()) throw std::runtime_error("invalid integer for " + name);
    return value;
}

inline int parse_integer(const std::string& token, const std::string& name, int minimum = 0) {
    const std::uint64_t value = parse_unsigned(token, name);
    if (value < static_cast<std::uint64_t>(minimum)
        || value > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error(name + " is outside the supported integer range");
    }
    return static_cast<int>(value);
}

inline Problem read_problem(std::istream& input) {
    const auto header = read_tokens(input, "header");
    if (header.size() != 4) {
        throw std::runtime_error("header must contain n, data_rows, data_cols, and output policy");
    }
    const int n = parse_integer(header[0], "n", 1);
    if (n > 64) throw std::runtime_error("n must be in [1,64]");
    const int data_rows = parse_integer(header[1], "data_rows", 2);
    const int data_cols = parse_integer(header[2], "data_cols", 2);
    if (static_cast<long long>(data_rows) * data_cols != n) {
        throw std::runtime_error("common layout requires a full data rectangle: n = data_rows * data_cols");
    }
    Problem request;
    request.mode = "vdp";
    if (header[3] != "free" && header[3] != "fixed")
        throw std::runtime_error("output policy must be free or fixed");
    request.fixed_output = header[3] == "fixed";

    request.layout = Layout::for_data_shape(n, data_rows, data_cols);
    const auto positions = read_tokens(input, "data positions");
    if (positions.size() != static_cast<std::size_t>(n)) {
        throw std::runtime_error("data positions must contain exactly n entries");
    }
    for (int wire = 0; wire < n; ++wire) {
        if (parse_integer(positions[static_cast<std::size_t>(wire)], "data position")
            != request.layout.data_pos[static_cast<std::size_t>(wire)]) {
            throw std::runtime_error("common layout requires unchanged full row-major odd/odd data positions");
        }
    }

    const auto matrix = read_tokens(input, "target rows");
    if (matrix.size() != static_cast<std::size_t>(n)) {
        throw std::runtime_error("target matrix must contain exactly n hexadecimal rows");
    }
    for (const std::string& token : matrix) {
        if (token.empty() || token.front() == '-' || token.front() == '+') {
            throw std::runtime_error("target row must be an unsigned hexadecimal integer");
        }
        std::size_t consumed = 0;
        const Row bits = std::stoull(token, &consumed, 16);
        if (consumed != token.size() || (bits & ~mask_for(n)) != 0) {
            throw std::runtime_error("target row is not an n-bit hexadecimal integer");
        }
        request.target.push_back(bits);
    }
    (void)invert_matrix(request.target);

    request.initial_permutation.resize(static_cast<std::size_t>(n));
    for (int wire = 0; wire < n; ++wire) request.initial_permutation[static_cast<std::size_t>(wire)] = wire;
    const auto permutation = read_tokens(input, "output permutation");
    if (request.fixed_output) {
        if (permutation.size() != static_cast<std::size_t>(n)) {
            throw std::runtime_error("fixed output requires exactly n permutation entries");
        }
        for (int wire = 0; wire < n; ++wire) {
            request.initial_permutation[static_cast<std::size_t>(wire)] =
                parse_integer(permutation[static_cast<std::size_t>(wire)], "output permutation entry");
        }
        validate_permutation_vector(request.initial_permutation, n);
    } else if (permutation != std::vector<std::string>{"-"}) {
        throw std::runtime_error("free output requires '-' instead of a fixed permutation");
    }
    std::string extra;
    if (input >> extra) throw std::runtime_error("extra token after matrix request: " + extra);
    return request;
}

inline void verify_result(const Problem& request, const Circuit& result) {
    if (result.n != request.layout.n || result.mode != request.mode
        || result.target_rows != request.target
        || result.layout.data_rows != request.layout.data_rows
        || result.layout.data_cols != request.layout.data_cols
        || result.layout.grid_rows != request.layout.grid_rows
        || result.layout.grid_cols != request.layout.grid_cols
        || result.layout.data_pos != request.layout.data_pos
        || result.layout.occupied_data != request.layout.occupied_data) {
        throw std::runtime_error("synthesis changed the requested matrix, mode or layout");
    }
    if (request.fixed_output && result.output_permutation != request.initial_permutation) {
        throw std::runtime_error("synthesis did not preserve the required output permutation");
    }
    if (!verify_geometry_and_stats(result)) {
        throw std::runtime_error("synthesis result failed matrix replay, route geometry or cost verification");
    }
}

inline void write_json_to_stream(const Circuit& result, std::ostream& out) {
    out << "{\n";
    out << "  \"n\": " << result.n << ",\n";
    out << "  \"mode\": \"" << result.mode << "\",\n";
    out << "  \"layout\": {\"data_rows\": " << result.layout.data_rows
        << ", \"data_cols\": " << result.layout.data_cols
        << ", \"grid_rows\": " << result.layout.grid_rows
        << ", \"grid_cols\": " << result.layout.grid_cols
        << ", \"data_pos\": [";
    for (std::size_t i = 0; i < result.layout.data_pos.size(); ++i) {
        if (i) out << ", ";
        out << "[" << result.layout.row(result.layout.data_pos[i])
            << ", " << result.layout.col(result.layout.data_pos[i]) << "]";
    }
    out << "]},\n";
    out << "  \"stats\": {\"cnots\": " << result.cnot_count()
        << ", \"layers\": " << result.layer_count()
        << ", \"cycles\": " << result.cycles();
    if (result.mode == "logical") {
        const LogicalSequenceBounds bounds = logical_sequence_bounds(flatten_layers(result.layers), result.n);
        out << ", \"logical_lower_bound\": " << bounds.lower_bound
            << ", \"logical_capacity_bound\": " << bounds.capacity_bound
            << ", \"logical_endpoint_bound\": " << bounds.endpoint_bound
            << ", \"logical_critical_bound\": " << bounds.critical_bound;
    }
    out << "},\n";
    out << "  \"output_permutation\": [";
    for (std::size_t i = 0; i < result.output_permutation.size(); ++i) {
        if (i) out << ", ";
        out << result.output_permutation[i];
    }
    out << "],\n";
    out << "  \"target_rows_hex\": [";
    for (std::size_t i = 0; i < result.target_rows.size(); ++i) {
        if (i) out << ", ";
        out << "\"0x" << std::hex << result.target_rows[i] << std::dec << "\"";
    }
    out << "],\n";
    out << "  \"layers\": [\n";
    for (std::size_t li = 0; li < result.layers.size(); ++li) {
        const auto& layer = result.layers[li];
        out << "    {\"mode\": \"" << layer.mode << "\", \"cycles\": "
            << layer.cycles << ", \"operations\": [";
        for (std::size_t oi = 0; oi < layer.ops.size(); ++oi) {
            const auto& op = layer.ops[oi];
            if (oi) out << ", ";
            out << "{\"control\": " << op.control << ", \"target\": " << op.target
                << ", \"path\": [";
            for (std::size_t pi = 0; pi < op.path.size(); ++pi) {
                if (pi) out << ", ";
                out << "[" << result.layout.row(op.path[pi]) << ", " << result.layout.col(op.path[pi]) << "]";
            }
            out << "]}";
        }
        out << "]}";
        if (li + 1 < result.layers.size()) out << ",";
        out << "\n";
    }
    out << "  ]\n";
    out << "}\n";
    if (!out) throw std::runtime_error("could not write JSON output");
}

} // namespace linear_surface
