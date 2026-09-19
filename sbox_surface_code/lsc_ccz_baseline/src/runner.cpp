#include "common_model.hpp"
#include <fstream>
#include <iostream>
int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::runtime_error("usage: lsc_ccz_runner input.json trace.jsonl");
        std::ifstream input(argv[1]);
        if (!input) throw std::runtime_error("cannot open input");
        nlohmann::json spec;
        input >> spec;
        std::ofstream trace(argv[2]);
        if (!trace) throw std::runtime_error("cannot open trace");
        auto result = lsqecc::compile_ccz_model(spec, trace);
        trace.flush();
        if (!trace) throw std::runtime_error("trace output failure");
        std::cout << result.dump(2) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
