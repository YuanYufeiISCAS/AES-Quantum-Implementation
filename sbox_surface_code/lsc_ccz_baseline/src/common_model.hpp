#pragma once
#include <nlohmann/json.hpp>
#include <ostream>
namespace lsqecc {
struct DenseSlice;
struct LSInstruction;
struct Router;
struct InstructionApplicationResult;
InstructionApplicationResult apply_common_model_instruction(DenseSlice&, LSInstruction&, Router&);
nlohmann::json compile_ccz_model(const nlohmann::json&, std::ostream&);
}
