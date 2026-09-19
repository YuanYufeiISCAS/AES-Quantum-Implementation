cmake_minimum_required(VERSION 3.16)
project(lsc_ccz_baseline LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(LSC ${CMAKE_CURRENT_SOURCE_DIR}/vendor/liblsqecc)
set(SOURCES
  scheduler/wave_scheduler.cpp dag/directed_graph.cpp dag/domain_dags.cpp
  gates/gates.cpp gates/parse_gates.cpp gates/gate_approximator.cpp
  gates/decompose_rotation_stream.cpp
  patches/dense_slice.cpp patches/dense_patch_computation.cpp
  patches/patches.cpp patches/sparse_slice.cpp
  layout/layout.cpp layout/router.cpp layout/ascii_layout_spec.cpp
  layout/graph_search/custom_graph_search.cpp
  layout/graph_search/boost_based_graph_search.cpp
  ls_instructions/ls_instructions.cpp ls_instructions/local_instructions.cpp
  ls_instructions/parse_utils.cpp ls_instructions/ls_instructions_parse.cpp
  ls_instructions/ls_instruction_stream.cpp ls_instructions/ls_instructions_from_gates.cpp)
list(TRANSFORM SOURCES PREPEND "${LSC}/src/")
add_library(lsc_ccz STATIC ${SOURCES} src/common_model.cpp)
target_include_directories(lsc_ccz PUBLIC src ${LSC}/include
  ${LSC}/external/include ${LSC}/external/ordered-map/include ${LSC}/external/json/include)
add_executable(lsc_ccz_runner src/runner.cpp)
target_link_libraries(lsc_ccz_runner PRIVATE lsc_ccz)
