cmake_minimum_required(VERSION 3.16)
project(sbox_surface_code LANGUAGES CXX)
if(NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$"
   OR CMAKE_CXX_SIMULATE_ID STREQUAL "MSVC")
  message(FATAL_ERROR
    "The C++ kernels require GCC or GNU-driver Clang; MSVC and clang-cl are unsupported. "
    "Reference toolchain: GCC 11.4/libstdc++ on Linux x86-64.")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_executable(rewrite sbox_compile/kernels/rewrite.cpp)
add_executable(fresh sbox_compile/kernels/fresh_main.cpp)
find_package(Threads REQUIRED)
target_link_libraries(fresh PRIVATE Threads::Threads)
foreach(target rewrite fresh)
  target_compile_options(${target} PRIVATE -O3 -Wall -Wextra)
endforeach()
