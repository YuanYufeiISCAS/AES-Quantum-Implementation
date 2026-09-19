cmake_minimum_required(VERSION 3.16)
project(linear_layer_surface_code LANGUAGES CXX)

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$"
   OR CMAKE_CXX_SIMULATE_ID STREQUAL "MSVC")
    message(FATAL_ERROR
        "The C++ backend requires GCC or GNU-driver Clang; MSVC and clang-cl are unsupported. "
        "Reference toolchain: GCC 11.4/libstdc++ on Linux x86-64.")
endif()

add_executable(linear_surface_backend cpp/main.cpp)
target_compile_features(linear_surface_backend PRIVATE cxx_std_17)
set_target_properties(linear_surface_backend PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(linear_surface_backend PRIVATE -Wall -Wextra -Wpedantic)

include(CTest)
if(BUILD_TESTING)
    add_executable(test_kernels tests/test_kernels.cpp)
    target_include_directories(test_kernels PRIVATE cpp)
    target_compile_features(test_kernels PRIVATE cxx_std_17)
    set_target_properties(test_kernels PROPERTIES CXX_EXTENSIONS OFF)
    add_test(NAME kernels COMMAND test_kernels)
    find_package(Python3 3.10 COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(NAME certificates COMMAND ${Python3_EXECUTABLE} -m unittest discover -s tests -v)
        set_tests_properties(certificates PROPERTIES
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            ENVIRONMENT "LINEAR_SURFACE_BACKEND=$<TARGET_FILE:linear_surface_backend>")
    endif()
endif()
