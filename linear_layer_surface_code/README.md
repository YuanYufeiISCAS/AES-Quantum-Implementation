# Linear-layer surface-code synthesis

AES MixColumns: **103 CNOTs, 10 VDP layers, 20 cycles** on a fixed 9 × 17 grid with 32 data sites and free output permutation. Each VDP layer costs two cycles. Wire indices and path coordinates are zero-based; the certificate implements `final[i] = target[output_permutation[i]]`.

The optimizer combines matrix beam synthesis, exact 2–4-wire rewrites, free-output GL(4,2) rewrites, commuting-DAG routing, and bounded stochastic restarts. Complete routed circuits are ranked by `(cycles, layers, CNOTs)`.

## Files

| File | Purpose |
| --- | --- |
| `cpp/algebra.hpp` | GF(2) operations, CNOT cancellation, precedence lower bounds, local-window helpers. |
| `cpp/circuit.hpp` | Fixed layout, circuit representation, C++ correctness checks. |
| `cpp/routing.hpp` | Directed paths, bounded joint repacking, commuting-DAG scheduling. |
| `cpp/synthesis.hpp` | Matrix beam search and diverse complete-candidate selection. |
| `cpp/gl4.hpp`, `cpp/free_output.hpp` | Exact four-wire synthesis and output-permutation rewrites. |
| `cpp/rewrite.hpp`, `cpp/policy.hpp` | Local rewriting, stochastic walk, fixed search/routing policies. |
| `cpp/protocol.hpp`, `cpp/main.cpp` | Standalone C++ backend and private input/output protocol. |
| `linear_surface/io.py`, `__main__.py`, `__init__.py` | JSON I/O and Python command-line package. |
| `linear_surface/reproduce.py`, `search.py` | Fixed-recipe reproduction and fresh population search. |
| `linear_surface/verify.py` | Independent Python matrix/path verifier. |
| `instances/aes.json` | AES matrix and physical layout. |
| `recipes/aes_sota.json` | Synthesis seed and six refinement stages with finite work budgets. |
| `circuits/aes.json` | Complete SOTA circuit: gates, routed paths, permutation, and costs. |
| `tests/test_kernels.cpp`, `tests/test_verify.py`, `tests/test_cli.py` | Kernel, corruption, and deterministic-search tests. |
| `CMakeLists.txt` | Build and test setup. |

## Build and reproduce

Requirements: CMake ≥ 3.16, Python ≥ 3.10, and GCC or GNU-driver Clang with
C++17 support. Python uses only the standard library. The tested reference
environment is Linux x86-64, GCC 11.4/libstdc++, CMake 3.22.1, and Python 3.13.9.
Clang and other operating systems are unverified; CMake rejects MSVC and
clang-cl. Exact seeded reproduction uses the reference compiler and standard
library because random/shuffle ordering is toolchain-dependent.

Run from this directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 2
ctest --test-dir build --output-on-failure
python -m linear_surface reproduce --output runs/reproduce
```

`reproduce` starts from the matrix and follows the fixed recipe; saved circuits are not inputs. The reference run takes about five minutes on one CPU core. It saves the final circuit to `runs/reproduce/aes.json`, intermediate circuits to `stage_00.json`–`stage_06.json`, and stage costs/work/time to `run.json` in the same directory. Choose a new output directory for each run.

## View and verify results

```sh
python -m linear_surface verify
python -m linear_surface inspect
python -m linear_surface inspect --paths
python -m linear_surface verify runs/reproduce/aes.json
```

The first three commands use `circuits/aes.json`; viewing and verification do not require a C++ build. `verify` checks the target matrix, output permutation, fixed layout, path orientation, data-site obstacles, vertex disjointness, and reported costs.

Fresh search uses the same core algorithm with new seeds and a pool of up to eight verified circuits:

```sh
python -m linear_surface search --seed 1 --restarts 8 --output runs/search
```

Its best circuit is `runs/search/best.json`; a fresh finite search can produce a different score from the fixed SOTA recipe.
