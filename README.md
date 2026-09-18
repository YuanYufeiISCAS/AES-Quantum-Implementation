# AES Quantum Implementation

This is the code repository for **EUROCRYPT 2027 submission #338**.

Circuit synthesis, surface-code compilation, and independent verification for
AES components and related binary linear layers. The repository contains
logical CNOT circuits, a routed linear-layer optimizer, and a joint AES S-box
compiler. Each module includes concrete circuit data and an independent verifier.

## Project structure

| Module | Contents | Main result or coverage |
| --- | --- | --- |
| [Logical linear layers](linear_layer_logical_depth/README.md) | Target matrices, layered CNOT circuits, output permutations, and a Python verifier | 30 linear maps on 16 or 32 wires; AES MixColumns uses 116 CNOTs at logical depth 9 |
| [Surface-code linear layers](linear_layer_surface_code/README.md) | C++ matrix synthesis, local circuit rewriting, routing, and Python verification | AES MixColumns: 103 CNOTs in 10 VDP layers, or 20 cycles, on a 9 × 17 grid |
| [Surface-code S-boxes](sbox_surface_code/README.md) | Joint optimization of linear blocks, nonlinear operations, placements, and schedules | In-place and Liao–Luo C* S-boxes on 1–5-row layouts; saved circuits, recipe replay, and fresh search |

The logical-depth module covers 12 cipher linear layers and 18 binary block
matrices. The two surface-code modules provide optimization algorithms as well
as saved results for direct inspection and verification.

## Requirements

- Python 3.10 or newer. Python code uses only the standard library.
- CMake 3.16 or newer and a C++17 compiler for the optimization backends.
- A POSIX shell for the command examples below.

Verification of bundled results runs directly in Python, without building the
C++ backends. Build the backends before starting matrix synthesis or fresh
search. Each module's README describes its toolchain and reproducibility settings.

## Verify the bundled results

Run these commands from the repository root. The subshells keep each module's
working directory separate:

```sh
python3 -B linear_layer_logical_depth/verify.py
(cd linear_layer_surface_code && python3 -B -m linear_surface verify)
(cd sbox_surface_code && python3 -B -m sbox_compile verify)
```

- The logical verifier checks all 30 target maps, output permutations, CNOT
  counts, disjoint-endpoint layers, and dependency-based depths.
- The routed linear-layer verifier checks the AES matrix, placement, directed
  paths, data-site obstacles, vertex disjointness, and cycle count.
- The S-box verifier checks complete physical circuits, footprints, causal
  dependencies, conditional corrections, resets, semantics, and reported costs.
  It checks ten main cases and the additional three-row placement-cover circuit.

Successful verification returns exit status 0; invalid circuits return a
nonzero status. The module READMEs describe commands for individual circuits,
machine-readable reports, and detailed inspection.

## Build the optimization backends

From the repository root:

```sh
cmake -S linear_layer_surface_code -B linear_layer_surface_code/build -DCMAKE_BUILD_TYPE=Release
cmake --build linear_layer_surface_code/build -j 2

cmake -S sbox_surface_code -B sbox_surface_code/build -DCMAKE_BUILD_TYPE=Release
cmake --build sbox_surface_code/build -j 2
```

### AES MixColumns synthesis and routing

Run from `linear_layer_surface_code/`:

```sh
python3 -B -m linear_surface reproduce --output runs/reproduce
python3 -B -m linear_surface verify runs/reproduce/aes.json
python3 -B -m linear_surface inspect --paths
python3 -B -m linear_surface search --seed 1 --restarts 8 --output runs/search
```

`reproduce` starts from the matrix and follows a fixed synthesis/refinement
recipe. `search` explores new candidates with the supplied seed and restart
budget. Both save complete routed circuits for independent verification.

### AES S-box compilation

Run from `sbox_surface_code/`:

```sh
python3 -B -m sbox_compile replay --case cstar_4row --output output/replay
python3 -B -m sbox_compile verify output/replay/circuit.json.gz
python3 -B -m sbox_compile search --case cstar_4row --seed 2026914746 --iterations 32 --budget 100 --output output/search
```

`replay` reconstructs the supplied transformations from a parent circuit.
`search` performs a new bounded optimization. Case names are `inplace_1row`
through `inplace_5row` and `cstar_1row` through `cstar_5row`.

Use a new output directory for every reproduction, replay, or search run.
Fresh searches can produce different scores from the saved results; verify
each generated circuit using its module's verifier.

## Tests

After building the backends, run from the repository root:

```sh
python3 -B linear_layer_logical_depth/verify.py --self-test
ctest --test-dir linear_layer_surface_code/build --output-on-failure
(cd sbox_surface_code && python3 -B -m unittest discover -s tests -v)
```

The suites include matrix and circuit checks, malformed-circuit rejection,
replay regressions, and bounded search tests. The logical-layer self-tests can
also run without a C++ build.

## Cost conventions

- **Logical depth** counts endpoint-disjoint CNOT layers before geometric
  routing. Output permutations are explicit in the circuit data.
- **VDP routed layers** use vertex-disjoint paths on the specified patch grid.
  In the linear-layer module, each VDP layer costs two cycles.
- **S-box latency** uses the strict native ledger documented in the S-box
  module, including linear routing, CCZ rounds, H batches, corrections,
  measurements, and resets. It is a complete-module cost, distinct from CNOT
  depth alone.

Compare results under the same target map, placement, output convention, and
cost model. A verified circuit establishes the reported cost for that circuit;
global optimality across different decompositions is a separate question.
