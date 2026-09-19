# Adapted LSC-CCZ Baseline

This directory reproduces the LSC-CCZ comparison using an adapted Lattice
Surgery Compiler backend. It is not a result from unmodified LSC: the adapter
adds the paper's CNOT, restoring H, Bell, CCZ-supply, correction, and reset
contracts to LSC's instruction dispatch. Routing uses LSC's Dijkstra search,
and ready instructions are selected by its wave scheduler. Independent
primitives may overlap within and between modules.

## Reproduce

Requirements: Python 3.10+, CMake 3.16+, a C++20 compiler, and the `openssl`
executable. Python uses only its standard library. The tested build uses
Linux x86-64, GCC 11.4/libstdc++, and CMake 3.22.1. All required C++ headers
are vendored; neither Boost, ProjectQ, nor a network download is needed.

Run from this directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
python3 -B -m unittest discover -s tests -v
python3 -B run.py reproduce --output output/reproduce
python3 -B run.py verify --output output/reproduce
```

`reproduce` lowers the bundled source gates, compiles both directions from
scratch, and independently checks the resulting traces. It does not load
precomputed routes or schedules. Use `--direction forward` or `--direction
inverse` to run one direction. Every reproduction requires a new output
directory; each direction writes about 100 MB of input and trace data.
`verify` rebuilds the input and rechecks an existing trace without compiling
again. Do not use Python's `-O` option.

| Direction | Logical cycles | CCZ states | Reserved patches |
|---|---:|---:|---:|
| AES evaluation | 7,159 | 10,080 | 56,918 |
| Uncomputation | 7,009 | 10,080 | 56,918 |
| Sequential total | 14,168 | 20,160 | 56,918 |

These totals cover AES evaluation and uncomputation, not the full Grover
oracle or diffusion operator. Both directions reserve the same rectangle;
their areas are not added. `expected.json` records the comparison targets,
including operation counts and patch-cycle volumes. It is used only after
compilation and independent verification.

## Source circuit and adaptation

The state and key S-boxes are the original 23-wire in-place and 28-wire C*
circuits of Liao and Luo, *Quantum Circuit Synthesis for AES with Low DW-cost*
(IACR ePrint 2025/1494). The AES assembly, fixed placements, module boundaries,
and explicit external output relocations are supplied by this experiment.
They are not presented as an AES construction supplied by LSC's authors.

MixColumns uses the published 105-CNOT network extracted from
`encryption_and_grover_oracle_ProjectQ/ours_Enc_AES128_IU.py` in the
`Minimal_T-depth_Width` source repository. Only that CNOT network is used;
ProjectQ does not compile or time this baseline. The extracted module gates
are bundled in `inputs/`, so the source repositories and ProjectQ are not
runtime dependencies. The extraction versions were
`16acfebc86b0ae78bfd3775600988eb6622b005f` for the S-box sources and
`1e7892f77d1acc94a8b1b4826af50dbfa6cf3a0c` for MixColumns.

`src/lower.py` emits the CCZ-consuming implementations and measured QAND
erasure. All conditional corrections retain their predicates and reserved
execution slots. Adjacent Hadamards are simplified only with the corresponding
physical Pauli conjugation. No optimized paper placements, paths, or timing
traces are imported. The adapter reserves complete CNOT paths for two cycles,
complete restoring H regions for three, Bell operations for three, and explicit
measurements/resets for one. Physical Pauli pulses and classical feedback use
the same leading-order timing convention as the paper.

Each accepted CCZ state must be available at its assigned sites before use.
Its sites remain occupied until consumption, and consumed resource sites and
injection ports are explicitly reset. Factory production, buffering, and
delivery are excluded under the paper's conditional supply assumption.

## Verification and files

`src/audit.py` reconstructs footprints and dependencies without importing the
router or scheduler. It checks full event coverage, fixed boundary directions,
spatial conflicts, feedback readiness, live resource states, resets, and the
reported makespan. `src/quantum.py` checks branch maps for each distinct gadget
and verifies every emitted occurrence against those templates; it also
recovers the supplied source gate sequence.

`src/source_verify.py` checks QAND preconditions, ciphertexts, final keys,
lookahead outputs, scratch cleanup, and inverse restoration on 259 inputs.
Ciphertexts are cross-checked with OpenSSL, including a standard known-answer
example. These are sampled full-AES functional tests, not exhaustive
verification over all plaintext/key pairs. They complement the local quantum
branch checks and the independent geometric replay.

| File | Purpose |
|---|---|
| `inputs/aes_forward.json.gz`, `aes_inverse.json.gz` | Complete source modules, data placement, and interfaces |
| `src/lower.py`, `inputs.py` | Source-preserving nonlinear lowering and deterministic resource-site allocation |
| `src/common_model.cpp`, `common_model.hpp`, `runner.cpp` | Primitive adapter and compiler entry point |
| `vendor/liblsqecc/` | Required LSC sources and bundled dependency headers |
| `src/audit.py`, `quantum.py`, `source_verify.py` | Timing, branch, and source-circuit checks |
| `run.py`, `expected.json` | Reproduction entry point and independently checked target costs |
| `tests/`, remaining `inputs/*.json` | Primitive tests and deliberately corrupted-input tests |

The vendored LSC version, modified dispatch files, and licenses are identified
in `THIRD_PARTY_NOTICES.md`.
