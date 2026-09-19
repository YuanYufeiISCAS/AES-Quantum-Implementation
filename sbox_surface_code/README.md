# AES S-box Compilation and Concurrent Scheduling

Standalone compilation of the in-place and Liao–Luo C* S-boxes on 1–5-row
surface-code layouts. Python 3.10+ is required; it has no third-party runtime
dependencies. Building the optimization kernels requires CMake 3.16+ and GCC
or GNU-driver Clang with C++17 support.

The tested reference environment is Linux x86-64, GCC 11.4/libstdc++,
CMake 3.22.1, and Python 3.13.9. Clang and other operating systems are unverified;
CMake rejects MSVC and clang-cl. Exact seeded search reproduction uses the
reference compiler and standard library. Verification and recipe replay run
directly in Python without building the kernels.

The full search combines three-block paid boundary changes
`M'_j = P_j M_j P_(j-1)^-1`, seeded/free-output GL(4,2) rewriting, independent
fresh synthesis, joint injection/correction portfolios, and a nonmonotone
four-state beam. Every admitted circuit is completely lowered and independently
checked. Restoring H footprints, measurement-conditioned corrections, and resets
remain explicit.

The final scheduler allows independent H, CNOT, measurement, and reset
operations to run concurrently. It preserves each data wire's operation order,
measurement dependencies, complete physical footprints, and live injection
resources. It changes neither gates nor routes. The original synthesis and
recipe replay remain available as the preceding compilation stage.

## Run

```sh
cmake -S . -B build
cmake --build build -j2
python3 -B -m sbox_compile verify
python3 -B -m sbox_compile verify --paper
python3 -B -m sbox_compile reproduce --output output/paper
python3 -B -m sbox_compile schedule --paper --case inplace_3row --starts 2048 --output output/state_schedule
python3 -B -m sbox_compile replay --case cstar_4row --output output/replay
python3 -B -m sbox_compile search --case cstar_4row --seed 2026914746 --iterations 32 --budget 100 --output output/search
python3 -B -m unittest discover -s tests -v
```

`verify --paper` checks the final schedules used in the paper. `reproduce`
reconstructs their timestamps from the recorded resource-use order and checks
them again; it does not rerun the search. `schedule` runs a fresh scheduling
portfolio, with 2,048 deterministic starts for each of two variants: preserving
linear CNOT layers or scheduling their CNOTs individually. It uses only the
source circuit, not the saved final timestamps. Both variants retain the
original directed routes and full H footprints.

Without `--paper`, `verify` checks the original synthesis witnesses, `replay`
applies a supplied recipe to its parent, and `search` performs a new synthesis
search without loading `results/` or `recipes/`. To schedule a new circuit, use
`schedule output/search/circuit.json.gz --output output/search_scheduled`.
All output directories must be new. Do not use Python's `-O` option: some
verification routines use assertions.
Search exposes only the seed, iteration count, and deterministic work multiplier.
Exact saved-result reconstruction uses `replay`; a fresh search uses finite work
budgets and can find a different score. The C* four-row replay trajectory is
`465 → 463 → 461 → 459`.

## Paper results

For a complete event trace, latency is `max(event.start + event.duration)`.
Adding durations of different operation types would count concurrent work more
than once. CNOTs take two cycles, restoring H gates and resource Bell operations
take three, and explicit measurement/reset operations take one. Physical Pauli
corrections and classical feed-forward retain the original leading-order
timing contract. A one-cycle preparation of injection ports is included at the
module boundary. Accepted CCZ states must be ready at their specified resource
sites; factories and delivery are outside the module cost.

`paper_results/` contains the ten integrated reference circuits and the
four-row C* refinement. Each artifact embeds its complete source circuit and
final trace. `N_H` counts gates, not batches; the `N_*` operation counts are
unchanged by scheduling.

| Data rows | In-place latency | C* latency |
|---|---:|---:|
| 1 | 495 | 459 |
| 2 | 478 | 440 |
| 3 | 472 | 447 |
| 4 | 512 | 428 |
| 5 | 520 | 439 |

The integrated state circuit has 431 linear CNOTs and 256 H gates; the
integrated four-row C* circuit has 356 linear CNOTs and 206 H gates.
`paper_results/cstar_4row_refinement.json.gz` instead has 349 linear CNOTs,
206 H gates, and latency 412. It is identified as a refinement, not substituted
for the integrated circuit in the table.

The `parents/`, `recipes/`, and `results/` directories preserve the synthesis
experiments and exact recipe replay, including the placement-cover witness.
Those circuits are not all the same as the integrated paper references.
Their stage-based costs remain synthesis metadata, not the final concurrent
latencies. Use `--paper` to select the integrated reference for a fresh schedule.

## Adapted LSC-CCZ baseline

`lsc_ccz_baseline/` contains the explicit nonlinear lowering, adapted LSC
backend, independent trace audit, branch checks, and source inputs. It
reproduces 7,159 cycles for AES evaluation and 7,009 for uncomputation.
See its [README](lsc_ccz_baseline/README.md) for build commands, the adaptation
boundary, and source-circuit provenance. This baseline additionally requires
C++20 and OpenSSL for the independent AES tests.

## Files

All `.json.gz` files are ordinary gzip-compressed JSON containing complete
physical witnesses. The package is self-contained.

Each bundled recipe contains a case identifier, seed, explicit transformation
steps, and expected costs. Replay uses the steps directly; fresh search takes
its seed, iteration count, and work budget from the command line.

| File(s) | Purpose |
|---|---|
| `CMakeLists.txt` | Build the two kernels. |
| `sbox_compile/__init__.py`, `__main__.py` | Package version and command-line entry points. |
| `sbox_compile/model.py` | Circuit I/O, exact state identity, and objective. |
| `sbox_compile/search.py` | Full outer Joint Rewrite algorithm and fixed policy constants. |
| `sbox_compile/windows.py`, `matrix.py` | Paid matrix boundaries, reconstruction, and binary algebra. |
| `sbox_compile/linear.py` | Kernel interfaces, matrix/route checks, and candidate portfolios. |
| `sbox_compile/kernels/rewrite.cpp`, `gl4.hpp` | Routed CNOT walk and exact 20,160-state GL(4,2) table. |
| `sbox_compile/kernels/fresh.cpp`, `fresh_main.cpp` | Independent matrix synthesis and fixed-profile driver. |
| `sbox_compile/nonlinear.py` | Regional joint search and complete-circuit combinations. |
| `sbox_compile/injection.py`, `repair.py` | Port assignment, directed rerouting, and changed-leg repair. |
| `sbox_compile/compiler.py`, `lower.py`, `native.py` | Explicit physical lowering, correction variants, and canonical native reconstruction. |
| `sbox_compile/layout.py`, `geometry.py` | Patch layout, routing, restoring H, and local lowering rules. |
| `sbox_compile/event.py`, `scheduler.py`, `schedule.py` | Feasible fallback, bounded state search, and H-run scheduling. |
| `sbox_compile/cost.py` | Gate counts and schedule-derived latency. |
| `sbox_compile/parallel.py`, `parallel_model.py`, `parallel_scheduler.py` | Final concurrent scheduling, dependency graph, and resource lifetimes. |
| `sbox_compile/parallel_verify.py` | Independent event coverage, occupancy, feedback, and makespan checks. |
| `sbox_compile/verify.py`, `resources.py`, `dependencies.py`, `timeline.py` | Independent interfaces, footprints, causal dependencies, and execution checks. |
| `sbox_compile/primitives.py`, `semantics.py` | Exact local branch identities and exhaustive coherent AES checks. |
| `sbox_compile/replay.py` | Reconstruct saved transformations without reading the answer. |
| `parents/*.json.gz` | Complete starting circuits for all ten cases. |
| `recipes/*.json.gz` | Case identifiers, seeds, transformation steps, schedules, and expected costs. |
| `results/*.json.gz` | Complete verified circuits and their cost index. |
| `paper_results/*.json.gz`, `index.json` | Integrated paper sources, concurrent traces, and costs; the C* refinement is separately identified. |
| `tests/test_submission.py` | Circuit/replay regressions, independent fresh synthesis, and corruption tests. |
| `tests/test_parallel.py` | Paper-schedule replay, fresh scheduling, and malformed-trace rejection. |
| `lsc_ccz_baseline/` | Standalone adapted LSC baseline, inputs, expected costs, and tests. |
