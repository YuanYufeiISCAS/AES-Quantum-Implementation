# AES S-box Joint Rewrite

Standalone compilation of the in-place and Liao–Luo C* S-boxes on 1–5-row
surface-code layouts. Python 3.10+, CMake 3.16+, and a C++17 compiler are required;
Python has no third-party runtime dependencies.
Tested with Python 3.13 and GCC 11.4; identical searches use the same toolchain.

The full search combines three-block paid boundary changes
`M'_j = P_j M_j P_(j-1)^-1`, seeded/free-output GL(4,2) rewriting, independent
fresh synthesis, joint injection/correction portfolios, and a nonmonotone
four-state beam. Every admitted circuit is completely lowered and independently
checked. Restoring H footprints, measurement-conditioned corrections, and resets
remain explicit.

## Run

```sh
cmake -S . -B build
cmake --build build -j2
python3 -B -m sbox_compile verify
python3 -B -m sbox_compile replay --case cstar_4row --output output/replay
python3 -B -m sbox_compile search --case cstar_4row --seed 2026914746 --iterations 32 --budget 100 --output output/search
python3 -B -m unittest discover -s tests -v
```

`verify` checks stored circuits; `replay` applies the supplied recipe to its
parent circuit; `search` runs a new search without loading `results/` or
`recipes/`. Outputs must use a new directory.
Search exposes only the seed, iteration count, and deterministic work multiplier.
Exact saved-result reconstruction uses `replay`; a fresh search uses finite work
budgets and can find a different score. The C* four-row replay trajectory is
`465 → 463 → 461 → 459`.

## Saved results

Strict native latency is
`D_CNOT + 5*R_CCZ + 3*H_batches + 2*correction_batches + measurement + reset + port_reset`.
`N_H` counts individual H gates, not H batches. The table below uses this strict
native ledger; comparison with paper-level totals requires matching the
accounting conventions.

| Data rows | In-place latency | C* latency |
|---|---:|---:|
| 1 | 576 | 506 |
| 2 | 539 | 475 |
| 3 | 505 | 471 |
| 4 | 546 | 459 |
| 5 | 546 | 473 |

The ten main results cover both S-box families on all five layouts. The
additional placement-cover circuit `results/inplace_3row_cover.json.gz` has
latency 505, 256 H gates, and 417 linear CNOTs; the main in-place three-row
circuit has 419 linear CNOTs. This additional circuit is checked by `verify`.

## Files

All `.json.gz` files are ordinary gzip-compressed JSON containing complete
physical witnesses. The package is self-contained.

| File(s) | Purpose |
|---|---|
| `CMakeLists.txt`, `.gitignore` | Build the two kernels; exclude generated files. |
| `sbox_compile/__init__.py`, `__main__.py` | Package version and three CLI commands. |
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
| `sbox_compile/verify.py`, `resources.py`, `dependencies.py`, `timeline.py` | Independent interfaces, footprints, causal dependencies, and execution checks. |
| `sbox_compile/primitives.py`, `semantics.py` | Exact local branch identities and exhaustive coherent AES checks. |
| `sbox_compile/replay.py` | Reconstruct saved transformations without reading the answer. |
| `parents/*.json.gz` | Complete starting circuits for all ten cases. |
| `recipes/*.json.gz` | Transformations, schedules, seeds, and parameters for replay. |
| `results/*.json.gz` | Complete verified circuits and their cost index. |
| `tests/test_submission.py` | Circuit/replay regressions, independent fresh synthesis, and corruption tests. |
