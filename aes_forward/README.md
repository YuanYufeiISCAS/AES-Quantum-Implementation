# Fixed Forward AES-128 Circuit

This directory reproduces the paper's adopted forward AES-128 implementation
from its selected circuits and routes. It does not rerun synthesis, placement
optimization, or routing search. The assembler instantiates the fixed modules,
connects their prescribed interfaces, and reconstructs the forward schedule.

## Run

Requirements: Python 3.10+ and the `openssl` executable. Python uses only its
standard library; no C++ build is needed. Keep this directory beside
`sbox_surface_code/`, whose verifier and selected S-box witnesses are reused.
No research-workspace paths or network downloads are required.

From this directory:

```sh
python3 -B reproduce.py --output output/paper
python3 -B reproduce.py --verify output/paper/circuit.json.gz
python3 -B -m unittest discover -s tests -v
```

The first command constructs and verifies the forward circuit. It creates:

- `circuit.json.gz`: both complete S-box witnesses, the fixed linear circuits
  and transfer routes, all module calls, their dependencies and timestamps,
  global patch reservations, and gate counts.
- `verification.json`: machine-readable timing, geometry, gate-count, and
  functional checks.

The second command independently replays that saved circuit. An existing
output directory is never overwritten. Running without `--output` performs
the same construction and checks without saving files. The default functional
check uses three fixed and 256 random plaintext/key pairs; `--random-trials`
changes only this test count, not the circuit. Do not use Python's `-O` option.

## Reproduced results

| Quantity | Forward AES-128 |
|---|---:|
| Latency, in logical cycles | 5,075 |
| Reserved patches, on the 149 × 382 grid | 56,918 |
| Reserved patch cycles | 288,858,850 |
| State S-box calls | 160 |
| Key S-box calls | 40 |
| CCZ states consumed | 10,080 |
| Linear CNOTs | 92,292 |
| CNOTs for CCZ injection | 30,240 |
| Reserved conditional CNOTs | 35,520 |
| Total explicit CNOTs | 158,052 |
| H gates | 49,200 |
| QAND / QAND erasure / ordinary Toffoli | 5,280 / 5,280 / 4,800 |

The state and key S-boxes take 472 and 428 cycles. Each of the four
MixColumns circuits uses 169 CNOTs in 14 routed layers, for 28 cycles;
the four blocks execute in parallel. Physical ShiftRows takes seven cycles
in Rounds 2–9, while the final state placement conversion takes three.
The initial assignment absorbs Round 1 ShiftRows, and the output placement
records Round 10 ShiftRows.

The module dependency graph gives

```text
initial AddRoundKey + Round 1 + Rounds 2–8 + Rounds 9–10
       4           +   504   +   7 × 511  +     990     = 5,075 cycles.
```

The key schedule overlaps the state computation. In particular, the final
key S-boxes use separate regions and may run during Round 9 AddRoundKey.
The global reservation check includes that overlap. The selected 428-cycle
key module is used throughout; the separate 412-cycle refinement is not
substituted into this construction.

Linear CNOTs outside the S-boxes comprise 6,084 in MixColumns, 576 in the
round-local key interfaces, 32 in the final key input copy, 992 in key-word
additions, and 1,408 in AddRoundKey. Movement primitives retain their scheduled
time and patch reservations but are not expanded into additional CNOT counts.
Conditional correction slots remain reserved for every measurement outcome.

The output includes the ciphertext, the final round key, and the copied
RotWord input and its S-box output needed by the final key update. Other
workspace is checked to be zero. The retained final-key work registers are
cleared by uncomputation, which is outside this forward-only reproduction.

## Fixed inputs and verification

`data/components.json.gz` contains the two distinct 52-wire MixColumns
circuits, the state and key placement transfers, the final state adapter,
and the key-interface, key-word, lookahead, and AddRoundKey routes.
Blocks 0 and 2 share one MixColumns circuit; blocks 1 and 3 share the other.
The S-box inputs are `inplace_3row.json.gz` and `cstar_4row.json.gz` under
`../sbox_surface_code/paper_results/`. Their complete contents are embedded
in each generated forward-circuit artifact.

The checks reconstruct each MixColumns target from AES arithmetic and the
prescribed byte permutations, then execute all CNOTs on the 32 input basis
vectors and verify cleanup of the remaining 20 outputs. Transfer checks
include the Bell matchings, boundary directions, retained outputs, and final
restore moves. The existing S-box verifier checks the native primitives,
measurement-dependent corrections, full H footprints, and concurrent traces.

The global checker reserves each S-box's entire module region during its
call and the complete paths of other routed operations during their batches.
It verifies that all reservations fit in the common grid without conflicting
in time. S-box operations within those regions are checked through their
embedded event traces.

A separate bit-sliced simulation executes the native S-box actions, every
linear CNOT, and each state transfer at its physical data sites. It checks
the intermediate round keys, final key, retained work registers, and scratch
cleanup. All 259 ciphertexts are compared with OpenSSL, including

```text
key        000102030405060708090a0b0c0d0e0f
plaintext  00112233445566778899aabbccddeeff
ciphertext 69c4e0d86a7b0430d8cdb78070b4c55a
```

The full-AES input checks are sampled functional tests; the local S-box
checks cover all 256 state inputs and all 65,536 key-module input pairs.

## Model boundary

The result assumes the paper's protected primitive contracts and accepted
CCZ states ready at their specified resource sites before consumption.
Physical Pauli corrections and classical feedback use the same leading-order
timing convention as the paper. Factory production, buffering, delivery,
uncomputation, ciphertext comparison, and diffusion are not included.
The patch count is a reserved logical layout, not a physical-qubit total.

## Code

| File | Purpose |
|---|---|
| `reproduce.py` | Construct or replay the forward artifact |
| `forward_replay/circuit.py`, `schedule.py` | Fixed module composition, timing, and gate accounting |
| `forward_replay/geometry.py` | Global reservations and conflict checks |
| `forward_replay/mixcolumns.py` | Matrix, interface, cleanup, and CNOT-route verification |
| `forward_replay/key_geometry.py`, `transport.py` | Key-register routes and teleportation checks |
| `forward_replay/semantics.py` | Physical-site functional simulation and OpenSSL comparison |
| `forward_replay/verify.py` | Complete verification and comparison with the paper's results |
| `tests/test_forward.py` | Reproduction and malformed-witness regressions |
