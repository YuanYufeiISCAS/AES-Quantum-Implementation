# Logical-depth witnesses for linear layers

This directory contains all 30 circuits in the **This work** column of
Appendix A.1, Table 6: 12 cipher linear layers and 18 binary block matrices.
Every circuit is an in-place CNOT network on 16 or 32 wires, with its explicit
target matrix, layer schedule, and output permutation.

## Run

Use Python 3.10 or newer; all dependencies are in the standard library.
From the repository root:

```bash
python3 linear_layer_logical_depth/verify.py
python3 linear_layer_logical_depth/verify.py --self-test
python3 linear_layer_logical_depth/verify.py --matrix aes --matrix anubis
python3 linear_layer_logical_depth/verify.py --json
```

The script locates data relative to its own file, so it also runs from other
working directories. The default command checks all 30 rows and prints
`PASS: 30 verified, 0 failed.` Exit status is 0 on success, 1 on verification
failure, and 2 for invalid command-line arguments. Checks also run under
`python3 -O`.

To inspect a modified witness without editing the bundled circuits:

```bash
python3 linear_layer_logical_depth/verify.py --circuit /path/to/candidate.json
```

This mode checks the candidate against the matching Table 6 ID, expected
wire/gate/depth counts, and bundled target matrix. It permits different JSON
formatting or a different valid circuit, so the candidate's file hash is omitted.
The bundled reference file's hash is still checked.

## Files and conventions

- `manifest.json`: the 30 table rows in paper order, expected counts/depths,
  citation numbers, group labels, and circuit-integrity hashes.
- `circuits/<id>.json`: the complete target matrix and logical circuit for one
  row, including every gate, its layer, and the output permutation.
- `verify.py`: independent matrix replay, bit-string simulation, depth
  accounting, input validation, and built-in self-tests.

Each circuit has these fields:

| Field | Meaning |
| --- | --- |
| `schema_version` | Format version, currently 1 |
| `id` | Stable matrix identifier from the manifest |
| `wires` | Number of in-place logical wires; indices are zero-based |
| `matrix_rows[i][j]` | Binary coefficient of input wire `j` in canonical output `i` |
| `output_permutation[i]` | Canonical output index carried by physical output wire `i` |
| `layers[k]` | The CNOT pairs `[control, target]` executed in layer `k` |

The **leftmost character** of each matrix row is column/input wire 0.
A gate `[c, t]` performs `wire[t] ^= wire[c]`. Layers execute in array order.
Every wire participates in at most one gate per layer, including gates that
share a control or a target. A CNOT takes one logical-depth unit; ancilla count
is zero. For AES, wire `8*b + k` is bit `k` of byte `b`, least-significant bit
first within each byte.

For column-vector input `x` and target matrix `M`, the physical output satisfies
`y[i] = (M*x)[output_permutation[i]]`. To recover canonical output, assign
`canonical[output_permutation[i]] = y[i]`. This is the paper's free output
relabeling convention; it adds no SWAP gates to the logical metric.

## What the verifier checks

1. Require exactly the 30 manifest IDs, in Table 6 order, with a complete circuit
   inventory. Check the saved circuit bytes against their manifest SHA-256.
2. Validate the JSON schema, binary matrix dimensions and rank, wire indices,
   output bijection, and disjoint endpoints in every nonempty layer.
3. Start from the identity matrix and replay every CNOT by XORing the control
   row into the target row. Require resulting row `i` to equal target row
   `output_permutation[i]`.
4. Independently simulate packed bit strings for zero and every input basis
   vector, undo the output relabeling, and compare with matrix multiplication.
   Since both maps are linear, equality on the basis establishes equality on
   all inputs. CNOT networks introduce no basis-dependent phases, so this also
   determines their action on arbitrary quantum states.
5. Recount gates and layers. Recompute earliest gate levels while preserving
   the order of gates on each wire: `level = 1 + max(last[c], last[t])`, followed
   by `last[c] = last[t] = level`. Both the saved schedule and this ASAP depth
   must equal the paper's depth.

For each row the explicit binary matrix defines the target map.
AES receives an additional check:
the verifier derives its 32-bit matrix from the standard MixColumns byte
coefficients over the polynomial `x^8 + x^4 + x^3 + x + 1` (`0x11b`).

The report includes the maximum wire load (endpoint lower bound) and, in JSON
mode, the gate-capacity bound `ceil(CNOTs / floor(wires/2))`. These are bounds
for scheduling the supplied gates. ASAP depth is exact for the supplied
per-wire gate order. The depth conclusions concern these circuits; global
optimality over alternative CNOT decompositions requires a separate argument.
These measurements precede geometric routing and surface-code cycle accounting.

Self-tests run the full dataset and deliberately corrupt gates, permutations,
matrices, layer endpoints, expected metrics, hashes, and package metadata.
A three-cycle test checks permutation direction; an altered-target test checks
that external witnesses remain tied to the bundled matrix. Corruption tests
call the mathematical verifier directly where appropriate, independently of
the file-hash check. Tests modify in-memory copies only.

## Table 6 coverage

Reference numbers identify the matrix constructions cited in Table 6.
For block-matrix IDs, the final three dimensions specify the block-array shape
and the size of each square binary block: `4x4-8`, for example, is a 4-by-4
array of 8-by-8 binary blocks, giving 32 wires. `-i-` denotes an involutory
variant. The `group` and `label` fields distinguish variants in the manifest.

| Row | Matrix ID | Reference | Wires | CNOTs | Logical depth |
| --- | --- | --- | ---: | ---: | ---: |
| 1 | `aes` | [43] | 32 | 116 | 9 |
| 2 | `anubis` | [5] | 32 | 121 | 11 |
| 3 | `clefia-m0` | [67] | 32 | 125 | 10 |
| 4 | `clefia-m1` | [67] | 32 | 129 | 10 |
| 5 | `fox-mu4` | [45] | 32 | 205 | 20 |
| 6 | `qarma128` | [2] | 32 | 48 | 3 |
| 7 | `twofish` | [65] | 32 | 167 | 14 |
| 8 | `whirlwind-m0` | [4] | 32 | 268 | 23 |
| 9 | `whirlwind-m1` | [4] | 32 | 225 | 21 |
| 10 | `joltik` | [41] | 16 | 48 | 7 |
| 11 | `midori` | [3] | 16 | 24 | 3 |
| 12 | `smallscale-aes` | [14] | 16 | 55 | 9 |
| 13 | `mds-skop15-4x4-4` | [70] | 16 | 60 | 10 |
| 14 | `mds-liusim16-4x4-4` | [57] | 16 | 60 | 10 |
| 15 | `mds-liwang16-4x4-4` | [51] | 16 | 64 | 10 |
| 16 | `mds-ctg16-4x4-4` | [13] | 16 | 54 | 9 |
| 17 | `mds-jpst17-4x4-4` | [42] | 16 | 52 | 8 |
| 18 | `mds-skop15-4x4-8` | [70] | 32 | 110 | 8 |
| 19 | `mds-liusim16-4x4-8` | [57] | 32 | 189 | 19 |
| 20 | `mds-liwang16-4x4-8` | [51] | 32 | 156 | 15 |
| 21 | `mds-ctg16-4x4-8` | [13] | 32 | 165 | 17 |
| 22 | `mds-tosc-sarsye16-4x4-8` | [63] | 32 | 179 | 19 |
| 23 | `mds-jpst17-4x4-8` | [42] | 32 | 107 | 9 |
| 24 | `mds-skop15-i-4x4-8` | [70] | 32 | 104 | 8 |
| 25 | `mds-liwang16-i-4x4-8` | [51] | 32 | 97 | 8 |
| 26 | `mds-tosc-sarsye16-i-4x4-8` | [63] | 32 | 109 | 10 |
| 27 | `mds-jpst17-i-4x4-8` | [42] | 32 | 100 | 9 |
| 28 | `mds-skop15-8x8-4` | [70] | 32 | 272 | 24 |
| 29 | `mds-ss17-8x8-4` | [64] | 32 | 281 | 27 |
| 30 | `mds-skop15-i-8x8-4` | [70] | 32 | 239 | 20 |
