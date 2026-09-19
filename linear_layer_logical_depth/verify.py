#!/usr/bin/env python3
"""Independently verify the 30 logical CNOT circuits in Appendix A.1, Table 6.

Python 3.10+; standard library only. Run from any working directory:
  python3 /path/to/linear_layer_logical_depth/verify.py
  python3 /path/to/linear_layer_logical_depth/verify.py --self-test
"""

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
import sys
import unittest


ROOT = Path(__file__).resolve().parent
TABLE6_IDS = """
aes anubis clefia-m0 clefia-m1 fox-mu4 qarma128 twofish whirlwind-m0
whirlwind-m1 joltik midori smallscale-aes
mds-skop15-4x4-4 mds-liusim16-4x4-4 mds-liwang16-4x4-4
mds-ctg16-4x4-4 mds-jpst17-4x4-4
mds-skop15-4x4-8 mds-liusim16-4x4-8 mds-liwang16-4x4-8
mds-ctg16-4x4-8 mds-tosc-sarsye16-4x4-8 mds-jpst17-4x4-8
mds-skop15-i-4x4-8 mds-liwang16-i-4x4-8
mds-tosc-sarsye16-i-4x4-8 mds-jpst17-i-4x4-8
mds-skop15-8x8-4 mds-ss17-8x8-4 mds-skop15-i-8x8-4
""".split()


class VerificationError(ValueError):
    """Invalid data, an incorrect circuit, or a mismatched paper result."""


def require(condition, message):
    # Explicit exceptions keep checks active under python -O.
    if not condition:
        raise VerificationError(message)


def integer(value, lower, upper):
    return type(value) is int and lower <= value <= upper


def keys(value, expected, context):
    require(isinstance(value, dict), f"{context}: expected an object")
    require(set(value) == set(expected), f"{context}: unexpected or missing fields")


def decode_json(raw):
    def unique_pairs(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def invalid_constant(value):
        raise VerificationError(f"nonstandard JSON number: {value}")

    return json.loads(raw, object_pairs_hook=unique_pairs,
                      parse_constant=invalid_constant)


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def validate_manifest(manifest):
    keys(manifest, ("schema_version", "paper", "circuits"), "manifest")
    require(integer(manifest["schema_version"], 2, 2), "manifest schema version: expected 2")
    paper = manifest["paper"]
    keys(paper, ("section", "table", "column"), "paper")
    require(paper["section"] == "Appendix A.1"
            and integer(paper["table"], 6, 6)
            and paper["column"] == "This work", "wrong paper/table/column")
    entries = manifest["circuits"]
    require(isinstance(entries, list) and len(entries) == len(TABLE6_IDS),
            "Table 6 must contain exactly 30 circuits")
    for row, (entry, identifier) in enumerate(zip(entries, TABLE6_IDS), 1):
        keys(entry, ("row", "id", "label", "reference", "group", "file",
                     "wires", "cnots", "logical_depth", "sha256"),
             f"manifest row {row}")
        require(integer(entry["row"], row, row) and entry["id"] == identifier,
                f"row {row}: expected {identifier}")
        require(entry["file"] == f"circuits/{identifier}.json",
                f"{identifier}: invalid circuit filename")
        for field, lower, upper in (("wires", 1, 64), ("cnots", 1, 100000),
                                    ("logical_depth", 1, 100000)):
            require(integer(entry[field], lower, upper),
                    f"{identifier}: invalid expected {field}")
        require(isinstance(entry["sha256"], str)
                and re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]),
                f"{identifier}: invalid circuit hash")
    return entries


def validate_inventory(entries, filenames):
    expected = {entry["file"] for entry in entries}
    actual = set(filenames)
    require(actual == expected,
            f"circuit inventory: missing={sorted(expected - actual)}, "
            f"unexpected={sorted(actual - expected)}")


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def aes_mixcolumns_rows():
    """Derive the 32-bit map from AES byte arithmetic, independently of JSON."""
    coefficients = ((2, 3, 1, 1), (1, 2, 3, 1),
                    (1, 1, 2, 3), (3, 1, 1, 2))

    def multiply(a, b):
        result = 0
        while b:
            if b & 1:
                result ^= a
            a <<= 1
            if a & 0x100:
                a ^= 0x11B  # x^8 + x^4 + x^3 + x + 1
            b >>= 1
        return result

    rows = [0] * 32
    for input_bit in range(32):
        inputs = [((1 << input_bit) >> (8 * byte)) & 255 for byte in range(4)]
        for byte, factors in enumerate(coefficients):
            value = 0
            for factor, x in zip(factors, inputs):
                value ^= multiply(factor, x)
            for bit in range(8):
                if (value >> bit) & 1:
                    rows[8 * byte + bit] |= 1 << input_bit
    return rows


def verify_circuit(data, expected):
    keys(data, ("schema_version", "id", "wires", "matrix_rows",
                "output_permutation", "layers"), "circuit")
    require(integer(data["schema_version"], 1, 1), "circuit schema version")
    require(data["id"] == expected["id"], "circuit/manifest ID mismatch")
    n = data["wires"]
    require(integer(n, 2, 64) and n == expected["wires"], "wrong wire count")
    strings = data["matrix_rows"]
    require(isinstance(strings, list) and len(strings) == n,
            "target matrix must have n rows")
    for i, row in enumerate(strings):
        require(isinstance(row, str) and len(row) == n and set(row) <= {"0", "1"},
                f"target matrix row {i}: expected exactly n binary digits")
    # Character j describes input wire j, including the leftmost character j=0.
    target = [int(row[::-1], 2) for row in strings]
    require(binary_rank(target) == n, "target matrix is singular")
    if data["id"] == "aes":
        require(target == aes_mixcolumns_rows(), "target differs from AES definition")

    permutation = data["output_permutation"]
    require(isinstance(permutation, list) and len(permutation) == n
            and all(integer(p, 0, n - 1) for p in permutation)
            and sorted(permutation) == list(range(n)),
            "output permutation must be a bijection on 0..n-1")
    layers = data["layers"]
    require(isinstance(layers, list) and layers, "expected nonempty layers")
    matrix = [1 << i for i in range(n)]
    last = [0] * n
    wire_loads = [0] * n
    gates = []
    for layer_index, layer in enumerate(layers, 1):
        require(isinstance(layer, list) and layer, f"layer {layer_index}: empty/invalid")
        endpoints = set()
        for gate_index, gate in enumerate(layer, 1):
            where = f"layer {layer_index}, gate {gate_index}"
            require(isinstance(gate, list) and len(gate) == 2,
                    f"{where}: expected [control, target]")
            control, target_wire = gate
            require(integer(control, 0, n - 1) and integer(target_wire, 0, n - 1),
                    f"{where}: invalid wire index")
            require(control != target_wire, f"{where}: control equals target")
            require(control not in endpoints and target_wire not in endpoints,
                    f"{where}: endpoints overlap within the layer")
            endpoints.update((control, target_wire))
            matrix[target_wire] ^= matrix[control]
            level = 1 + max(last[control], last[target_wire])
            last[control] = last[target_wire] = level
            wire_loads[control] += 1
            wire_loads[target_wire] += 1
            gates.append((control, target_wire))

    cnots, depth, asap = len(gates), len(layers), max(last)
    require(cnots == expected["cnots"],
            f"CNOT count {cnots} != paper {expected['cnots']}")
    require(depth == expected["logical_depth"],
            f"schedule depth {depth} != paper {expected['logical_depth']}")
    require(asap == expected["logical_depth"],
            f"fixed-wire-order ASAP depth {asap} != paper {expected['logical_depth']}")
    for i, p in enumerate(permutation):
        require(matrix[i] == target[p],
                f"matrix mismatch at output wire {i}, canonical bit {p}: "
                f"computed 0x{matrix[i]:x}, expected 0x{target[p]:x}")

    # An independent packed-bit simulator checks zero and every basis vector.
    # Linearity makes the basis comparison exhaustive for these CNOT-only maps.
    for x in [0] + [1 << j for j in range(n)]:
        state = x
        for control, target_wire in gates:
            if (state >> control) & 1:
                state ^= 1 << target_wire
        canonical = 0
        for i, p in enumerate(permutation):
            canonical |= ((state >> i) & 1) << p
        reference = 0
        for i, row in enumerate(target):
            reference |= ((row & x).bit_count() % 2) << i
        require(canonical == reference, f"basis replay failed at input 0x{x:x}")

    return {
        "id": data["id"], "wires": n, "cnots": cnots,
        "logical_depth": depth, "asap_depth": asap,
        "endpoint_lower_bound": max(wire_loads),
        "capacity_lower_bound": (cnots + n // 2 - 1) // (n // 2),
        "layer_gate_counts": [len(layer) for layer in layers],
        "wire_loads": wire_loads, "matrix_rank": n,
        "output_permutation_is_identity": permutation == list(range(n)),
        "matrix_equal": True, "basis_vectors_checked": n,
        "zero_input_checked": True, "status": "PASS",
    }


def check_bytes(raw, entry, check_hash=True):
    if check_hash:
        require(sha256(raw) == entry["sha256"], f"{entry['id']}: circuit SHA-256 mismatch")
    return verify_circuit(decode_json(raw), entry)


def check_external_bytes(raw, canonical_raw, entry):
    require(sha256(canonical_raw) == entry["sha256"],
            f"{entry['id']}: bundled circuit SHA-256 mismatch")
    result = check_bytes(raw, entry, check_hash=False)
    require(decode_json(raw)["matrix_rows"] == decode_json(canonical_raw)["matrix_rows"],
            "external target matrix differs from the bundled Table 6 target")
    return result


def load_package():
    entries = validate_manifest(decode_json((ROOT / "manifest.json").read_bytes()))
    validate_inventory(entries, [p.relative_to(ROOT).as_posix()
                                for p in (ROOT / "circuits").rglob("*.json")])
    return entries


def run_self_tests():
    entries = load_package()
    entry = entries[0]
    data = decode_json((ROOT / entry["file"]).read_bytes())
    manifest = decode_json((ROOT / "manifest.json").read_bytes())

    class Tests(unittest.TestCase):
        def setUp(self):
            self.data = copy.deepcopy(data)
            self.expected = copy.deepcopy(entry)

        def reject(self, message):
            with self.assertRaisesRegex(VerificationError, message):
                verify_circuit(self.data, self.expected)

        def test_all_30_original_circuits(self):
            for record in entries:
                with self.subTest(matrix=record["id"]):
                    result = check_bytes((ROOT / record["file"]).read_bytes(), record)
                    self.assertEqual(result["status"], "PASS")

        def test_deleted_gate(self):
            self.data["layers"][0].pop()
            self.reject("CNOT count")

        def test_reversed_gate(self):
            self.data["layers"][0][0].reverse()
            self.reject("matrix mismatch")

        def test_corrupted_permutation(self):
            p = self.data["output_permutation"]
            p[0], p[1] = p[1], p[0]
            self.reject("matrix mismatch")

        def test_nonbijective_permutation(self):
            self.data["output_permutation"][0] = self.data["output_permutation"][1]
            self.reject("bijection")

        def test_boolean_permutation_index(self):
            self.data["output_permutation"][0] = False
            self.reject("bijection")

        def test_endpoint_conflicts(self):
            for pair in ([17, 3], [3, 25], [25, 3], [3, 17]):
                with self.subTest(pair=pair):
                    self.data = copy.deepcopy(data)
                    self.data["layers"][0][1] = pair
                    self.reject("endpoints overlap")

        def test_invalid_wire_indices(self):
            for index in (-1, 32, True, 1.0, "1"):
                with self.subTest(index=index):
                    self.data["layers"][0][0][0] = index
                    self.reject("invalid wire index")

        def test_equal_endpoints(self):
            self.data["layers"][0][0] = [17, 17]
            self.reject("control equals target")

        def test_malformed_gate(self):
            self.data["layers"][0][0] = [17, 25, 1]
            self.reject("expected ")

        def test_empty_layer(self):
            self.data["layers"].append([])
            self.reject("empty/invalid")

        def test_missing_layers(self):
            self.data["layers"] = []
            self.reject("nonempty layers")

        def test_wrong_wire_count(self):
            self.data["wires"] = True
            self.reject("wire count")

        def test_wrong_id(self):
            self.data["id"] = "anubis"
            self.reject("ID mismatch")

        def test_wrong_matrix_dimension(self):
            self.data["matrix_rows"].pop()
            self.reject("must have n rows")

        def test_nonbinary_matrix(self):
            self.data["matrix_rows"][0] = "2" + self.data["matrix_rows"][0][1:]
            self.reject("binary digits")

        def test_singular_matrix(self):
            self.data["matrix_rows"][0] = self.data["matrix_rows"][1]
            self.reject("singular")

        def test_wrong_aes_target(self):
            rows = self.data["matrix_rows"]
            rows[0], rows[1] = rows[1], rows[0]
            self.reject("AES definition")

        def test_wrong_paper_count(self):
            self.expected["cnots"] += 1
            self.reject("CNOT count")

        def test_wrong_paper_depth(self):
            self.expected["logical_depth"] -= 1
            self.reject("schedule depth")

        def test_unexpected_fields(self):
            self.data["verified"] = True
            self.reject("unexpected or missing fields")

        def test_invalid_schema_version(self):
            self.data["schema_version"] = True
            self.reject("schema version")

        def test_hash_mismatch(self):
            raw = (ROOT / entry["file"]).read_bytes() + b" "
            with self.assertRaisesRegex(VerificationError, "SHA-256"):
                check_bytes(raw, entry)

        def test_external_circuit_uses_bundled_target(self):
            record = entries[1]  # ANUBIS: also cover a target other than AES.
            canonical_raw = (ROOT / record["file"]).read_bytes()
            damaged = decode_json(canonical_raw)
            rows = damaged["matrix_rows"]
            rows[0], rows[1] = rows[1], rows[0]
            damaged["output_permutation"] = [
                1 if p == 0 else 0 if p == 1 else p
                for p in damaged["output_permutation"]]
            # The gates still match this altered matrix and altered relabeling.
            self.assertEqual(verify_circuit(damaged, record)["status"], "PASS")
            with self.assertRaisesRegex(VerificationError, "bundled Table 6 target"):
                check_external_bytes(json.dumps(damaged), canonical_raw, record)

        def test_duplicate_json_keys(self):
            with self.assertRaisesRegex(VerificationError, "duplicate JSON key"):
                decode_json('{"id":"aes","id":"anubis"}')

        def test_nonstandard_json_numbers(self):
            for number in ("NaN", "Infinity", "-Infinity"):
                with self.subTest(number=number):
                    with self.assertRaisesRegex(VerificationError, "nonstandard JSON"):
                        decode_json('{"value":' + number + '}')

        def test_manifest_schema_version(self):
            for version in (1, 3, True, 2.0, "2"):
                with self.subTest(version=version):
                    damaged = copy.deepcopy(manifest)
                    damaged["schema_version"] = version
                    with self.assertRaisesRegex(VerificationError, "manifest schema version"):
                        validate_manifest(damaged)

        def test_manifest_field_sets(self):
            for section, fields in (("paper", ("section", "table", "column")),
                                    ("row", ("row", "id", "label", "reference", "group",
                                             "file", "wires", "cnots", "logical_depth", "sha256"))):
                for field in (*fields, "extra"):
                    with self.subTest(section=section, field=field):
                        damaged = copy.deepcopy(manifest)
                        target = damaged["paper"] if section == "paper" else damaged["circuits"][0]
                        if field == "extra":
                            target[field] = "unused"
                        else:
                            del target[field]
                        with self.assertRaisesRegex(VerificationError, "unexpected or missing fields"):
                            validate_manifest(damaged)

        def test_manifest_circuit_hash(self):
            for value in (None, "", "0" * 63, "g" * 64):
                with self.subTest(value=value):
                    damaged = copy.deepcopy(manifest)
                    damaged["circuits"][0]["sha256"] = value
                    with self.assertRaisesRegex(VerificationError, "invalid circuit hash"):
                        validate_manifest(damaged)

        def test_missing_manifest_row(self):
            damaged = copy.deepcopy(manifest)
            damaged["circuits"].pop()
            with self.assertRaisesRegex(VerificationError, "exactly 30"):
                validate_manifest(damaged)

        def test_duplicate_manifest_row(self):
            damaged = copy.deepcopy(manifest)
            damaged["circuits"][1] = damaged["circuits"][0]
            with self.assertRaisesRegex(VerificationError, "row 2"):
                validate_manifest(damaged)

        def test_manifest_path_traversal(self):
            damaged = copy.deepcopy(manifest)
            damaged["circuits"][0]["file"] = "../aes.json"
            with self.assertRaisesRegex(VerificationError, "filename"):
                validate_manifest(damaged)

        def test_missing_or_extra_circuit(self):
            files = [record["file"] for record in entries]
            for damaged in (files[:-1], files + ["circuits/extra.json"]):
                with self.assertRaisesRegex(VerificationError, "inventory"):
                    validate_inventory(entries, damaged)

        def test_three_cycle_permutation_direction(self):
            # Two SWAP decompositions give physical output (x1,x2,x0).
            # A 3-cycle distinguishes a permutation from its inverse.
            toy = {"schema_version": 1, "id": "three-cycle", "wires": 3,
                   "matrix_rows": ["100", "010", "001"],
                   "output_permutation": [1, 2, 0],
                   "layers": [[[0, 1]], [[1, 0]], [[0, 1]],
                              [[1, 2]], [[2, 1]], [[1, 2]]]}
            claim = {"id": "three-cycle", "wires": 3, "cnots": 6, "logical_depth": 6}
            self.assertEqual(verify_circuit(toy, claim)["status"], "PASS")
            toy["output_permutation"] = [2, 0, 1]
            with self.assertRaisesRegex(VerificationError, "matrix mismatch"):
                verify_circuit(toy, claim)

        def test_correct_map_with_redundant_schedule_wait(self):
            toy = {"schema_version": 1, "id": "idle", "wires": 4,
                   "matrix_rows": ["1000", "1100", "0010", "0011"],
                   "output_permutation": [0, 1, 2, 3],
                   "layers": [[[0, 1]], [[2, 3]]]}
            claim = {"id": "idle", "wires": 4, "cnots": 2, "logical_depth": 2}
            with self.assertRaisesRegex(VerificationError, "ASAP depth"):
                verify_circuit(toy, claim)

    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
    return 0 if result.wasSuccessful() else 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--matrix", action="append", choices=TABLE6_IDS,
                           help="verify selected matrix (repeatable); default: all 30")
    selection.add_argument("--circuit", type=Path,
                           help="verify an external circuit against its Table 6 row; omit file hash check")
    selection.add_argument("--self-test", action="store_true",
                           help="run all valid circuits and deliberate corruption tests")
    parser.add_argument("--json", action="store_true", help="print a machine-readable report")
    args = parser.parse_args(argv)
    if args.self_test and args.json:
        parser.error("--json and --self-test are separate output modes")
    results, failures = [], []
    try:
        if args.self_test:
            return run_self_tests()
        entries = load_package()
        if args.circuit:
            raw = args.circuit.read_bytes()
            data = decode_json(raw)
            require(isinstance(data, dict), "external circuit must be an object")
            selected = [e for e in entries if e["id"] == data.get("id")]
            require(len(selected) == 1, "external circuit has no matching Table 6 ID")
            result = check_external_bytes(
                raw, (ROOT / selected[0]["file"]).read_bytes(), selected[0])
            result["artifact_hash_checked"] = False
            results.append(result)
        else:
            selected = [e for e in entries if not args.matrix or e["id"] in args.matrix]
            for entry in selected:
                try:
                    result = check_bytes((ROOT / entry["file"]).read_bytes(), entry)
                    result["artifact_hash_checked"] = True
                    results.append(result)
                except (OSError, ValueError, UnicodeError) as error:
                    failures.append({"id": entry["id"], "error": str(error)})
    except (OSError, ValueError, UnicodeError) as error:
        failures.append({"id": "package", "error": str(error)})

    if args.json:
        print(json.dumps({"status": "FAIL" if failures else "PASS",
                          "verified": len(results), "results": results,
                          "errors": failures}, indent=2))
    else:
        print(f"{'Matrix':<30} {'Wires':>5} {'CNOTs':>6} {'Depth':>5} {'ASAP':>5} {'Endpoint LB':>11}")
        for result in results:
            print(f"{result['id']:<30} {result['wires']:>5} {result['cnots']:>6} "
                  f"{result['logical_depth']:>5} {result['asap_depth']:>5} "
                  f"{result['endpoint_lower_bound']:>11}  PASS")
        for failure in failures:
            print(f"FAIL {failure['id']}: {failure['error']}", file=sys.stderr)
        print(f"{'FAIL' if failures else 'PASS'}: {len(results)} verified, {len(failures)} failed.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
