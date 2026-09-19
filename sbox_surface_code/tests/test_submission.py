"""Regression, adversarial, and algorithm-structure tests; standard library only."""

from __future__ import annotations

import copy
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

from sbox_compile import linear, nonlinear, windows
from sbox_compile.compiler import compile_circuit
from sbox_compile.matrix import replay_rows
from sbox_compile.model import CASES, ROOT, canonical, fixture, objective, read, write
from sbox_compile.replay import replay
from sbox_compile.search import State, select_beam
from sbox_compile.verify import require_verified, verify

EXPECTED = {
    "inplace_1row": (576, 254, 330),
    "inplace_2row": (539, 254, 314),
    "inplace_3row": (505, 256, 419),
    "inplace_4row": (546, 256, 336),
    "inplace_5row": (546, 254, 421),
    "cstar_1row": (506, 208, 332),
    "cstar_2row": (475, 204, 356),
    "cstar_3row": (471, 206, 367),
    "cstar_4row": (459, 206, 349),
    "cstar_5row": (473, 210, 370),
}


class Fixtures(unittest.TestCase):
    def test_all_complete_circuits(self):
        for case in CASES:
            with self.subTest(case=case):
                circuit = fixture(case)
                report = require_verified(circuit)
                self.assertEqual(objective(circuit), EXPECTED[case])
                self.assertEqual(
                    report["checked_inputs"], 65536 if case.startswith("cstar") else 256
                )
                self.assertTrue(report["scratch_zero"])

    def test_result_index(self):
        index = read(ROOT / "results/index.json.gz")
        expected_files = {f"{case}.json.gz" for case in CASES} | {"inplace_3row_cover.json.gz"}
        self.assertEqual(len(index), len(expected_files))
        self.assertEqual({entry["file"] for entry in index}, expected_files)
        for entry in index:
            with self.subTest(file=entry["file"]):
                self.assertEqual(set(entry), {"case", "file", "cost", "status"})
                circuit = read(ROOT / "results" / entry["file"])
                self.assertEqual(entry["case"], circuit["case"])
                self.assertEqual(entry["cost"], circuit["cost"])
                self.assertEqual(entry["status"], "verified_search_incumbent")

    def test_logical_segment_fields(self):
        required = {"index", "kind", "component", "gates", "cnot_count", "ccz_count"}
        allowed = required | {"target_rows_hex"}
        for directory in ("parents", "results"):
            for path in sorted((ROOT / directory).glob("*.json.gz")):
                if path.name == "index.json.gz":
                    continue
                with self.subTest(directory=directory, case=path.stem):
                    for segment in read(path)["logical_source"]["segments"]:
                        self.assertLessEqual(required, set(segment))
                        self.assertLessEqual(set(segment), allowed)

    def test_recipe_fields(self):
        for case in CASES:
            with self.subTest(case=case):
                recipe = read(ROOT / "recipes" / f"{case}.json.gz")
                self.assertEqual(set(recipe), {"case", "seed", "steps", "expected_cost"})
                self.assertEqual(recipe["case"], case)
                self.assertIs(type(recipe["seed"]), int)
                for step in recipe["steps"]:
                    self.assertEqual(
                        set(step), {"start", "stop", "boundaries", "window", "nonlinear"}
                    )

    def test_all_recipe_replays(self):
        for case in CASES:
            with self.subTest(case=case):
                parent = fixture(case, parent=True)
                snapshot = canonical(parent)
                recipe = read(ROOT / "recipes" / f"{case}.json.gz")
                # The replay engine is not allowed to open any result file.
                with patch(
                    "pathlib.Path.read_bytes", side_effect=AssertionError("unexpected file read")
                ):
                    rebuilt = replay(parent, recipe)
                self.assertEqual(parent, read(ROOT / "parents" / f"{case}.json.gz"))
                self.assertEqual(canonical(parent), snapshot)
                self.assertEqual(rebuilt, fixture(case))

    def test_recompile_physical_inputs(self):
        circuit = fixture("cstar_4row")
        rebuilt = compile_circuit(circuit, stages=circuit["schedule"]["stages"])
        self.assertEqual(rebuilt, circuit)

    def test_later_placement_cover_result(self):
        circuit = read(ROOT / "results/inplace_3row_cover.json.gz")
        require_verified(circuit)
        self.assertEqual(objective(circuit), (505, 256, 417))

    def test_exclusive_deterministic_io(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "data.json.gz"
            write(path, {"b": 1, "a": [2]})
            raw = path.read_bytes()
            self.assertEqual(read(path), {"a": [2], "b": 1})
            with self.assertRaises(FileExistsError):
                write(path, {})
            second = Path(directory) / "second.json.gz"
            write(second, {"a": [2], "b": 1})
            self.assertEqual(raw, second.read_bytes())


class Corruption(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = fixture("cstar_4row")

    def reject(self, mutate):
        circuit = copy.deepcopy(self.circuit)
        mutate(circuit)
        self.assertFalse(verify(circuit)["passed"])

    def test_cost(self):
        self.reject(lambda c: c["cost"].update(latency=458))

    def test_route(self):
        def mutate(c):
            op = next(op for op in c["operations"] if op["kind"] == "linear_cnot_layer")
            op["operations"][0]["path"][1] = [999, 999]

        self.reject(mutate)

    def test_full_h_footprint(self):
        def mutate(c):
            op = next(op for op in c["operations"] if op["kind"] == "h")
            op["footprint"]["reserved_patch_coords"].pop()

        self.reject(mutate)

    def test_conditional_correction(self):
        def mutate(c):
            op = next(op for op in c["operations"] if op["kind"] == "conditional_cnot")
            op["condition"] = "unproduced"

        self.reject(mutate)

    def test_missing_reset(self):
        def mutate(c):
            i = next(i for i, op in enumerate(c["operations"]) if op["kind"] == "qand_reset")
            c["operations"].pop(i)

        self.reject(mutate)

    def test_primitive_branch(self):
        def mutate(c):
            op = next(op for op in c["operations"] if op["kind"] == "linear_cnot_layer")
            op["native_cnot_contracts"][0]["explicit_paulis"].pop()

        self.reject(mutate)

    def test_unpaid_placement(self):
        def mutate(c):
            p = c["linear_candidates"][1]["choice"]["input_placement"]
            p[0], p[1] = p[1], p[0]

        self.reject(mutate)

    def test_duplicate_stage(self):
        self.reject(lambda c: c["schedule"]["stages"].append(c["schedule"]["stages"][0]))

    def test_recipe_boundary(self):
        recipe = read(ROOT / "recipes/cstar_4row.json.gz")
        recipe["steps"][1]["boundaries"][1] = list(range(28))
        with self.assertRaises(ValueError):
            replay(fixture("cstar_4row", parent=True), recipe)


class Algorithm(unittest.TestCase):
    def test_paid_independent_boundaries(self):
        circuit = fixture("cstar_4row", parent=True)
        rng = random.Random(47)
        for index in (0, 6, 14, 16):
            pin, pout = list(range(28)), list(range(28))
            rng.shuffle(pin)
            rng.shuffle(pout)
            task = windows.task(circuit, index, pin, pout)
            old = replay_rows(circuit["linear_candidates"][index]["candidate"])
            expected = [0] * 28
            for row in range(28):
                expected[pout[row]] = sum(((old[row] >> col) & 1) << pin[col] for col in range(28))
            self.assertEqual(task["rows"], expected)
            self.assertEqual(replay_rows(task["reference"]), expected)

    def test_independent_multiregion_combinations(self):
        regions = [
            {"name": str(i), "entries": [{"key": f"{i}:{j}"} for j in range(3)]} for i in range(3)
        ]
        combinations = nonlinear.combinations(regions)
        self.assertEqual(len(combinations[0]), 3)
        self.assertTrue(any(len(c) == 2 for c in combinations))
        self.assertTrue(any(len(c) == 1 for c in combinations))
        self.assertLessEqual(len(combinations), 24)

    def test_uphill_beam_and_monotone_best(self):
        def state(cost, placement, key):
            circuit = {
                "cost": dict(latency=cost, N_H=1, N_CNOT_linear=1),
                "linear_candidates": [{"choice": {"output_placement": placement}}],
            }
            return State(circuit, [], key)

        candidates = [
            state(100, [0, 1], "best"),
            state(105, [1, 0], "uphill"),
            state(133, [2, 1], "outside"),
        ]
        selected = select_beam(candidates)
        self.assertEqual(selected[0].key, "best")
        self.assertIn("uphill", [s.key for s in selected])
        self.assertNotIn("outside", [s.key for s in selected])

    @unittest.skipUnless((ROOT / "build/rewrite").is_file(), "build C++ kernels first")
    def test_deterministic_rewrite_and_free_output(self):
        circuit = fixture("cstar_4row", parent=True)
        geom = windows.geometry(circuit)
        identity = list(range(28))
        task = windows.task(circuit, 6, identity, identity)
        for free in (False, True):
            a = linear.rewrite(
                task["rows"], geom, seed=19, work=32, warm=task["reference"], free=free
            )
            b = linear.rewrite(
                task["rows"], geom, seed=19, work=32, warm=task["reference"], free=free
            )
            self.assertEqual(a, b)
            self.assertTrue(a)
            for candidate in a:
                linear.validate(candidate, task["rows"], geom, free=free)

    @unittest.skipUnless((ROOT / "build/fresh").is_file(), "build C++ kernels first")
    def test_independent_fresh_rediscovers_winning_matrix(self):
        circuit = fixture("cstar_4row", parent=True)
        geom = windows.geometry(circuit)
        identity = list(range(28))
        task = windows.task(circuit, 6, identity, identity)
        actual = linear.fresh(
            task["rows"], geom, seed=553628090, work=300000, depth_goal=5, restarts=1
        )
        # Reference is read only after the independent search has finished.
        recipe = read(ROOT / "recipes/cstar_4row.json.gz")
        winner = recipe["steps"][0]["window"]["linear_candidates"][6]["candidate"]

        def physical(candidate):
            return [
                [(op["control"], op["target"], op["path"]) for op in layer["operations"]]
                for layer in candidate["layers"]
            ]

        self.assertEqual(actual[0]["stats"], {"surface_depth": 10, "layers": 5, "cnots": 22})
        self.assertEqual(physical(actual[0]), physical(winner))


if __name__ == "__main__":
    unittest.main()
