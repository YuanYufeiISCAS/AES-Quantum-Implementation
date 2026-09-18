"""Independent certificate checks, including malformed/corrupted witnesses."""

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

from linear_surface.io import ROOT, load
from linear_surface.reproduce import check_recipe
from linear_surface.search import Population, diversity_key
from linear_surface.verify import check_instance, verify


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.instance = load(ROOT / "instances/aes.json")
        cls.circuit = load(ROOT / "circuits/aes.json")

    def reject(self, change, message=None):
        damaged = deepcopy(self.circuit)
        change(damaged)
        with self.assertRaises(ValueError) as context:
            verify(damaged, self.instance)
        if message:
            self.assertIn(message, str(context.exception))

    def test_saved_witness(self):
        self.assertEqual(verify(self.circuit, self.instance),
                         {"cnots": 103, "layers": 10, "cycles": 20})

    def test_aes_matrix_from_field_definition(self):
        # Independently expand [02 03 01 01] and its rotations over
        # GF(2^8)/(x^8+x^4+x^3+x+1), with LSB-first bits within each byte.
        coefficients = ((2, 3, 1, 1), (1, 2, 3, 1), (1, 1, 2, 3), (3, 1, 1, 2))
        matrix = [0] * 32
        for column in range(32):
            input_byte, bit = divmod(column, 8)
            value = 1 << bit
            times_two = ((value << 1) ^ (0x11B if value & 0x80 else 0)) & 255
            for output_byte in range(4):
                multiplier = coefficients[output_byte][input_byte]
                image = value if multiplier == 1 else times_two if multiplier == 2 else times_two ^ value
                for output_bit in range(8):
                    matrix[8 * output_byte + output_bit] |= ((image >> output_bit) & 1) << column
        self.assertEqual(check_instance(self.instance), matrix)

    def test_wrong_cost(self):
        self.reject(lambda c: c["stats"].update(cycles=18), "reported cost")

    def test_wrong_layer_cost(self):
        self.reject(lambda c: c["layers"][0].update(cycles=1), "two cycles")

    def test_repeated_endpoint(self):
        self.reject(lambda c: c["layers"][0]["operations"].append(
            deepcopy(c["layers"][0]["operations"][0])), "repeated endpoint")

    def test_reversed_cnot(self):
        def change(c):
            op = c["layers"][0]["operations"][0]
            op["control"], op["target"] = op["target"], op["control"]
        self.reject(change, "endpoints")

    def test_path_out_of_grid(self):
        self.reject(lambda c: c["layers"][0]["operations"][0]["path"].__setitem__(1, [-1, 15]), "path row")

    def test_nonadjacent_path(self):
        self.reject(lambda c: c["layers"][0]["operations"][0]["path"].__setitem__(1, [0, 0]), "nonadjacent")

    def test_repeated_vertex(self):
        self.reject(lambda c: c["layers"][0]["operations"][0]["path"].insert(1, [5, 15]), "self-intersecting")

    def test_data_obstacle(self):
        self.reject(lambda c: c["layers"][0]["operations"][0]["path"].__setitem__(1, [5, 13]), "data site")

    def test_wrong_output_permutation(self):
        def change(c):
            p = c["output_permutation"]
            p[0], p[1] = p[1], p[0]
        self.reject(change, "CNOT replay")

    def test_wrong_target_is_not_self_certified(self):
        self.reject(lambda c: c["target_rows_hex"].__setitem__(0, hex(
            int(c["target_rows_hex"][0], 16) ^ int(c["target_rows_hex"][1], 16))), "changed the target")

    def test_layout_is_fixed(self):
        self.reject(lambda c: c["layout"]["data_pos"].__setitem__(0, [0, 0]), "physical layout")

    def test_boolean_is_not_wire_index(self):
        self.reject(lambda c: c["layers"][0]["operations"][0].update(control=True), "control")

    def test_recipe_is_parameter_only(self):
        recipe = load(ROOT / "recipes/aes_sota.json")
        check_recipe(recipe)
        recipe["stages"][0]["gates"] = [[0, 1]]
        with self.assertRaisesRegex(ValueError, "unexpected stage fields"):
            check_recipe(recipe)

    def test_duplicate_json_keys(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.json"
            path.write_text('{"n": 32, "n": 64}', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate JSON field"):
                load(path)

    def test_population_is_validated_and_deduplicated(self):
        pool = Population(self.instance)
        pool.admit(self.circuit)
        pool.admit(deepcopy(self.circuit))
        self.assertEqual(len(pool.entries), 1)
        choice = pool.select()
        self.assertEqual(diversity_key(choice), diversity_key(self.circuit))
        choice["stats"]["cycles"] = 999
        self.assertEqual(pool.best["stats"]["cycles"], 20)
        with self.assertRaises(ValueError):
            pool.admit(choice)


class GeometryTests(unittest.TestCase):
    def witness(self, gates):
        layout = {"data_rows": 2, "data_cols": 2, "grid_rows": 5, "grid_cols": 5,
                  "data_pos": [[1, 1], [1, 3], [3, 1], [3, 3]]}
        state = [1, 2, 4, 8]
        for c, t, _ in gates:
            state[t] ^= state[c]
        instance = {"n": 4, "layout": layout, "target_rows_hex": list(map(hex, state))}
        circuit = {**instance, "mode": "vdp", "output_permutation": [0, 1, 2, 3],
                   "stats": {"cnots": len(gates), "layers": 1, "cycles": 2},
                   "layers": [{"cycles": 2, "operations": [
                       {"control": c, "target": t, "path": p} for c, t, p in gates]}]}
        return circuit, instance

    def test_valid_path(self):
        c, i = self.witness([(0, 1, [[1, 1], [0, 1], [0, 2], [1, 2], [1, 3]])])
        verify(c, i)

    def test_control_port(self):
        c, i = self.witness([(0, 1, [[1, 1], [1, 2], [1, 3]])])
        with self.assertRaisesRegex(ValueError, "orientation"):
            verify(c, i)

    def test_target_port(self):
        c, i = self.witness([(0, 1, [[1, 1], [0, 1], [0, 2], [0, 3], [1, 3]])])
        with self.assertRaisesRegex(ValueError, "orientation"):
            verify(c, i)

    def test_vertex_sharing_with_distinct_edges(self):
        c, i = self.witness([
            (0, 3, [[1, 1], [2, 1], [2, 2], [2, 3], [2, 4], [3, 4], [3, 3]]),
            (1, 2, [[1, 3], [0, 3], [0, 2], [1, 2], [2, 2], [3, 2], [3, 1]]),
        ])
        with self.assertRaisesRegex(ValueError, "share a vertex"):
            verify(c, i)


if __name__ == "__main__":
    unittest.main()
