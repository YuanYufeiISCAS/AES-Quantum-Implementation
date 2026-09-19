from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from audit import audit
from quantum import verify_templates, verify_bell_primitive


class TemplateTests(unittest.TestCase):
    def fixture(self, name):
        return json.loads((ROOT / "inputs" / f"{name}.json").read_text())

    def compile(self, spec, valid=True):
        with tempfile.TemporaryDirectory(prefix="ccz_template_test_") as tmp:
            path = Path(tmp)
            (path / "input.json").write_text(json.dumps(spec))
            result = subprocess.run([str(ROOT / "build/lsc_ccz_runner"), str(path / "input.json"), str(path / "trace.jsonl")],
                                    capture_output=True, text=True, timeout=40)
            if not valid:
                self.assertNotEqual(result.returncode, 0)
                return result.stderr
            self.assertEqual(result.returncode, 0, result.stderr)
            obj = json.loads(result.stdout)
            check = audit(spec, path / "trace.jsonl", obj)
            trace = [json.loads(line) for line in (path / "trace.jsonl").read_text().splitlines()]
            return obj, check, trace

    def test_bell_fine_grained_measurements(self):
        self.assertTrue(verify_bell_primitive()["passed"])

    def test_all_ccz_branches(self):
        self.assertTrue(verify_templates(self.fixture("ccz"))["passed"])

    def test_toffoli_branches(self):
        self.assertTrue(verify_templates(self.fixture("toffoli"))["passed"])

    def test_matched_and_erase(self):
        self.assertTrue(verify_templates(self.fixture("qand_pair"))["passed"])

    def test_pauli_conjugation_templates(self):
        self.assertTrue(verify_templates(self.fixture("toffoli_pauli"))["passed"])
        self.assertTrue(verify_templates(self.fixture("qand_pair_pauli"))["passed"])

    def test_wrong_conjugated_pauli_rejected(self):
        spec = self.fixture("toffoli_pauli")
        next(o for o in spec["operations"] if o["kind"] == "x_cond")["kind"] = "z_cond"
        with self.assertRaises(AssertionError): verify_templates(spec)

    def test_qand_dagger_timing(self):
        result, checked, trace = self.compile(self.fixture("qand_dagger"))
        self.assertEqual(result["latency_logical_cycles"], 10)
        self.assertEqual(checked["CCZ_states"], 0)

    def test_ccz_consumption(self):
        _, checked, _ = self.compile(self.fixture("ccz"))
        self.assertEqual(checked["CCZ_states"], 1)
        self.assertEqual(checked["primitive_counts"]["bell"], 3)
        self.assertEqual(checked["primitive_counts"]["reset"], 6)

    def test_matched_port_reset_contract(self):
        spec = self.fixture("toffoli_matched")
        _, checked, _ = self.compile(spec)
        self.assertEqual(checked["primitive_counts"]["reset"], 3)
        self.assertTrue(verify_templates(spec)["passed"])

    def test_consumed_resource_requires_declared_replacement(self):
        spec = self.fixture("toffoli_matched")
        spec.pop("resource_reuse")
        self.assertIn("unclean final resource", self.compile(spec, False))

    def test_ideal_supply_never_removes_port_reset(self):
        spec = self.fixture("toffoli_matched")
        next(o for o in spec["operations"] if o["kind"] == "reset")["kind"] = "barrier"
        self.assertIn("unclean final injection", self.compile(spec, False))

    def test_deleted_cz_predicate_fails_semantics(self):
        spec = self.fixture("ccz")
        next(o for o in spec["operations"] if "guard" in o).pop("guard")
        with self.assertRaises(AssertionError): verify_templates(spec)

    def test_wrong_cz_predicate_fails_semantics(self):
        spec = self.fixture("ccz")
        next(o for o in spec["operations"] if "guard" in o)["guard"]["equals"] = 0
        with self.assertRaises(AssertionError): verify_templates(spec)

    def test_wrong_pauli_feedback_fails_semantics(self):
        spec = self.fixture("ccz")
        op = next(o for o in spec["operations"] if o["kind"] == "z_cond")
        op["controls"][0]["component"] = "k"
        with self.assertRaises(AssertionError): verify_templates(spec)

    def test_deleted_basis_change_fails_semantics(self):
        spec = self.fixture("toffoli")
        next(o for o in spec["operations"] if o["kind"] == "h")["kind"] = "barrier"
        with self.assertRaises(AssertionError): verify_templates(spec)

    def test_dirty_resource_rejected(self):
        spec = self.fixture("ccz")
        op = next(o for o in spec["operations"] if o["kind"] == "reset")
        op["kind"] = "barrier"
        self.assertIn("unclean", self.compile(spec, False))

    def test_nonadjacent_bell_rejected(self):
        spec = self.fixture("ccz")
        q = next(p for p in spec["layout"]["patches"] if p["role"] == "CCZ_resource")
        q["cell"][0] += 2
        self.assertIn("adjacent", self.compile(spec, False))

    def test_wrong_bell_orientation_rejected(self):
        spec = self.fixture("ccz")
        next(o for o in spec["operations"] if o["kind"] == "bell")["parity"] = "XX"
        self.assertIn("Bell parity", self.compile(spec, False))

    def test_arbitrary_duration_rejected(self):
        spec = self.fixture("ccz")
        next(o for o in spec["operations"] if o["kind"] == "bell")["duration"] = 0
        self.assertIn("duration", self.compile(spec, False))

    def test_h_footprint_obstacle(self):
        spec = self.fixture("toffoli")
        r, c = spec["layout"]["patches"][2]["cell"]
        spec["layout"]["patches"].append({"id": len(spec["layout"]["patches"]), "cell": [r + 1, c + 1], "role": "data"})
        self.assertIn("deadlock", self.compile(spec, False))

    def test_feedback_delay_is_charged(self):
        spec = self.fixture("qand_dagger")
        spec["feedback_latency"] = 7
        result, _, trace = self.compile(spec)
        measure = next(e for e in trace if e["kind"] == "mx")
        correction = next(e for e in trace if "guard" in e)
        self.assertGreaterEqual(correction["start"], measure["end"] + 7)
        self.assertGreater(result["latency_logical_cycles"], 10)

    def test_module_stage_policy(self):
        spec = self.fixture("qand_dagger")
        spec["scheduling_policy"] = "dependency_overlap"
        result, _, _ = self.compile(spec)
        self.assertEqual(result["latency_logical_cycles"], 8)

    def test_trace_corruptions(self):
        spec = self.fixture("ccz")
        result, _, trace = self.compile(spec)
        variants = []
        bad = deepcopy(trace)
        next(e for e in bad if e["kind"] == "h")["footprint"].pop()
        variants.append(bad)
        bad = deepcopy(trace)
        next(e for e in bad if e["kind"] == "bell")["end"] -= 1
        variants.append(bad)
        bad = deepcopy(trace)
        next(e for e in bad if "guard" in e)["feedback_ready"] += 1
        variants.append(bad)
        bad = deepcopy(trace)
        next(e for e in bad if e["kind"] == "cx")["footprint"][1] = [-1, -1]
        variants.append(bad)
        bad = deepcopy(trace)
        bad.pop()
        variants.append(bad)
        for bad in variants:
            with tempfile.TemporaryDirectory(prefix="ccz_negative_trace_") as tmp:
                path = Path(tmp) / "trace.jsonl"
                path.write_text("\n".join(json.dumps(e) for e in bad))
                with self.assertRaises(ValueError): audit(spec, path, result)


if __name__ == "__main__":
    unittest.main()
