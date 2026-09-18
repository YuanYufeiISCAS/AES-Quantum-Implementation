"""Small matrix-to-circuit integration tests, suitable for the default suite."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from linear_surface.io import ROOT, load, save
from linear_surface.verify import verify

BACKEND = Path(os.environ.get("LINEAR_SURFACE_BACKEND", ROOT / "build/linear_surface_backend"))


@unittest.skipUnless(BACKEND.is_file(), "build the C++ backend first")
class CommandTests(unittest.TestCase):
    def test_small_fresh_search(self):
        # This independently constructed instance uses no AES result/recipe.
        instance = {
            "n": 4,
            "layout": {"data_rows": 2, "data_cols": 2, "grid_rows": 5, "grid_cols": 5,
                       "data_pos": [[1, 1], [1, 3], [3, 1], [3, 3]]},
            "target_rows_hex": ["0x5", "0x3", "0x4", "0xb"],
        }
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "instance.json"
            save(source, instance)
            command = [sys.executable, "-m", "linear_surface", "search", "--instance", str(source),
                       "--backend", str(BACKEND), "--seed", "3", "--restarts", "2"]
            first = directory / "first"
            completed = subprocess.run([*command, "--output", str(first)], cwd=ROOT,
                                       text=True, capture_output=True, timeout=30)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            best = load(first / "best.json")
            verify(best, instance)
            self.assertEqual(len(load(first / "run.json")["attempts"]), 3)
            second = directory / "second"
            repeated = subprocess.run([*command, "--output", str(second)], cwd=ROOT,
                                      text=True, capture_output=True, timeout=30)
            self.assertEqual(repeated.returncode, 0, repeated.stderr)
            self.assertEqual(load(second / "best.json"), best)

            # Existing output directories are never overwritten by a run.
            rejected = subprocess.run([*command, "--output", str(first)], cwd=ROOT,
                                      text=True, capture_output=True, timeout=10)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertEqual(load(first / "best.json"), best)

    def test_backend_rejects_oversized_budget(self):
        result = subprocess.run([str(BACKEND), "refine", "1", "999999", "1", "1", "0", "0", "0"],
                                input="", text=True, capture_output=True, timeout=5)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("exceeds supported bounds", result.stderr)


if __name__ == "__main__":
    unittest.main()
