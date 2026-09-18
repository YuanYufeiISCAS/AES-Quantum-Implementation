"""Portable circuit I/O and exact, non-cryptographic state comparisons."""

from __future__ import annotations

import gzip
import json
import zlib
from pathlib import Path

SCHEMA = "sbox-surface-code-v1"
ROOT = Path(__file__).resolve().parent.parent
CASES = tuple(f"{family}_{row}row" for family in ("inplace", "cstar") for row in range(1, 6))


def canonical(value) -> str:
    """An exact canonical representation, used for equality and stable ordering."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def read(path: str | Path):
    path = Path(path)
    raw = path.read_bytes()
    if path.suffix == ".gz":
        raw = gzip.decompress(raw)
    return json.loads(raw)


def write(path: str | Path, value) -> None:
    """Create a new result; never replace an existing circuit."""
    path = Path(path)
    raw = (canonical(value) + "\n").encode()
    if path.suffix == ".gz":
        raw = gzip.compress(raw, compresslevel=9, mtime=0)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write(raw)


def fixture(case: str, *, parent: bool = False) -> dict:
    if case not in CASES:
        raise ValueError(f"unknown case: {case}")
    return read(ROOT / ("parents" if parent else "results") / f"{case}.json.gz")


def objective(circuit: dict) -> tuple[int, int, int]:
    return tuple(circuit["cost"][key] for key in ("latency", "N_H", "N_CNOT_linear"))


def physical_key(circuit: dict) -> bytes:
    # Lossless compression bounds beam/tabu memory while retaining exact equality.
    raw = canonical([circuit["layout"], circuit["operations"], circuit["schedule"]["stages"]])
    return zlib.compress(raw.encode(), level=1)
