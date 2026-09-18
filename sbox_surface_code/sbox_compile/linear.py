"""Seeded/free-output rewrite and independent fresh synthesis of paid matrices."""

from __future__ import annotations

import copy
import json
import subprocess

from .geometry import _layer_geometry, _require
from .matrix import replay_rows
from .model import ROOT, canonical
from .native import NativeGeometry, color_path_conflicts


def word(candidate: dict) -> list[list[int]]:
    return [
        [op["control"], op["target"]] for layer in candidate["layers"] for op in layer["operations"]
    ]


def gaussian(rows: list[int]) -> list[list[int]]:
    """Construct a fixed-output reference using only row additions."""
    reduced = list(rows)
    gates = []
    for i in range(len(rows)):
        if not (reduced[i] >> i) & 1:
            j = next((j for j in range(i + 1, len(rows)) if (reduced[j] >> i) & 1), None)
            _require(j is not None, "singular matrix")
            reduced[i] ^= reduced[j]
            gates.append([j, i])
        for j in range(len(rows)):
            if j != i and (reduced[j] >> i) & 1:
                reduced[j] ^= reduced[i]
                gates.append([i, j])
    return list(reversed(gates))


def routed_reference(rows, geom, gates=None) -> dict:
    gates = gaussian(rows) if gates is None else gates
    layers = [
        {
            "mode": "vdp",
            "logical_depth": 2,
            "operations": [{"control": c, "target": t, "path": geom.cnot(c, t)}],
        }
        for c, t in gates
    ]
    return {
        "n": len(rows),
        "layout": geom.layout.to_dict(),
        "verified": True,
        "output_permutation": list(range(len(rows))),
        "layers": layers,
        "stats": {"surface_depth": 2 * len(layers), "layers": len(layers), "cnots": len(gates)},
    }


def validate(
    candidate: dict, rows: list[int], geom: NativeGeometry, *, free=False
) -> tuple[int, int]:
    n = len(rows)
    q = candidate["output_permutation"]
    _require(
        candidate["n"] == n and sorted(q) == list(range(n)), "invalid matrix dimensions/permutation"
    )
    _require(free or q == list(range(n)), "unpaid free-output permutation")
    _require(replay_rows(candidate) == [rows[i] for i in q], "matrix replay failed")
    _require(
        candidate["layout"]["data_rows"] == geom.layout.data_rows
        and candidate["layout"]["data_cols"] == geom.layout.data_cols,
        "matrix layout mismatch",
    )
    batches = count = reported = 0
    offsets = {
        (
            op["path"][0][0] - geom.coord(op["control"])[0],
            op["path"][0][1] - geom.coord(op["control"])[1],
        )
        for layer in candidate["layers"]
        for op in layer["operations"]
    }
    _require(len(offsets) <= 1 and offsets <= {(0, 0), (1, 1)}, "mixed route coordinates")
    shifted = next(iter(offsets), (1, 1)) == (1, 1)
    for layer in candidate["layers"]:
        ops, duration = _layer_geometry(geom, layer, source_shift=shifted)
        batches += len(color_path_conflicts([op["path"] for op in ops])["batches"])
        reported += duration
        count += len(ops)
    _require(
        candidate["stats"]
        == {"surface_depth": reported, "layers": len(candidate["layers"]), "cnots": count},
        "matrix statistics mismatch",
    )
    return 2 * batches, count


def _run(name: str, request: str, args=()) -> dict:
    binary = ROOT / "build" / name
    if not binary.is_file():
        raise FileNotFoundError(
            "Build the kernels first: cmake -S . -B build && cmake --build build -j2"
        )
    result = subprocess.run(
        [str(binary), *args], input=request, text=True, capture_output=True, check=True
    )
    return json.loads(result.stdout)


def rewrite(
    rows: list[int],
    geom: NativeGeometry,
    *,
    seed: int,
    work: int,
    warm: dict | None = None,
    free: bool = False,
) -> list[dict]:
    layout = geom.layout
    values = (
        len(rows),
        layout.data_rows,
        layout.data_cols,
        layout.grid_rows,
        layout.grid_cols,
        seed,
        int(free),
        3,
        64 * work,
        work,
        2048,
    )
    lines = [
        "JOINT_LINEAR_V1",
        " ".join(map(str, values)),
        " ".join(map(str, rows)),
        "1" if warm is not None else "0",
    ]
    if warm is not None:
        _require(replay_rows(warm) == rows, "warm network does not realize paid matrix")
        gates = word(warm)
        lines.extend((" ".join(map(str, range(len(rows)))), str(len(gates))))
        lines.extend(f"{c} {t}" for c, t in gates)
    report = _run("rewrite", "\n".join(lines) + "\n", ("--stdin-v1",))
    candidates = []
    for item in report["candidates"]:
        c = item["candidate"]
        c["verified"] = True
        validate(c, rows, geom, free=free)
        candidates.append(c)
    return candidates


def fresh(
    rows: list[int],
    geom: NativeGeometry,
    *,
    seed: int,
    work: int,
    depth_goal: int,
    restarts: int = 8,
) -> list[dict]:
    """No incumbent network or answer is passed to the independent engine."""
    candidates = []
    for ordinal in range(restarts):
        inverse = ordinal % 2 == 1
        target = invert_matrix(rows) if inverse else rows
        profile = min(2, ordinal // 2)
        values = (
            len(rows),
            geom.layout.data_rows,
            geom.layout.data_cols,
            (seed + 104729 * ordinal) % (2**63),
            depth_goal,
            profile,
            work,
        )
        report = _run(
            "fresh", " ".join(map(str, values)) + "\n" + " ".join(map(str, target)) + "\n"
        )
        for c in report["candidates"]:
            if inverse:
                c = copy.deepcopy(c)
                c["layers"].reverse()
            validate(c, rows, geom)
            candidates.append(c)
    return portfolio(candidates, rows, geom)


def invert_matrix(rows: list[int]) -> list[int]:
    n = len(rows)
    left, right = list(rows), [1 << i for i in range(n)]
    for i in range(n):
        pivot = next((j for j in range(i, n) if (left[j] >> i) & 1), None)
        _require(pivot is not None, "singular matrix")
        left[i], left[pivot] = left[pivot], left[i]
        right[i], right[pivot] = right[pivot], right[i]
        for j in range(n):
            if j != i and (left[j] >> i) & 1:
                left[j] ^= left[i]
                right[j] ^= right[i]
    return right


def portfolio(candidates, rows, geom, *, limit=3, fresh_candidates=()) -> list[dict]:
    unique = {canonical(c["layers"]): c for c in candidates}
    ordered = sorted(
        unique.values(), key=lambda c: (*validate(c, rows, geom), canonical(c["layers"]))
    )
    selected = ordered[:limit]
    if fresh_candidates and limit > 1 and not any(c in fresh_candidates for c in selected):
        selected[-1:] = [
            min(fresh_candidates, key=lambda c: (*validate(c, rows, geom), canonical(c["layers"])))
        ]
    return selected
