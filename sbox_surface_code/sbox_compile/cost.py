"""cost."""

# Adapted from the authors' strict_sbox/scheduled_compiler.py; campaign infrastructure removed.
from __future__ import annotations
from collections import Counter
from . import native as native_compiler
from .geometry import _require


def scheduled_cost(operations, stages, *, parallel_h: bool = False) -> dict[str, int]:
    """Counts come from gates; latency comes from the materialized stages."""
    cost = native_compiler._native_cost(operations)
    by_id = {op["id"]: op for op in operations}
    cycles: Counter[str] = Counter()
    batches: Counter[str] = Counter()
    now = 0
    for stage in stages:
        ids = stage["operation_ids"]
        _require(bool(ids), "empty execution stage")
        selected = [by_id[identifier] for identifier in ids]
        _require(len({op["kind"] for op in selected}) == 1, "mixed-kind execution stage")
        _require(len({op["category"] for op in selected}) == 1, "mixed-category stage")
        duration = stage["end"] - stage["start"]
        _require(type(duration) is int and duration >= 0, "invalid stage duration")
        _require(stage["start"] == now, "schedule has an overlap or uncharged gap")
        _require(all((op["duration"] == duration for op in selected)), "stage duration mismatch")
        now = stage["end"]
        cycles[selected[0]["category"]] += duration
        batches[selected[0]["kind"]] += 1
    cost.update(
        latency=now,
        D_CNOT=cycles["linear"],
        H_cycles=cycles["hadamard"],
        CCZ_injection_cycles=cycles["ccz_injection"],
        correction_CNOT_cycles=cycles["correction"],
        QAND_measurement_cycles=cycles["qand_measurement"],
        QAND_reset_cycles=cycles["qand_reset"],
        port_reset_cycles=cycles["port_reset"],
        conditional_CNOT_batches=batches["conditional_cnot"],
        QAND_measurement_batches=batches["qand_measure_x"],
        QAND_reset_batches=batches["qand_reset"],
    )
    if parallel_h:
        cost["H_batches"] = batches["h"]
        cost["H_parallel_saved_cycles"] = 3 * (cost["N_H"] - cost["H_batches"])
    expected = (
        cost["D_CNOT"]
        + 5 * cost["R_CCZ"]
        + 3 * (cost["H_batches"] if parallel_h else cost["N_H"])
        + 2 * cost["conditional_CNOT_batches"]
        + cost["QAND_measurement_cycles"]
        + cost["QAND_reset_cycles"]
        + cost["port_reset_cycles"]
    )
    _require(now == expected, "event-scheduled primitive cost identity failed")
    return cost
