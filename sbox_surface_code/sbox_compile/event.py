"""event."""

# Adapted from the authors' strict_sbox/event_scheduler.py; campaign infrastructure removed.
from __future__ import annotations
from collections import Counter
from typing import Any, Mapping, Sequence

SCHEMA = "strict-sbox-fixed-native-event-schedule-v1"
PARALLEL_H_SCHEMA = "strict-sbox-fixed-native-event-schedule-v2"
_CCZ_KINDS = frozenset({"conditional_z", "h", "conditional_cnot"})
_QAND_KINDS = frozenset({"qand_measure_x", "h", "conditional_cnot", "qand_reset"})
_PARALLEL_KINDS = frozenset({"conditional_cnot", "qand_measure_x", "qand_reset"})
_DURATIONS = {
    "h": 3,
    "conditional_z": 0,
    "conditional_cnot": 2,
    "qand_measure_x": 1,
    "qand_reset": 1,
}


class _BudgetExceeded(Exception):
    pass


def _outputs(operation: Mapping[str, Any]) -> set[str]:
    if operation["kind"] == "native_ccz_bell":
        return {
            value for gadget in operation.get("gadgets", ()) for value in gadget.get("outcomes", ())
        }
    if operation["kind"] == "qand_measure_x":
        return {operation["outcome"]}
    return set()


def _h_patch_footprint(operation: Mapping[str, Any]) -> set[tuple[int, int]] | None:
    """Require the complete declared 2x2 macro, not just its data patch."""
    try:
        if (
            operation["kind"] != "h"
            or operation["duration"] != 3
            or type(operation["wire"]) is not int
            or (operation["data_qubits"] != [operation["wire"]])
        ):
            return None
        footprint = operation["footprint"]
        raw = footprint["reserved_patch_coords"]
        if len(raw) != 4 or any(
            (len(point) != 2 or any((type(value) is not int for value in point)) for point in raw)
        ):
            return None
        patches = set(map(tuple, raw))
        rows, cols = ({point[0] for point in patches}, {point[1] for point in patches})
        if (
            len(rows) != 2
            or len(cols) != 2
            or max(rows) - min(rows) != 1
            or (max(cols) - min(cols) != 1)
            or (patches != {(row, col) for row in rows for col in cols})
        ):
            return None
        data = footprint["data_coord"]
        auxiliary = footprint["auxiliary_coords"]
        if (
            len(data) != 2
            or any((type(value) is not int for value in data))
            or tuple(data) not in patches
            or (len(auxiliary) != 3)
            or any(
                (
                    len(point) != 2 or any((type(value) is not int for value in point))
                    for point in auxiliary
                )
            )
            or (set(map(tuple, auxiliary)) != patches - {tuple(data)})
        ):
            return None
        return patches
    except (KeyError, TypeError, ValueError, IndexError):
        return None


def _window_reason(
    pipelines: Sequence[Sequence[Mapping[str, Any]]],
    source_ids: Sequence[str],
    records: Mapping[str, Mapping[str, Any]],
    ready: set[str],
    kind: str,
    *,
    parallel_h: bool = False,
) -> str | None:
    """Conservative admission; unrecognized inputs remain serial."""
    if not 1 <= len(pipelines) <= 5:
        return "unsupported_pipeline_count"
    used: set[int] = set()
    for source_id, pipeline in zip(source_ids, pipelines, strict=True):
        record = records.get(source_id)
        if record is None or not pipeline or len(pipeline) > 16:
            return "unsupported_or_missing_pipeline"
        physical = record.get("physical_indices", ())
        if len(physical) != 3 or len(set(physical)) != 3 or used & set(physical):
            return "non_disjoint_gadget_supports"
        used.update(physical)
        roles = {"qand", "toffoli"} if kind == "ccz_corrections" else {"qand_dagger"}
        if record.get("role") not in roles:
            return "unsupported_gadget_role"
        available = set(ready)
        for operation in pipeline:
            operation_kind = operation["kind"]
            if operation.get("duration") != _DURATIONS.get(operation_kind):
                return "unsupported_primitive_duration"
            if parallel_h and operation_kind == "h" and (_h_patch_footprint(operation) is None):
                return "missing_or_invalid_complete_H_footprint"
            support = operation.get("data_qubits", ())
            if not support or not set(support) <= set(physical):
                return "nonlocal_or_missing_operation_support"
            if operation_kind == "conditional_cnot":
                if support != [operation.get("control"), operation.get("target")]:
                    return "cnot_support_mismatch"
                path = operation.get("path", ())
                if len(path) < 2 or any((len(point) != 2 for point in path)):
                    return "missing_complete_cnot_path"
                if operation.get("condition") not in available:
                    return "unavailable_or_cross_pipeline_outcome"
            else:
                if support != [operation.get("wire")]:
                    return "single_wire_support_mismatch"
                if operation_kind == "conditional_z":
                    if (
                        "condition_variables" not in operation
                        or not set(operation["condition_variables"]) <= available
                    ):
                        return "unavailable_or_undeclared_z_outcome"
                elif operation_kind in {"qand_measure_x", "qand_reset"}:
                    if operation["wire"] != physical[2]:
                        return "qand_target_mismatch"
            available.update(_outputs(operation))
    return None


def _solve_window(
    pipelines: Sequence[Sequence[Mapping[str, Any]]],
    budget: dict[str, int],
    *,
    parallel_h: bool = False,
    dependencies: Mapping[str, Sequence[str]] | None = None,
) -> list[list[str]]:
    """Exact cursor DP, with an eager singleton rule for ready zero-time Zs.

    Advancing a ready zero-time Z cannot delay another pipeline or consume a
    resource interval.  Disjoint gadget supports and the absence of cross-
    pipeline outcome dependencies make that canonicalization cost preserving.
    """
    ends = tuple(map(len, pipelines))
    footprints = [
        [
            (
                set(map(tuple, operation["path"]))
                if operation["kind"] == "conditional_cnot"
                else (
                    _h_patch_footprint(operation)
                    if parallel_h and operation["kind"] == "h"
                    else set(operation["data_qubits"])
                )
            )
            for operation in pipeline
        ]
        for pipeline in pipelines
    ]
    positions = {
        operation["id"]: (index, position)
        for index, pipeline in enumerate(pipelines)
        for position, operation in enumerate(pipeline)
    }
    predecessors = {
        identifier: [positions[previous] for previous in previous_ids]
        for identifier, previous_ids in (dependencies or {}).items()
    }
    costs: dict[tuple[int, ...], int] = {}
    choices: dict[tuple[int, ...], tuple[int, ...]] = {}

    def advance(state: tuple[int, ...], indices: tuple[int, ...]) -> tuple[int, ...]:
        result = list(state)
        for index in indices:
            result[index] += 1
        return tuple(result)

    def moves(state: tuple[int, ...]):
        heads = [
            (index, pipeline[state[index]])
            for index, pipeline in enumerate(pipelines)
            if state[index] < len(pipeline)
            and all(
                (
                    state[other] > position
                    for other, position in predecessors.get(pipeline[state[index]]["id"], ())
                )
            )
        ]
        zero = next(
            (index for index, operation in heads if operation["kind"] == "conditional_z"), None
        )
        if zero is not None:
            yield (0, (zero,))
            return
        if not parallel_h:
            for index, operation in heads:
                if operation["kind"] == "h":
                    yield (operation["duration"], (index,))
        for kind in sorted(_PARALLEL_KINDS | {"h"} if parallel_h else _PARALLEL_KINDS):
            ready = [index for index, operation in heads if operation["kind"] == kind]
            for mask in range(1, 1 << len(ready)):
                indices = tuple((ready[bit] for bit in range(len(ready)) if mask & 1 << bit))
                occupied: set[Any] = set()
                data_used: set[int] = set()
                for index in indices:
                    support = footprints[index][state[index]]
                    data = set(pipelines[index][state[index]]["data_qubits"])
                    if occupied & support or data_used & data:
                        break
                    occupied.update(support)
                    data_used.update(data)
                else:
                    yield (pipelines[indices[0]][state[indices[0]]]["duration"], indices)

    def solve(state: tuple[int, ...]) -> int:
        if state in costs:
            return costs[state]
        if budget["remaining"] == 0:
            raise _BudgetExceeded
        budget["remaining"] -= 1
        budget["visited"] += 1
        if state == ends:
            costs[state] = 0
            return 0
        best: tuple[int, tuple[int, ...]] | None = None
        for duration, indices in moves(state):
            candidate = (duration + solve(advance(state, indices)), indices)
            if best is None or candidate < best:
                best = candidate
        if best is None:
            raise ValueError("admitted local pipeline has no legal ready operation")
        costs[state], choices[state] = best
        return best[0]

    state = (0,) * len(pipelines)
    solve(state)
    batches = []
    while state != ends:
        indices = choices[state]
        batches.append([pipelines[index][state[index]]["id"] for index in indices])
        state = advance(state, indices)
    return batches


def schedule_operations(
    operations: Sequence[Mapping[str, Any]],
    logical_gadgets: Sequence[Mapping[str, Any]],
    *,
    node_budget: int = 250000,
    parallel_h: bool = False,
) -> dict[str, Any]:
    """Return stages/evidence without mutating the canonical native operations.

    ``node_budget`` bounds newly visited DP states across the complete call,
    not separately for every window.  A failed window is emitted serially;
    all later windows can still be identified and safely emitted without
    exceeding that budget.  Existing operation timestamps are ignored.

    ``parallel_h=False`` retains the v1 serialized-H behavior.  Opting in
    selects v2: only equal-duration H operations with disjoint complete 2x2
    footprints can run together.  Standalone H runs also preserve same-wire
    and declared same-gadget order and never cross even a zero-time barrier.
    """
    if type(node_budget) is not int or node_budget < 0:
        raise ValueError("node_budget must be a nonnegative integer")
    if type(parallel_h) is not bool:
        raise ValueError("parallel_h must be a boolean")
    by_id = {}
    for operation in operations:
        identifier = operation.get("id")
        if not isinstance(identifier, str) or not identifier or identifier in by_id:
            raise ValueError("native operation IDs must be nonempty and unique")
        if not isinstance(operation.get("kind"), str):
            raise ValueError("native operation kind must be declared")
        if type(operation.get("duration")) is not int or operation["duration"] < 0:
            raise ValueError("native durations must be nonnegative integers")
        by_id[identifier] = operation
    records = {record["source_id"]: record for record in logical_gadgets}
    if len(records) != len(logical_gadgets):
        raise ValueError("logical gadget source IDs must be unique")
    prefixes = sorted(records, key=lambda source_id: (-len(source_id), source_id))
    source_of = {
        identifier: next(
            (source_id for source_id in prefixes if identifier.startswith(source_id + ":")), None
        )
        for identifier in by_id
    }
    stages: list[dict[str, Any]] = []
    windows: list[dict[str, Any]] = []
    ready: set[str] = set()
    budget = {"remaining": node_budget, "visited": 0}
    now = 0

    def emit(batch: list[str]) -> None:
        nonlocal now
        duration = by_id[batch[0]]["duration"]
        stages.append({"start": now, "end": now + duration, "operation_ids": batch})
        now += duration
        for identifier in batch:
            ready.update(_outputs(by_id[identifier]))

    def window(start: int, end: int, source_ids: list[str], kind: str, identifier: str) -> None:
        original = operations[start:end]
        pipelines = [
            [operation for operation in original if source_of[operation["id"]] == source_id]
            for source_id in source_ids
        ]
        stage_start, time_start, nodes_before = (len(stages), now, budget["visited"])
        reason = _window_reason(pipelines, source_ids, records, ready, kind, parallel_h=parallel_h)
        budget_unknown = False
        batches = None
        if reason is None:
            try:
                batches = _solve_window(pipelines, budget, parallel_h=parallel_h)
            except _BudgetExceeded:
                reason, budget_unknown = ("node_budget_exhausted", True)
        if batches is None:
            batches = [[operation["id"]] for operation in original]
        for batch in batches:
            emit(batch)
        serial_latency = sum((operation["duration"] for operation in original))
        windows.append(
            {
                "id": identifier,
                "kind": kind,
                "source_ids": source_ids,
                "operation_ids": [operation["id"] for operation in original],
                "pipelines": [
                    {
                        "source_id": source_id,
                        "operation_ids": [operation["id"] for operation in pipeline],
                    }
                    for source_id, pipeline in zip(source_ids, pipelines, strict=True)
                ],
                "stage_range": [stage_start, len(stages)],
                "start": time_start,
                "end": now,
                "serial_latency": serial_latency,
                "latency": now - time_start,
                "saved_cycles": serial_latency - (now - time_start),
                "status": "EXACT_LOCAL" if reason is None else "UNKNOWN",
                "local_minimum_proved": reason is None,
                "budget_unknown": budget_unknown,
                "reason": reason,
                "nodes": budget["visited"] - nodes_before,
                "proof_scope": (
                    "fixed window, full per-gadget order, stored paths, disjoint complete H footprints, singleton explicit Z and no mixed kinds; no global optimality claim"
                    if parallel_h
                    else "fixed window, full per-gadget order, stored paths, serial singleton H and explicit Z; no global optimality claim"
                ),
            }
        )

    def hadamard_run(start: int, end: int) -> None:
        original = operations[start:end]
        wires = list(dict.fromkeys((operation.get("wire") for operation in original)))
        pipelines = [
            [operation for operation in original if operation.get("wire") == wire] for wire in wires
        ]
        dependencies: dict[str, list[str]] = {}
        last_owner = {}
        for operation in original:
            owner = source_of[operation["id"]]
            if owner is not None:
                if owner in last_owner:
                    dependencies[operation["id"]] = [last_owner[owner]]
                last_owner[owner] = operation["id"]
        stage_start, time_start, nodes_before = (len(stages), now, budget["visited"])
        reason = None
        if not 1 <= len(pipelines) <= 5 or any((len(pipeline) > 16 for pipeline in pipelines)):
            reason = "unsupported_H_run_pipeline_count"
        elif any((_h_patch_footprint(operation) is None for operation in original)):
            reason = "missing_or_invalid_complete_H_footprint"
        batches, budget_unknown = (None, False)
        if reason is None:
            try:
                batches = _solve_window(
                    pipelines, budget, parallel_h=True, dependencies=dependencies
                )
            except _BudgetExceeded:
                reason, budget_unknown = ("node_budget_exhausted", True)
        if batches is None:
            batches = [[operation["id"]] for operation in original]
        for batch in batches:
            emit(batch)
        serial_latency = sum((operation["duration"] for operation in original))
        windows.append(
            {
                "id": "Hrun:" + original[0]["id"],
                "kind": "hadamard_run",
                "source_ids": list(
                    dict.fromkeys(
                        (
                            source_of[operation["id"]]
                            for operation in original
                            if source_of[operation["id"]] is not None
                        )
                    )
                ),
                "operation_ids": [operation["id"] for operation in original],
                "pipelines": [
                    {"wire": wire, "operation_ids": [operation["id"] for operation in pipeline]}
                    for wire, pipeline in zip(wires, pipelines, strict=True)
                ],
                "owner_dependencies": dependencies,
                "stage_range": [stage_start, len(stages)],
                "start": time_start,
                "end": now,
                "serial_latency": serial_latency,
                "latency": now - time_start,
                "saved_cycles": serial_latency - (now - time_start),
                "status": "EXACT_LOCAL" if reason is None else "UNKNOWN",
                "local_minimum_proved": reason is None,
                "budget_unknown": budget_unknown,
                "reason": reason,
                "nodes": budget["visited"] - nodes_before,
                "proof_scope": "fixed maximal consecutive H run, original same-wire/same-owner order, complete 2x2 patch conflicts; no barrier crossing or global optimum claim",
            }
        )

    index = 0
    while index < len(operations):
        operation = operations[index]
        if operation["kind"] == "native_ccz_bell":
            emit([operation["id"]])
            index += 1
            source_ids = list(
                dict.fromkeys((gadget["source_id"] for gadget in operation.get("gadgets", ())))
            )
            allowed = set(source_ids)
            end = index
            while (
                end < len(operations)
                and operations[end]["kind"] in _CCZ_KINDS
                and (source_of[operations[end]["id"]] in allowed)
            ):
                end += 1
            if end > index:
                window(index, end, source_ids, "ccz_corrections", "ccz:" + operation["id"])
                index = end
        elif operation["kind"] == "qand_measure_x" and source_of[operation["id"]] is not None:
            source_id = source_of[operation["id"]]
            segment = source_id.rsplit(":", 1)[0]
            end, source_ids = (index, [])
            while end < len(operations):
                candidate = operations[end]
                candidate_source = source_of[candidate["id"]]
                if (
                    candidate["kind"] not in _QAND_KINDS
                    or candidate_source is None
                    or candidate_source.rsplit(":", 1)[0] != segment
                    or (records[candidate_source].get("role") != "qand_dagger")
                ):
                    break
                if candidate_source not in source_ids:
                    source_ids.append(candidate_source)
                end += 1
            if end > index:
                window(index, end, source_ids, "qand_daggers", "qand:" + operation["id"])
                index = end
            else:
                emit([operation["id"]])
                index += 1
        elif parallel_h and operation["kind"] == "h":
            end = index + 1
            while end < len(operations) and operations[end]["kind"] == "h":
                end += 1
            hadamard_run(index, end)
            index = end
        else:
            emit([operation["id"]])
            index += 1
    emitted = [identifier for stage in stages for identifier in stage["operation_ids"]]
    if Counter(emitted) != Counter(by_id.keys()):
        raise ValueError("event schedule did not preserve every operation exactly once")
    serial_latency = sum((operation["duration"] for operation in operations))
    result = {
        "schema": PARALLEL_H_SCHEMA if parallel_h else SCHEMA,
        "stages": stages,
        "latency": now,
        "windows": windows,
        "stats": {
            "serial_latency": serial_latency,
            "saved_cycles": serial_latency - now,
            "operations": len(operations),
            "stages": len(stages),
            "parallel_stages": sum((len(stage["operation_ids"]) > 1 for stage in stages)),
            "H_singleton_stages": sum(
                (
                    by_id[stage["operation_ids"][0]]["kind"] == "h"
                    and len(stage["operation_ids"]) == 1
                    for stage in stages
                )
            ),
            "windows": len(windows),
            "exact_windows": sum((item["local_minimum_proved"] for item in windows)),
            "unknown_windows": sum((item["status"] == "UNKNOWN" for item in windows)),
            "node_budget": node_budget,
            "dp_nodes": budget["visited"],
        },
        "scope": "fixed canonical native operations; local CCZ correction/QAND erasure retiming only; not a global routing or scheduling optimum",
    }
    if parallel_h:
        h_stages = [stage for stage in stages if by_id[stage["operation_ids"][0]]["kind"] == "h"]
        result["parallel_h"] = True
        result["H_policy"] = (
            "same-kind complete-2x2-footprint VDP; no H/non-H overlap; canonical same-wire/per-gadget order"
        )
        result["stats"].update(
            H_stages=len(h_stages),
            H_cycles=sum((stage["end"] - stage["start"] for stage in h_stages)),
            H_parallel_stages=sum((len(stage["operation_ids"]) > 1 for stage in h_stages)),
            H_operations=sum((len(stage["operation_ids"]) for stage in h_stages)),
        )
        result["scope"] = (
            "fixed canonical native operations; local CCZ correction/QAND erasure and maximal contiguous H-run retiming with full native patch footprints; not a global routing or scheduling optimum"
        )
    return result
