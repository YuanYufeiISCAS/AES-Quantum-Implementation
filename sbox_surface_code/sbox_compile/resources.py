"""resources."""

# Adapted from the authors' strict_sbox/schedule_verify.py; campaign infrastructure removed.
from __future__ import annotations
from collections.abc import Mapping, Sequence
from typing import Any

_PRIMITIVES = {
    "linear_cnot_layer": (2, "linear"),
    "h": (3, "hadamard"),
    "conditional_cnot": (2, "correction"),
    "conditional_z": (0, "pauli"),
    "x": (0, "pauli"),
    "relabel": (0, "relabel"),
    "port_reset": (1, "port_reset"),
    "qand_measure_x": (1, "qand_measurement"),
    "qand_reset": (1, "qand_reset"),
    "native_ccz_data_to_port": (2, "ccz_injection"),
    "native_ccz_bell": (3, "ccz_injection"),
}
_BATCHABLE = {"conditional_cnot", "qand_measure_x", "qand_reset"}
_CCZ_CORRECTIONS = {"conditional_z", "h", "conditional_cnot"}
_DAGGER_OPERATIONS = {"qand_measure_x", "h", "conditional_cnot", "qand_reset"}
_SCOPE = "independent event-timeline, full-path VDP, explicit-outcome and serial-H check; requires separate canonical native semantic/geometry replay"
_PARALLEL_SCOPE = "independent event-timeline, full-path VDP, explicit-outcome and disjoint complete-2x2-H-footprint check; requires separate canonical native semantic/geometry replay to bind footprints to the actual layout"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and (not isinstance(value, (str, bytes)))


def _wire(value: Any) -> int:
    _require(type(value) is int and value >= 0, "physical data wire must be a nonnegative integer")
    return value


def _wire_list(value: Any) -> list[int]:
    _require(_sequence(value), "physical data support must be a sequence")
    result = [_wire(wire) for wire in value]
    _require(len(result) == len(set(result)), "duplicate physical data support")
    return result


def _path_vertices(raw: Any) -> set[tuple[int, int]]:
    _require(_sequence(raw) and len(raw) >= 2, "complete CNOT path is missing")
    path = []
    for coord in raw:
        _require(
            _sequence(coord) and len(coord) == 2 and all((type(v) is int for v in coord)),
            "CNOT path coordinate is malformed",
        )
        path.append((coord[0], coord[1]))
    _require(len(set(path)) == len(path), "CNOT path is not simple")
    _require(
        all((abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1 for a, b in zip(path, path[1:]))),
        "CNOT path is not grid-connected",
    )
    return set(path)


def _disjoint_paths(paths: Sequence[Any]) -> int:
    used: set[tuple[int, int]] = set()
    for raw in paths:
        vertices = _path_vertices(raw)
        _require(not used & vertices, "simultaneous CNOT paths are not full-path vertex-disjoint")
        used.update(vertices)
    return len(paths)


def _h_footprint(operation: Mapping[str, Any]) -> tuple[set[tuple[int, int]], tuple[int, int]]:
    """Check the complete declared H allocation, not a scheduler resource mask.

    Native replay separately checks this declaration against geometry.h_footprint;
    this interface has no layout bounds and does not pretend to supply them.
    """
    footprint = operation.get("footprint")
    _require(isinstance(footprint, Mapping), "H footprint declaration is missing")

    def coordinate(raw: Any) -> tuple[int, int]:
        _require(
            _sequence(raw) and len(raw) == 2 and all((type(value) is int for value in raw)),
            "H footprint coordinates must be integer grid pairs",
        )
        return (raw[0], raw[1])

    raw = footprint.get("reserved_patch_coords")
    _require(
        _sequence(raw) and len(raw) == 4,
        "H footprint must reserve all four patches of a complete 2x2 square",
    )
    vertices = {coordinate(point) for point in raw}
    rows, cols = ({point[0] for point in vertices}, {point[1] for point in vertices})
    _require(
        len(vertices) == 4
        and len(rows) == len(cols) == 2
        and (max(rows) - min(rows) == max(cols) - min(cols) == 1)
        and (vertices == {(row, col) for row in rows for col in cols}),
        "H footprint is not a complete adjacent 2x2 square",
    )
    data = coordinate(footprint.get("data_coord"))
    _require(data in vertices, "H footprint does not include its declared data patch")
    auxiliary = footprint.get("auxiliary_coords")
    _require(
        _sequence(auxiliary)
        and len(auxiliary) == 3
        and ({coordinate(point) for point in auxiliary} == vertices - {data}),
        "H footprint auxiliary allocation must be exactly the other three patches",
    )
    return (vertices, data)


def _condition_variables(expression: Any) -> set[str]:
    if isinstance(expression, str):
        _require(bool(expression), "condition outcome name is empty")
        return {expression}
    _require(
        isinstance(expression, Mapping) and len(expression) == 1,
        "condition expression is malformed",
    )
    key = next(iter(expression))
    _require(
        key in {"xor", "and"} and _sequence(expression[key]) and (len(expression[key]) >= 2),
        "condition expression must use explicit xor/and operands",
    )
    result: set[str] = set()
    for operand in expression[key]:
        result.update(_condition_variables(operand))
    return result


def _owner(identifier: str, records: Mapping[str, Any]) -> str | None:
    matches = [source_id for source_id in records if identifier.startswith(source_id + ":")]
    _require(len(matches) <= 1, f"ambiguous operation source ID: {identifier}")
    return matches[0] if matches else None


def _support(operation: Mapping[str, Any], records: Mapping[str, Any]) -> set[int]:
    kind = operation["kind"]
    if kind in {"h", "x", "conditional_z", "qand_measure_x", "qand_reset"}:
        expected = {_wire(operation["wire"])}
    elif kind == "conditional_cnot":
        expected = {_wire(operation["control"]), _wire(operation["target"])}
        _require(len(expected) == 2, "conditional CNOT repeats a data endpoint")
    elif kind == "linear_cnot_layer":
        children = operation["operations"]
        _require(_sequence(children) and bool(children), "empty native linear layer")
        expected = {_wire(child[key]) for child in children for key in ("control", "target")}
    elif kind in {"native_ccz_data_to_port", "native_ccz_bell"}:
        gadgets = operation["gadgets"]
        _require(_sequence(gadgets) and bool(gadgets), "empty native CCZ group")
        source_ids = [gadget["source_id"] for gadget in gadgets]
        _require(
            all((isinstance(sid, str) and sid in records for sid in source_ids)),
            "native CCZ group contains an undeclared gadget",
        )
        _require(len(source_ids) == len(set(source_ids)), "native CCZ group repeats a gadget")
        triples = [wire for sid in source_ids for wire in records[sid]["physical_indices"]]
        _require(len(triples) == len(set(triples)), "native CCZ group has overlapping data triples")
        expected = set(triples)
    else:
        expected = set()
    actual = set(_wire_list(operation["data_qubits"]))
    _require(actual == expected, f"operation data support does not replay: {operation['id']}")
    return expected


def _derive_regions(
    operations: Sequence[Mapping[str, Any]],
    records: Mapping[str, Any],
    owners: Mapping[str, str | None],
    *,
    parallel_h: bool,
) -> list[dict[str, Any]]:
    """Derive allowed regions from the immutable canonical stream, not stages."""
    regions: list[dict[str, Any]] = []

    def append(kind: str, begin: int, end: int) -> None:
        selected = operations[begin:end]
        pipelines: dict[str, list[str]] = {}
        if kind == "h_run":
            for operation in selected:
                pipelines.setdefault(f"wire:{operation['wire']}", []).append(operation["id"])
                owner = owners[operation["id"]]
                if owner is not None:
                    pipelines.setdefault(f"owner:{owner}", []).append(operation["id"])
        elif kind != "barrier":
            for operation in selected:
                owner = owners[operation["id"]]
                _require(owner is not None, "reorderable operation has no declared gadget")
                pipelines.setdefault(owner, []).append(operation["id"])
            used: set[int] = set()
            for source_id in pipelines:
                triple = set(records[source_id]["physical_indices"])
                _require(
                    not used & triple, "concurrent pipelines have overlapping physical data triples"
                )
                used.update(triple)
        regions.append({"kind": kind, "ids": [op["id"] for op in selected], "pipelines": pipelines})

    index = 0
    while index < len(operations):
        operation = operations[index]
        if operation["kind"] == "native_ccz_bell":
            group = {gadget["source_id"] for gadget in operation["gadgets"]}
            _require(
                all((records[sid]["role"] != "qand_dagger" for sid in group)),
                "QAND dagger cannot appear in a CCZ Bell group",
            )
            append("barrier", index, index + 1)
            index += 1
            begin = index
            while index < len(operations):
                candidate = operations[index]
                if (
                    candidate["kind"] not in _CCZ_CORRECTIONS
                    or owners[candidate["id"]] not in group
                ):
                    break
                _require(
                    not candidate["id"].endswith(":entry_H"),
                    "CCZ entry H occurs inside post-Bell corrections",
                )
                index += 1
            _require(index > begin, "native Bell group has no explicit correction pipeline")
            _require(
                {owners[op["id"]] for op in operations[begin:index]} == group,
                "native Bell group is missing a correction pipeline",
            )
            append("ccz_correction", begin, index)
        elif operation["kind"] == "qand_measure_x":
            owner = owners[operation["id"]]
            _require(
                owner is not None and records[owner]["role"] == "qand_dagger",
                "QAND measurement has no matching dagger declaration",
            )
            segment = owner.rsplit(":", 1)[0]
            begin = index
            while index < len(operations):
                candidate = operations[index]
                source_id = owners[candidate["id"]]
                if (
                    candidate["kind"] not in _DAGGER_OPERATIONS
                    or source_id is None
                    or records[source_id]["role"] != "qand_dagger"
                    or (source_id.rsplit(":", 1)[0] != segment)
                ):
                    break
                index += 1
            append("qand_dagger", begin, index)
        elif parallel_h and operation["kind"] == "h":
            begin = index
            while index < len(operations) and operations[index]["kind"] == "h":
                index += 1
            append("h_run", begin, index)
        else:
            _require(
                operation["kind"] not in {"conditional_z", "conditional_cnot", "qand_reset"},
                "correction/reset occurs outside its canonical allowed window",
            )
            append("barrier", index, index + 1)
            index += 1
    return regions
