"""dependencies."""

# Adapted from the authors' strict_sbox/segment_verify.py; campaign infrastructure removed.
from __future__ import annotations
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from typing import Any
from . import resources as primitive_checks

SCHEMA = "strict-source-segment-constraints-v1"
BATCHABLE_KINDS = frozenset({"h", "conditional_cnot", "qand_measure_x", "qand_reset"})
_INJECTION = {"native_ccz_data_to_port", "native_ccz_bell"}
_BARRIERS = {"linear_cnot_layer", "x", "relabel"}
_PIPELINE = {"h", "conditional_z", "conditional_cnot", "qand_measure_x", "qand_reset"}
_SCOPE = "independently derived same-source-segment gadget/data/outcome DAG, canonical singleton injection order, Bell-END port/resource release with explicit reset-before-reuse, complete H/CNOT footprints; requires separate canonical native branch/geometry/AES replay"


def _require(value, message):
    if not value:
        raise ValueError(message)


def _coord(value):
    _require(
        isinstance(value, (list, tuple))
        and len(value) == 2
        and all((type(component) is int for component in value)),
        "invalid patch coordinate",
    )
    return tuple(value)


def _coords(value, count=None):
    _require(isinstance(value, (list, tuple)), "patch allocation must be a sequence")
    result = [_coord(item) for item in value]
    _require(
        len(set(result)) == len(result) and (count is None or len(result) == count),
        "duplicate or incomplete patch allocation",
    )
    return result


def derive_segment_constraints(
    operations: Sequence[Mapping[str, Any]], logical_gadgets: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """Derive JSON-serializable constraints; ignore operation timestamps.

    Dependencies can span regions (for example an earlier segment's Bell and a
    later segment's port reset). A region-wise scheduler treats predecessors in
    already completed regions as satisfied. Returned masks are conveniences for
    search; verify_segment_timeline always rederives them from the operations.
    """
    _require(
        primitive_checks._sequence(operations) and primitive_checks._sequence(logical_gadgets),
        "canonical operations and gadgets must be sequences",
    )
    records, segment_records = ({}, defaultdict(list))
    for record in logical_gadgets:
        _require(isinstance(record, Mapping), "gadget declaration must be an object")
        sid = record["source_id"]
        _require(
            isinstance(sid, str) and ":" in sid and (sid not in records),
            "duplicate or malformed gadget source ID",
        )
        _require(record["role"] in {"qand", "qand_dagger", "toffoli"}, "unknown gadget role")
        support = primitive_checks._wire_list(record["physical_indices"])
        _require(len(support) == 3, "gadget must have three distinct physical wires")
        records[sid] = record
        segment_records[sid.rsplit(":", 1)[0]].append(sid)
    for segment, ids in segment_records.items():
        support = [q for sid in ids for q in records[sid]["physical_indices"]]
        _require(
            len(support) == len(set(support)),
            f"source nonlinear segment has overlapping gadget supports: {segment}",
        )
    by_id, positions, supports, owners, footprints = ({}, {}, {}, {}, {})
    h_data, h_wires = ({}, {})
    for position, operation in enumerate(operations):
        _require(isinstance(operation, Mapping), "canonical operation must be an object")
        identifier, kind = (operation["id"], operation["kind"])
        _require(
            isinstance(identifier, str) and identifier and (identifier not in by_id),
            "duplicate/malformed canonical operation ID",
        )
        _require(kind in primitive_checks._PRIMITIVES, f"unknown native primitive: {kind}")
        duration, category = primitive_checks._PRIMITIVES[kind]
        _require(
            type(operation["duration"]) is int
            and operation["duration"] == duration
            and (operation["category"] == category),
            "native primitive duration/category mismatch",
        )
        by_id[identifier], positions[identifier] = (operation, position)
        supports[identifier] = primitive_checks._support(operation, records)
        owner = primitive_checks._owner(identifier, records)
        owners[identifier] = [owner] if owner is not None else []
        footprint = set()
        if kind == "h":
            footprint, data = primitive_checks._h_footprint(operation)
            wire = operation["wire"]
            _require(
                h_data.setdefault(wire, data) == data and h_wires.setdefault(data, wire) == wire,
                "H data-patch/wire anchoring is inconsistent",
            )
        elif kind == "conditional_cnot":
            footprint = primitive_checks._path_vertices(operation["path"])
        elif kind == "linear_cnot_layer":
            primitive_checks._disjoint_paths([child["path"] for child in operation["operations"]])
            footprint = {
                point
                for child in operation["operations"]
                for point in primitive_checks._path_vertices(child["path"])
            }
        elif kind == "port_reset":
            footprint = set(_coords(operation["ports"]))
            _require(
                bool(footprint) and operation["prepared_state"] == "|0>",
                "port reset has no real zero-state allocation",
            )
        if kind in _PIPELINE:
            _require(
                owner is not None,
                f"nonlinear pipeline operation has no declared owner: {identifier}",
            )
            _require(
                supports[identifier] <= set(records[owner]["physical_indices"]),
                "pipeline touches another gadget's data",
            )
        footprints[identifier] = footprint
    dependencies = set()

    def edge(before, after):
        _require(
            before in by_id and after in by_id and (positions[before] < positions[after]),
            f"invalid/reversed canonical dependency: {before} -> {after}",
        )
        dependencies.add((before, after))

    groups, group_for_reset = ([], {})
    used_ports, last_port_bell, last_resource_bell = (set(), {}, {})
    injected_sources = Counter()
    for operation in operations:
        if operation["kind"] != "native_ccz_data_to_port":
            continue
        acquire = operation["id"]
        _require(acquire.endswith(":data_to_port"), "native acquire ID has no group suffix")
        group = acquire.removesuffix(":data_to_port")
        bell_id, reset_id = (group + ":Bell", group + ":port_reset")
        _require(
            bell_id in by_id and by_id[bell_id]["kind"] == "native_ccz_bell",
            "injection group has no matching Bell",
        )
        bell = by_id[bell_id]
        source_ids = [item["source_id"] for item in operation["gadgets"]]
        _require(
            source_ids == [item["source_id"] for item in bell["gadgets"]],
            "acquire/Bell gadget cover differs",
        )
        _require(
            len({sid.rsplit(":", 1)[0] for sid in source_ids}) == 1,
            "injection group crosses source segments",
        )
        _require(
            all((records[sid]["role"] != "qand_dagger" for sid in source_ids)),
            "QAND dagger consumes a CCZ group",
        )
        injected_sources.update(source_ids)
        owners[acquire] = owners[bell_id] = source_ids
        ports, resources, paths = ([], [], [])
        for item in operation["gadgets"]:
            sid = item["source_id"]
            _require(
                item["physical_indices"] == records[sid]["physical_indices"],
                "acquire changes gadget physical operands",
            )
            local_ports, local_resources = (
                _coords(item["port_coords"], 3),
                _coords(item["resource_coords"], 3),
            )
            local_paths = item["data_to_port_cnot_paths"]
            _require(
                isinstance(local_paths, (list, tuple)) and len(local_paths) == 3,
                "native acquire needs three complete paths",
            )
            for path, port in zip(local_paths, local_ports, strict=True):
                primitive_checks._path_vertices(path)
                _require(_coord(path[-1]) == port, "acquire path does not end at its declared port")
            ports.extend(local_ports)
            resources.extend(local_resources)
            paths.extend(local_paths)
        _require(
            len(set(ports)) == len(ports)
            and len(set(resources)) == len(resources)
            and (not set(ports) & set(resources)),
            "injection group duplicates or aliases port/resource allocations",
        )
        primitive_checks._disjoint_paths(paths)
        path_tiles = {point for path in paths for point in primitive_checks._path_vertices(path)}
        _require(not path_tiles & set(resources), "data-to-port path crosses stored CCZ resources")
        footprints[acquire] = path_tiles | set(resources)
        footprints[bell_id] = set(ports) | set(resources)
        edge(acquire, bell_id)
        reused = set(ports) & used_ports
        reset = None
        if reused:
            _require(
                reset_id in by_id and by_id[reset_id]["kind"] == "port_reset",
                "port reuse is missing an explicit reset",
            )
            _require(
                set(_coords(by_id[reset_id]["ports"])) == reused,
                "reset does not cover exactly the reused ports",
            )
            reset = reset_id
            group_for_reset[reset_id] = source_ids
            for port in reused:
                edge(last_port_bell[port], reset_id)
            edge(reset_id, acquire)
        else:
            _require(reset_id not in by_id, "first-use ports have an unexplained reset")
        for resource in resources:
            if resource in last_resource_bell:
                edge(last_resource_bell[resource], acquire)
            last_resource_bell[resource] = bell_id
        for port in ports:
            last_port_bell[port] = bell_id
        used_ports.update(ports)
        groups.append(
            {
                "group": group,
                "acquire": acquire,
                "bell": bell_id,
                "ports": [list(p) for p in ports],
                "resources": [list(p) for p in resources],
                "reset": reset,
                "source_ids": source_ids,
            }
        )
    _require(
        {op["id"] for op in operations if op["kind"] == "native_ccz_bell"}
        == {group["bell"] for group in groups},
        "unpaired or repeated native Bell group",
    )
    for operation in operations:
        if operation["kind"] == "port_reset":
            _require(
                operation["id"] in group_for_reset,
                "port reset is not bound to the next native acquire",
            )
            owners[operation["id"]] = group_for_reset[operation["id"]]
    regions, operation_region = ([], {})
    for operation in operations:
        identifier, kind = (operation["id"], operation["kind"])
        source_ids = owners[identifier]
        if kind in _BARRIERS:
            regions.append(
                {
                    "kind": "barrier",
                    "segment": None,
                    "operation_ids": [identifier],
                    "source_ids": [],
                }
            )
        else:
            segments = {sid.rsplit(":", 1)[0] for sid in source_ids}
            _require(len(segments) == 1, "nonlinear operation is not bound to one source segment")
            segment = next(iter(segments))
            if (
                not regions
                or regions[-1]["kind"] != "nonlinear_segment"
                or regions[-1]["segment"] != segment
            ):
                regions.append(
                    {
                        "kind": "nonlinear_segment",
                        "segment": segment,
                        "operation_ids": [],
                        "source_ids": [],
                    }
                )
            regions[-1]["operation_ids"].append(identifier)
            for sid in source_ids:
                if sid not in regions[-1]["source_ids"]:
                    regions[-1]["source_ids"].append(sid)
        operation_region[identifier] = len(regions) - 1
    last_gadget, last_wire = ({}, {})
    for operation in operations:
        identifier = operation["id"]
        if operation["kind"] != "port_reset":
            for sid in owners[identifier]:
                if sid in last_gadget:
                    edge(last_gadget[sid], identifier)
                last_gadget[sid] = identifier
        for wire in supports[identifier]:
            if wire in last_wire:
                edge(last_wire[wire], identifier)
            last_wire[wire] = identifier
    injection_order = [op["id"] for op in operations if op["kind"] in _INJECTION]
    for before, after in zip(injection_order, injection_order[1:]):
        edge(before, after)
    producers, consumers = ({}, {})
    for operation in operations:
        identifier, kind = (operation["id"], operation["kind"])
        if kind == "native_ccz_bell":
            _require(
                operation["aggregate_outputs_ready_relative_cycle"] == 2
                and operation["corrections_consume_outputs_after_reserved_cycle"] == 3,
                "Bell outcome readiness contract changed",
            )
            outcomes = []
            for item in operation["gadgets"]:
                sid = item["source_id"]
                expected = {f"{sid}:{letter}{i}" for letter in ("k", "r") for i in range(3)}
                _require(
                    set(item["outcomes"]) == expected and len(item["outcomes"]) == 6,
                    "incomplete Bell outcomes",
                )
                outcomes.extend(item["outcomes"])
            relative = 2
        elif kind == "qand_measure_x":
            owner = owners[identifier][0]
            _require(
                records[owner]["role"] == "qand_dagger"
                and operation["wire"] == records[owner]["physical_indices"][2],
                "QAND measurement target/role mismatch",
            )
            _require(operation["outcome"] == owner + ":m", "QAND outcome name mismatch")
            outcomes, relative = ([operation["outcome"]], 1)
        else:
            continue
        for outcome in outcomes:
            _require(
                isinstance(outcome, str) and outcome not in producers, "duplicate outcome producer"
            )
            producers[outcome] = {
                "operation_id": identifier,
                "ready_relative_cycle": relative,
                "usable_after_relative_cycle": operation["duration"],
            }
    for operation in operations:
        identifier, kind = (operation["id"], operation["kind"])
        if kind == "conditional_cnot":
            _require(
                isinstance(operation["condition"], str)
                and operation["all_branch_reserved"] is True,
                "conditional CNOT lacks explicit all-branch reservation",
            )
            variables = {operation["condition"]}
        elif kind == "conditional_z":
            variables = primitive_checks._condition_variables(operation["condition"])
            declared = operation["condition_variables"]
            _require(
                primitive_checks._sequence(declared)
                and len(declared) == len(set(declared))
                and (set(declared) == variables),
                "Z outcome declaration differs from actual condition",
            )
            _require(
                operation["execution"] == "explicit_physical_Pauli_not_frame",
                "conditional Z was replaced by a deferred frame",
            )
        else:
            continue
        owner = owners[identifier][0]
        for variable in variables:
            _require(
                variable in producers and variable.startswith(owner + ":"),
                "conditional operation has no matching gadget outcome",
            )
            edge(producers[variable]["operation_id"], identifier)
        consumers[identifier] = sorted(variables)
    for sid, record in records.items():
        if record["role"] == "qand_dagger":
            _require(injected_sources[sid] == 0, "QAND dagger consumes a CCZ resource")
            for suffix, kind in (("measureX", "qand_measure_x"), ("reset", "qand_reset")):
                identifier = sid + ":" + suffix
                _require(
                    identifier in by_id and by_id[identifier]["kind"] == kind,
                    "missing explicit QAND measurement/reset",
                )
                _require(
                    by_id[identifier]["wire"] == record["physical_indices"][2],
                    "QAND reset/measurement target mismatch",
                )
        else:
            _require(
                injected_sources[sid] == 1, "gadget does not have exactly one native injection"
            )
            for i, wire in enumerate(record["physical_indices"]):
                identifier = f"{sid}:Z{i}"
                _require(
                    identifier in by_id
                    and by_id[identifier]["kind"] == "conditional_z"
                    and (by_id[identifier]["wire"] == wire),
                    "missing explicit zero-cycle CCZ Z correction",
                )
    return {
        "schema": SCHEMA,
        "regions": regions,
        "dependencies": [
            list(pair)
            for pair in sorted(
                dependencies, key=lambda pair: (positions[pair[0]], positions[pair[1]])
            )
        ],
        "operation_region": operation_region,
        "operation_supports": {identifier: sorted(value) for identifier, value in supports.items()},
        "operation_footprints": {
            identifier: [list(point) for point in sorted(value)]
            for identifier, value in footprints.items()
        },
        "source_ids_by_operation": owners,
        "batchable_kinds": sorted(BATCHABLE_KINDS),
        "port_lifetimes": groups,
        "producer_by_outcome": producers,
        "outcomes_by_consumer": consumers,
        "fixed_injection_order": injection_order,
        "scope": _SCOPE,
    }
