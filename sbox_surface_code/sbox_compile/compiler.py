"""Complete physical reconstruction from matrices, routes, and corrections."""

from __future__ import annotations

import copy

from . import geometry as front, native
from .cost import scheduled_cost
from .event import schedule_operations
from .lower import _dagger_operations, _group_operations, routes_from_artifact
from .model import SCHEMA
from .schedule import schedule_causal
from .verify import check_interfaces, require_verified


def set_schedule(circuit: dict, stages: list[dict]) -> dict:
    result = copy.deepcopy(circuit)
    by_id = {op["id"]: op for op in result["operations"]}
    for stage in stages:
        for sid in stage["operation_ids"]:
            by_id[sid]["start"], by_id[sid]["end"] = stage["start"], stage["end"]
    result["schedule"] = {"stages": copy.deepcopy(stages)}
    result["cost"] = scheduled_cost(result["operations"], stages, parallel_h=True)
    result["space_time"]["reserved_patch_cycles"] = (
        result["space_time"]["reserved_module_patches"] * result["cost"]["latency"]
    )
    return result


def compile_circuit(
    source: dict,
    *,
    routes=None,
    corrections=None,
    stages=None,
    schedule_nodes: int = 25000,
    seed: int = 0,
) -> dict:
    """Lower the entire logical source; unchanged linear blocks are also replayed.

    ``stages`` is used only for deterministic witness replay. Search constructs
    fresh stages with bounded dynamic programming and the causal H-run solver.
    """
    check_interfaces(source)
    geom = native.NativeGeometry(
        source["row"], source["layout"]["data_cols"], source["logical_width"]
    )
    linears, partitions = native._native_linear_map(source, geom, optimize=True)
    records = {r["source_id"]: r for r in source["logical_gadgets"]}
    proofs = {r["source_id"]: copy.deepcopy(r) for r in source["search"]["correction_windows"]}
    front._require(set(corrections or {}) <= set(proofs), "unknown correction override")
    proofs.update(copy.deepcopy(corrections or {}))
    route_map = routes_from_artifact(source) if routes is None else copy.deepcopy(routes)
    operations, windows, used_ports = [], [], set()
    for segment in source["logical_schedule"]:
        name, kind = segment["segment"], segment["kind"]
        if kind == "linear":
            operations.extend(linears[name])
        elif kind == "relabel":
            operations.append(
                {
                    "id": f"{name}:relabel",
                    "kind": "relabel",
                    "duration": 0,
                    "category": "relabel",
                    "data_qubits": [],
                    "input_placement": segment["input_placement"],
                    "output_placement": segment["output_placement"],
                    "meaning": "logical wire renaming, not a Pauli or Clifford frame",
                }
            )
        elif kind == "affine_x":
            for i, gate in enumerate(segment["mapped_gates"]):
                physical = gate["physical_indices"]
                operations.append(
                    {
                        "id": f"{name}:X:{i}",
                        "kind": "x",
                        "wire": physical[0],
                        "data_qubits": physical,
                        "duration": 0,
                        "category": "pauli",
                        "execution": "explicit_physical_Pauli_not_frame",
                    }
                )
        elif kind == "ccz":
            expected = {
                f"{name}:{i}"
                for i in range(len(segment["mapped_gates"]))
                if records[f"{name}:{i}"]["role"] != "qand_dagger"
            }
            groups = route_map.get(name, [])
            actual = [item["source_id"] for group in groups for item in group]
            front._require(
                len(actual) == len(set(actual)) and set(actual) == expected,
                "incomplete or duplicate CCZ cover",
            )
            for group_index, group in enumerate(groups):
                for item in group:
                    record = records[item["source_id"]]
                    front._require(
                        all(
                            item[key] == record[key]
                            for key in ("logical_indices", "physical_indices")
                        ),
                        "route support changed",
                    )
                operations.extend(
                    _group_operations(
                        geom,
                        records,
                        group,
                        proofs,
                        used_ports,
                        f"{name}:native_round:{group_index}",
                    )
                )
                windows.extend(proofs[item["source_id"]] for item in group)
            for i in range(len(segment["mapped_gates"])):
                sid = f"{name}:{i}"
                if records[sid]["role"] == "qand_dagger":
                    operations.extend(_dagger_operations(records[sid], proofs[sid], geom))
                    windows.append(proofs[sid])
        else:
            raise ValueError(f"unsupported source segment: {kind}")
    operations, cancellations = front.cancel_disjoint_h_pairs(operations)
    native._decorate_conditional_cnots(operations, geom)
    front._timeline(operations)
    result = copy.deepcopy(source)
    result.update(schema=SCHEMA, operations=operations)
    result["search"] = {
        "optimized": True,
        "linear_VDP_partitions": partitions,
        "correction_windows": windows,
        "final_H_cancellations": cancellations,
    }
    if stages is None:
        fallback = schedule_operations(
            operations, source["logical_gadgets"], parallel_h=True, node_budget=schedule_nodes
        )
        stages = schedule_causal(
            operations,
            source["logical_gadgets"],
            fallback["stages"],
            node_budget=schedule_nodes,
            frontier_limit=16000,
            greedy_starts=12,
            seed=seed,
            cpu_seconds=None,
        )["stages"]
    result = set_schedule(result, stages)
    require_verified(result)
    return result
