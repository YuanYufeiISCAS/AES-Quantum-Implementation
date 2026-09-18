"""Independent internal placements paid by complete neighboring matrices."""

from __future__ import annotations

import copy

from . import geometry as front, linear, native
from .compiler import compile_circuit
from .lower import routes_from_artifact
from .matrix import (
    invert_perm,
    linear_target_rows_from_span,
    replay_rows,
    update_placement_for_relabel,
)
from .repair import _repair_injection_group
from .semantics import segment_name


def geometry(circuit):
    return native.NativeGeometry(
        circuit["row"], circuit["layout"]["data_cols"], circuit["logical_width"]
    )


def move_placement(placement, old_to_new):
    front._require(
        sorted(old_to_new) == list(range(len(placement))), "invalid boundary permutation"
    )
    result = [None] * len(placement)
    for old, new in enumerate(old_to_new):
        result[new] = placement[old]
    return result


def task(circuit, index, pin, pout):
    """Construct M' = Pout M Pin^-1 and a paid constructive warm start."""
    source = circuit["linear_candidates"][index]
    old = source["choice"]
    choice = copy.deepcopy(old)
    choice["input_placement"] = move_placement(old["input_placement"], pin)
    choice["output_placement"] = move_placement(old["output_placement"], pout)
    target, middle = linear_target_rows_from_span(
        circuit["logical_source"]["segments"],
        old["start_pos"],
        old["end_pos"],
        tuple(choice["input_placement"]),
    )
    inverse = invert_perm(middle)
    choice["output_permutation"] = [inverse[q] for q in choice["output_placement"]]
    rows = [target[i] for i in choice["output_permutation"]]
    n = circuit["logical_width"]
    # Constructive row permutations are explicitly charged; they are only a seed.
    gates = (
        linear.gaussian([1 << pin[i] for i in range(n)])
        + linear.word(source["candidate"])
        + linear.gaussian([1 << i for i in invert_perm(pout)])
    )
    reference = linear.routed_reference(rows, geometry(circuit), gates)
    front._require(replay_rows(reference) == rows, "paid reference matrix mismatch")
    return {
        "index": index,
        "choice": choice,
        "rows": rows,
        "reference": reference,
        "pin": list(pin),
        "pout": list(pout),
    }


def bind(candidate, matrix_task, geom):
    """A free-Q offer is admitted only after its physical matrix matches the task."""
    result = copy.deepcopy(candidate)
    result["output_permutation"] = list(range(len(matrix_task["rows"])))
    linear.validate(result, matrix_task["rows"], geom)
    return result


def remap_correction(old_record, record, proof, geom):
    sid, physical = record["source_id"], record["physical_indices"]
    dagger = record["role"] == "qand_dagger"
    specs = {s["id"]: s for s in native._cz_specs(sid, physical, dagger)}
    mapping = dict(zip(old_record["physical_indices"], physical, strict=True))
    chosen, operations = [], []
    for old in proof["chosen_czs"]:
        edge = {**specs[old["id"]], "target": mapping[old["target"]]}
        if edge["pair"] == old["pair"] and edge["target"] == old["target"] and "path" in old:
            edge["path"] = copy.deepcopy(old["path"])
        expansion = front._cz_expansion(
            tuple(edge["pair"]),
            edge["condition"],
            edge["id"],
            edge["target"],
            geom,
            path=edge.get("path"),
        )
        edge["path"] = next(op["path"] for op in expansion if op["kind"] == "conditional_cnot")
        chosen.append(edge)
        operations.extend(expansion)
    if not dagger:
        operations.append(front._h(physical[2], f"{sid}:trailing_H:{physical[2]}", geom))
    _, cancellations = front.cancel_disjoint_h_pairs(operations)
    return {"source_id": sid, "chosen_czs": chosen, "cancellations": cancellations}


def rebuild_schedule(circuit):
    choices = {r["choice"]["start_pos"]: r["choice"] for r in circuit["linear_candidates"]}
    source = circuit["logical_source"]["segments"]
    placement = list(circuit["initial_placement"])
    schedule, position = [], 0
    while position < len(source):
        segment = source[position]
        name = segment_name(segment)
        if position in choices:
            choice = choices[position]
            front._require(choice["input_placement"] == placement, "unpaid boundary discontinuity")
            schedule.append(
                {
                    "kind": "linear",
                    "segment": name,
                    "candidate": choice["candidate"],
                    "cost": choice["cost"],
                    "covered_segments": choice["covered_segments"],
                    "span": [choice["start_pos"], choice["end_pos"]],
                    **{
                        k: choice[k]
                        for k in ("input_placement", "output_placement", "output_permutation")
                    },
                }
            )
            placement = list(choice["output_placement"])
            position = choice["end_pos"] + 1
            continue
        if segment["kind"] == "relabel":
            new = list(update_placement_for_relabel(segment, tuple(placement)))
            schedule.append(
                {
                    "kind": "relabel",
                    "segment": name,
                    "input_placement": placement,
                    "output_placement": new,
                }
            )
            placement = new
        else:
            inverse = invert_perm(placement)
            schedule.append(
                {
                    "kind": segment["kind"],
                    "segment": name,
                    "placement": placement,
                    "mapped_gates": [
                        {
                            "kind": gate["kind"],
                            "logical_indices": gate["indices"],
                            "physical_indices": [inverse[q] for q in gate["indices"]],
                        }
                        for gate in segment["gates"]
                    ],
                }
            )
        position += 1
    front._require(placement == circuit["final_placement"], "external placement changed")
    circuit["logical_schedule"] = schedule


def materialize(circuit, start, boundaries, candidates, *, schedule_nodes=25000):
    n = circuit["logical_width"]
    identity = list(range(n))
    front._require(
        boundaries[0] == boundaries[-1] == identity, "window external boundaries must be fixed"
    )
    front._require(len(boundaries) == len(candidates) + 1, "incomplete matrix cover")
    result = copy.deepcopy(circuit)
    geom = geometry(circuit)
    for offset, candidate in enumerate(candidates):
        index = start + offset
        matrix_task = task(circuit, index, boundaries[offset], boundaries[offset + 1])
        c = bind(candidate, matrix_task, geom)
        choice = matrix_task["choice"]
        c["output_permutation"] = choice["output_permutation"]
        c["verified"] = True
        choice["cost"] = [c["stats"][k] for k in ("surface_depth", "layers", "cnots")]
        proof = front.validate_candidate(c, choice, result["logical_source"]["segments"], geom)
        result["linear_candidates"][index] = {
            "candidate": c,
            "choice": choice,
            "choice_index": index,
            "matrix_geometry_replay": {
                "passed": True,
                **{k: proof[k] for k in ("cost", "target_rows_hex", "source_coordinate_offset")},
            },
        }
    rebuild_schedule(result)
    mapped = {
        f"{s['segment']}:{i}": gate
        for s in result["logical_schedule"]
        if s["kind"] == "ccz"
        for i, gate in enumerate(s["mapped_gates"])
    }
    records = {r["source_id"]: r for r in result["logical_gadgets"]}
    old_records = {r["source_id"]: r for r in circuit["logical_gadgets"]}
    proofs = {r["source_id"]: r for r in circuit["search"]["correction_windows"]}
    corrections, changed = {}, set()
    for sid, record in records.items():
        record["physical_indices"] = mapped[sid]["physical_indices"]
        if record["physical_indices"] != old_records[sid]["physical_indices"]:
            changed.add(sid.rsplit(":", 1)[0])
            corrections[sid] = remap_correction(old_records[sid], record, proofs[sid], geom)
    routes = routes_from_artifact(circuit)
    for name in changed & set(routes):
        groups = []
        for group in routes[name]:
            modified = copy.deepcopy(group)
            for option in modified:
                option["physical_indices"] = records[option["source_id"]]["physical_indices"]
                option["data_coords"] = [list(geom.coord(q)) for q in option["physical_indices"]]
            repaired, _ = _repair_injection_group(geom, group, modified)
            groups.extend(repaired)
        routes[name] = groups
    return compile_circuit(
        result, routes=routes, corrections=corrections, schedule_nodes=schedule_nodes
    )
