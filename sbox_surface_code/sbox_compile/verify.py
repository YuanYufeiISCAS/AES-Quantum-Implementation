"""Independent reconstruction of native primitives, execution, and AES action."""

from __future__ import annotations

import copy

from . import geometry, native
from .cost import scheduled_cost
from .matrix import invert_perm, update_placement_for_relabel
from .model import SCHEMA
from .semantics import segment_name
from .timeline import verify_causal_timeline


def check_interfaces(circuit: dict) -> None:
    """Bind every paid linear span and nonlinear operand to the inline source."""
    require = geometry._require
    width = circuit["logical_width"]
    placement = list(circuit["initial_placement"])
    require(sorted(placement) == list(range(width)), "invalid input placement")
    source = circuit["logical_source"]["segments"]
    choices = circuit["linear_candidates"]
    records = {item["source_id"]: item for item in circuit["logical_gadgets"]}
    require(len(records) == len(circuit["logical_gadgets"]), "duplicate logical gadget")
    position = 0
    linear_index = 0
    seen = set()
    for item in circuit["logical_schedule"]:
        require(position < len(source), "schedule exceeds source")
        if item["kind"] == "linear":
            record = choices[linear_index]
            choice = record["choice"]
            require(choice["start_pos"] == position, "linear span has a gap or overlap")
            require(choice["input_placement"] == placement, "unpaid input placement change")
            require(
                item["candidate"] == choice["candidate"]
                and item["segment"] == segment_name(source[position]),
                "linear schedule/candidate mismatch",
            )
            for key in ("input_placement", "output_placement", "output_permutation"):
                require(item[key] == choice[key], "linear schedule boundary mismatch")
            placement = list(choice["output_placement"])
            position = choice["end_pos"] + 1
            linear_index += 1
            continue
        segment = source[position]
        name = segment_name(segment)
        require(
            item["segment"] == name and item["kind"] == segment["kind"], "source segment mismatch"
        )
        if item["kind"] == "relabel":
            require(item["input_placement"] == placement, "relabel input mismatch")
            placement = list(update_placement_for_relabel(segment, tuple(placement)))
            require(item["output_placement"] == placement, "relabel output mismatch")
        else:
            require(item["placement"] == placement, "unpaid nonlinear placement change")
            inverse = invert_perm(placement)
            expected = [
                {
                    "kind": gate["kind"],
                    "logical_indices": gate["indices"],
                    "physical_indices": [inverse[q] for q in gate["indices"]],
                }
                for gate in segment["gates"]
            ]
            require(item["mapped_gates"] == expected, "mapped logical operations mismatch")
            if item["kind"] == "ccz":
                for i, gate in enumerate(expected):
                    sid = f"{name}:{i}"
                    require(
                        records[sid]["logical_indices"] == gate["logical_indices"]
                        and records[sid]["physical_indices"] == gate["physical_indices"],
                        "logical gadget support mismatch",
                    )
                    seen.add(sid)
        position += 1
    require(position == len(source) and linear_index == len(choices), "incomplete source cover")
    require(seen == set(records), "incomplete gadget cover")
    require(placement == circuit["final_placement"], "external output placement changed")


def verify(circuit: dict) -> dict:
    """Recompute all proofs; saved success flags have no authority."""
    try:
        if circuit.get("schema") == "sbox-surface-code-parallel-v1":
            from .parallel import verify as verify_parallel

            return verify_parallel(circuit)
        geometry._require(circuit["schema"] == SCHEMA, "unsupported submission circuit")
        check_interfaces(circuit)
        serial = copy.deepcopy(circuit)
        serial.update(schema=native.SCHEMA, model=native.MODEL)
        geometry._timeline(serial["operations"])
        serial["cost"] = native._native_cost(serial["operations"])
        serial["space_time"]["reserved_patch_cycles"] = (
            serial["space_time"]["reserved_module_patches"] * serial["cost"]["latency"]
        )
        physical = native.verify_native(serial)
        geometry._require(physical["passed"], f"native reconstruction: {physical.get('errors')}")
        temporal = verify_causal_timeline(
            circuit["operations"], circuit["logical_gadgets"], circuit["schedule"]["stages"]
        )
        geometry._require(temporal["passed"], f"causal timeline: {temporal.get('errors')}")
        cost = scheduled_cost(circuit["operations"], circuit["schedule"]["stages"], parallel_h=True)
        geometry._require(cost == circuit["cost"], "cost differs from reconstructed schedule")
        geometry._require(
            circuit["space_time"]["reserved_patch_cycles"]
            == circuit["space_time"]["reserved_module_patches"] * cost["latency"],
            "space-time cost mismatch",
        )
        aes = physical["lowering"]["actual_physical_AES"]
        return {
            "passed": True,
            "case": circuit["case"],
            "cost": cost,
            "checked_inputs": aes["checked_inputs"],
            "scratch_zero": aes["scratch_zero"],
            "checked_operations": physical["checked_operations"],
            "checked_stages": temporal["checked_stages"],
            "errors": [],
        }
    except (AssertionError, ValueError, KeyError, TypeError, IndexError, AttributeError, RuntimeError) as error:
        return {"passed": False, "errors": [str(error)]}


def require_verified(circuit: dict) -> dict:
    report = verify(circuit)
    if not report["passed"]:
        raise ValueError("; ".join(report["errors"]))
    return report
