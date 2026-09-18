"""timeline."""

# Adapted from the authors' strict_sbox/causal_verify.py; campaign infrastructure removed.
from __future__ import annotations
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
import json
from typing import Any
from . import resources as primitive_checks
from .dependencies import derive_segment_constraints

SCHEMA = "strict-causal-h-run-constraints-v1"
DAG_POLICY = "canonical-gadget-h-runs-disjoint-support-v1"
BATCHABLE_KINDS = frozenset({"h", "conditional_cnot", "qand_measure_x", "qand_reset"})
_SCOPE = "independent same-source-segment temporal replay with canonical-to-actual adjacent disjoint-support/full-footprint H swaps only inside gadget H runs; all non-H gadget boundaries, data/outcome dependencies, canonical singleton injection order, Bell-END resource release and explicit resets retained; requires separate canonical native branch/geometry/restoring-H/AES replay"


def _require(value, message):
    if not value:
        raise ValueError(message)


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def derive_causal_constraints(
    operations: Sequence[Mapping[str, Any]],
    logical_gadgets: Sequence[Mapping[str, Any]],
    *,
    dag_policy: str = DAG_POLICY,
) -> dict[str, Any]:
    """Parse canonical contracts and independently derive the H-run DAG.

    The returned legacy-compatible maps are conveniences for search, never
    trusted timing evidence.  ``verify_causal_timeline`` rederives all of them.
    Operation timestamps and scheduler-supplied dependency fields are ignored.
    """
    _require(type(dag_policy) is str and dag_policy == DAG_POLICY, "unknown causal DAG policy")
    parsed = derive_segment_constraints(operations, logical_gadgets)
    by_id = {op["id"]: op for op in operations}
    positions = {op["id"]: i for i, op in enumerate(operations)}
    supports = {sid: set(wires) for sid, wires in parsed["operation_supports"].items()}
    footprints = {
        sid: set(map(tuple, points)) for sid, points in parsed["operation_footprints"].items()
    }
    owners = parsed["source_ids_by_operation"]
    region_of = parsed["operation_region"]
    reasons = defaultdict(set)

    def edge(before, after, reason):
        _require(
            before in by_id and after in by_id and (positions[before] < positions[after]),
            f"invalid/reversed canonical causal dependency: {before} -> {after}",
        )
        reasons[before, after].add(reason)

    blocks = {record["source_id"]: [] for record in logical_gadgets}
    for op in operations:
        sid, kind = (op["id"], op["kind"])
        if kind == "port_reset":
            continue
        for owner in owners[sid]:
            pipeline = blocks[owner]
            if (
                kind == "h"
                and pipeline
                and (pipeline[-1]["kind"] == "h_run")
                and (pipeline[-1]["region"] == region_of[sid])
            ):
                pipeline[-1]["operation_ids"].append(sid)
            else:
                pipeline.append(
                    {
                        "kind": "h_run" if kind == "h" else "singleton",
                        "region": region_of[sid],
                        "operation_ids": [sid],
                    }
                )
    h_runs, relaxed_pairs = ([], [])
    for owner, pipeline in blocks.items():
        for previous, following in zip(pipeline, pipeline[1:]):
            for before in previous["operation_ids"]:
                for after in following["operation_ids"]:
                    edge(before, after, "non_H_gadget_boundary")
        for block in pipeline:
            if block["kind"] != "h_run":
                continue
            ids = block["operation_ids"]
            local_relaxed = []
            for j, after in enumerate(ids):
                for before in ids[:j]:
                    if supports[before] & supports[after]:
                        edge(before, after, "same_wire_H_order")
                    elif footprints[before] & footprints[after]:
                        edge(before, after, "overlapping_footprint_H_order")
                    else:
                        local_relaxed.append([before, after])
            if len(ids) > 1:
                h_runs.append(
                    {
                        "source_id": owner,
                        "region": block["region"],
                        "operation_ids": list(ids),
                        "independent_pairs": local_relaxed,
                    }
                )
                relaxed_pairs.extend(local_relaxed)
    last_wire = {}
    for op in operations:
        sid = op["id"]
        for wire in supports[sid]:
            if wire in last_wire:
                edge(last_wire[wire], sid, "data_wire")
            last_wire[wire] = sid
        if op["kind"] == "qand_reset":
            _require(
                op.get("prepared_state") == "|0>",
                "QAND reset does not explicitly prepare zero state",
            )
    last_port_bell, last_resource_bell = ({}, {})
    for group in parsed["port_lifetimes"]:
        acquire, bell, reset = (group["acquire"], group["bell"], group["reset"])
        edge(acquire, bell, "native_injection_body")
        if reset is not None:
            for point in by_id[reset]["ports"]:
                edge(last_port_bell[tuple(point)], reset, "Bell_END_before_port_reset")
            edge(reset, acquire, "explicit_reset_before_port_reuse")
        for point in map(tuple, group["resources"]):
            if point in last_resource_bell:
                edge(last_resource_bell[point], acquire, "Bell_END_before_resource_reuse")
            last_resource_bell[point] = bell
        for point in map(tuple, group["ports"]):
            last_port_bell[point] = bell
    injection_order = parsed["fixed_injection_order"]
    for before, after in zip(injection_order, injection_order[1:]):
        edge(before, after, "canonical_singleton_injection_order")
    for consumer, outcomes in parsed["outcomes_by_consumer"].items():
        for outcome in outcomes:
            edge(
                parsed["producer_by_outcome"][outcome]["operation_id"],
                consumer,
                "classical_outcome",
            )
    ordered = sorted(reasons, key=lambda pair: (positions[pair[0]], positions[pair[1]]))
    result = {
        key: value
        for key, value in parsed.items()
        if key not in {"dependencies", "schema", "scope"}
    }
    result.update(
        schema=SCHEMA,
        dag_policy=DAG_POLICY,
        scope=_SCOPE,
        proof_scope=_SCOPE,
        dependencies=[list(pair) for pair in ordered],
        dependency_reasons=[
            {"before": before, "after": after, "reasons": sorted(reasons[before, after])}
            for before, after in ordered
        ],
        batchable_kinds=sorted(BATCHABLE_KINDS),
        gadget_blocks=blocks,
        h_runs=h_runs,
        relaxed_H_pairs=relaxed_pairs,
    )
    return result


def _h_exchange_proof(operations, constraints, stage_of):
    """Construct adjacent-swap witnesses, independently of scheduler claims.

    Distinct-gadget interleaving retains the previously admitted disjoint-data
    contract.  The only *new* permutation within a gadget is proved below by
    explicit adjacent H exchanges.  Simultaneous H pairs receive separate
    full-support/full-footprint certificates even when no swap is necessary.
    """
    by_id = {op["id"]: op for op in operations}
    positions = {op["id"]: i for i, op in enumerate(operations)}
    supports = {sid: set(wires) for sid, wires in constraints["operation_supports"].items()}
    footprints = {
        sid: set(map(tuple, points)) for sid, points in constraints["operation_footprints"].items()
    }
    swaps, parallel = ([], [])
    projections = defaultdict(list)
    for op in operations:
        if op["kind"] != "port_reset":
            for owner in constraints["source_ids_by_operation"][op["id"]]:
                projections[owner].append(op["id"])

    def independent(left, right):
        _require(
            by_id[left]["kind"] == by_id[right]["kind"] == "h"
            and (not supports[left] & supports[right])
            and (not footprints[left] & footprints[right])
            and (constraints["operation_region"][left] == constraints["operation_region"][right]),
            "H exchange is not on independent complete native footprints",
        )

    for owner, ids in projections.items():
        target = sorted(ids, key=lambda sid: (stage_of[sid], positions[sid]))
        working, local_swaps = (list(ids), [])
        for index, sid in enumerate(target):
            found = working.index(sid, index)
            while found > index:
                left, right = (working[found - 1], working[found])
                independent(left, right)
                local_swaps.append([left, right])
                working[found - 1], working[found] = (right, left)
                found -= 1
        _require(working == target, "canonical-to-actual H exchange witness failed")
        if local_swaps:
            swaps.append(
                {
                    "source_id": owner,
                    "canonical_order": list(ids),
                    "actual_order": target,
                    "adjacent_swaps": local_swaps,
                }
            )
    for run in constraints["h_runs"]:
        ids = run["operation_ids"]
        for j, right in enumerate(ids):
            for left in ids[:j]:
                if stage_of[left] == stage_of[right]:
                    independent(left, right)
                    parallel.append([left, right])
    return {
        "policy": DAG_POLICY,
        "method": "independently replayed complete gadget projections by adjacent disjoint-H exchanges plus complete-footprint batch checks",
        "checked_complete_gadget_projections": len(projections),
        "checked_H_runs": len(constraints["h_runs"]),
        "checked_H_adjacent_swaps": sum((len(item["adjacent_swaps"]) for item in swaps)),
        "checked_same_gadget_parallel_H_pairs": len(parallel),
        "adjacent_H_swap_witness": swaps,
        "parallel_H_witness": parallel,
        "scope": "new within-gadget H-run permutations only; other permutations retain the source-segment disjoint-gadget contract",
    }


def verify_causal_timeline(
    operations: Sequence[Mapping[str, Any]],
    logical_gadgets: Sequence[Mapping[str, Any]],
    stages: Sequence[Mapping[str, Any]],
    *,
    dag_policy: str = DAG_POLICY,
    claimed_constraints: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Independently replay an exact-cover v9 timeline; malformed input fails.

    Full-circuit callers must bind ``dag_policy`` to their distinct circuit and
    schedule schemas.  No caller-supplied DAG, costs, masks or proof is trusted.
    If a constraint certificate is retained, it must exactly match rederivation.
    """
    try:
        constraints = derive_causal_constraints(operations, logical_gadgets, dag_policy=dag_policy)
        if claimed_constraints is not None:
            _require(
                isinstance(claimed_constraints, Mapping)
                and _json(claimed_constraints) == _json(constraints),
                "claimed causal DAG certificate does not exactly replay",
            )
        _require(primitive_checks._sequence(stages), "stages must be a sequence")
        by_id = {op["id"]: op for op in operations}
        for op in operations:
            _require(
                all((type(op[key]) is int and op[key] >= 0 for key in ("start", "end")))
                and op["end"] - op["start"] == op["duration"],
                "operation timestamp/duration mismatch",
            )
        stage_of, cycles, batches = ({}, Counter(), Counter())
        now, previous_region = (0, -1)
        for index, stage in enumerate(stages):
            _require(isinstance(stage, Mapping), "stage must be an object")
            ids = stage["operation_ids"]
            _require(
                primitive_checks._sequence(ids)
                and bool(ids)
                and all((isinstance(sid, str) and sid in by_id for sid in ids)),
                "empty stage or unknown operation ID",
            )
            _require(
                len(set(ids)) == len(ids) and (not set(ids) & stage_of.keys()),
                "duplicate operation in stages",
            )
            selected = [by_id[sid] for sid in ids]
            kinds = {op["kind"] for op in selected}
            _require("h" not in kinds or kinds == {"h"}, "H cannot overlap a non-H body")
            _require(len(kinds) == 1, "stage mixes primitive kinds")
            kind = selected[0]["kind"]
            duration, category = primitive_checks._PRIMITIVES[kind]
            _require(
                duration != 0 or len(ids) == 1,
                "zero-duration explicit pulse must have a singleton stage",
            )
            _require(
                len(ids) == 1 or kind in BATCHABLE_KINDS,
                "injection/reset/barrier macro must remain singleton",
            )
            _require(
                type(stage["start"]) is int
                and type(stage["end"]) is int
                and (stage["start"] == now)
                and (stage["end"] == now + duration),
                "stage overlaps, has a gap or changes primitive duration",
            )
            _require(
                all(
                    (op["start"] == stage["start"] and op["end"] == stage["end"] for op in selected)
                ),
                "operation timestamps differ from exact-cover stage",
            )
            regions = {constraints["operation_region"][sid] for sid in ids}
            _require(len(regions) == 1, "stage crosses a source segment or hard barrier")
            region = next(iter(regions))
            _require(
                region >= previous_region,
                "operation crosses a canonical source segment or hard barrier",
            )
            previous_region = region
            used_wires, used_patches = (set(), set())
            for sid in ids:
                wires = set(constraints["operation_supports"][sid])
                patches = set(map(tuple, constraints["operation_footprints"][sid]))
                _require(not wires & used_wires, "parallel operations share physical data")
                _require(
                    not patches & used_patches,
                    "parallel complete operation footprints are not vertex-disjoint",
                )
                used_wires.update(wires)
                used_patches.update(patches)
                stage_of[sid] = index
            cycles[category] += duration
            batches[kind] += 1
            now = stage["end"]
        _require(
            set(stage_of) == set(by_id),
            "stages omit canonical operations or explicit zero-cycle pulses",
        )
        for before, after in constraints["dependencies"]:
            _require(
                stage_of[before] < stage_of[after]
                and by_id[before]["end"] <= by_id[after]["start"],
                f"canonical causal/data/injection/reset/outcome dependency reversed: {before} -> {after}",
            )
        for sid, outcomes in constraints["outcomes_by_consumer"].items():
            for outcome in outcomes:
                producer = constraints["producer_by_outcome"][outcome]
                producing = by_id[producer["operation_id"]]
                _require(
                    producing["start"] + producer["usable_after_relative_cycle"]
                    <= by_id[sid]["start"],
                    "correction precedes outcome readiness or Bell reserved completion",
                )
        lifetimes = constraints["port_lifetimes"]
        for index, lifetime in enumerate(lifetimes):
            acquire, bell = (by_id[lifetime["acquire"]], by_id[lifetime["bell"]])
            allocated = set(map(tuple, lifetime["ports"] + lifetime["resources"]))
            begin, end = (acquire["start"], bell["end"])
            _require(begin < end and acquire["end"] <= bell["start"], "invalid injection lifetime")
            for other in lifetimes[index + 1 :]:
                other_begin = by_id[other["acquire"]]["start"]
                other_end = by_id[other["bell"]]["end"]
                if max(begin, other_begin) < min(end, other_end):
                    _require(
                        not allocated & set(map(tuple, other["ports"] + other["resources"])),
                        "simultaneous CCZ lifetimes reuse occupied port/resource patches",
                    )
            for op in operations:
                if op["id"] in {lifetime["acquire"], lifetime["bell"]} or op["duration"] == 0:
                    continue
                if max(begin, op["start"]) < min(end, op["end"]):
                    _require(
                        not allocated
                        & set(map(tuple, constraints["operation_footprints"][op["id"]])),
                        "operation touches port/resource patches still occupied by a live CCZ group",
                    )
        proof = _h_exchange_proof(operations, constraints, stage_of)
        _require(now == sum(cycles.values()), "category cost sum differs from latency")
        return {
            "schema": SCHEMA + ":replay",
            "dag_policy": DAG_POLICY,
            "passed": True,
            "errors": [],
            "latency": now,
            "cost_by_category": dict(sorted(cycles.items())),
            "batch_counts": dict(sorted(batches.items())),
            "operation_counts": dict(sorted(Counter((op["kind"] for op in operations)).items())),
            "checked_operations": len(operations),
            "checked_stages": len(stages),
            "checked_dependencies": len(constraints["dependencies"]),
            "checked_outcomes": len(constraints["producer_by_outcome"]),
            "checked_port_lifetimes": len(lifetimes),
            "checked_H_footprints": sum((op["kind"] == "h" for op in operations)),
            "region_counts": dict(
                sorted(Counter((region["kind"] for region in constraints["regions"])).items())
            ),
            "causal_reordering_proof": proof,
            "proof_scope": _SCOPE,
            "scope": _SCOPE,
        }
    except (
        ValueError,
        KeyError,
        TypeError,
        IndexError,
        AttributeError,
        RuntimeError,
        OverflowError,
    ) as error:
        return {
            "schema": SCHEMA + ":replay",
            "dag_policy": DAG_POLICY,
            "passed": False,
            "errors": [str(error)],
            "latency": None,
            "cost_by_category": {},
            "batch_counts": {},
            "proof_scope": _SCOPE,
            "scope": _SCOPE,
        }
