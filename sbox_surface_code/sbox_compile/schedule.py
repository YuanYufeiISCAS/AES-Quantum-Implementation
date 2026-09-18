"""schedule."""

# Adapted from the authors' strict_sbox/causal_scheduler.py; campaign infrastructure removed.
from __future__ import annotations
from collections import defaultdict
import time
from .scheduler import RegionProblem, solve_region
from .timeline import DAG_POLICY, derive_causal_constraints

SCHEMA = "strict-causal-h-run-schedule-v1"


def solve_local(
    operations,
    records,
    *,
    node_budget=1200,
    frontier_limit=1500,
    greedy_starts=1,
    cpu_seconds=0.2,
    seed=0,
    fallback=None,
):
    """Score one complete nonlinear region; explicit singleton fallback exists."""
    constraints = derive_causal_constraints(operations, records)
    if len(constraints["regions"]) != 1 or not constraints["regions"][0].get("source_ids"):
        raise ValueError("local causal search requires exactly one nonlinear region")
    problem = RegionProblem(operations, constraints, constraints["regions"][0])
    if fallback is None:
        fallback = [[op["id"]] for op in operations]
    result = solve_region(
        problem,
        fallback,
        node_budget=node_budget,
        frontier_limit=frontier_limit,
        greedy_starts=greedy_starts,
        cpu_seconds=cpu_seconds,
        seed=seed,
    )
    return dict(
        result,
        dag_policy=DAG_POLICY,
        proof_scope="fixed explicit operations/routes and canonical H-run DAG only",
    )


def schedule_causal(
    operations,
    logical_gadgets,
    fallback_stages,
    *,
    node_budget=25000,
    frontier_limit=16000,
    greedy_starts=12,
    seed=0,
    cpu_seconds=None,
):
    """Preserve every source/linear boundary and optimize its actual H-run DAG."""
    constraints = derive_causal_constraints(operations, logical_gadgets)
    by_id = {op["id"]: op for op in operations}
    fallback = defaultdict(list)
    for stage in fallback_stages:
        region_ids = {constraints["operation_region"][sid] for sid in stage["operation_ids"]}
        if len(region_ids) != 1:
            raise ValueError("causal fallback crosses a source-region boundary")
        fallback[next(iter(region_ids))].append(stage["operation_ids"])
    stages, reports, now = ([], [], 0)
    deadline = None if cpu_seconds is None else time.process_time() + max(0.0, cpu_seconds)
    for index, region in enumerate(constraints["regions"]):
        remaining = None if deadline is None else max(0.0, deadline - time.process_time())
        problem = RegionProblem(operations, constraints, region)
        result = solve_region(
            problem,
            fallback[index],
            node_budget=node_budget if region.get("source_ids") else 0,
            frontier_limit=frontier_limit,
            greedy_starts=greedy_starts,
            seed=seed + 1009 * index,
            cpu_seconds=remaining,
        )
        reports.append(
            {
                "region": index,
                "segment": region["segment"],
                "kind": region["kind"],
                "operations": len(problem.ids),
                **{key: value for key, value in result.items() if key != "batches"},
            }
        )
        for batch in result["batches"]:
            duration = by_id[batch[0]]["duration"]
            stages.append({"operation_ids": batch, "start": now, "end": now + duration})
            now += duration
    return {
        "schema": SCHEMA,
        "dag_policy": DAG_POLICY,
        "stages": stages,
        "latency": now,
        "regions": reports,
        "scope": "fixed native interface; independently justified disjoint H-run concurrency; unchanged non-H/source barriers and homogeneous primitives",
        "stats": {
            "regions": len(reports),
            "improved_regions": sum(
                (item["latency"] < item["initial_latency"] for item in reports)
            ),
            "expanded": sum((item["stats"]["expanded"] for item in reports)),
            "exact_nonlinear_regions": sum(
                (
                    item["status"] == "EXACT_FIXED_SEGMENT"
                    and bool(constraints["regions"][item["region"]].get("source_ids"))
                    for item in reports
                )
            ),
        },
    }
