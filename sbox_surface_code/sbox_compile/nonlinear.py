"""Joint port-cover, correction-direction, routing, and H-run portfolios."""

from __future__ import annotations

import copy
import itertools
import random

from . import geometry as front, native
from .compiler import compile_circuit
from .injection import _ripup_group, _route_group
from .lower import _dagger_operations, _group_operations, correction_variants, routes_from_artifact
from .model import canonical, objective, physical_key
from .schedule import solve_local
from .windows import geometry


def region_names(circuit, start, stop, *, downstream=False):
    linear_index = -1
    names, tail = [], []
    for segment in circuit["logical_schedule"]:
        if segment["kind"] == "linear":
            linear_index += 1
        elif segment["kind"] == "ccz":
            if start <= linear_index < stop:
                names.append(segment["segment"])
            elif linear_index >= stop:
                tail.append(segment["segment"])
    return names + tail[: int(downstream)]


def _context(circuit, name):
    routes = routes_from_artifact(circuit)
    used = set()
    for segment in circuit["logical_schedule"]:
        if segment["segment"] == name:
            break
        for group in routes.get(segment["segment"], []):
            used.update(tuple(p) for item in group for p in item["port_coords"])
    records = [r for r in circuit["logical_gadgets"] if r["source_id"].rsplit(":", 1)[0] == name]
    proofs = {
        p["source_id"]: p
        for p in circuit["search"]["correction_windows"]
        if p["source_id"].rsplit(":", 1)[0] == name
    }
    return {
        "name": name,
        "records": records,
        "seed_groups": routes.get(name, []),
        "used_ports_before": [list(p) for p in sorted(used)],
        "proofs": proofs,
    }


def _score(context, groups, proofs, geom, *, nodes, seed):
    records = {r["source_id"]: r for r in context["records"]}
    injected = {sid for sid, r in records.items() if r["role"] != "qand_dagger"}
    actual = [g["source_id"] for group in groups for g in group]
    front._require(
        set(actual) == injected and len(actual) == len(injected), "invalid nonlinear source cover"
    )
    operations, used = [], set()
    incoming = set(map(tuple, context["used_ports_before"]))
    external_resets = 0
    for i, group in enumerate(groups):
        ports = {tuple(p) for item in group for p in item["port_coords"]}
        external_resets += int(bool(ports & (used | incoming))) - int(bool(ports & used))
        operations.extend(
            _group_operations(
                geom, records, group, proofs, used, f"{context['name']}:native_round:{i}"
            )
        )
    for record in context["records"]:
        if record["role"] == "qand_dagger":
            operations.extend(_dagger_operations(record, proofs[record["source_id"]], geom))
    operations, _ = front.cancel_disjoint_h_pairs(operations)
    native._decorate_conditional_cnots(operations, geom)
    score = solve_local(
        operations,
        context["records"],
        node_budget=nodes,
        frontier_limit=1500,
        greedy_starts=1,
        cpu_seconds=None,
        seed=seed,
    )
    return score["latency"] + external_resets, sum(op["kind"] == "h" for op in operations)


def _correction_route(proof, geom, rng):
    result = copy.deepcopy(proof)
    edge = rng.choice(result["chosen_czs"])
    target = edge["target"]
    control = next(q for q in edge["pair"] if q != target)
    old = edge.get("path", geom.cnot(control, target))
    interior = [
        tuple(p) for p in old[1:-1] if tuple(p) not in geom.occupied | geom.ports | geom.resources
    ]
    blocked = [rng.choice(interior)] if interior and rng.randrange(2) else []
    path = geom.bfs(
        geom.coord(control), geom.coord(target), used=blocked, neighbor_order=rng.randrange(4)
    )
    if path is not None:
        edge["path"] = path
    return result


def regional_portfolio(circuit, name, *, seed, work, schedule_nodes):
    """Complete regional candidates, including non-improving structural changes."""
    rng, geom = random.Random(seed), geometry(circuit)
    context = _context(circuit, name)
    variants = {r["source_id"]: correction_variants(r, geom) for r in context["records"]}
    baseline = {"groups": context["seed_groups"], "proofs": context["proofs"]}
    baseline["key"] = canonical([baseline["groups"], baseline["proofs"]])
    baseline["score"] = _score(
        context,
        baseline["groups"],
        baseline["proofs"],
        geom,
        nodes=min(1200, schedule_nodes),
        seed=seed,
    )
    retained = {baseline["key"]: baseline}
    originals = [g for group in context["seed_groups"] for g in group]
    count = len(originals)
    masks = (
        ([(1 << count) - 1] + sorted(range(1, (1 << count) - 1), key=lambda m: (m.bit_count(), m)))
        if count
        else []
    )
    failures = 0
    for iteration in range(work):
        pool = sorted(retained.values(), key=lambda e: (e["score"], e["key"]))
        base = pool[0] if iteration % 3 == 0 else rng.choice(pool)
        groups, proofs = copy.deepcopy(base["groups"]), copy.deepcopy(base["proofs"])
        try:
            if groups and iteration % 3 != 1:
                if iteration % 4 == 0:
                    # Route every bounded subset over time, then complete its cover.
                    mask = masks[(iteration // 4) % len(masks)]
                    selected = [g for i, g in enumerate(originals) if mask & (1 << i)]
                    routed = _route_group(
                        context, selected, geom, rng, starts=4, fresh_only=iteration % 8 == 0
                    )
                    if routed is not None:
                        ids = {g["source_id"] for g in routed}
                        groups = [routed] + [
                            [g for g in group if g["source_id"] not in ids]
                            for group in context["seed_groups"]
                        ]
                        groups = [g for g in groups if g]
                else:
                    i = rng.randrange(len(groups))
                    routed = (
                        _ripup_group(context, groups[i], geom, rng, starts=6)
                        if iteration % 2
                        else _route_group(context, groups[i], geom, rng, starts=4)
                    )
                    if routed is not None:
                        groups[i] = routed
                    elif len(groups[i]) > 1:
                        group = groups.pop(i)
                        groups[i:i] = [[g] for g in group]
                    if len(groups) > 1 and iteration % 5 == 0:
                        rng.shuffle(groups)
            # Direction/order enumeration and detours are both paid explicitly.
            records = context["records"] if iteration % 4 == 0 else [rng.choice(context["records"])]
            for record in records:
                sid = record["source_id"]
                choices = variants[sid]
                proof = copy.deepcopy(
                    choices[(iteration + rng.randrange(len(choices))) % len(choices)]["proof"]
                )
                if iteration % 2:
                    proof = _correction_route(proof, geom, rng)
                proofs[sid] = proof
            key = canonical([groups, proofs])
            if key in retained:
                continue
            score = _score(
                context,
                groups,
                proofs,
                geom,
                nodes=min(1200, schedule_nodes),
                seed=seed + iteration,
            )
            entry = {"groups": groups, "proofs": proofs, "key": key, "score": score}
            retained[key] = entry
            ranked = sorted(retained.values(), key=lambda e: (e["score"], e["key"]))
            chosen = ranked[:3]
            varied = next((e for e in ranked[3:] if e["groups"] != ranked[0]["groups"]), None)
            if varied is not None:
                chosen.append(varied)
            retained = {e["key"]: e for e in chosen}
        except ValueError:
            failures += 1
    entries = [
        e
        for e in sorted(retained.values(), key=lambda e: (e["score"], e["key"]))
        if e["key"] != baseline["key"]
    ]
    return {"name": name, "entries": entries, "failed_proposals": failures}


def combinations(regions, limit=24):
    """All-region waves, mixed waves, singles, prefixes, and pairs."""
    lanes = [
        [(region["name"], entry) for entry in region["entries"]]
        for region in regions
        if region["entries"]
    ]
    selected, seen = [], set()

    def admit(entries):
        key = tuple((name, e["key"]) for name, e in entries)
        if entries and key not in seen and len(selected) < limit:
            seen.add(key)
            selected.append(entries)

    for wave in range(min(3, max(map(len, lanes), default=0))):
        admit([lane[min(wave, len(lane) - 1)] for lane in lanes])
    if len(lanes) > 2:
        admit([lane[i % len(lane)] for i, lane in enumerate(lanes)])
    for depth in range(4):
        for lane in lanes:
            if depth < len(lane):
                admit([lane[depth]])
        if depth == 0:
            for length in range(2, len(lanes)):
                admit([lane[0] for lane in lanes[:length]])
            for left, right in itertools.combinations(lanes, 2):
                admit([left[0], right[0]])
    return selected


def optimize(circuit, names, *, seed, work, schedule_nodes=25000, comparisons=24):
    """Local scores order proposals; only whole-circuit replays admit candidates."""
    regions = [
        regional_portfolio(
            circuit, name, seed=seed + 1009 * i, work=work, schedule_nodes=schedule_nodes
        )
        for i, name in enumerate(names)
    ]
    retained = {physical_key(circuit): circuit}
    full_comparisons = 0
    for selected in combinations(regions, limit=comparisons):
        routes = routes_from_artifact(circuit)
        proofs = {}
        for name, entry in selected:
            routes[name] = entry["groups"]
            proofs.update(entry["proofs"])
        result = compile_circuit(
            circuit, routes=routes, corrections=proofs, schedule_nodes=schedule_nodes
        )
        retained[physical_key(result)] = result
        full_comparisons += 1
        ranked = sorted(retained.items(), key=lambda item: (objective(item[1]), item[0]))
        retained = dict(ranked[:3])
    return list(retained.values()), {
        "regions": len(regions),
        "full_comparisons": full_comparisons,
        "failed_route_proposals": sum(r["failed_proposals"] for r in regions),
    }
