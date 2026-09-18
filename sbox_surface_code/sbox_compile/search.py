"""Full nonmonotone Joint Rewrite search with a monotone verified incumbent."""

from __future__ import annotations

from collections import Counter
import copy
from dataclasses import dataclass
import itertools
import random

from . import linear, nonlinear, windows
from .geometry import _require
from .matrix import invert_perm
from .model import canonical, objective, physical_key
from .replay import replay, transition
from .verify import require_verified

BEAM_WIDTH = 4
UPHILL_SLACK = 32
WINDOW_BLOCKS = 3
MATRIX_CANDIDATES = 3
WINDOW_COMBINATIONS = 2
MAX_CHAIN = 24


@dataclass
class State:
    circuit: dict
    steps: list
    key: bytes


def select_beam(states):
    distinct = {}
    for state in states:
        old = distinct.get(state.key)
        if old is None or (objective(state.circuit), len(state.steps)) < (
            objective(old.circuit),
            len(old.steps),
        ):
            distinct[state.key] = state
    ranked = sorted(distinct.values(), key=lambda s: (objective(s.circuit), len(s.steps), s.key))
    eligible = [
        s
        for s in ranked
        if objective(s.circuit)[0] <= objective(ranked[0].circuit)[0] + UPHILL_SLACK
    ]
    selected, seen = [], set()
    for state in eligible:
        placements = canonical(
            [r["choice"]["output_placement"] for r in state.circuit["linear_candidates"]]
        )
        if placements not in seen:
            seen.add(placements)
            selected.append(state)
        if len(selected) == BEAM_WIDTH:
            return selected
    selected_keys = {s.key for s in selected}
    return (selected + [s for s in eligible if s.key not in selected_keys])[:BEAM_WIDTH]


def window_order(circuit, seed):
    geom = windows.geometry(circuit)
    weights = []
    for record in circuit["linear_candidates"]:
        c = copy.deepcopy(record["candidate"])
        rows = linear.replay_rows(c)
        c["output_permutation"] = list(range(len(rows)))
        weights.append(linear.validate(c, rows, geom)[0])
    rng = random.Random(seed)
    pairs = [(i, i + WINDOW_BLOCKS - 1) for i in range(len(weights) - WINDOW_BLOCKS + 1)]
    ranked = sorted(pairs, key=lambda w: (-sum(weights[w[0] : w[1] + 1]), rng.random()))
    rest = list(pairs)
    rng.shuffle(rest)
    order = []
    for pair in itertools.zip_longest(ranked, rest):
        for window in pair:
            if window is not None and window not in order:
                order.append(window)
    return order


def _mutate(width, rng):
    permutation = list(range(width))
    cycle = rng.sample(range(width), 3 if rng.random() < 0.25 else 2)
    for old, new in zip(cycle, cycle[1:] + cycle[:1]):
        permutation[old] = new
    return permutation


def _combinations(portfolios, tasks, geom):
    costs = {
        id(c): linear.validate(c, t["rows"], geom) for p, t in zip(portfolios, tasks) for c in p
    }

    def key(entries):
        return (
            sum(costs[id(c)][0] for c in entries),
            sum(costs[id(c)][1] for c in entries),
            tuple(canonical(c["layers"]) for c in entries),
        )

    unique = {}
    for entries in itertools.product(*portfolios):
        unique.setdefault(tuple(canonical(c["layers"]) for c in entries), entries)
    ranked = sorted(unique.values(), key=key)
    selected = ranked[: WINDOW_COMBINATIONS - 1]
    if len(ranked) > len(selected):
        best = ranked[0]
        selected.append(
            min(
                ranked[len(selected) :],
                key=lambda es: (
                    -sum(c["layers"] != b["layers"] for c, b in zip(es, best)),
                    key(es),
                ),
            )
        )
    return selected


def optimize(parent, *, seed=0, iterations=32, budget=100, progress=None):
    """Search uses only its parent, seed, and work budget; no result/recipe I/O.

    A budget unit allocates 32 routed rewrite words, 2000 independent fresh
    proposals, and two regional joint proposals. Scheduling uses up to 25,000
    states per region. Fixed work counts replace machine-dependent CPU cutoffs.
    """
    _require(type(seed) is int and 0 <= seed < 2**63, "invalid search seed")
    _require(type(iterations) is int and 1 <= iterations <= 128, "iterations must be in [1,128]")
    _require(type(budget) is int and 1 <= budget <= 1000, "budget must be in [1,1000]")
    require_verified(parent)
    parent_snapshot = canonical(parent)
    rng = random.Random(seed)
    incumbent = State(copy.deepcopy(parent), [], physical_key(parent))
    beam, seen = [incumbent], {incumbent.key: 0}
    counters, trace = Counter(), []
    order = window_order(parent, seed)
    nodes = min(25000, 1000 * budget)
    for iteration in range(iterations):
        source_state = beam[0] if iteration % 3 == 0 else rng.choice(beam)
        if len(source_state.steps) >= MAX_CHAIN:
            continue
        source = source_state.circuit
        start, stop = order[(iteration + seed % len(order)) % len(order)]
        n, geom = source["logical_width"], windows.geometry(source)
        identity = list(range(n))
        offers = {}
        lane = (
            "identity"
            if iteration % 4 == 1
            else "free_output" if iteration % 2 == 0 else "independent_cycles"
        )
        if progress:
            progress(
                {
                    "iteration": iteration,
                    "phase": "matrix_search",
                    "window": [start, stop],
                    "lane": lane,
                    "best": objective(incumbent.circuit),
                }
            )
        boundaries = [identity[:]]
        for offset in range(WINDOW_BLOCKS - 1):
            if lane == "free_output":
                t = windows.task(source, start + offset, boundaries[-1], identity)
                proposals = linear.rewrite(
                    t["rows"],
                    geom,
                    seed=rng.randrange(2**31),
                    work=32 * budget,
                    warm=t["reference"],
                    free=True,
                )
                varied = [p for p in proposals if p["output_permutation"] != identity]
                pool = varied if varied and iteration % 4 == 0 else proposals
                if pool:
                    offer = pool[(iteration // 2) % len(pool)]
                    offers[offset] = offer
                    boundaries.append(invert_perm(offer["output_permutation"]))
                else:
                    boundaries.append(_mutate(n, rng))
                counters["free_output_calls"] += 1
            else:
                boundaries.append(identity[:] if lane == "identity" else _mutate(n, rng))
        boundaries.append(identity[:])
        tasks = [
            windows.task(source, start + offset, boundaries[offset], boundaries[offset + 1])
            for offset in range(WINDOW_BLOCKS)
        ]
        portfolios = []
        for offset, t in enumerate(tasks):
            unchanged = t["pin"] == t["pout"] == identity
            candidates = []
            if unchanged:
                retained = copy.deepcopy(source["linear_candidates"][t["index"]]["candidate"])
                retained["output_permutation"] = identity[:]
                candidates.append(retained)
            if offset in offers:
                candidates.append(windows.bind(offers[offset], t, geom))
            candidates.extend(
                linear.rewrite(
                    t["rows"],
                    geom,
                    seed=rng.randrange(2**31),
                    work=32 * budget,
                    warm=t["reference"],
                )
            )
            independent = []
            if offset == iteration % WINDOW_BLOCKS or not candidates:
                independent = linear.fresh(
                    t["rows"],
                    geom,
                    seed=rng.randrange(2**31),
                    work=2000 * budget,
                    depth_goal=max(
                        1,
                        source["linear_candidates"][t["index"]]["candidate"]["stats"][
                            "surface_depth"
                        ]
                        // 2,
                    ),
                )
                candidates.extend(independent)
                counters["independent_fresh_calls"] += 1
                counters["independent_fresh_candidates"] += len(independent)
            portfolios.append(
                linear.portfolio(candidates, t["rows"], geom, fresh_candidates=independent)
            )
        if not all(portfolios):
            counters["incomplete_covers"] += 1
            continue
        costs = []
        for entries in _combinations(portfolios, tasks, geom):
            try:
                window = windows.materialize(
                    source, start, boundaries, entries, schedule_nodes=nodes
                )
            except ValueError as error:
                counters["failed_window_reconstructions"] += 1
                if progress:
                    progress(
                        {"iteration": iteration, "phase": "window_rejected", "error": str(error)}
                    )
                continue
            if progress:
                progress(
                    {
                        "iteration": iteration,
                        "phase": "nonlinear_search",
                        "window_cost": objective(window),
                    }
                )
            names = nonlinear.region_names(window, start, stop, downstream=iteration % 4 == 0)
            variants, report = nonlinear.optimize(
                window,
                names,
                seed=rng.randrange(2**31),
                work=2 * budget,
                schedule_nodes=nodes,
                comparisons=min(24, 2 * budget),
            )
            counters.update(report)
            candidates = [window] + variants
            children = []
            for circuit in candidates:
                key = physical_key(circuit)
                counters["complete_comparisons"] += 1
                costs.append(objective(circuit))
                if key in seen and seen[key] <= len(source_state.steps) + 1:
                    continue
                step = {
                    "start": start,
                    "stop": stop,
                    "boundaries": boundaries,
                    "window": transition(source, window),
                    "nonlinear": None if circuit == window else transition(window, circuit),
                }
                child = State(circuit, source_state.steps + [step], key)
                seen[key] = len(child.steps)
                children.append(child)
                if objective(circuit) < objective(incumbent.circuit):
                    incumbent = child
                    counters["incumbent_improvements"] += 1
            beam = select_beam([*beam, incumbent, *children])
            child_keys = {s.key for s in children}
            counters["uphill_admissions"] += sum(
                s.key in child_keys and objective(s.circuit)[0] > objective(source)[0] for s in beam
            )
        record = {
            "iteration": iteration,
            "window": [start, stop],
            "lane": lane,
            "best": objective(incumbent.circuit),
            "comparisons": costs,
            "beam": [objective(s.circuit) for s in beam],
        }
        trace.append(record)
        if progress:
            progress(record)
    _require(canonical(parent) == parent_snapshot, "search mutated its parent")
    recipe = {
        "case": parent["case"],
        "seed": seed,
        "iterations": iterations,
        "budget": budget,
        "steps": incumbent.steps,
        "expected_cost": incumbent.circuit["cost"],
    }
    rebuilt = replay(parent, recipe)
    _require(rebuilt == incumbent.circuit, "independent final replay changed the circuit")
    return {
        "circuit": rebuilt,
        "recipe": recipe,
        "report": {
            "status": (
                "verified_search_incumbent" if incumbent.steps else "verified_incumbent_fallback"
            ),
            "seed": seed,
            "iterations": iterations,
            "budget": budget,
            "counters": dict(counters),
            "trace": trace,
            "parent": objective(parent),
            "result": objective(rebuilt),
            "termination": "finite_work_budget",
            "optimality": "unknown",
        },
    }
