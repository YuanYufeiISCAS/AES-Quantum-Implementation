"""scheduler."""

# Adapted from the authors' strict_sbox/segment_scheduler.py; campaign infrastructure removed.
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
import heapq
import random
import time

SCHEMA = "strict-whole-nonlinear-segment-schedule-v1"


def _bits(mask):
    while mask:
        bit = mask & -mask
        yield (bit.bit_length() - 1)
        mask ^= bit


@dataclass(slots=True)
class _Node:
    state: int
    cost: int
    parent: "_Node | None"
    batches: tuple[int, ...]


class RegionProblem:
    """Compact bitset instance; inputs are a separately derived contract."""

    def __init__(self, operations, constraints, region):
        by_id = {op["id"]: op for op in operations}
        self.ids = list(region["operation_ids"])
        self.ops = [by_id[sid] for sid in self.ids]
        self.index = {sid: i for i, sid in enumerate(self.ids)}
        self.n = len(self.ids)
        self.full = (1 << self.n) - 1
        self.pred = [0] * self.n
        self.succ = [0] * self.n
        for before, after in constraints["dependencies"]:
            if before in self.index and after in self.index:
                i, j = (self.index[before], self.index[after])
                self.pred[j] |= 1 << i
                self.succ[i] |= 1 << j
        self.durations = [op["duration"] for op in self.ops]
        self.kinds = [op["kind"] for op in self.ops]
        self.batchable = set(constraints["batchable_kinds"])
        self.footprints = [
            set(map(tuple, constraints["operation_footprints"][sid])) for sid in self.ids
        ]
        self.supports = [set(constraints["operation_supports"][sid]) for sid in self.ids]
        self.conflict = [0] * self.n
        for i in range(self.n):
            for j in range(i):
                if self.footprints[i] & self.footprints[j] or self.supports[i] & self.supports[j]:
                    self.conflict[i] |= 1 << j
                    self.conflict[j] |= 1 << i
        self.by_kind = defaultdict(int)
        for i, kind in enumerate(self.kinds):
            self.by_kind[kind, self.durations[i], self.ops[i]["category"]] |= 1 << i
        self.congestion = []
        self.singleton = []
        for (kind, duration, _), mask in self.by_kind.items():
            if not duration:
                continue
            if kind not in self.batchable:
                self.singleton.append((duration, mask))
                continue
            groups = defaultdict(int)
            for i in _bits(mask):
                for point in self.footprints[i]:
                    groups["patch", point] |= 1 << i
                for wire in self.supports[i]:
                    groups["wire", wire] |= 1 << i
            masks = set(groups.values()) | {1 << i for i in _bits(mask)}
            ordered = sorted(masks, key=int.bit_count, reverse=True)
            maximal = []
            for item in ordered:
                if not any((item & other == item for other in maximal)):
                    maximal.append(item)
            self.congestion.append((duration, tuple(maximal)))
        if any((self.pred[i] >> i for i in range(self.n))):
            raise ValueError("region dependency contradicts canonical order")
        self.tail = [0] * self.n
        for i in reversed(range(self.n)):
            self.tail[i] = self.durations[i] + max(
                (self.tail[j] for j in _bits(self.succ[i])), default=0
            )

    def ready(self, state):
        return sum(
            (1 << i for i in _bits(self.full ^ state) if self.pred[i] & state == self.pred[i])
        )

    def close_zero(self, state):
        """Execute ready zero-time operations as explicit singleton stages."""
        batches = []
        while True:
            i = next((i for i in _bits(self.ready(state)) if self.durations[i] == 0), None)
            if i is None:
                return (state, tuple(batches))
            state |= 1 << i
            batches.append(1 << i)

    def _maximal_sets(self, ready):
        """All maximal compatible subsets of a same-kind ready set.

        Adding another already-ready compatible operation cannot hurt future
        feasibility in the homogeneous, non-overlapping-stage model.  Thus
        non-maximal subsets are dominated; no mixed-kind overlap is assumed.
        """
        if not any((self.conflict[i] & ready for i in _bits(ready))):
            yield ready
            return

        def visit(chosen, possible, excluded):
            if not possible and (not excluded):
                yield chosen
                return
            union = possible | excluded
            pivot = max(
                _bits(union), key=lambda i: (possible & ~self.conflict[i] & ~(1 << i)).bit_count()
            )
            candidates = possible & (self.conflict[pivot] | 1 << pivot)
            for i in tuple(_bits(candidates)):
                bit = 1 << i
                neighbors = ready & ~self.conflict[i] & ~bit
                yield from visit(chosen | bit, possible & neighbors, excluded & neighbors)
                possible &= ~bit
                excluded |= bit

        yield from visit(0, ready, 0)

    def moves(self, state):
        ready = self.ready(state)
        for (kind, duration, _), mask in self.by_kind.items():
            available = mask & ready
            if not available:
                continue
            if kind not in self.batchable:
                for i in _bits(available):
                    yield (duration, 1 << i)
            else:
                for batch in self._maximal_sets(available):
                    yield (duration, batch)

    def lower_bound(self, state):
        remaining = self.full ^ state
        if not remaining:
            return 0
        congestion = sum(
            (duration * (mask & remaining).bit_count() for duration, mask in self.singleton)
        )
        congestion += sum(
            (
                duration * max(((mask & remaining).bit_count() for mask in masks))
                for duration, masks in self.congestion
            )
        )
        critical = max((self.tail[i] for i in _bits(remaining)))
        return max(congestion, critical)

    def batches_to_ids(self, batches):
        return [[self.ids[i] for i in _bits(batch)] for batch in batches]

    def ids_to_batches(self, batches):
        return [sum((1 << self.index[sid] for sid in batch)) for batch in batches]

    def check_batches(self, batches):
        state, cost = (0, 0)
        for batch in batches:
            if not batch or batch & state or batch & ~self.full:
                raise ValueError("duplicate, empty or foreign fallback batch")
            indices = list(_bits(batch))
            if (
                len({(self.kinds[i], self.durations[i], self.ops[i]["category"]) for i in indices})
                != 1
            ):
                raise ValueError("inhomogeneous fallback batch")
            if len(indices) > 1 and self.kinds[indices[0]] not in self.batchable:
                raise ValueError("nonbatchable fallback stage")
            for i in indices:
                if self.pred[i] & state != self.pred[i] or self.conflict[i] & batch:
                    raise ValueError("fallback violates dependencies/resources")
            cost += self.durations[indices[0]]
            state |= batch
        if state != self.full:
            raise ValueError("incomplete fallback region")
        return cost


def _path(node):
    parts = []
    while node is not None:
        parts.append(node.batches)
        node = node.parent
    return [batch for part in reversed(parts) for batch in part]


def solve_region(
    problem,
    fallback,
    *,
    node_budget=20000,
    frontier_limit=12000,
    greedy_starts=8,
    seed=0,
    cpu_seconds=None,
    checkpoint=None,
):
    """Anytime feasible incumbent plus a honestly bounded A* search."""
    if (
        type(node_budget) is not int
        or node_budget < 0
        or type(frontier_limit) is not int
        or (frontier_limit < 1)
        or (type(greedy_starts) is not int)
        or (greedy_starts < 0)
    ):
        raise ValueError("invalid region search bounds")
    started = time.process_time()
    deadline = None if cpu_seconds is None else started + max(0.0, cpu_seconds)
    best = problem.ids_to_batches(fallback)
    upper = problem.check_batches(best)
    old_upper = upper
    root_state, zero = problem.close_zero(0)
    root = _Node(root_state, 0, None, zero)
    lower = problem.lower_bound(root_state)
    rng = random.Random(seed)
    updates = []

    def admit(batches, cost, method):
        nonlocal best, upper
        if cost < upper:
            if problem.check_batches(batches) != cost:
                raise ValueError("optimizer's reconstructed cost mismatch")
            best, upper = (batches, cost)
            updates.append(
                {"latency": cost, "method": method, "cpu_seconds": time.process_time() - started}
            )
            if checkpoint:
                checkpoint(problem.batches_to_ids(best), upper)

    for start in range(greedy_starts if node_budget else 0):
        if deadline is not None and time.process_time() >= deadline:
            break
        state, batches, cost = (root_state, list(zero), 0)
        while state != problem.full:
            options = []
            for duration, batch in problem.moves(state):
                next_state, closure = problem.close_zero(state | batch)
                heuristic = problem.lower_bound(next_state)
                jitter = rng.uniform(-4.0, 4.0) * min(start, 3) if start else 0.0
                options.append(
                    (
                        (
                            duration + heuristic + jitter,
                            -batch.bit_count(),
                            -sum((problem.tail[i] for i in _bits(batch))),
                            batch,
                        ),
                        duration,
                        batch,
                        next_state,
                        closure,
                    )
                )
            if not options:
                raise ValueError("canonical region contains a deadlock")
            _, duration, batch, state, closure = min(options)
            cost += duration
            batches.extend((batch, *closure))
        admit(batches, cost, "greedy" if not start else "randomized_dispatch")
    frontier = [(lower, 0, -root_state.bit_count(), 0, root)]
    seen = {root_state: 0}
    counter, expanded, discarded = (1, 0, 0)
    timed_out, bounded = (False, False)
    while frontier and expanded < node_budget:
        if expanded % 64 == 0 and deadline is not None and (time.process_time() >= deadline):
            timed_out = True
            break
        f, _, _, _, node = heapq.heappop(frontier)
        if f >= upper or seen.get(node.state) != node.cost:
            continue
        expanded += 1
        for duration, batch in problem.moves(node.state):
            cost = node.cost + duration
            if cost >= upper:
                continue
            state, closure = problem.close_zero(node.state | batch)
            if cost >= seen.get(state, 10**20):
                continue
            estimate = cost + problem.lower_bound(state)
            if estimate >= upper:
                continue
            child = _Node(state, cost, node, (batch, *closure))
            seen[state] = cost
            if state == problem.full:
                admit(_path(child), cost, "bounded_A_star")
                continue
            heapq.heappush(frontier, (estimate, -cost, -state.bit_count(), counter, child))
            counter += 1
        if len(frontier) > 2 * frontier_limit:
            retained = heapq.nsmallest(frontier_limit, frontier)
            discarded += len(frontier) - len(retained)
            bounded = True
            frontier = retained
            heapq.heapify(frontier)
    unfinished = any(
        (item[0] < upper and seen.get(item[4].state) == item[4].cost for item in frontier)
    )
    exact = upper == lower or (
        not bounded
        and (not unfinished)
        and (not timed_out)
        and (expanded < node_budget or not frontier)
    )
    return {
        "batches": problem.batches_to_ids(best),
        "latency": upper,
        "initial_latency": old_upper,
        "lower_bound": upper if exact else lower,
        "status": "EXACT_FIXED_SEGMENT" if exact else "FEASIBLE_SEARCH_UNKNOWN",
        "stats": {
            "expanded": expanded,
            "seen": len(seen),
            "frontier_discarded": discarded,
            "node_budget": node_budget,
            "frontier_limit": frontier_limit,
            "greedy_starts": greedy_starts,
            "timed_out": timed_out,
            "cpu_seconds": time.process_time() - started,
        },
        "improvements": updates,
    }
