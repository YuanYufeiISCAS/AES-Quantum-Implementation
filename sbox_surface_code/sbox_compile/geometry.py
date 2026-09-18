"""geometry."""

# Adapted from the authors' strict_sbox/compiler.py; campaign infrastructure removed.
from __future__ import annotations
import copy
import itertools
from collections import deque
from functools import lru_cache
from typing import Any, Iterable, Mapping, Sequence
from .layout import FoundryLayout, resource_coord_for_port
from .matrix import apply_output_permutation, linear_target_rows_from_span, replay_rows

Coord = tuple[int, int]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _coord(value: Sequence[int]) -> Coord:
    _require(len(value) == 2, "coordinate needs two entries")
    return (int(value[0]), int(value[1]))


def _edge(a: Coord, b: Coord) -> tuple[Coord, Coord]:
    return (min(a, b), max(a, b))


class Geometry:
    """Full occupied rectangle, directed terminals, and explicit scratch."""

    def __init__(self, rows: int, cols: int, width: int):
        self.layout = FoundryLayout(rows, cols)
        self.width = width
        _require(0 < width <= rows * cols, "width exceeds the data rectangle")
        self.occupied = {self.layout.data_coord(q) for q in range(width)}
        self.ports = set(self.layout.port_coords())
        self.resources = {resource_coord_for_port(self.layout, p) for p in self.ports}
        self.vertices = self.layout.all_vertices()
        self.route_cache: dict[tuple[int, int], list[list[int]]] = {}

    def coord(self, wire: int) -> Coord:
        _require(0 <= wire < self.width, f"wire {wire} out of range")
        return self.layout.data_coord(wire)

    def port_axis(self, port: Coord) -> str:
        _require(port in self.ports, f"not a port: {port}")
        return "vertical" if port[0] in {self.layout.min_row, self.layout.max_row} else "horizontal"

    @staticmethod
    def _axis_ok(a: Coord, b: Coord, axis: str) -> bool:
        return a[1] == b[1] if axis == "vertical" else a[0] == b[0]

    def validate_path(
        self,
        raw: Sequence[Sequence[int]],
        start: Coord,
        goal: Coord,
        *,
        start_axis: str = "vertical",
        goal_axis: str = "horizontal",
        extra_banned: Iterable[Coord] = (),
    ) -> list[Coord]:
        path = [_coord(p) for p in raw]
        _require(
            len(path) >= 2 and path[0] == start and (path[-1] == goal),
            "route endpoints do not match",
        )
        _require(len(set(path)) == len(path), "route is not simple")
        banned = self.occupied | self.ports | self.resources | set(extra_banned)
        banned -= {start, goal}
        _require(
            not set(path) & banned, "route crosses reserved data, port, resource or another route"
        )
        _require(all((p in self.vertices for p in path)), "route leaves the foundry")
        _require(
            all((abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1 for a, b in zip(path, path[1:]))),
            "non-neighbor route step",
        )
        _require(self._axis_ok(path[0], path[1], start_axis), "route leaves wrong control boundary")
        _require(self._axis_ok(path[-2], path[-1], goal_axis), "route enters wrong target boundary")
        return path

    def bfs(
        self,
        start: Coord,
        goal: Coord,
        *,
        start_axis: str = "vertical",
        goal_axis: str = "horizontal",
        used: Iterable[Coord] = (),
        neighbor_order: int = 0,
    ) -> list[list[int]] | None:
        banned = self.occupied | self.ports | self.resources | set(used)
        banned -= {start, goal}
        if start in set(used) or goal in set(used):
            return None
        queue: deque[Coord] = deque([start])
        parent: dict[Coord, Coord | None] = {start: None}
        offsets = ((-1, 0), (0, -1), (0, 1), (1, 0))
        offsets = offsets[neighbor_order % 4 :] + offsets[: neighbor_order % 4]
        while queue:
            current = queue.popleft()
            if current == goal:
                path: list[list[int]] = []
                while current is not None:
                    path.append(list(current))
                    current = parent[current]
                path.reverse()
                self.validate_path(
                    path, start, goal, start_axis=start_axis, goal_axis=goal_axis, extra_banned=used
                )
                return path
            neighbors = [(current[0] + dr, current[1] + dc) for dr, dc in offsets]
            neighbors.sort(key=lambda p: abs(p[0] - goal[0]) + abs(p[1] - goal[1]))
            for nxt in neighbors:
                if nxt not in self.vertices or nxt in banned or nxt in parent:
                    continue
                if current == start and (not self._axis_ok(current, nxt, start_axis)):
                    continue
                if nxt == goal and (not self._axis_ok(current, nxt, goal_axis)):
                    continue
                parent[nxt] = current
                queue.append(nxt)
        return None

    def cnot(self, control: int, target: int) -> list[list[int]]:
        key = (control, target)
        if key not in self.route_cache:
            path = self.bfs(self.coord(control), self.coord(target))
            _require(path is not None, f"no directed serial CNOT route for {control}->{target}")
            self.route_cache[key] = path
        return copy.deepcopy(self.route_cache[key])

    def h_footprint(self, wire: int) -> dict[str, Any]:
        r, c = self.coord(wire)
        for dr, dc in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
            auxiliary = [(r + dr, c), (r, c + dc), (r + dr, c + dc)]
            if all(
                (
                    p in self.vertices and p not in self.occupied | self.ports | self.resources
                    for p in auxiliary
                )
            ):
                footprint = [list((r, c)), *(list(p) for p in auxiliary)]
                return {
                    "template": "Beverland-AppA-three-cycle-H-with-square-scratch-contract",
                    "data_coord": [r, c],
                    "auxiliary_coords": [list(p) for p in auxiliary],
                    "reserved_patch_coords": footprint,
                    "phases": [
                        {
                            "relative_cycle": 0,
                            "action": "transversal_H_and_expand",
                            "reserved_patch_coords": footprint,
                        },
                        {
                            "relative_cycle": 1,
                            "action": "deform_and_rotate_boundaries",
                            "reserved_patch_coords": footprint,
                        },
                        {
                            "relative_cycle": 2,
                            "action": "contract_and_restore_native_orientation",
                            "reserved_patch_coords": footprint,
                        },
                    ],
                    "native_orientation_before": "control-vertical,target-horizontal",
                    "native_orientation_after": "control-vertical,target-horizontal",
                    "certification_scope": "allocation_and_serialization_under_declared_native_H_primitive; not_stabilizer_simulation",
                }
        raise ValueError(f"no reserved square scratch for H on {wire}")


def _layer_geometry(
    geometry: Geometry, layer: Mapping[str, Any], *, source_shift: bool
) -> tuple[list[dict[str, Any]], int]:
    mode = str(layer.get("mode", "")).lower()
    _require(mode in {"vdp", "edp"}, "a physical candidate must use VDP or EDP layers")
    used_data: set[int] = set()
    used_vertices: set[Coord] = set()
    used_edges: set[tuple[Coord, Coord]] = set()
    vertex_disjoint = True
    ops: list[dict[str, Any]] = []
    for original in layer["operations"]:
        control, target = (int(original["control"]), int(original["target"]))
        _require(
            control != target and (not {control, target} & used_data),
            "linear layer reuses a data endpoint",
        )
        used_data.update((control, target))
        raw = original["path"]
        if source_shift:
            raw = [[int(r) - 1, int(c) - 1] for r, c in raw]
        path = geometry.validate_path(raw, geometry.coord(control), geometry.coord(target))
        vertex_disjoint &= not bool(used_vertices & set(path))
        edges = {_edge(a, b) for a, b in zip(path, path[1:])}
        _require(not used_edges & edges, "linear layer repeats a routing edge")
        used_edges.update(edges)
        used_vertices.update(path)
        ops.append({"control": control, "target": target, "path": [list(p) for p in path]})
    _require(mode != "vdp" or vertex_disjoint, "VDP layer contains a vertex collision")
    duration = 2 if vertex_disjoint else 4
    _require(
        int(layer.get("logical_depth", duration)) == duration,
        "linear layer duration disagrees with replayed VDP/EDP geometry",
    )
    return (ops, duration)


def validate_candidate(
    data: Mapping[str, Any],
    choice: Mapping[str, Any],
    segments: Sequence[dict[str, Any]],
    geometry: Geometry,
) -> dict[str, Any]:
    """Replay matrix, full placement, directed geometry, and cost in memory."""
    _require(data.get("verified") is True, "candidate does not declare completed verification")
    _require(data.get("n") == geometry.width, "candidate logical width mismatch")
    layout = data.get("layout", {})
    _require(
        layout.get("data_rows") == geometry.layout.data_rows
        and layout.get("data_cols") == geometry.layout.data_cols,
        "candidate layout mismatch",
    )
    target, middle = linear_target_rows_from_span(
        segments, int(choice["start_pos"]), int(choice["end_pos"]), tuple(choice["input_placement"])
    )
    permutation = [int(q) for q in data["output_permutation"]]
    _require(sorted(permutation) == list(range(geometry.width)), "candidate permutation malformed")
    _require(permutation == choice["output_permutation"], "override changes the boundary placement")
    _require(
        replay_rows(dict(data)) == [target[q] for q in permutation],
        "candidate matrix replay failed",
    )
    _require(
        list(apply_output_permutation(middle, permutation)) == choice["output_placement"],
        "candidate output placement mismatch",
    )
    offsets = set()
    for layer in data["layers"]:
        for operation in layer["operations"]:
            _require(bool(operation.get("path")), "candidate lacks routed path")
            for raw, wire in (
                (operation["path"][0], int(operation["control"])),
                (operation["path"][-1], int(operation["target"])),
            ):
                endpoint = geometry.coord(wire)
                offsets.add((int(raw[0]) - endpoint[0], int(raw[1]) - endpoint[1]))
    _require(
        len(offsets) <= 1 and offsets <= {(0, 0), (1, 1)},
        "candidate mixes native and C++ route coordinate conventions",
    )
    coordinate_offset = next(iter(offsets), (1, 1))
    layers = []
    for index, layer in enumerate(data["layers"]):
        ops, duration = _layer_geometry(geometry, layer, source_shift=coordinate_offset == (1, 1))
        layers.append(
            {
                "kind": "linear_cnot_layer",
                "duration": duration,
                "category": "linear",
                "mode": layer["mode"],
                "operations": ops,
                "source_layer": index,
            }
        )
    cost = [
        sum((l["duration"] for l in layers)),
        len(layers),
        sum((len(l["operations"]) for l in layers)),
    ]
    stats = data["stats"]
    _require(
        cost == [int(stats["surface_depth"]), int(stats["layers"]), int(stats["cnots"])],
        "candidate statistics do not match the complete route replay",
    )
    return {
        "layers": layers,
        "cost": cost,
        "target_rows_hex": [hex(r) for r in target],
        "output_permutation": permutation,
        "source_coordinate_offset": list(coordinate_offset),
    }


def _touched(operation: Mapping[str, Any]) -> set[int]:
    return set((int(q) for q in operation.get("data_qubits", [])))


def cancel_disjoint_h_pairs(
    operations: Sequence[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Only H_q H_q=I; commute past gates with disjoint logical support."""
    result: list[dict[str, Any]] = []
    proofs: list[dict[str, Any]] = []
    for operation in operations:
        if operation["kind"] != "h":
            result.append(operation)
            continue
        wire = int(operation["wire"])
        previous = len(result) - 1
        while previous >= 0 and wire not in _touched(result[previous]):
            previous -= 1
        if (
            previous >= 0
            and result[previous]["kind"] == "h"
            and (int(result[previous]["wire"]) == wire)
        ):
            proofs.append(
                {
                    "rule": "H_q_H_q_identity_through_disjoint_support",
                    "wire": wire,
                    "removed_ids": [result[previous]["id"], operation["id"]],
                    "crossed_ids": [op["id"] for op in result[previous + 1 :]],
                }
            )
            del result[previous]
        else:
            result.append(operation)
    return (result, proofs)


def _h(wire: int, identifier: str, geometry: Geometry) -> dict[str, Any]:
    return {
        "id": identifier,
        "kind": "h",
        "wire": wire,
        "data_qubits": [wire],
        "duration": 3,
        "category": "hadamard",
        "footprint": geometry.h_footprint(wire),
    }


def _cz_expansion(
    pair: tuple[int, int],
    condition: str,
    identifier: str,
    target: int,
    geometry: Geometry,
    *,
    path: Sequence[Sequence[int]] | None = None,
) -> list[dict[str, Any]]:
    """Expand a CZ using an independently checked, fixed-native CNOT route.

    An omitted path retains the original deterministic routing choice.  A
    supplied path is a geometric witness, not a cost or certification claim;
    replay always checks its endpoints, native axes and reserved obstacles.
    No H basis or Pauli/Clifford frame is carried across this expansion.
    """
    _require(
        len(pair) == 2 and all((type(wire) is int for wire in pair)) and (pair[0] != pair[1]),
        "CZ support must contain two distinct integer wires",
    )
    _require(type(target) is int and target in pair, "CZ target must be an endpoint")
    control = pair[1] if target == pair[0] else pair[0]
    raw = geometry.cnot(control, target) if path is None else path
    _require(
        isinstance(raw, (list, tuple))
        and all(
            (
                isinstance(point, (list, tuple))
                and len(point) == 2
                and all((type(value) is int for value in point))
                for point in raw
            )
        ),
        "correction path must be a sequence of integer coordinate pairs",
    )
    routed = geometry.validate_path(raw, geometry.coord(control), geometry.coord(target))
    return [
        _h(target, f"{identifier}:H0", geometry),
        {
            "id": f"{identifier}:CX",
            "kind": "conditional_cnot",
            "control": control,
            "target": target,
            "data_qubits": [control, target],
            "path": [list(point) for point in routed],
            "condition": condition,
            "duration": 2,
            "category": "correction",
            "all_branch_reserved": True,
            "source_cz": identifier,
        },
        _h(target, f"{identifier}:H1", geometry),
    ]


def lower_cz_window(
    specifications: Sequence[dict[str, Any]],
    trailing_h: Sequence[int],
    geometry: Geometry,
    *,
    optimize: bool,
    identifier: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Search explicit CZ order/direction, not a deferred correction frame."""
    _require(len(specifications) <= 3, "local exact CZ window limited to three edges")
    orders = (
        itertools.permutations(range(len(specifications)))
        if optimize
        else [tuple(range(len(specifications)))]
    )
    best: tuple[Any, list[dict[str, Any]], dict[str, Any]] | None = None
    considered = 0
    for order in orders:
        directions = (
            itertools.product((0, 1), repeat=len(specifications))
            if optimize
            else [(1,) * len(specifications)]
        )
        for bits in directions:
            ops: list[dict[str, Any]] = []
            chosen = []
            for index in order:
                specification = specifications[index]
                pair = tuple((int(q) for q in specification["pair"]))
                target = pair[bits[index]]
                ops.extend(
                    _cz_expansion(
                        pair,
                        str(specification["condition"]),
                        str(specification["id"]),
                        target,
                        geometry,
                    )
                )
                chosen.append(
                    {
                        "id": specification["id"],
                        "pair": list(pair),
                        "target": target,
                        "condition": specification["condition"],
                    }
                )
            ops.extend(
                (_h(wire, f"{identifier}:trailing_H:{wire}", geometry) for wire in trailing_h)
            )
            simplified, cancellations = cancel_disjoint_h_pairs(ops) if optimize else (ops, [])
            considered += 1
            key = (
                sum((op["duration"] for op in simplified)),
                sum((len(op.get("path", [])) for op in simplified)),
                tuple(order),
                tuple(bits),
            )
            proof = {
                "rule": "commuting_diagonal_CZs_then_CZ_equals_H_CNOT_H",
                "chosen_czs": chosen,
                "cancellations": cancellations,
            }
            if best is None or key < best[0]:
                best = (key, simplified, proof)
    _require(best is not None, "empty correction search")
    proof = best[2]
    proof["variants_considered"] = considered
    proof["scope"] = "exact_for_this_three-edge_order-and-direction-window_only"
    return (best[1], proof)


def _validate_injection_group(geometry: Geometry, options: Sequence[dict[str, Any]]) -> None:
    used_vertices: set[Coord] = set()
    used_data: set[int] = set()
    used_ports: set[Coord] = set()
    used_resources: set[Coord] = set()
    for option in options:
        physical = [int(q) for q in option["physical_indices"]]
        ports = [_coord(p) for p in option["port_coords"]]
        _require(
            len(physical) == len(ports) == len(option["data_to_port_cnot_paths"]) == 3,
            "CCZ gadget needs three legs",
        )
        _require(
            len(set(physical)) == 3 and (not set(physical) & used_data), "CCZ group reuses data"
        )
        used_data.update(physical)
        triple = int(option["port_triple_index"])
        _require(
            0 <= triple < len(geometry.layout.port_triples())
            and set(ports) == set(geometry.layout.port_triples()[triple]),
            "CCZ ports are not one prepared resource triple",
        )
        resources = [_coord(p) for p in option["resource_coords"]]
        _require(
            resources == [resource_coord_for_port(geometry.layout, p) for p in ports],
            "CCZ resource is not at declared delivery site",
        )
        _require(
            not set(ports) & used_ports and (not set(resources) & used_resources),
            "CCZ resource or port reuse",
        )
        used_ports.update(ports)
        used_resources.update(resources)
        for q, port, raw in zip(physical, ports, option["data_to_port_cnot_paths"], strict=True):
            path = geometry.validate_path(
                raw,
                geometry.coord(q),
                port,
                goal_axis=geometry.port_axis(port),
                extra_banned=used_vertices,
            )
            used_vertices.update(path)


def _route_fixed_group(
    geometry: Geometry, options: Sequence[dict[str, Any]]
) -> list[dict[str, Any]] | None:
    """Deterministic multi-order directed VDP; failure means UNKNOWN."""
    legs = [
        (i, j, geometry.coord(int(option["physical_indices"][j])), _coord(option["port_coords"][j]))
        for i, option in enumerate(options)
        for j in range(3)
    ]
    natural = list(range(len(legs)))
    distances = [abs(a[0] - b[0]) + abs(a[1] - b[1]) for _, _, a, b in legs]
    orders = [
        natural,
        natural[::-1],
        sorted(natural, key=lambda j: distances[j]),
        sorted(natural, key=lambda j: (-distances[j], j)),
    ]
    if len(legs) == 3:
        orders = list(itertools.permutations(natural))
    for order in orders:
        for neighbor_order in range(2):
            used: set[Coord] = set()
            routes: dict[tuple[int, int], list[list[int]]] = {}
            for leg in order:
                i, j, start, goal = legs[leg]
                path = geometry.bfs(
                    start,
                    goal,
                    goal_axis=geometry.port_axis(goal),
                    used=used,
                    neighbor_order=neighbor_order,
                )
                if path is None:
                    break
                routes[i, j] = path
                used.update((_coord(p) for p in path))
            if len(routes) == len(legs):
                result = copy.deepcopy(list(options))
                for i, option in enumerate(result):
                    option["data_to_port_cnot_paths"] = [routes[i, j] for j in range(3)]
                    option["data_coords"] = [
                        list(geometry.coord(int(q))) for q in option["physical_indices"]
                    ]
                    option["resource_coords"] = [
                        list(resource_coord_for_port(geometry.layout, _coord(p)))
                        for p in option["port_coords"]
                    ]
                    option["route_length"] = sum(
                        (len(path) - 1 for path in option["data_to_port_cnot_paths"])
                    )
                    option["manhattan_score"] = sum(
                        (
                            abs(a[0] - b[0]) + abs(a[1] - b[1])
                            for a, b in zip(option["data_coords"], option["port_coords"])
                        )
                    )
                _validate_injection_group(geometry, result)
                return result
    return None


def route_injection_group(
    geometry: Geometry, source: Sequence[dict[str, Any]], *, optimize: bool = True
) -> tuple[list[list[dict[str, Any]]], dict[str, int]]:
    if not source:
        return ([], {"retained_groups": 0, "rerouted_groups": 0, "split_groups": 0})
    try:
        _validate_injection_group(geometry, source)
        return (
            [copy.deepcopy(list(source))],
            {"retained_groups": 1, "rerouted_groups": 0, "split_groups": 0},
        )
    except ValueError:
        pass
    routed = _route_fixed_group(geometry, source)
    if routed is not None:
        return ([routed], {"retained_groups": 0, "rerouted_groups": 1, "split_groups": 0})
    feasible: dict[int, list[dict[str, Any]]] = {}
    unknown_subsets = 0
    if optimize:
        for mask in range(1, (1 << len(source)) - 1):
            if mask.bit_count() < 2:
                continue
            subset = [option for i, option in enumerate(source) if mask & 1 << i]
            subset_routed = _route_fixed_group(geometry, subset)
            if subset_routed is not None:
                feasible[mask] = subset_routed
            else:
                unknown_subsets += 1
    for option_index, option in enumerate(source):
        routed = _route_fixed_group(geometry, [option])
        if routed is None:
            assignments = []
            physical = option["physical_indices"]
            for triple_index, triple in enumerate(geometry.layout.port_triples()):
                for permutation in itertools.permutations(range(3)):
                    ports = [triple[j] for j in permutation]
                    score = sum(
                        (
                            abs(geometry.coord(q)[0] - p[0]) + abs(geometry.coord(q)[1] - p[1])
                            for q, p in zip(physical, ports)
                        )
                    )
                    assignments.append((score, triple_index, permutation, ports))
            for _, triple_index, permutation, ports in sorted(assignments):
                alternative = copy.deepcopy(option)
                alternative.update(
                    port_triple_index=triple_index,
                    port_permutation=list(permutation),
                    port_coords=[list(p) for p in ports],
                )
                routed = _route_fixed_group(geometry, [alternative])
                if routed is not None:
                    break
        _require(
            routed is not None,
            f"directed three-leg CCZ routing unresolved for gate {option.get('gate_index')}; no timing result emitted",
        )
        feasible[1 << option_index] = routed
    full_mask = (1 << len(source)) - 1

    @lru_cache(maxsize=None)
    def cover(mask: int) -> tuple[int, ...]:
        if not mask:
            return ()
        first = mask & -mask
        candidates = []
        for subset, witness in feasible.items():
            if subset & first and subset & mask == subset:
                tail = cover(mask ^ subset)
                route_length = sum(
                    (
                        len(path) - 1
                        for m in (subset, *tail)
                        for option in feasible[m]
                        for path in option["data_to_port_cnot_paths"]
                    )
                )
                candidates.append((1 + len(tail), route_length, (subset, *tail)))
        _require(bool(candidates), "no materialized subset cover")
        return min(candidates)[2]

    groups = [feasible[mask] for mask in cover(full_mask)]
    return (
        groups,
        {
            "retained_groups": 0,
            "rerouted_groups": len(groups),
            "split_groups": 1,
            "materialized_subsets": len(feasible),
            "heuristic_unknown_subsets": unknown_subsets,
        },
    )


def _timeline(operations: list[dict[str, Any]]) -> None:
    now = 0
    for operation in operations:
        operation["start"] = now
        now += int(operation["duration"])
        operation["end"] = now


def _without_time(operation: Mapping[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in operation.items() if key not in {"start", "end"}}


def _chosen_cz_window(
    specifications: Sequence[dict[str, Any]],
    trailing_h: Sequence[int],
    rewrite: Mapping[str, Any],
    geometry: Geometry,
    *,
    optimized: bool,
    identifier: str,
) -> list[dict[str, Any]]:
    """Replay a chosen rewrite; do not trust its claimed cost or proof status."""
    expected = {item["id"]: item for item in specifications}
    chosen = rewrite["chosen_czs"]
    _require(
        len(chosen) == len(expected) and {item["id"] for item in chosen} == set(expected),
        "correction window drops or duplicates CZs",
    )
    operations = []
    for item in chosen:
        specification = expected[item["id"]]
        _require(
            item["pair"] == specification["pair"]
            and item["condition"] == specification["condition"],
            "chosen CZ changed its branch or support",
        )
        _require(
            "path" not in item or item["path"] is not None,
            "explicit correction path must not be null",
        )
        operations.extend(
            _cz_expansion(
                tuple(item["pair"]),
                item["condition"],
                item["id"],
                item["target"],
                geometry,
                path=item.get("path"),
            )
        )
    operations.extend(
        (_h(wire, f"{identifier}:trailing_H:{wire}", geometry) for wire in trailing_h)
    )
    if optimized:
        operations, cancellations = cancel_disjoint_h_pairs(operations)
    else:
        cancellations = []
    _require(
        cancellations == rewrite["cancellations"], "local HH cancellation witness does not replay"
    )
    return operations
