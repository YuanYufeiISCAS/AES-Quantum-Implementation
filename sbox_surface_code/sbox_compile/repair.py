"""repair."""

# Adapted from the authors' strict_sbox/optimization_v15/windows.py; campaign infrastructure removed.
from __future__ import annotations
import copy
import itertools
from . import geometry as front
from .geometry import _require

MAX_REPAIR_ROUTE_ATTEMPTS = 48


def _repair_injection_group(geometry, original, changed):
    """Pin valid unmoved legs; try bounded repairs before full group rerouting."""
    front._validate_injection_group(geometry, original)
    _require(
        len(original) == len(changed)
        and [option["source_id"] for option in original]
        == [option["source_id"] for option in changed],
        "injection repair changed the source group cover/order",
    )
    try:
        front._validate_injection_group(geometry, changed)
        return (
            [copy.deepcopy(list(changed))],
            {
                "retained_groups": 1,
                "rerouted_groups": 0,
                "split_groups": 0,
                "repair_attempts": 0,
                "preserved_unchanged_legs": 3 * len(changed),
                "repaired_moved_legs": 0,
                "repair_fallback_groups": 0,
            },
        )
    except ValueError:
        pass
    paths, fixed_vertices, moved = ({}, set(), [])
    for index, (before, after) in enumerate(zip(original, changed, strict=True)):
        _require(
            len(after["physical_indices"]) == len(after["port_coords"]) == 3,
            "injection repair requires exactly three data/port legs",
        )
        for leg in range(3):
            start = geometry.coord(after["physical_indices"][leg])
            goal = tuple(after["port_coords"][leg])
            if (
                before["physical_indices"][leg] == after["physical_indices"][leg]
                and before["port_coords"][leg] == after["port_coords"][leg]
            ):
                path = geometry.validate_path(
                    before["data_to_port_cnot_paths"][leg],
                    start,
                    goal,
                    goal_axis=geometry.port_axis(goal),
                    extra_banned=fixed_vertices,
                )
                paths[index, leg] = [list(point) for point in path]
                fixed_vertices.update(path)
            else:
                moved.append((index, leg, start, goal))
    attempts = 0
    if paths and moved:
        natural = tuple(range(len(moved)))
        distances = [
            abs(start[0] - goal[0]) + abs(start[1] - goal[1]) for _, _, start, goal in moved
        ]
        if len(moved) <= 3:
            orders = list(itertools.permutations(natural))
        else:
            orders = list(
                dict.fromkeys(
                    (
                        natural,
                        natural[::-1],
                        tuple(sorted(natural, key=lambda j: (distances[j], j))),
                        tuple(sorted(natural, key=lambda j: (-distances[j], j))),
                    )
                )
            )
        for order in orders:
            for neighbor_order in range(4):
                if attempts >= MAX_REPAIR_ROUTE_ATTEMPTS:
                    break
                attempts += 1
                routed, used = (copy.deepcopy(paths), set(fixed_vertices))
                for slot in order:
                    index, leg, start, goal = moved[slot]
                    path = geometry.bfs(
                        start,
                        goal,
                        goal_axis=geometry.port_axis(goal),
                        used=used,
                        neighbor_order=neighbor_order,
                    )
                    if path is None:
                        break
                    routed[index, leg] = path
                    used.update(map(tuple, path))
                if len(routed) != 3 * len(changed):
                    continue
                result = copy.deepcopy(list(changed))
                for index, option in enumerate(result):
                    option["data_to_port_cnot_paths"] = [
                        copy.deepcopy(routed[index, j]) for j in range(3)
                    ]
                    option["data_coords"] = [
                        list(geometry.coord(wire)) for wire in option["physical_indices"]
                    ]
                    option["resource_coords"] = [
                        list(front.resource_coord_for_port(geometry.layout, tuple(port)))
                        for port in option["port_coords"]
                    ]
                    option["route_length"] = sum(
                        (len(path) - 1 for path in option["data_to_port_cnot_paths"])
                    )
                    option["manhattan_score"] = sum(
                        (
                            abs(a[0] - b[0]) + abs(a[1] - b[1])
                            for a, b in zip(
                                option["data_coords"], option["port_coords"], strict=True
                            )
                        )
                    )
                front._validate_injection_group(geometry, result)
                return (
                    [result],
                    {
                        "retained_groups": 0,
                        "rerouted_groups": 1,
                        "split_groups": 0,
                        "repair_attempts": attempts,
                        "preserved_unchanged_legs": len(paths),
                        "repaired_moved_legs": len(moved),
                        "repair_fallback_groups": 0,
                    },
                )
    groups, stats = front.route_injection_group(geometry, changed, optimize=True)
    _require(
        sorted((option["source_id"] for group in groups for option in group))
        == sorted((option["source_id"] for option in changed)),
        "fallback injection repair changed the exact gadget cover",
    )
    for group in groups:
        front._validate_injection_group(geometry, group)
    return (
        groups,
        {
            **stats,
            "repair_attempts": attempts,
            "preserved_unchanged_legs": 0,
            "repaired_moved_legs": 0,
            "repair_fallback_groups": 1,
        },
    )
