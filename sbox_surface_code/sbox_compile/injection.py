"""injection."""

# Adapted from the authors' strict_sbox/optimization_v13/injection.py; campaign infrastructure removed.
from __future__ import annotations
import copy
from .geometry import resource_coord_for_port, _validate_injection_group
from .model import canonical


def _route_group(context, options, geometry, rng, *, starts=4, fresh_only=False):
    """Construct a genuine new matching/path witness; bounded failure is UNKNOWN.

    Paths are routed jointly in randomized leg orders on the unchanged native
    geometry.  At least one gate's port matching is perturbed.  Optional banned
    scratch vertices diversify paths; they never relax a physical constraint.
    """
    changed = copy.deepcopy(options)
    triples = geometry.layout.port_triples()
    if len(changed) > len(triples):
        return None
    selected = (
        set(range(len(changed)))
        if fresh_only
        else set(rng.sample(range(len(changed)), 1 + rng.randrange(min(2, len(changed)))))
    )
    seen = set()
    for index, option in enumerate(changed):
        triple = option["port_triple_index"]
        if triple in seen:
            selected.add(index)
        else:
            seen.add(triple)
    reserved = {
        option["port_triple_index"] for index, option in enumerate(changed) if index not in selected
    }
    initial_ports = set(map(tuple, context["used_ports_before"]))
    for index in sorted(selected):
        option = changed[index]
        choices = [triple for triple in range(len(triples)) if triple not in reserved]
        if fresh_only:
            choices = [triple for triple in choices if not set(triples[triple]) & initial_ports]
        if not choices:
            return None
        rng.shuffle(choices)
        prefer_fresh = rng.randrange(2) == 0
        choices.sort(
            key=lambda triple: bool(set(triples[triple]) & initial_ports) != (not prefer_fresh)
        )
        shortlist = choices[: min(4, len(choices))]
        triple = rng.choice(shortlist)
        permutation = list(range(3))
        rng.shuffle(permutation)
        if triple == option["port_triple_index"] and permutation == option["port_permutation"]:
            permutation = permutation[1:] + permutation[:1]
        option.update(
            port_triple_index=triple,
            port_permutation=permutation,
            port_coords=[list(triples[triple][j]) for j in permutation],
        )
        reserved.add(triple)
    legs = [
        (index, j, geometry.coord(option["physical_indices"][j]), tuple(option["port_coords"][j]))
        for index, option in enumerate(changed)
        for j in range(3)
    ]
    endpoints = {point for _, _, start, goal in legs for point in (start, goal)}
    old_interior = sorted(
        {
            tuple(point)
            for option in options
            for path in option["data_to_port_cnot_paths"]
            for point in path[1:-1]
        }
        - endpoints
        - geometry.occupied
        - geometry.ports
        - geometry.resources
    )
    for attempt in range(starts):
        order = list(range(len(legs)))
        rng.shuffle(order)
        if attempt == 0:
            order.sort(
                key=lambda j: -(
                    abs(legs[j][2][0] - legs[j][3][0]) + abs(legs[j][2][1] - legs[j][3][1])
                )
            )
        detours = set(rng.sample(old_interior, min(len(old_interior), attempt % 3)))
        used = set(detours)
        paths = {}
        for leg in order:
            index, j, start, goal = legs[leg]
            path = geometry.bfs(
                start, goal, goal_axis="horizontal", used=used, neighbor_order=rng.randrange(4)
            )
            if path is None:
                break
            paths[index, j] = path
            used.update(map(tuple, path))
        if len(paths) != len(legs):
            continue
        for index, option in enumerate(changed):
            option["data_to_port_cnot_paths"] = [paths[index, j] for j in range(3)]
            option["resource_coords"] = [
                list(resource_coord_for_port(geometry.layout, tuple(port)))
                for port in option["port_coords"]
            ]
            option["route_length"] = sum(
                (len(path) - 1 for path in option["data_to_port_cnot_paths"])
            )
            option["manhattan_score"] = sum(
                (
                    abs(a[0] - b[0]) + abs(a[1] - b[1])
                    for a, b in zip(option["data_coords"], option["port_coords"], strict=True)
                )
            )
        _validate_injection_group(geometry, changed)
        return changed
    return None


def _ripup_group(context, options, geometry, rng, *, starts=6):
    """Reroute one gate while retaining every other gate's actual paths.

    The first trials keep the existing terminals and change path geometry;
    subsequent trials also change the selected gate's prepared-port matching.
    Thus a successful mutation never requires splitting the original group.
    """
    index = rng.randrange(len(options))
    original = options[index]
    original_hash = canonical(options)
    fixed = {
        tuple(point)
        for i, option in enumerate(options)
        if i != index
        for path in option["data_to_port_cnot_paths"]
        for point in path
    }
    reserved = {option["port_triple_index"] for i, option in enumerate(options) if i != index}
    triples = geometry.layout.port_triples()
    available = [i for i in range(len(triples)) if i not in reserved]
    interior = sorted(
        {tuple(point) for path in original["data_to_port_cnot_paths"] for point in path[1:-1]}
        - geometry.occupied
        - geometry.ports
        - geometry.resources
    )
    for attempt in range(starts):
        candidate = copy.deepcopy(original)
        if attempt >= 2:
            triple = rng.choice(available)
            permutation = list(range(3))
            rng.shuffle(permutation)
            candidate.update(
                port_triple_index=triple,
                port_permutation=permutation,
                port_coords=[list(triples[triple][j]) for j in permutation],
            )
        used = set(fixed)
        if attempt % 2 and interior:
            used.add(rng.choice(interior))
        order = list(range(3))
        rng.shuffle(order)
        paths = {}
        for leg in order:
            path = geometry.bfs(
                tuple(candidate["data_coords"][leg]),
                tuple(candidate["port_coords"][leg]),
                goal_axis="horizontal",
                used=used,
                neighbor_order=rng.randrange(4),
            )
            if path is None:
                break
            paths[leg] = path
            used.update(map(tuple, path))
        if len(paths) != 3:
            continue
        candidate["data_to_port_cnot_paths"] = [paths[j] for j in range(3)]
        candidate["resource_coords"] = [
            list(resource_coord_for_port(geometry.layout, tuple(port)))
            for port in candidate["port_coords"]
        ]
        candidate["route_length"] = sum(
            (len(path) - 1 for path in candidate["data_to_port_cnot_paths"])
        )
        candidate["manhattan_score"] = sum(
            (
                abs(a[0] - b[0]) + abs(a[1] - b[1])
                for a, b in zip(candidate["data_coords"], candidate["port_coords"], strict=True)
            )
        )
        changed = copy.deepcopy(options)
        changed[index] = candidate
        _validate_injection_group(geometry, changed)
        checked = changed
        if canonical(checked) != original_hash:
            return checked
    return None
