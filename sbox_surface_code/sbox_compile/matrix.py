"""matrix."""

# Adapted from the authors' sbox/selector/permutation_aware_multiblock.py; campaign infrastructure removed.
from __future__ import annotations
from typing import Dict, List, Optional, Sequence, Tuple

Placement = tuple[int, ...]


def identity_rows(n: int) -> List[int]:
    return [1 << i for i in range(n)]


def invert_perm(perm: Sequence[int]) -> List[int]:
    inv = [0 for _ in perm]
    for idx, value in enumerate(perm):
        inv[int(value)] = idx
    return inv


def apply_output_permutation(placement: Placement, output_perm: Sequence[int]) -> Placement:
    return tuple((placement[int(idx)] for idx in output_perm))


def swap_placement_values(placement: Placement, a: int, b: int) -> Placement:
    out = list(placement)
    for idx, value in enumerate(out):
        if value == a:
            out[idx] = b
        elif value == b:
            out[idx] = a
    return tuple(out)


def apply_cx_to_rows(rows: List[int], control: int, target: int) -> None:
    rows[target] ^= rows[control]


def replay_rows(data: Dict[str, object]) -> Optional[List[int]]:
    n = data.get("n")
    layers = data.get("layers")
    if not isinstance(n, int) or not isinstance(layers, list):
        return None
    rows = identity_rows(n)
    for layer in layers:
        if not isinstance(layer, dict):
            return None
        ops = layer.get("operations")
        if not isinstance(ops, list):
            return None
        for op in ops:
            if not isinstance(op, dict):
                return None
            control = int(op["control"])
            target = int(op["target"])
            if control == target or control < 0 or target < 0 or (control >= n) or (target >= n):
                return None
            apply_cx_to_rows(rows, control, target)
    return rows


def apply_linear_gates_to_rows(
    rows: List[int], seg: Dict[str, object], placement: Placement
) -> None:
    logical_to_physical = invert_perm(placement)
    gates = seg.get("gates", [])
    if not isinstance(gates, list):
        gates = []
    for gate in gates:
        if not isinstance(gate, dict):
            raise ValueError("malformed linear gate")
        if gate.get("kind") != "cx":
            raise ValueError(f"linear segment contains non-CX gate: {gate.get('kind')}")
        indices = gate.get("indices", [])
        if not isinstance(indices, list) or len(indices) != 2:
            raise ValueError("CX gate does not have two indices")
        control = logical_to_physical[int(indices[0])]
        target = logical_to_physical[int(indices[1])]
        apply_cx_to_rows(rows, control, target)


def linear_target_rows_from_span(
    segments: Sequence[Dict[str, object]], start: int, end: int, placement: Placement
) -> Tuple[Tuple[int, ...], Placement]:
    """Build the physical-row target for a linear/relabel span.

    The returned placement includes relabels inside the span.  Candidate
    output_permutation is applied by the caller after this placement, matching
    the end-of-span output frame.
    """
    rows = identity_rows(len(placement))
    out_placement = placement
    for pos in range(start, end + 1):
        seg = segments[pos]
        kind = str(seg.get("kind", ""))
        if kind == "linear":
            apply_linear_gates_to_rows(rows, seg, out_placement)
        elif kind == "relabel":
            out_placement = update_placement_for_relabel(seg, out_placement)
        else:
            raise ValueError(f"cannot merge segment kind into linear span: {kind}")
    return (tuple(rows), out_placement)


def update_placement_for_relabel(seg: Dict[str, object], placement: Placement) -> Placement:
    out = placement
    gates = seg.get("gates", [])
    if not isinstance(gates, list):
        gates = []
    for gate in gates:
        if not isinstance(gate, dict):
            raise ValueError("malformed relabel gate")
        indices = gate.get("indices", [])
        if not isinstance(indices, list) or len(indices) != 2:
            raise ValueError("relabel gate does not have two indices")
        out = swap_placement_values(out, int(indices[0]), int(indices[1]))
    return out
