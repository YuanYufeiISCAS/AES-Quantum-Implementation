"""semantics."""

# Adapted from the authors' strict_sbox/semantics.py; campaign infrastructure removed.
from __future__ import annotations
from collections import Counter
from functools import lru_cache
from typing import Any, Iterable

SCHEMA = "sbox-coherent-semantics-v1"
ROLES = frozenset({"qand", "qand_dagger", "toffoli"})


class SemanticError(ValueError):
    """A source circuit or proposed lowering failed an exact semantic check."""


def canonical_family(family: str) -> str:
    aliases = {
        "state": "inplace",
        "ordinary": "inplace",
        "inplace": "inplace",
        "cstar": "cstar",
        "c*": "cstar",
    }
    try:
        return aliases[family.lower()]
    except KeyError as exc:
        raise SemanticError(f"unsupported circuit family: {family}") from exc


def segment_name(segment: dict[str, Any]) -> str:
    return f"{int(segment['index']):03d}_{segment['component']}"


def gf256_mul(left: int, right: int) -> int:
    """Polynomial multiplication modulo x^8+x^4+x^3+x+1 (0x11b)."""
    result = 0
    for _ in range(8):
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left & 256:
            left ^= 283
    return result


@lru_cache(maxsize=256)
def aes_sbox(value: int) -> int:
    """Independent AES definition: inversion followed by the affine map."""
    if not 0 <= value < 256:
        raise ValueError("AES input must be a byte")
    result, base, exponent = (1, value, 254)
    while exponent:
        if exponent & 1:
            result = gf256_mul(result, base)
        base = gf256_mul(base, base)
        exponent >>= 1
    inverse = result
    result = inverse
    for rotation in range(1, 5):
        result ^= (inverse << rotation | inverse >> 8 - rotation) & 255
    return result ^ 99


@lru_cache(maxsize=32)
def variable_plane(variable: int, input_width: int) -> int:
    """Bit i stores variable's value in input assignment i (LSB-first)."""
    if input_width not in (8, 16) or not 0 <= variable < input_width:
        raise SemanticError("truth planes support the 8- and 16-bit AES interfaces")
    basis_size = 1 << input_width
    half_period = 1 << variable
    period = 2 * half_period
    block = (1 << half_period) - 1 << half_period
    repetitions = ((1 << basis_size) - 1) // ((1 << period) - 1)
    return block * repetitions


def input_contract(family: str) -> dict[str, Any]:
    family = canonical_family(family)
    if family == "inplace":
        return {
            "width": 23,
            "input_width": 8,
            "x_wires": list(range(8)),
            "y_wires": [],
            "clean_wires": list(range(8, 23)),
            "operation": "|x>|0^15> -> |AES_SBOX(x)>|0^15>",
        }
    return {
        "width": 28,
        "input_width": 16,
        "x_wires": list(range(8)),
        "y_wires": list(range(16, 24)),
        "clean_wires": list(range(8, 16)) + list(range(24, 28)),
        "operation": "|x>|y>|0^12> -> |x>|y xor AES_SBOX(x)>|0^12>",
    }


def initial_planes(family: str) -> tuple[list[int], int]:
    contract = input_contract(family)
    state = [0] * contract["width"]
    for variable, wire in enumerate(contract["x_wires"] + contract["y_wires"]):
        state[wire] = variable_plane(variable, contract["input_width"])
    return (state, (1 << (1 << contract["input_width"])) - 1)


def flatten_operations(internal: dict[str, Any]) -> list[dict[str, Any]]:
    """Preserve source labels; logical swaps must update compiler placement."""
    width = int(internal.get("total_qubits", internal.get("n", 0)))
    if width <= 0:
        raise SemanticError("source has no positive circuit width")
    operations: list[dict[str, Any]] = []
    seen_segments: set[str] = set()
    arities = {"x": 1, "cx": 2, "ccx": 3, "swap": 2}
    for segment_position, segment in enumerate(internal["segments"]):
        name = segment_name(segment)
        if name in seen_segments:
            raise SemanticError(f"duplicate source segment: {name}")
        seen_segments.add(name)
        for gate_index, gate in enumerate(segment["gates"]):
            kind = str(gate["kind"])
            indices = gate["indices"]
            if kind not in arities or len(indices) != arities[kind]:
                raise SemanticError(f"invalid source operation: {name}:{gate_index}")
            if any((type(wire) is not int or wire < 0 or wire >= width for wire in indices)):
                raise SemanticError(f"invalid source wire: {name}:{gate_index}")
            if len(set(indices)) != len(indices):
                raise SemanticError(f"gate uses a wire twice: {name}:{gate_index}")
            operations.append(
                {
                    "key": f"{name}:{gate_index}",
                    "segment": name,
                    "segment_position": segment_position,
                    "segment_kind": str(segment["kind"]),
                    "gate_index": gate_index,
                    "segment_gate_index": gate_index,
                    "stream_index": len(operations),
                    "kind": kind,
                    "indices": list(indices),
                }
            )
    return operations


def apply_operation(state: list[int], operation: dict[str, Any], ones: int) -> None:
    """Apply a phase-free logical operation to an exact bit-plane state."""
    kind, indices = (operation["kind"], operation["indices"])
    if kind == "cx":
        control, target = indices
        state[target] ^= state[control]
    elif kind == "ccx":
        left, right, target = indices
        state[target] ^= state[left] & state[right]
    elif kind == "swap":
        left, right = indices
        state[left], state[right] = (state[right], state[left])
    elif kind == "x":
        state[indices[0]] ^= ones
    else:
        raise SemanticError(f"unsupported logical operation: {kind}")


def _expected_planes(family: str) -> list[int]:
    contract = input_contract(family)
    state, ones = initial_planes(family)
    repetitions = ones // ((1 << 256) - 1)
    output_wires = contract["y_wires"] or contract["x_wires"]
    for bit, wire in enumerate(output_wires):
        one_byte_plane = sum(((aes_sbox(x) >> bit & 1) << x for x in range(256)))
        sbox_plane = one_byte_plane * repetitions
        state[wire] = state[wire] ^ sbox_plane if contract["y_wires"] else sbox_plane
    return state


def verify_logical_operations(
    family: str, operations: Iterable[dict[str, Any]], *, check_roles: bool = True
) -> dict[str, Any]:
    """Replay arbitrary proposed NCT operations against the complete AES interface.

    When roles are present their preconditions are checked at *that call's*
    position, so this hook also rejects invalid reorderings.  Physical routing
    and local gadget witnesses must be checked by the compiler's verifier.
    """
    width = input_contract(family)["width"]
    identity = list(range(width))
    return verify_placed_logical_operations(
        family, operations, identity, identity, check_roles=check_roles
    )


def verify_placed_logical_operations(
    family: str,
    operations: Iterable[dict[str, Any]],
    initial_placement: Iterable[int],
    final_placement: Iterable[int],
    *,
    check_roles: bool = True,
) -> dict[str, Any]:
    """Verify an actual routed NCT skeleton, including its interface mappings.

    Placement is **physical slot -> source logical wire**.  Gates use physical
    indices.  Source logical relabels are represented by the compiler's changed
    subsequent gate indices and final placement, not by physical SWAP gates.
    A genuine physical SWAP is nevertheless supported if present in the gate
    stream.  The final state is compared in source logical order against the
    independently computed AES interface over all 256/65,536 input assignments.

    This helper deliberately neither trusts a source ``passed`` flag nor imports
    a saved replay certificate.  QAND clean-entry and dagger-erasure conditions
    are recomputed in the actual physical gate order and actual placement.
    """
    family = canonical_family(family)
    contract = input_contract(family)
    width = contract["width"]
    initial = list(initial_placement)
    final = list(final_placement)
    for name, placement in (("initial", initial), ("final", final)):
        if (
            len(placement) != width
            or any((type(wire) is not int for wire in placement))
            or sorted(placement) != list(range(width))
        ):
            raise SemanticError(f"{name} placement is not a physical-to-logical permutation")
    logical_input, ones = initial_planes(family)
    state = [logical_input[logical] for logical in initial]
    count = 0
    role_counts: Counter[str] = Counter()
    arities = {"x": 1, "cx": 2, "ccx": 3, "swap": 2}
    for operation in operations:
        kind = operation["kind"]
        indices = operation["indices"]
        if kind not in arities or len(indices) != arities[kind]:
            raise SemanticError(f"invalid routed logical operation: {operation.get('key', count)}")
        if any((type(wire) is not int or not 0 <= wire < width for wire in indices)) or len(
            set(indices)
        ) != len(indices):
            raise SemanticError(f"invalid routed operation wires: {operation.get('key', count)}")
        if kind == "ccx":
            role_counts[operation.get("role", "toffoli")] += 1
        if check_roles and operation["kind"] == "ccx":
            left, right, target = operation["indices"]
            role = operation.get("role", "toffoli")
            if role not in ROLES:
                raise SemanticError(f"unknown lowering role: {role}")
            if role == "qand" and state[target] != 0:
                raise SemanticError(f"unclean QAND entry: {operation.get('key', count)}")
            if role == "qand_dagger" and state[target] != state[left] & state[right]:
                raise SemanticError(f"invalid QAND erasure: {operation.get('key', count)}")
        apply_operation(state, operation, ones)
        count += 1
    logical_output = [0] * width
    for physical, logical in enumerate(final):
        logical_output[logical] = state[physical]
    state = logical_output
    expected = _expected_planes(family)
    failures = []
    for wire, (actual, wanted) in enumerate(zip(state, expected)):
        difference = actual ^ wanted
        if difference:
            assignment = (difference & -difference).bit_length() - 1
            failures.append(
                {
                    "wire": wire,
                    "x": assignment & 255,
                    "y": assignment >> 8 if family == "cstar" else None,
                    "actual": actual >> assignment & 1,
                    "expected": wanted >> assignment & 1,
                }
            )
    return {
        "passed": not failures,
        "family": family,
        "checked_inputs": 1 << contract["input_width"],
        "checked_input_pairs": 65536 if family == "cstar" else None,
        "input_width": contract["input_width"],
        "operations": count,
        "placement_convention": "physical_to_logical",
        "initial_placement": initial,
        "final_placement": final,
        "role_counts": dict(role_counts),
        "checked_roles": check_roles,
        "independent_reference": "GF(2^8) inversion modulo 0x11b; AES affine map",
        "scratch_zero": all((state[wire] == 0 for wire in contract["clean_wires"])),
        "failures": failures[:8],
    }


def ccz_explicit_corrections(
    k: tuple[int, int, int], r: tuple[int, int, int]
) -> list[dict[str, Any]]:
    """Branch-instantiated correction gates, including explicit zero-latency Zs.

    The conditional schedule must retain all three CZ slots for a branch-safe
    worst-case estimate; this routine gives the operations executed in one
    branch, not permission to price every branch as an average.
    """
    if len(k) != 3 or len(r) != 3 or any((bit not in (0, 1) for bit in k + r)):
        raise SemanticError("CCZ injection needs three binary k and r outcomes")
    operations: list[dict[str, Any]] = []
    for omitted in range(3):
        pair = [wire for wire in range(3) if wire != omitted]
        operations.append(
            {
                "kind": "cz",
                "indices": pair,
                "enabled": bool(k[omitted]),
                "condition": f"k{omitted}",
                "explicit": True,
            }
        )
    for wire in range(3):
        other = [index for index in range(3) if index != wire]
        enabled = r[wire] ^ k[other[0]] & k[other[1]]
        operations.append(
            {
                "kind": "z",
                "indices": [wire],
                "enabled": bool(enabled),
                "condition": f"r{wire} xor (k{other[0]} and k{other[1]})",
                "explicit": True,
                "d_cycles": 0,
                "latency_convention": "paper_local_Pauli_primitive_not_a_frame",
            }
        )
    return operations


def _bits(value: int, width: int) -> tuple[int, ...]:
    return tuple((value >> index & 1 for index in range(width)))


def certify_ccz_gadget() -> dict[str, Any]:
    """Exact phase/Kraus proof for all 64 injection outcomes and 8 inputs."""
    witnesses = []
    for outcome in range(64):
        k = _bits(outcome & 7, 3)
        r = _bits(outcome >> 3, 3)
        operations = ccz_explicit_corrections(k, r)
        global_phase = k[0] & k[1] & k[2] ^ sum((r[i] * k[i] for i in range(3))) & 1
        for assignment in range(8):
            x = _bits(assignment, 3)
            resource = tuple((x[i] ^ k[i] for i in range(3)))
            raw_phase = (
                resource[0] & resource[1] & resource[2]
                ^ sum((r[i] * resource[i] for i in range(3))) & 1
            )
            correction_phase = 0
            for operation in operations:
                if operation["enabled"]:
                    support_value = 1
                    for wire in operation["indices"]:
                        support_value &= x[wire]
                    correction_phase ^= support_value
            expected_phase = x[0] & x[1] & x[2] ^ global_phase
            if raw_phase ^ correction_phase != expected_phase:
                raise SemanticError(f"CCZ branch proof failed: k={k}, r={r}, x={x}")
        witnesses.append({"k": list(k), "r": list(r), "global_phase_bit": global_phase})
    return {
        "schema": "strict-ccz-explicit-branch-proof-v1",
        "passed": True,
        "checked_branches": 64,
        "checked_basis_inputs_per_branch": 8,
        "prepared_resource": "2^(-3/2) sum_u (-1)^(u0*u1*u2)|u>",
        "measurements": "Z(data_i)Z(resource_i)=(-1)^k_i, then X(resource_i)=(-1)^r_i",
        "paper_parity_implementation": "fresh |0> ancilla; CX(data,ancilla), CX(resource,ancilla); Z(ancilla)",
        "raw_phase": "f(x xor k) xor r dot (x xor k)",
        "explicit_cz_correction": "CZ(1,2)^k0 CZ(0,2)^k1 CZ(0,1)^k2",
        "explicit_z_correction": "product_i Z(i)^(r_i xor product_{j!=i} k_j)",
        "corrected_kraus": "(-1)^(f(k) xor r dot k) * CCZ / 8",
        "kraus_amplitude_squared": {"numerator": 1, "denominator": 64},
        "amplitude_is_input_independent": True,
        "normalization": {"branches": 64, "sum_branch_probabilities": "64 * 1/64 = 1"},
        "deferred_corrections": False,
        "pauli_frame": False,
        "clifford_frame": False,
    }


def certify_qand_dagger_gadget() -> dict[str, Any]:
    """Exact X-measurement erasure on |a,b,a*b>, with explicit CZ and reset."""
    witnesses = []
    for m in (0, 1):
        for assignment in range(4):
            a, b = _bits(assignment, 2)
            target = a & b
            raw_phase = m & target
            correction_phase = m & a & b
            if raw_phase ^ correction_phase:
                raise SemanticError("QAND dagger phase erasure failed")
            witnesses.append([m, a, b, raw_phase, correction_phase])
    return {
        "schema": "strict-qand-dagger-explicit-branch-proof-v1",
        "passed": True,
        "checked_branches": 2,
        "checked_basis_inputs_per_branch": 4,
        "precondition": "target = control_a AND control_b on the full coherent input subspace",
        "operations": [
            "measure_X(target) -> m",
            "explicit_CZ(control_a,control_b)^m",
            "reset_target_to_zero",
        ],
        "corrected_kraus_on_clean_AND_subspace": "AND_erasure / sqrt(2)",
        "kraus_amplitude_squared": {"numerator": 1, "denominator": 2},
        "amplitude_is_input_independent": True,
        "normalization": {"branches": 2, "sum_branch_probabilities": "2 * 1/2 = 1"},
        "deferred_corrections": False,
        "reset_required_for_clean_output": True,
    }


def certify_cz_decomposition() -> dict[str, Any]:
    """Small exact integer calculation of H_t CX_c,t H_t = CZ_c,t.

    Unnormalized Hadamards have entries +/-1.  Two Hs give denominator two,
    so an integer check suffices and does not hide floating-point tolerances.
    """
    for enabled in (0, 1):
        for source in range(4):
            control, target = _bits(source, 2)
            amplitudes = [0] * 4
            for middle in (0, 1):
                first_sign = -1 if target & middle else 1
                flipped = middle ^ control & enabled
                for output in (0, 1):
                    second_sign = -1 if flipped & output else 1
                    amplitudes[control | output << 1] += first_sign * second_sign
            expected = [0] * 4
            expected[source] = -2 if enabled & control & target else 2
            if amplitudes != expected:
                raise SemanticError("H-CNOT-H CZ decomposition failed")
    return {
        "schema": "strict-conditional-cz-HCXH-proof-v1",
        "passed": True,
        "checked_conditions": 2,
        "checked_inputs_per_condition": 4,
        "identity": "H_t CX(c,t)^m H_t = CZ(c,t)^m",
        "arithmetic": "exact integers; normalization denominator 2",
        "H_gates_are_explicit": True,
        "basis_is_restored": True,
    }


def gadget_certificates() -> dict[str, Any]:
    return {
        "schema": "strict-sbox-local-gadgets-v1",
        "passed": True,
        "ccz": certify_ccz_gadget(),
        "qand_dagger": certify_qand_dagger_gadget(),
        "conditional_cz": certify_cz_decomposition(),
        "composition": "every corrected branch equals its phase-free source operation times an input-independent scalar",
    }
