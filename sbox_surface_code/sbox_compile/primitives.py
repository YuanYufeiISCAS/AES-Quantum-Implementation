"""primitives."""

# Adapted from the authors' strict_sbox/primitive_contracts.py; campaign infrastructure removed.
from __future__ import annotations
from copy import deepcopy
from functools import lru_cache
import itertools
from typing import Any, Sequence

SCHEMA = "native-sbox-primitive-contracts-v1"
ORIENTATION = "standard_Z_vertical_X_horizontal"
Coord = tuple[int, int]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _coord(value: Sequence[int]) -> Coord:
    _require(
        len(value) == 2 and all((type(v) is int for v in value)),
        "coordinate must contain two integers",
    )
    return (value[0], value[1])


def _xor(*expressions: Any) -> dict[str, Any]:
    """Keep classical parity expressions explicit, including empty parity zero."""
    return {"xor": list(expressions)}


def evaluate_expression(expression: Any, outcomes: dict[str, int]) -> int:
    if isinstance(expression, str):
        value = outcomes[expression]
        _require(type(value) is int and value in (0, 1), "nonbinary measurement outcome")
        return value
    if type(expression) is int:
        _require(expression in (0, 1), "nonbinary classical constant")
        return expression
    _require(
        isinstance(expression, dict) and set(expression) == {"xor"},
        "unsupported classical expression",
    )
    result = 0
    for part in expression["xor"]:
        result ^= evaluate_expression(part, outcomes)
    return result


def _axis(left: Coord, right: Coord) -> str:
    _require(
        abs(left[0] - right[0]) + abs(left[1] - right[1]) == 1,
        "native joint parity requires adjacent patches",
    )
    return "vertical" if left[1] == right[1] else "horizontal"


def _pauli(
    patch: Coord, kind: str, expression: Any, identifier: str, at_cycle: int, condition: str | None
) -> dict[str, Any]:
    return {
        "id": identifier,
        "kind": f"explicit_conditional_{kind}",
        "patch": list(patch),
        "condition_expression": expression,
        "execute_when": condition,
        "relative_cycle": at_cycle,
        "duration": 0,
        "execution": "physical_Pauli_pulse_before_next_dependent_stage",
        "d_scaled_latency": 0,
        "physical_pulse_latency_is_zero": False,
        "deferred_frame": False,
    }


def _bell_measurement(
    left: Coord,
    right: Coord,
    identifier: str,
    k_name: str,
    r_name: str,
    *,
    relative_cycle: int,
    readout_cycle: int,
    ready_cycle: int,
    condition: str | None,
) -> dict[str, Any]:
    """The local Bell instrument, with all three fine-grained output bits."""
    direction = _axis(left, right)
    joint = "ZZ" if direction == "vertical" else "XX"
    readout = "X" if joint == "ZZ" else "Z"
    joint_name = f"{identifier}:joint"
    left_readout = f"{identifier}:readout_left"
    right_readout = f"{identifier}:readout_right"
    k_expression = joint_name if joint == "ZZ" else _xor(left_readout, right_readout)
    r_expression = joint_name if joint == "XX" else _xor(left_readout, right_readout)
    return {
        "id": identifier,
        "kind": "destructive_native_Bell_measurement",
        "patches": [list(left), list(right)],
        "native_orientation": ORIENTATION,
        "native_axis": direction,
        "execute_when": condition,
        "joint": {
            "kind": f"measure_{joint}",
            "patches": [list(left), list(right)],
            "outcome": joint_name,
            "relative_cycle": relative_cycle,
            "duration": 1,
        },
        "readouts": [
            {
                "kind": f"measure_{readout}",
                "patch": list(left),
                "outcome": left_readout,
                "relative_cycle": readout_cycle,
            },
            {
                "kind": f"measure_{readout}",
                "patch": list(right),
                "outcome": right_readout,
                "relative_cycle": readout_cycle,
            },
        ],
        "logical_outputs": [
            {
                "name": k_name,
                "meaning": "ZZ_Bell_parity",
                "expression": k_expression,
                "ready_relative_cycle": ready_cycle,
            },
            {
                "name": r_name,
                "meaning": "XX_Bell_parity",
                "expression": r_expression,
                "ready_relative_cycle": ready_cycle,
            },
        ],
        "raw_outcomes": [joint_name, left_readout, right_readout],
        "destructive_readout_folded_into_joint_window": readout_cycle == relative_cycle,
        "coherent_instrument": "each fine branch is a Bell bra times a known input-independent phase divided by sqrt(2)",
        "deferred_frame": False,
    }


def cnot_contract(
    control: Sequence[int],
    target: Sequence[int],
    path: Sequence[Sequence[int]],
    id: str,
    condition: str | None = None,
) -> dict[str, Any]:
    """Materialize one uniform-orientation, depth-two CNOT corridor.

    Endpoints use coordinates, not logical-wire indices.  The caller checks
    global occupancy and vertex-disjoint parallel corridors.  This contract
    checks every internal location and native interface.  EDP crossings are
    not certified by putting individual contracts in the same time slot.

    For m=2p+1 internal patches: stage 1 measures ZZ(c,a1) and prepares p
    corrected Bell pairs (a2,a3), ... .  Stage 2 Bell-measures (a1,a2), ...,
    measures XX(am,t), and reads am in Z.  All endpoint Pauli corrections are
    physically executed before the next dependent primitive.
    """
    control = _coord(control)
    target = _coord(target)
    route = [_coord(p) for p in path]
    _require(isinstance(id, str) and bool(id), "CNOT needs a stable identifier")
    _require(
        condition is None or (isinstance(condition, str) and bool(condition)),
        "invalid execution condition",
    )
    _require(
        len(route) >= 3 and route[0] == control and (route[-1] == target),
        "CNOT corridor requires its two endpoints and at least one actual ancilla",
    )
    _require(len(set(route)) == len(route), "CNOT corridor must be simple")
    for left, right in zip(route, route[1:]):
        _axis(left, right)
    _require(_axis(route[0], route[1]) == "vertical", "CNOT control must expose native vertical ZZ")
    _require(
        _axis(route[-2], route[-1]) == "horizontal", "CNOT target must expose native horizontal XX"
    )
    ancillas = route[1:-1]
    _require(
        len(ancillas) % 2 == 1,
        "two-stage Bell-chain template needs an odd number of internal ancillas",
    )
    pair_count = (len(ancillas) - 1) // 2
    u, v, w = (f"{id}:{label}" for label in ("u", "v", "w"))
    first_stage = [
        {
            "id": f"{id}:prepare_a1",
            "kind": "prepare_plus",
            "patch": list(ancillas[0]),
            "relative_cycle": 0,
            "duration": 0,
            "protection": "first_parity_syndrome_window",
            "execute_when": condition,
        },
        {
            "id": f"{id}:control_parity",
            "kind": "measure_ZZ",
            "patches": [list(control), list(ancillas[0])],
            "outcome": u,
            "relative_cycle": 0,
            "duration": 1,
            "execute_when": condition,
        },
    ]
    second_stage = []
    explicit_paulis = []
    outcome_names = [u, v, w]
    teleport_k, teleport_r = ([], [])
    for pair_index in range(pair_count):
        left, right = ancillas[2 * pair_index + 1 : 2 * pair_index + 3]
        direction = _axis(left, right)
        parity = "ZZ" if direction == "vertical" else "XX"
        basis = "+" if parity == "ZZ" else "0"
        b = f"{id}:bell_prepare_{pair_index}:b"
        first_stage.append(
            {
                "id": f"{id}:bell_prepare_{pair_index}",
                "kind": "prepare_native_Bell_pair",
                "patches": [list(left), list(right)],
                "initial_states": [basis, basis],
                "native_axis": direction,
                "joint_parity": parity,
                "outcome": b,
                "relative_cycle": 0,
                "duration": 1,
                "execute_when": condition,
                "protection": "Bell_preparation_parity_syndrome_window",
            }
        )
        explicit_paulis.append(
            _pauli(
                right,
                "X" if parity == "ZZ" else "Z",
                b,
                f"{id}:bell_prepare_{pair_index}:physical_correction",
                1,
                condition,
            )
        )
        outcome_names.append(b)
        measure_left, measure_right = ancillas[2 * pair_index : 2 * pair_index + 2]
        k, r = (f"{id}:teleport_{pair_index}:k", f"{id}:teleport_{pair_index}:r")
        measurement = _bell_measurement(
            measure_left,
            measure_right,
            f"{id}:teleport_{pair_index}",
            k,
            r,
            relative_cycle=1,
            readout_cycle=1,
            ready_cycle=2,
            condition=condition,
        )
        second_stage.append(measurement)
        outcome_names.extend(measurement["raw_outcomes"])
        teleport_k.append(k)
        teleport_r.append(r)
    second_stage.extend(
        [
            {
                "id": f"{id}:target_parity",
                "kind": "measure_XX",
                "patches": [list(ancillas[-1]), list(target)],
                "outcome": v,
                "relative_cycle": 1,
                "duration": 1,
                "execute_when": condition,
            },
            {
                "id": f"{id}:last_ancilla_readout",
                "kind": "measure_Z",
                "patch": list(ancillas[-1]),
                "outcome": w,
                "relative_cycle": 1,
                "folded_into_final_parity_window": True,
                "execute_when": condition,
            },
        ]
    )
    z_expression = _xor(v, *teleport_r)
    x_expression = _xor(u, w, *teleport_k)
    explicit_paulis.extend(
        [
            _pauli(control, "Z", z_expression, f"{id}:physical_Z_control", 2, condition),
            _pauli(target, "X", x_expression, f"{id}:physical_X_target", 2, condition),
        ]
    )
    return {
        "schema": "native-cnot-bell-chain-contract-v1",
        "id": id,
        "control": list(control),
        "target": list(target),
        "path": [list(p) for p in route],
        "condition": condition,
        "native_orientation": ORIENTATION,
        "reserved_cycles": 2,
        "outputs_ready_relative_cycle": 2,
        "footprint": [list(p) for p in route],
        "ancilla_coords": [list(p) for p in ancillas],
        "ancilla_count": len(ancillas),
        "bell_pair_count": pair_count,
        "paper_conservative_ancilla_upper_bound": len(route) - 1,
        "extra_endpoint_stitching_patch_required": False,
        "stitching_justification": "explicit odd-internal-vertex Bell-chain template includes every endpoint ancilla; m=1 is the one-ancilla local CNOT",
        "stages": [
            {"relative_cycle": 0, "duration": 1, "operations": first_stage},
            {"relative_cycle": 1, "duration": 1, "operations": second_stage},
        ],
        "explicit_paulis": explicit_paulis,
        "raw_outcomes": outcome_names,
        "endpoint_corrections": {"Z_control": z_expression, "X_target": x_expression},
        "disabled_branch": (
            {
                "condition_value": 0,
                "data_action": "identity",
                "execute_no_quantum_suboperations": True,
                "raw_outcome_convention": "all_zero_when_not_measured",
                "same_two_cycle_slot_reserved": True,
            }
            if condition is not None
            else None
        ),
        "measurement_branch_count_when_enabled": 1 << 3 + 4 * pair_count,
        "corrected_branch_probability": f"1/{1 << 3 + 4 * pair_count}",
        "branch_identity": "each corrected enabled Kraus branch equals input-independent scalar times CNOT",
        "no_pauli_frame": True,
        "no_clifford_frame": True,
        "readout_convention": "internal terminal readout is folded into the final protected CNOT window as in the paper depth-two primitive",
        "reuse_convention": "corridor physical initialization is part of the first Bell/parity syndrome window; data and port standalone resets remain separate",
    }


def bell_contract(
    resource: Sequence[int], port: Sequence[int], id: str, k_name: str, r_name: str
) -> dict[str, Any]:
    """Replace resource-CNOT plus terminal measurements by an actual Bell macro.

    The macro reserves three cycles: native joint parity, individual destructive
    readout, and one idle cycle.  k/r are ready at relative time 2; a caller may
    conservatively wait to relative time 3.  Incoming resource and port patches
    must BOTH have the common standard orientation.  No extra patch is needed.
    """
    resource, port = (_coord(resource), _coord(port))
    _require(isinstance(id, str) and bool(id), "Bell macro needs an identifier")
    _require(
        all((isinstance(v, str) and bool(v) for v in (k_name, r_name))) and k_name != r_name,
        "Bell outputs must have distinct stable names",
    )
    measurement = _bell_measurement(
        resource,
        port,
        id,
        k_name,
        r_name,
        relative_cycle=0,
        readout_cycle=1,
        ready_cycle=2,
        condition=None,
    )
    return {
        "schema": "native-resource-port-Bell-contract-v1",
        "id": id,
        "resource": list(resource),
        "port": list(port),
        "k_name": k_name,
        "r_name": r_name,
        "native_orientation": ORIENTATION,
        "supply_native_orientation": ORIENTATION,
        "reserved_cycles": 3,
        "outputs_ready_relative_cycle": 2,
        "footprint": [list(resource), list(port)],
        "extra_ancilla_coords": [],
        "native_joint": measurement["joint"],
        "readouts": measurement["readouts"],
        "logical_outputs": measurement["logical_outputs"],
        "raw_outcomes": measurement["raw_outcomes"],
        "measurement": measurement,
        "stages": [
            {"relative_cycle": 0, "duration": 1, "kind": "native_joint_parity"},
            {"relative_cycle": 1, "duration": 1, "kind": "individual_destructive_readout"},
            {"relative_cycle": 2, "duration": 1, "kind": "reserved_idle"},
        ],
        "refinement_of": "CNOT(resource,port) followed by Z(port),X(resource)",
        "refinement_is_destructive_instrument_not_standalone_CNOT": True,
        "fine_branch_count": 8,
        "coarse_Bell_outcome_count": 4,
        "each_fine_instrument": "Bell_bra(k,r) * input-independent_phase / sqrt(2)",
        "no_pauli_frame": True,
        "no_clifford_frame": True,
    }


def _bits(value: int, width: int) -> tuple[int, ...]:
    return tuple((value >> bit & 1 for bit in range(width)))


def _certify_local_cnot() -> dict[str, Any]:
    witnesses = []
    for outcome in range(8):
        u, v, w = _bits(outcome, 3)
        phases = []
        for control, target in itertools.product((0, 1), repeat=2):
            ancilla = control ^ u
            flip = w ^ ancilla
            raw_target = target ^ flip
            raw_phase = v & flip
            corrected_target = raw_target ^ u ^ w
            corrected_phase = raw_phase ^ v & control
            _require(corrected_target == target ^ control, "local CNOT output mismatch")
            _require(corrected_phase == v & (u ^ w), "local CNOT input-dependent residual phase")
            phases.append(corrected_phase)
        witnesses.append([u, v, w, phases[0]])
    return {
        "passed": True,
        "checked_branches": 8,
        "basis_inputs_per_branch": 4,
        "raw_byproduct": "Z_control^v X_target^(u xor w)",
        "explicit_corrections": ["physical_Z_control^v", "physical_X_target^(u xor w)"],
        "corrected_kraus": "(-1)^(v*(u xor w)) CNOT / sqrt(8)",
        "normalization": "8 * 1/8 = 1",
        "input_independent_scalar": True,
    }


def _certify_bell_preparation() -> dict[str, Any]:
    witnesses = []
    for axis in ("horizontal", "vertical"):
        for outcome in (0, 1):
            corrected = [0] * 4
            for left in (0, 1):
                right = left ^ (outcome if axis == "vertical" else 0)
                sign = -1 if axis == "horizontal" and outcome & left else 1
                if axis == "vertical":
                    right ^= outcome
                elif outcome & right:
                    sign *= -1
                corrected[left | right << 1] += sign
            _require(corrected == [1, 0, 0, 1], "physical Bell preparation correction failed")
            witnesses.append([axis, outcome, corrected])
    return {
        "passed": True,
        "checked_native_axes": 2,
        "checked_outcomes_per_axis": 2,
        "horizontal": "initialize 00; measure XX=b; physically apply Z_right^b",
        "vertical": "initialize ++; measure ZZ=b; physically apply X_right^b",
        "output": "(|00>+|11>)/sqrt(2)",
        "deferred_frame": False,
    }


def _native_bell_row(
    axis: str, joint: int, left_readout: int, right_readout: int
) -> tuple[list[int], int, int, int]:
    """Integer row coefficients, normalized by 1/2, and k,r,global phase."""
    row = [0] * 4
    if axis == "vertical":
        k, r = (joint, left_readout ^ right_readout)
        global_phase = right_readout & k
        for left in (0, 1):
            right = left ^ k
            row[left | right << 1] = -1 if left_readout & left ^ right_readout & right else 1
    else:
        k, r = (left_readout ^ right_readout, joint)
        global_phase = r & left_readout
        row[left_readout | right_readout << 1] = 1
        row[left_readout ^ 1 | (right_readout ^ 1) << 1] = -1 if r else 1
    return (row, k, r, global_phase)


def _certify_native_bell() -> dict[str, Any]:
    witnesses = []
    for axis in ("vertical", "horizontal"):
        gram = [[0] * 4 for _ in range(4)]
        for joint, left_readout, right_readout in itertools.product((0, 1), repeat=3):
            row, k, r, global_phase = _native_bell_row(axis, joint, left_readout, right_readout)
            expected = [0] * 4
            for left in (0, 1):
                expected[left | (left ^ k) << 1] = -1 if r & left ^ global_phase else 1
            _require(row == expected, "native Bell instrument phase mismatch")
            for i, j in itertools.product(range(4), repeat=2):
                gram[i][j] += row[i] * row[j]
            witnesses.append([axis, joint, left_readout, right_readout, k, r, global_phase])
        _require(
            gram == [[4 if i == j else 0 for j in range(4)] for i in range(4)],
            "native Bell fine-grained instrument is not normalized",
        )
    return {
        "passed": True,
        "checked_native_axes": 2,
        "fine_branches_per_axis": 8,
        "basis_inputs_per_branch": 4,
        "vertical": "joint ZZ=k; X_left=p,X_right=q; r=p xor q; extra global phase=q*k",
        "horizontal": "joint XX=r; Z_left=p,Z_right=q; k=p xor q; extra global phase=r*p",
        "fine_kraus": "(-1)^global_phase Bell_bra(k,r)/sqrt(2)",
        "normalization": "sum_fine K_dagger K = I_4, exact integer Gram matrix / 4",
    }


def _teleportation_matrix(k: int, r: int, bell_pair: list[int] | None = None) -> list[list[int]]:
    """Contract Phi+(left,output) with Bell_bra(source,left), omitting 1/2.

    The optional integer Bell-pair coefficients support a direct corruption
    regression.  Coefficient index is left + 2*output throughout.
    """
    pair = [1, 0, 0, 1] if bell_pair is None else bell_pair
    _require(
        len(pair) == 4 and all((type(value) is int for value in pair)),
        "invalid Bell-pair coefficients",
    )
    bra = [0] * 4
    for source in (0, 1):
        bra[source | (source ^ k) << 1] = -1 if r & source else 1
    return [
        [
            sum((pair[left | output << 1] * bra[source | left << 1] for left in (0, 1)))
            for source in (0, 1)
        ]
        for output in (0, 1)
    ]


def _certify_teleportation_and_remote_parities() -> dict[str, Any]:
    witnesses = []
    for k, r in itertools.product((0, 1), repeat=2):
        matrix = _teleportation_matrix(k, r)
        for source in (0, 1):
            actual = [matrix[output][source] for output in (0, 1)]
            expected = [0, 0]
            expected[source ^ k] = -1 if r & source else 1
            _require(actual == expected, "Bell teleportation tensor contraction failed")
            witnesses.append([k, r, source, actual])
    for u, v, w, K, R in itertools.product((0, 1), repeat=5):
        expected_global = R & u ^ v & (u ^ K ^ w)
        for control, target in itertools.product((0, 1), repeat=2):
            ancilla = control ^ u
            after_chain = ancilla ^ K
            chain_phase = R & ancilla
            target_flip = w ^ after_chain
            raw_target = target ^ target_flip
            raw_phase = chain_phase ^ v & target_flip
            corrected_target = raw_target ^ u ^ K ^ w
            corrected_phase = raw_phase ^ (v ^ R) & control
            _require(corrected_target == target ^ control, "remote CNOT target correction failed")
            _require(corrected_phase == expected_global, "remote CNOT phase depends on input")
    return {
        "passed": True,
        "teleportation_branches": 4,
        "teleportation_inputs_per_branch": 2,
        "teleportation_kraus": "X^k Z^r / 2, up to outcome-only global phase",
        "chain_rule": "successive Bell contractions XOR their k bits and XOR their r bits; anticommutation only adds an outcome-only global phase",
        "checked_aggregate_outcomes": 32,
        "basis_inputs_per_aggregate": 4,
        "endpoint_corrections": "physical Z_control^(v xor XOR_j r_j), physical X_target^(u xor w xor XOR_j k_j)",
        "no_intermediate_frame_state": True,
    }


def _pauli_projector_integer(state: list[int], x_mask: int, z_mask: int, outcome: int) -> list[int]:
    """Apply I+(-1)^outcome P; callers keep the power-of-two normalization."""
    result = list(state)
    for index, amplitude in enumerate(state):
        phase = outcome ^ (index & z_mask).bit_count() & 1
        result[index ^ x_mask] += -amplitude if phase else amplitude
    return result


def _certify_one_pair_native_cnot() -> dict[str, Any]:
    """Direct 5-qubit contraction of every native stage, not an ideal-CNOT stub.

    Both Bell-preparation directions and both Bell-measurement directions are
    checked.  The seven raw measurement bits give 128 equal-probability branches.
    This supplements the arbitrary-chain compositional algebra above.
    """
    phases = []
    for prep_axis, bell_axis in itertools.product(("vertical", "horizontal"), repeat=2):
        initial_plus = [2, 3, 4] if prep_axis == "vertical" else [2]
        x_readouts = 2 if bell_axis == "vertical" else 0
        expected_integer_magnitude = 1 << 1 + (len(initial_plus) - 1) // 2 + x_readouts // 2
        for packed in range(128):
            u, b, joint, readout_a1, readout_a2, v, w = _bits(packed, 7)
            k = joint if bell_axis == "vertical" else readout_a1 ^ readout_a2
            r = joint if bell_axis == "horizontal" else readout_a1 ^ readout_a2
            branch_sign = None
            for control, target in itertools.product((0, 1), repeat=2):
                state = [0] * 32
                for plus_assignment in range(1 << len(initial_plus)):
                    basis = control | target << 1
                    for bit, wire in enumerate(initial_plus):
                        basis |= (plus_assignment >> bit & 1) << wire
                    state[basis] = 1
                state = _pauli_projector_integer(state, 0, 1 << 0 | 1 << 2, u)
                if prep_axis == "vertical":
                    state = _pauli_projector_integer(state, 0, 1 << 3 | 1 << 4, b)
                    if b:
                        state = [state[index ^ 1 << 4] for index in range(32)]
                else:
                    state = _pauli_projector_integer(state, 1 << 3 | 1 << 4, 0, b)
                    if b:
                        state = [
                            -amplitude if index & 1 << 4 else amplitude
                            for index, amplitude in enumerate(state)
                        ]
                state = _pauli_projector_integer(
                    state,
                    1 << 2 | 1 << 3 if bell_axis == "horizontal" else 0,
                    1 << 2 | 1 << 3 if bell_axis == "vertical" else 0,
                    joint,
                )
                state = _pauli_projector_integer(state, 1 << 4 | 1 << 1, 0, v)
                output = [0] * 4
                for index, amplitude in enumerate(state):
                    if index >> 4 & 1 != w:
                        continue
                    a1, a2 = (index >> 2 & 1, index >> 3 & 1)
                    if bell_axis == "horizontal":
                        if a1 != readout_a1 or a2 != readout_a2:
                            continue
                        readout_phase = 0
                    else:
                        readout_phase = a1 & readout_a1 ^ a2 & readout_a2
                    corrected_data = index & 3 ^ (u ^ k ^ w) << 1
                    correction_phase = (v ^ r) & (index & 1)
                    output[corrected_data] += (
                        -amplitude if readout_phase ^ correction_phase else amplitude
                    )
                destination = control | (target ^ control) << 1
                _require(
                    all((value == 0 for index, value in enumerate(output) if index != destination)),
                    "native five-qubit CNOT has a wrong output basis component",
                )
                amplitude = output[destination]
                _require(
                    abs(amplitude) == expected_integer_magnitude,
                    "native five-qubit CNOT has wrong branch normalization",
                )
                sign = amplitude < 0
                if branch_sign is None:
                    branch_sign = sign
                _require(
                    branch_sign == sign, "native five-qubit CNOT has input-dependent branch phase"
                )
            phases.append([prep_axis, bell_axis, packed, branch_sign])
    return {
        "passed": True,
        "native_axis_combinations": 4,
        "fine_branches_per_combination": 128,
        "basis_inputs_per_branch": 4,
        "simulated_qubits": 5,
        "arithmetic": "integer amplitudes with exact analytic powers-of-two normalization",
        "corrected_kraus": "input-independent sign * CNOT / sqrt(128)",
        "includes": [
            "physical Bell preparation correction",
            "native parity projectors",
            "all three destructive ancilla readouts",
            "physical endpoint Pauli corrections",
        ],
    }


def _certify_ccz_native_bells() -> dict[str, Any]:
    checks = 0
    for axes in itertools.product(("vertical", "horizontal"), repeat=3):
        for packed in range(512):
            triples = [_bits(packed >> 3 * i & 7, 3) for i in range(3)]
            legs = [_native_bell_row(axis, *raw) for axis, raw in zip(axes, triples)]
            k = tuple((leg[1] for leg in legs))
            r = tuple((leg[2] for leg in legs))
            extra = legs[0][3] ^ legs[1][3] ^ legs[2][3]
            global_phase = k[0] & k[1] & k[2] ^ sum((k[i] * r[i] for i in range(3))) & 1 ^ extra
            for assignment in range(8):
                x = _bits(assignment, 3)
                resource = tuple((x[i] ^ k[i] for i in range(3)))
                raw_phase = resource[0] & resource[1] & resource[2]
                for i in range(3):
                    coefficient = legs[i][0][resource[i] | x[i] << 1]
                    _require(
                        coefficient in (-1, 1), "native Bell branch annihilates a data basis input"
                    )
                    raw_phase ^= coefficient == -1
                correction_phase = 0
                for i in range(3):
                    j, l = [bit for bit in range(3) if bit != i]
                    correction_phase ^= k[i] & x[j] & x[l]
                    correction_phase ^= (r[i] ^ k[j] & k[l]) & x[i]
                _require(
                    raw_phase ^ correction_phase == x[0] & x[1] & x[2] ^ global_phase,
                    "native-Bell CCZ corrected branch failed",
                )
                checks += 1
    return {
        "passed": True,
        "checked_orientation_patterns": 8,
        "fine_branches_per_pattern": 512,
        "basis_inputs_per_branch": 8,
        "total_basis_branch_checks": checks,
        "corrected_kraus": "input-independent sign * CCZ / sqrt(512)",
        "normalization": "512 * 1/512 = 1 per native orientation pattern",
        "explicit_data_corrections": "CZ23^k1 CZ13^k2 CZ12^k3 and Z_i^(r_i xor k_j*k_l), all physically applied",
    }


@lru_cache(maxsize=1)
def _cached_certificates() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "passed": True,
        "native_orientation": ORIENTATION,
        "local_cnot": _certify_local_cnot(),
        "bell_preparation": _certify_bell_preparation(),
        "native_bell_measurement": _certify_native_bell(),
        "remote_cnot_parities": _certify_teleportation_and_remote_parities(),
        "native_one_pair_cnot": _certify_one_pair_native_cnot(),
        "ccz_native_bells": _certify_ccz_native_bells(),
        "no_pauli_frame": True,
        "no_clifford_frame": True,
        "proof_arithmetic": "exact integer amplitudes and binary phase polynomials; no floating point",
        "scope": "coherent protected-primitive refinement with explicit native ancillary templates; not a stabilizer noise simulation or factory proof",
        "timing_contracts": [
            "All logical Pauli corrections are materialized physical pulses; their O(1) physical execution time is excluded only from leading d-scaled cycles.",
            "Internal CNOT terminal readouts are folded into the second parity window as specified by the paper's depth-two CNOT primitive.",
            "Corridor initialization is protected by its first Bell/parity window; standalone data and port reinitializations remain separately charged.",
            "Controller-derived parities and physical correction pulses are available before the next dependent protected stage.",
            "Only vertex-disjoint complete corridors share a two-cycle CNOT batch; crossing EDP layers require a separate proven expansion or explicit VDP serialization.",
            "All incoming CCZ-state qubits are provided in the common standard native orientation, not the frozen compiler's radial-control convention.",
        ],
        "references": [
            "Beverland, Kliuchnikov, Schoute, arXiv:2110.11493, Figs.3,5,18 and Appendix A",
            "AES_estimator/main.tex, Local and Long-Range CNOT Gates; Logical Hadamard; CCZ-State Consumption",
        ],
    }


def primitive_certificates() -> dict[str, Any]:
    """Return immutable-by-copy exact branch certificates; no source-file I/O."""
    return deepcopy(_cached_certificates())


def verify_cnot_contract(contract: dict[str, Any]) -> dict[str, Any]:
    try:
        expected = cnot_contract(
            contract["control"],
            contract["target"],
            contract["path"],
            contract["id"],
            contract["condition"],
        )
        _require(
            contract == expected,
            "materialized CNOT contract differs from exact native reconstruction",
        )
        return {
            "passed": True,
            "id": contract["id"],
            "checked_ancillas": contract["ancilla_count"],
            "checked_explicit_paulis": len(contract["explicit_paulis"]),
            "errors": [],
        }
    except (KeyError, TypeError, ValueError) as exc:
        return {"passed": False, "errors": [str(exc)]}


def verify_bell_contract(contract: dict[str, Any]) -> dict[str, Any]:
    try:
        expected = bell_contract(
            contract["resource"],
            contract["port"],
            contract["id"],
            contract["k_name"],
            contract["r_name"],
        )
        _require(
            contract == expected,
            "materialized Bell contract differs from exact native reconstruction",
        )
        return {
            "passed": True,
            "id": contract["id"],
            "fine_branches": 8,
            "outputs_ready_relative_cycle": 2,
            "errors": [],
        }
    except (KeyError, TypeError, ValueError) as exc:
        return {"passed": False, "errors": [str(exc)]}
