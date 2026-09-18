"""native."""

# Adapted from the authors' strict_sbox/native_compiler.py; campaign infrastructure removed.
from __future__ import annotations
from collections import Counter
import copy
from typing import Any, Mapping, Sequence
from . import geometry as frontend_compiler
from .semantics import (
    flatten_operations,
    gadget_certificates,
    verify_logical_operations,
    verify_placed_logical_operations,
)

SCHEMA = "strict-sbox-native-materialized-v1"
MODEL = "uniform-native-ancilla-CNOT-explicit-Bell-CCZ-serial-H-v1"
_require = frontend_compiler._require


class NativeGeometry(frontend_compiler.Geometry):

    def port_axis(self, port: tuple[int, int]) -> str:
        _require(port in self.ports, f"not a native port: {port}")
        return "horizontal"


def _primitive_certificates() -> dict[str, Any]:
    from .primitives import primitive_certificates

    return primitive_certificates()


def _cnot_contract(control, target, path, identifier: str, condition=None) -> dict[str, Any]:
    from .primitives import cnot_contract

    return cnot_contract(control, target, path, identifier, condition=condition)


def _bell_contract(resource, port, identifier: str, k_name: str, r_name: str) -> dict[str, Any]:
    from .primitives import bell_contract

    return bell_contract(resource, port, identifier, k_name, r_name)


def color_path_conflicts(
    paths: Sequence[Sequence[Sequence[int]]], *, optimize: bool = True, node_budget: int = 10000
) -> dict[str, Any]:
    """Materialize a complete VDP partition; exact color search is bounded."""
    vertices = [set(map(tuple, path)) for path in paths]
    size = len(vertices)
    _require(size > 0, "empty source CNOT layer")
    adjacency = [
        {j for j in range(size) if j != i and vertices[i] & vertices[j]} for i in range(size)
    ]

    def next_vertex(colors: list[int]) -> int:
        return max(
            (i for i in range(size) if colors[i] < 0),
            key=lambda i: (
                len({colors[j] for j in adjacency[i] if colors[j] >= 0}),
                len(adjacency[i]),
                -i,
            ),
        )

    greedy = [-1] * size
    for _ in range(size):
        index = next_vertex(greedy)
        banned = {greedy[j] for j in adjacency[index] if greedy[j] >= 0}
        color = 0
        while color in banned:
            color += 1
        greedy[index] = color
    selected = greedy
    upper = 1 + max(greedy)
    lower = 2 if any(adjacency) else 1
    proved_minimum = lower == upper
    nodes = 0
    unknown = False
    if optimize:
        for limit in range(lower, upper):
            colors = [-1] * size

            def search(remaining: int) -> bool:
                nonlocal nodes
                nodes += 1
                if nodes > node_budget:
                    raise TimeoutError("local coloring node budget")
                if remaining == 0:
                    return True
                index = next_vertex(colors)
                banned = {colors[j] for j in adjacency[index] if colors[j] >= 0}
                used_colors = 1 + max(colors)
                for color in range(min(limit, used_colors + 1)):
                    if color in banned:
                        continue
                    colors[index] = color
                    if search(remaining - 1):
                        return True
                    colors[index] = -1
                return False

            try:
                feasible = search(size)
            except TimeoutError:
                unknown = True
                break
            if feasible:
                selected = colors
                proved_minimum = True
                break
        else:
            proved_minimum = True
    batches = [
        [i for i, value in enumerate(selected) if value == color]
        for color in range(1 + max(selected))
    ]
    _require(
        sorted((i for batch in batches for i in batch)) == list(range(size)),
        "coloring misses an operation",
    )
    for batch in batches:
        used: set[tuple[int, int]] = set()
        for index in batch:
            _require(not used & vertices[index], "color class is not vertex-disjoint")
            used.update(vertices[index])
    return {
        "batches": batches,
        "colors": len(batches),
        "native_cycles": 2 * len(batches),
        "local_minimum_proved": proved_minimum and (not unknown),
        "search_nodes": nodes,
        "node_budget": node_budget,
        "budget_unknown": unknown,
        "scope": "fixed-source-layer path-conflict coloring only; not global scheduling optimality",
    }


def _native_linear_map(
    frontend: Mapping[str, Any], geometry: NativeGeometry, *, optimize: bool
) -> tuple[dict[str, list[dict[str, Any]]], list[dict[str, Any]]]:
    candidates = {record["choice"]["candidate"]: record for record in frontend["linear_candidates"]}
    by_segment = {}
    proofs = []
    for segment in frontend["logical_schedule"]:
        if segment["kind"] != "linear":
            continue
        record = candidates[segment["candidate"]]
        candidate_proof = frontend_compiler.validate_candidate(
            record["candidate"], record["choice"], frontend["logical_source"]["segments"], geometry
        )
        _require(
            candidate_proof["target_rows_hex"]
            == record["matrix_geometry_replay"]["target_rows_hex"],
            "native retained linear target differs from inline source semantics",
        )
        source_offset = candidate_proof["source_coordinate_offset"]
        _require(
            source_offset == record["matrix_geometry_replay"]["source_coordinate_offset"],
            "native source coordinate convention does not replay",
        )
        _require(source_offset in ([0, 0], [1, 1]), "source coordinate offset malformed")
        operations = []
        for layer_index, source_layer in enumerate(record["candidate"]["layers"]):
            routed, original_duration = frontend_compiler._layer_geometry(
                geometry, source_layer, source_shift=source_offset == [1, 1]
            )
            coloring = color_path_conflicts([item["path"] for item in routed], optimize=optimize)
            source_id = f"{segment['segment']}:linear:{layer_index}"
            proofs.append(
                {
                    "source_linear_id": source_id,
                    "source_mode": source_layer["mode"],
                    "source_reported_cycles": original_duration,
                    **coloring,
                }
            )
            for batch_index, indices in enumerate(coloring["batches"]):
                selected = [copy.deepcopy(routed[index]) for index in indices]
                identifier = f"{source_id}:vdp:{batch_index}"
                operation = {
                    "id": identifier,
                    "kind": "linear_cnot_layer",
                    "duration": 2,
                    "category": "linear",
                    "mode": "vdp",
                    "operations": selected,
                    "source_layer": layer_index,
                    "source_linear_id": source_id,
                    "source_operation_indices": indices,
                    "batch_index": batch_index,
                    "data_qubits": sorted(
                        {q for item in selected for q in (item["control"], item["target"])}
                    ),
                }
                operation["native_cnot_contracts"] = [
                    _cnot_contract(
                        geometry.coord(item["control"]),
                        geometry.coord(item["target"]),
                        item["path"],
                        f"{identifier}:CNOT:{index}",
                    )
                    for index, item in enumerate(selected)
                ]
                operations.append(operation)
        by_segment[segment["segment"]] = operations
    return (by_segment, proofs)


def _native_injection(
    geometry: NativeGeometry,
    group: Sequence[dict[str, Any]],
    identifier: str,
    used_ports: set[tuple[int, int]],
) -> list[dict[str, Any]]:
    frontend_compiler._validate_injection_group(geometry, group)
    ports = [tuple(port) for option in group for port in option["port_coords"]]
    data = [int(wire) for option in group for wire in option["physical_indices"]]
    reused = sorted(set(ports) & used_ports)
    used_ports.update(ports)
    operations = []
    if reused:
        operations.append(
            {
                "id": f"{identifier}:port_reset",
                "kind": "port_reset",
                "duration": 1,
                "category": "port_reset",
                "data_qubits": [],
                "ports": [list(p) for p in reused],
                "prepared_state": "|0>",
                "reason": "previous Bell readout was destructive; native port reinitialization is explicit",
            }
        )
    cnot_contracts = []
    bell_contracts = []
    for option in group:
        source_id = option["source_id"]
        for index, (wire, port, resource, path) in enumerate(
            zip(
                option["physical_indices"],
                option["port_coords"],
                option["resource_coords"],
                option["data_to_port_cnot_paths"],
                strict=True,
            )
        ):
            cnot_contracts.append(
                _cnot_contract(
                    geometry.coord(int(wire)), port, path, f"{source_id}:native_data_CNOT:{index}"
                )
            )
            bell_contracts.append(
                _bell_contract(
                    resource,
                    port,
                    f"{source_id}:native_Bell:{index}",
                    f"{source_id}:k{index}",
                    f"{source_id}:r{index}",
                )
            )
    operations.extend(
        [
            {
                "id": f"{identifier}:data_to_port",
                "kind": "native_ccz_data_to_port",
                "duration": 2,
                "category": "ccz_injection",
                "data_qubits": data,
                "gadgets": copy.deepcopy(list(group)),
                "mode": "vdp",
                "native_cnot_contracts": cnot_contracts,
                "control_axis": "vertical",
                "target_axis": "horizontal",
            },
            {
                "id": f"{identifier}:Bell",
                "kind": "native_ccz_bell",
                "duration": 3,
                "category": "ccz_injection",
                "data_qubits": data,
                "gadgets": [
                    {
                        "source_id": option["source_id"],
                        "outcomes": [f"{option['source_id']}:k{j}" for j in range(3)]
                        + [f"{option['source_id']}:r{j}" for j in range(3)],
                    }
                    for option in group
                ],
                "bell_contracts": bell_contracts,
                "reserved_schedule": [
                    {"relative_start": 0, "relative_end": 1, "action": "native_joint_parity"},
                    {
                        "relative_start": 1,
                        "relative_end": 2,
                        "action": "individual_destructive_readouts",
                    },
                    {"relative_start": 2, "relative_end": 3, "action": "reserved_idle"},
                ],
                "aggregate_outputs_ready_relative_cycle": 2,
                "corrections_consume_outputs_after_reserved_cycle": 3,
            },
        ]
    )
    return operations


def _z_corrections(source_id: str, physical: Sequence[int]) -> list[dict[str, Any]]:
    result = []
    for index in range(3):
        others = [j for j in range(3) if j != index]
        variables = [f"{source_id}:r{index}", *(f"{source_id}:k{j}" for j in others)]
        result.append(
            {
                "id": f"{source_id}:Z{index}",
                "kind": "conditional_z",
                "wire": int(physical[index]),
                "data_qubits": [int(physical[index])],
                "duration": 0,
                "category": "pauli",
                "condition": {"xor": [variables[0], {"and": variables[1:]}]},
                "condition_variables": variables,
                "execution": "explicit_physical_Pauli_not_frame",
            }
        )
    return result


def _cz_specs(source_id: str, physical: Sequence[int], dagger: bool) -> list[dict[str, Any]]:
    if dagger:
        return [
            {
                "id": f"{source_id}:eraseCZ",
                "pair": list(physical[:2]),
                "condition": f"{source_id}:m",
            }
        ]
    return [
        {
            "id": f"{source_id}:CZ{index}",
            "pair": [physical[j] for j in range(3) if j != index],
            "condition": f"{source_id}:k{index}",
        }
        for index in range(3)
    ]


def _decorate_conditional_cnots(operations: list[dict[str, Any]], geometry: NativeGeometry) -> None:
    for operation in operations:
        if operation["kind"] == "conditional_cnot":
            operation["native_contract"] = _cnot_contract(
                geometry.coord(int(operation["control"])),
                geometry.coord(int(operation["target"])),
                operation["path"],
                operation["id"],
                condition=operation["condition"],
            )


def _native_cost(operations: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    kinds = Counter((op["kind"] for op in operations))
    cycles: Counter[str] = Counter()
    for operation in operations:
        cycles[operation["category"]] += int(operation["duration"])
    nonlinear = sum((len(op["gadgets"]) for op in operations if op["kind"] == "native_ccz_bell"))
    linear_cnots = sum(
        (len(op["operations"]) for op in operations if op["kind"] == "linear_cnot_layer")
    )
    native_cnots = linear_cnots + kinds["conditional_cnot"] + 3 * nonlinear
    contracts = [contract for op in operations for contract in op.get("native_cnot_contracts", [])]
    contracts.extend(
        (op["native_contract"] for op in operations if op["kind"] == "conditional_cnot")
    )
    return {
        "latency": sum((int(op["duration"]) for op in operations)),
        "D_CNOT": cycles["linear"],
        "N_CNOT_linear": linear_cnots,
        "linear_VDP_batches": kinds["linear_cnot_layer"],
        "N_CCZ": nonlinear,
        "R_CCZ": kinds["native_ccz_bell"],
        "N_H": kinds["h"],
        "H_cycles": cycles["hadamard"],
        "CCZ_injection_cycles": cycles["ccz_injection"],
        "N_conditional_CNOT": kinds["conditional_cnot"],
        "correction_CNOT_cycles": cycles["correction"],
        "N_QAND_dagger": kinds["qand_measure_x"],
        "QAND_measurement_cycles": cycles["qand_measurement"],
        "QAND_reset_cycles": cycles["qand_reset"],
        "port_reset_cycles": cycles["port_reset"],
        "N_explicit_conditional_Z": kinds["conditional_z"],
        "N_explicit_affine_X": kinds["x"],
        "N_native_ancilla_CNOTs": native_cnots,
        "N_CNOT_endpoint_Pauli_slots": 2 * native_cnots,
        "N_CNOT_explicit_Pauli_slots": sum(
            (len(contract["explicit_paulis"]) for contract in contracts)
        ),
        "N_CNOT_internal_measurement_outcomes": sum(
            (len(contract["raw_outcomes"]) for contract in contracts)
        ),
        "N_resource_Bell_joint_measurements": 3 * nonlinear,
        "N_resource_Bell_individual_readouts": 6 * nonlinear,
        "Bell_reserved_idle_cycles": kinds["native_ccz_bell"],
    }


def _without_time(operation: Mapping[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in operation.items() if key not in {"start", "end"}}


def _reconstruct_and_replay(
    artifact: Mapping[str, Any], geometry: NativeGeometry
) -> dict[str, Any]:
    records = {record["source_id"]: record for record in artifact["logical_gadgets"]}
    windows = {record["source_id"]: record for record in artifact["search"]["correction_windows"]}
    _require(
        len(records) == len(artifact["logical_gadgets"]) and set(records) == set(windows),
        "native gadget/window binding mismatch",
    )
    optimized = bool(artifact["search"]["optimized"])
    expected_linear_map, linear_proofs = _native_linear_map(artifact, geometry, optimize=optimized)
    _require(
        linear_proofs == artifact["search"]["linear_VDP_partitions"],
        "native linear partition proof does not replay",
    )
    expected_linears = {op["id"]: op for ops in expected_linear_map.values() for op in ops}
    canonical = []
    logical = []
    seen_gadgets: set[str] = set()
    seen_linears: set[str] = set()
    used_ports: set[tuple[int, int]] = set()
    for actual in artifact["operations"]:
        kind = actual["kind"]
        if kind == "linear_cnot_layer":
            identifier = actual["id"]
            _require(
                identifier in expected_linears and identifier not in seen_linears,
                "unknown or duplicate native linear layer",
            )
            _require(
                _without_time(actual) == expected_linears[identifier],
                "native linear operations differ from retained candidate and VDP partition",
            )
            seen_linears.add(identifier)
            canonical.append(expected_linears[identifier])
            logical.extend(
                (
                    {"kind": "cx", "indices": [item["control"], item["target"]]}
                    for item in actual["operations"]
                )
            )
        elif kind == "relabel":
            canonical.append(_without_time(actual))
        elif kind == "x":
            canonical.append(_without_time(actual))
            logical.append({"kind": "x", "indices": [int(actual["wire"])]})
        elif kind == "native_ccz_data_to_port":
            identifier = actual["id"].removesuffix(":data_to_port")
            group = actual["gadgets"]
            for option in group:
                source_id = option["source_id"]
                _require(
                    source_id in records and source_id not in seen_gadgets,
                    "undeclared or duplicate native CCZ",
                )
                seen_gadgets.add(source_id)
                record = records[source_id]
                physical = [int(q) for q in option["physical_indices"]]
                _require(
                    record["role"] in {"qand", "toffoli"}
                    and record["physical_indices"] == physical,
                    "native CCZ role/support mismatch",
                )
                canonical.append(
                    frontend_compiler._h(physical[2], f"{source_id}:entry_H", geometry)
                )
                logical.append(
                    {"kind": "ccx", "indices": physical, "role": record["role"], "key": source_id}
                )
            canonical.extend(_native_injection(geometry, group, identifier, used_ports))
            for option in group:
                source_id, physical = (option["source_id"], option["physical_indices"])
                canonical.extend(_z_corrections(source_id, physical))
                canonical.extend(
                    frontend_compiler._chosen_cz_window(
                        _cz_specs(source_id, physical, False),
                        [physical[2]],
                        windows[source_id],
                        geometry,
                        optimized=optimized,
                        identifier=source_id,
                    )
                )
        elif kind == "qand_measure_x":
            source_id = actual["id"].removesuffix(":measureX")
            _require(
                source_id in records and source_id not in seen_gadgets,
                "undeclared or duplicate native QAND erasure",
            )
            seen_gadgets.add(source_id)
            record = records[source_id]
            physical = record["physical_indices"]
            _require(
                record["role"] == "qand_dagger" and record["reset_required"] is True,
                "native erasure role/clean-output contract mismatch",
            )
            logical.append(
                {"kind": "ccx", "indices": physical, "role": "qand_dagger", "key": source_id}
            )
            canonical.append(
                {
                    "id": f"{source_id}:measureX",
                    "kind": "qand_measure_x",
                    "wire": physical[2],
                    "data_qubits": [physical[2]],
                    "duration": 1,
                    "category": "qand_measurement",
                    "outcome": f"{source_id}:m",
                    "precondition": "target=control0 AND control1 on the certified input subspace",
                }
            )
            canonical.extend(
                frontend_compiler._chosen_cz_window(
                    _cz_specs(source_id, physical, True),
                    [],
                    windows[source_id],
                    geometry,
                    optimized=optimized,
                    identifier=source_id,
                )
            )
            canonical.append(
                {
                    "id": f"{source_id}:reset",
                    "kind": "qand_reset",
                    "wire": physical[2],
                    "data_qubits": [physical[2]],
                    "duration": 1,
                    "category": "qand_reset",
                    "prepared_state": "|0>",
                    "reason": "source clean-scratch interface and/or subsequent reuse",
                }
            )
    _require(
        seen_gadgets == set(records) and seen_linears == set(expected_linears),
        "native physical timeline is incomplete",
    )
    if optimized:
        canonical, cancellations = frontend_compiler.cancel_disjoint_h_pairs(canonical)
    else:
        cancellations = []
    _require(
        cancellations == artifact["search"]["final_H_cancellations"],
        "native HH cancellation proof does not replay",
    )
    _decorate_conditional_cnots(canonical, geometry)
    _require(
        canonical == [_without_time(op) for op in artifact["operations"]],
        "native physical timeline differs from independently reconstructed explicit lowering",
    )
    proof = verify_placed_logical_operations(
        artifact["family"], logical, artifact["initial_placement"], artifact["final_placement"]
    )
    _require(
        proof["passed"] is True,
        f"native actual physical AES replay failed: {proof.get('failures')}",
    )
    return {
        "passed": True,
        "actual_physical_AES": proof,
        "native_gadgets_reconstructed": len(seen_gadgets),
        "native_linear_batches_reconstructed": len(seen_linears),
        "composition": "explicit-CNOT-byproduct and native-Bell identities, corrected CCZ/QAND Kraus proofs, exact timeline reconstruction, full-domain actual physical NCT replay",
    }


def verify_native(artifact: Mapping[str, Any]) -> dict[str, Any]:
    """Replay native geometry, complete microcontracts, phases and actual AES."""
    try:
        _require(
            artifact["schema"] == SCHEMA and artifact["model"] == MODEL,
            "not a native-uniform artifact",
        )
        _require(artifact["no_deferred_frames"] is True, "native artifact defers correction frames")
        geometry = NativeGeometry(
            int(artifact["row"]),
            int(artifact["layout"]["data_cols"]),
            int(artifact["logical_width"]),
        )
        _require(
            artifact["case"] == f"{artifact['family']}_{artifact['row']}row",
            "native case label does not match family/layout",
        )
        for key, value in geometry.layout.to_dict().items():
            _require(
                artifact["layout"].get(key) == value, f"native layout metadata mismatch: {key}"
            )
        _require(
            artifact["layout"]["occupied_data_coords"]
            == [list(coord) for coord in sorted(geometry.occupied)],
            "native occupied-data metadata mismatch",
        )
        _require(
            artifact["layout"]["native_orientation"]
            == "all_data_ports_resources_Z_vertical_X_horizontal",
            "nonuniform native orientation",
        )
        _require(
            artifact["layout"]["port_target_axis"] == "horizontal",
            "native port target metadata is not horizontal",
        )
        _require(
            artifact["primitive_certificates"] == _primitive_certificates(),
            "native primitive proof differs from recomputation",
        )
        _require(
            artifact["gadget_certificates"] == gadget_certificates(),
            "corrected CCZ/QAND proof differs from recomputation",
        )
        logical_records = {record["source_id"]: record for record in artifact["logical_gadgets"]}
        _require(
            int(
                artifact["logical_source"].get(
                    "total_qubits", artifact["logical_source"].get("n", 0)
                )
            )
            == geometry.width,
            "inline logical source width mismatch",
        )
        logical_source_ops = flatten_operations(artifact["logical_source"])
        source_keys = set()
        for operation in logical_source_ops:
            if operation["kind"] == "ccx":
                key = operation["key"]
                _require(
                    key in logical_records
                    and logical_records[key]["logical_indices"] == operation["indices"],
                    "inline logical source and native nonlinear declarations differ",
                )
                operation["role"] = logical_records[key]["role"]
                source_keys.add(key)
        _require(
            source_keys == set(logical_records),
            "inline logical source does not contain every native nonlinear declaration",
        )
        logical_source_proof = verify_logical_operations(artifact["family"], logical_source_ops)
        _require(
            logical_source_proof["passed"] is True,
            "inline native logical source fails independent AES replay",
        )
        categories = {
            "linear_cnot_layer": "linear",
            "h": "hadamard",
            "conditional_cnot": "correction",
            "conditional_z": "pauli",
            "x": "pauli",
            "relabel": "relabel",
            "port_reset": "port_reset",
            "qand_measure_x": "qand_measurement",
            "qand_reset": "qand_reset",
            "native_ccz_data_to_port": "ccz_injection",
            "native_ccz_bell": "ccz_injection",
        }
        now = 0
        identifiers = set()
        ready = set()
        routed_cnots = 0
        for operation in artifact["operations"]:
            kind, identifier = (operation["kind"], operation["id"])
            _require(
                all((type(operation[key]) is int for key in ("duration", "start", "end")))
                and operation["duration"] >= 0,
                "native timeline must use nonnegative integer protected cycles",
            )
            _require(identifier not in identifiers, "duplicate native operation identifier")
            identifiers.add(identifier)
            _require(
                kind in categories and operation["category"] == categories[kind],
                "native cost category/kind mismatch",
            )
            _require(
                operation["start"] == now, "native operations overlap or have an uncharged gap"
            )
            now += int(operation["duration"])
            _require(operation["end"] == now, "native operation duration mismatch")
            if kind == "linear_cnot_layer":
                support = sorted(
                    {
                        int(q)
                        for item in operation["operations"]
                        for q in (item["control"], item["target"])
                    }
                )
                _require(
                    operation["data_qubits"] == support
                    and operation["mode"] == "vdp"
                    and (operation["duration"] == 2),
                    "native linear support/VDP timing mismatch",
                )
                frontend_compiler._layer_geometry(
                    geometry,
                    {"mode": "vdp", "operations": operation["operations"], "logical_depth": 2},
                    source_shift=False,
                )
                expected = [
                    _cnot_contract(
                        geometry.coord(item["control"]),
                        geometry.coord(item["target"]),
                        item["path"],
                        f"{identifier}:CNOT:{index}",
                    )
                    for index, item in enumerate(operation["operations"])
                ]
                _require(
                    expected == operation["native_cnot_contracts"],
                    "native linear ancilla/byproduct contract mismatch",
                )
                routed_cnots += len(expected)
            elif kind == "conditional_cnot":
                _require(
                    operation["data_qubits"] == [operation["control"], operation["target"]]
                    and operation["duration"] == 2,
                    "native conditional CNOT support/time mismatch",
                )
                _require(
                    operation["condition"] in ready and operation["all_branch_reserved"] is True,
                    "native CNOT correction precedes outcome or is not reserved",
                )
                geometry.validate_path(
                    operation["path"],
                    geometry.coord(operation["control"]),
                    geometry.coord(operation["target"]),
                )
                expected = _cnot_contract(
                    geometry.coord(operation["control"]),
                    geometry.coord(operation["target"]),
                    operation["path"],
                    identifier,
                    condition=operation["condition"],
                )
                _require(
                    operation["native_contract"] == expected,
                    "native correction ancilla/byproduct contract mismatch",
                )
                routed_cnots += 1
            elif kind == "native_ccz_data_to_port":
                frontend_compiler._validate_injection_group(geometry, operation["gadgets"])
                _require(
                    operation["duration"] == 2
                    and operation["control_axis"] == "vertical"
                    and (operation["target_axis"] == "horizontal"),
                    "native data-to-port orientation/timing mismatch",
                )
                _require(
                    operation["data_qubits"]
                    == [q for item in operation["gadgets"] for q in item["physical_indices"]],
                    "native CCZ support metadata mismatch",
                )
                routed_cnots += 3 * len(operation["gadgets"])
            elif kind == "native_ccz_bell":
                _require(
                    operation["duration"] == 3
                    and operation["aggregate_outputs_ready_relative_cycle"] == 2,
                    "native Bell macro timing mismatch",
                )
                ready.update(
                    (outcome for gadget in operation["gadgets"] for outcome in gadget["outcomes"])
                )
            elif kind == "h":
                _require(
                    operation["duration"] == 3 and operation["data_qubits"] == [operation["wire"]],
                    "H is not individually serialized",
                )
                _require(
                    operation["footprint"] == geometry.h_footprint(operation["wire"]),
                    "native H allocation/return-orientation mismatch",
                )
            elif kind == "qand_measure_x":
                _require(
                    operation["duration"] == 1 and operation["data_qubits"] == [operation["wire"]],
                    "QAND native X readout mismatch",
                )
                ready.add(operation["outcome"])
            elif kind == "conditional_z":
                _require(
                    operation["duration"] == 0 and operation["data_qubits"] == [operation["wire"]],
                    "explicit native Z convention mismatch",
                )
                _require(
                    all((item in ready for item in operation["condition_variables"])),
                    "native Z correction precedes Bell outcomes",
                )
            elif kind in {"port_reset", "qand_reset"}:
                _require(
                    operation["duration"] == 1
                    and operation["data_qubits"]
                    == ([operation["wire"]] if kind == "qand_reset" else []),
                    "native protected reset mismatch",
                )
            elif kind in {"x", "relabel"}:
                _require(
                    operation["duration"] == 0
                    and operation["data_qubits"] == ([operation["wire"]] if kind == "x" else []),
                    "native X/relabel mismatch",
                )
                if kind == "x":
                    _require(
                        operation.get("execution") == "explicit_physical_Pauli_not_frame",
                        "native affine X must be an explicit physical Pauli, not a deferred frame",
                    )
        cost = _native_cost(artifact["operations"])
        _require(cost == artifact["cost"], "native costs do not replay")
        reserved = len(geometry.vertices | geometry.resources)
        _require(
            artifact["space_time"]["reserved_module_patches"] == reserved
            and artifact["space_time"]["reserved_patch_cycles"] == reserved * now,
            "native reserved space-time mismatch",
        )
        _require(
            artifact["space_time"]["data_wires"] == geometry.width
            and artifact["space_time"]["port_registers"] == len(geometry.ports)
            and (
                artifact["space_time"]["reserved_resource_storage_sites"] == len(geometry.resources)
            ),
            "native space-time detail metadata mismatch",
        )
        lowering = _reconstruct_and_replay(artifact, geometry)
        return {
            "passed": True,
            "status": "validated_under_declared_native_primitive_timing_model",
            "checked_operations": len(artifact["operations"]),
            "checked_native_ancilla_CNOTs": routed_cnots,
            "checked_individual_H": cost["N_H"],
            "latency_recomputed": now,
            "lowering": lowering,
            "inline_logical_source": logical_source_proof,
            "scope": "logical branch identities, native endpoint/ancilla allocation, VDP scheduling and complete AES; not a stabilizer-level or circuit-level noise simulation",
            "errors": [],
        }
    except (ValueError, KeyError, TypeError, IndexError, RuntimeError) as error:
        return {"passed": False, "errors": [str(error)]}
