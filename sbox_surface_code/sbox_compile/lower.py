"""lower."""

# Adapted from the authors' strict_sbox/nonlinear_search.py; campaign infrastructure removed.
from __future__ import annotations
from collections import defaultdict
import copy
import itertools
from . import geometry as front, native
from .model import canonical


def routes_from_artifact(seed):
    routes = defaultdict(list)
    for operation in seed["operations"]:
        if operation["kind"] == "native_ccz_data_to_port":
            group = operation["gadgets"]
            segments = {item["source_id"].rsplit(":", 1)[0] for item in group}
            front._require(len(segments) == 1, "native group crosses source segments")
            routes[next(iter(segments))].append(copy.deepcopy(group))
    return dict(routes)


def correction_variants(record, geometry):
    """Enumerate all <=48 legal local rewrites, including longer proposals."""
    source_id, physical = (record["source_id"], record["physical_indices"])
    dagger = record["role"] == "qand_dagger"
    specs = native._cz_specs(source_id, physical, dagger)
    trailing = [] if dagger else [physical[2]]
    unique = {}
    for order in itertools.permutations(range(len(specs))):
        for bits in itertools.product((0, 1), repeat=len(specs)):
            chosen, operations = ([], [])
            for index in order:
                item = specs[index]
                target = item["pair"][bits[index]]
                chosen.append({**item, "target": target})
                operations.extend(
                    front._cz_expansion(
                        tuple(item["pair"]), item["condition"], item["id"], target, geometry
                    )
                )
            operations.extend(
                (front._h(wire, f"{source_id}:trailing_H:{wire}", geometry) for wire in trailing)
            )
            operations, cancellations = front.cancel_disjoint_h_pairs(operations)
            proof = {
                "source_id": source_id,
                "rule": "commuting_diagonal_CZs_then_CZ_equals_H_CNOT_H",
                "chosen_czs": chosen,
                "cancellations": cancellations,
                "scope": "selected joint-search rewrite; no local or global optimality claim",
            }
            rebuilt = front._chosen_cz_window(
                specs, trailing, proof, geometry, optimized=True, identifier=source_id
            )
            front._require(rebuilt == operations, "local correction proof does not replay")
            key = canonical(operations)
            unique[key] = {
                "operations": operations,
                "proof": proof,
                "key": key,
                "serial_latency": sum((op["duration"] for op in operations)),
            }
    return sorted(unique.values(), key=lambda item: (item["serial_latency"], item["key"]))


def _dagger_operations(record, proof, geometry):
    source_id, physical = (record["source_id"], record["physical_indices"])
    front._require(record["reset_required"] is True, "erasure requires explicit reset")
    operations = [
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
    ]
    operations.extend(
        front._chosen_cz_window(
            native._cz_specs(source_id, physical, True),
            [],
            proof,
            geometry,
            optimized=True,
            identifier=source_id,
        )
    )
    operations.append(
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
    return operations


def _group_operations(
    geometry, records, group, proofs, used_ports=None, identifier="local:native_round:0"
):
    operations = [
        front._h(item["physical_indices"][2], f"{item['source_id']}:entry_H", geometry)
        for item in group
    ]
    operations.extend(
        native._native_injection(
            geometry, group, identifier, set() if used_ports is None else used_ports
        )
    )
    for item in group:
        sid, physical = (item["source_id"], item["physical_indices"])
        front._require(records[sid]["physical_indices"] == physical, "group changes bound support")
        operations.extend(native._z_corrections(sid, physical))
        operations.extend(
            front._chosen_cz_window(
                native._cz_specs(sid, physical, False),
                [physical[2]],
                proofs[sid],
                geometry,
                optimized=True,
                identifier=sid,
            )
        )
    return operations
