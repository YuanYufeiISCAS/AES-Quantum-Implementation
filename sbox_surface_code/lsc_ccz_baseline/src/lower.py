"""Explicit CCZ, Toffoli, QAND, and measured-erasure lowering."""
from copy import deepcopy


def add_resource_pairs(patches, targets, rows, cols):
    """Deterministic nearest free vertical pairs; no latency-guided placement.

    Protect every data H square. Each port has a free horizontal approach;
    vertical adjacency fixes Bell parity to ZZ. Old T/Y sites are not retained.
    """
    patches = deepcopy(patches)
    forbidden = set()
    for p in patches:
        r, c = p["cell"]
        dr, dc = p["h_corner"]
        forbidden.update(((r, c), (r + dr, c), (r, c + dc), (r + dr, c + dc)))
    candidates = [(r, c) for r in range(2, rows - 2, 4) for c in range(2, cols - 2, 4)]
    ports = {}
    for q in sorted(targets):
        qr, qc = patches[q]["cell"]
        ordered = sorted((p for p in candidates if cols != 382 or (p[1] < 193) == (qc < 193)),
                         key=lambda p: (abs(p[0] - qr) + abs(p[1] - qc), p))
        for r, c in ordered:
            # Keep approaches on both sides clear of other allocated pairs.
            halo = {(rr, cc) for rr in (r, r + 1) for cc in (c - 1, c, c + 1)}
            if halo & forbidden:
                continue
            p, m = len(patches), len(patches) + 1
            patches.extend([{"id": p, "cell": [r, c], "role": "injection_port", "owner": q},
                            {"id": m, "cell": [r + 1, c], "role": "CCZ_resource", "owner": q}])
            forbidden.update(halo)
            ports[q] = (p, m)
            break
        else:
            raise ValueError(f"No declared CCZ pair for wire {q}")
    return patches, ports


class Lowerer:
    def __init__(self, ports, h_policy="grouped", reset_resources=True):
        self.ports, self.h_policy = ports, h_policy
        self.reset_resources = reset_resources
        self.ops, self.modules, self.gadgets = [], [], []

    def emit(self, kind, qubits, module, **fields):
        i = len(self.ops)
        op = {"id": i, "kind": kind, "qubits": list(qubits), "module": module, **fields}
        if kind in ("mx", "mz", "bell"):
            op["result"] = i
        self.ops.append(op)
        return i

    def gate(self, gate, module):
        kind, *q = gate
        if kind == "swap":
            a, b = q
            for c, t in ((a, b), (b, a), (a, b)):
                self.emit("cx", [c, t], module)
            return
        if kind in ("cx", "x"):
            self.emit(kind, q, module)
            return
        assert kind in ("qand", "qand_dagger", "toffoli", "ccz")
        a, b, t = q
        start = len(self.ops)
        if kind == "qand_dagger":
            bit = self.emit("mx", [t], module)
            self.emit("h", [b], module)
            self.emit("cx", [a, b], module, guard={"bit": bit, "component": "m", "equals": 1})
            finish = self.emit("h", [b], module)
            self.emit("reset", [t], module)
            self.emit("barrier", q, module, after=[finish])
            resources = []
        else:
            resources = [self.ports[w][1] for w in q]
            ports = [self.ports[w][0] for w in q]
            supply = self.emit("ccz_supply", resources, module)
            if kind != "ccz":
                self.emit("h", [t], module)
            for data, port in zip(q, ports):
                self.emit("cx", [data, port], module, after=[supply])
            bits = [self.emit("bell", [p, m], module, parity="ZZ") for p, m in zip(ports, resources)]
            # CZ01^(k2), CZ02^(k1), CZ12^(k0). The grouped variant uses
            # the elementary common-target H cancellation, with no route search.
            self.emit("h", [b], module)
            self.emit("cx", [a, b], module, guard={"bit": bits[2], "component": "k", "equals": 1})
            self.emit("h", [b], module)
            self.emit("h", [t], module)
            self.emit("cx", [a, t], module, guard={"bit": bits[1], "component": "k", "equals": 1})
            if self.h_policy == "ungrouped":
                self.emit("h", [t], module)
                self.emit("h", [t], module)
            self.emit("cx", [b, t], module, guard={"bit": bits[0], "component": "k", "equals": 1})
            if self.h_policy != "grouped_pauli" or kind == "ccz":
                self.emit("h", [t], module)
            z_ops = []
            for j in range(3):
                other = [k for k in range(3) if k != j]
                correction = "x_cond" if self.h_policy == "grouped_pauli" and kind != "ccz" and j == 2 else "z_cond"
                z_ops.append(self.emit(correction, [q[j]], module,
                    controls=[{"event": bits[j], "component": "r"},
                              {"event": bits[other[0]], "component": "k"},
                              {"event": bits[other[1]], "component": "k"}]))
            # The grouped_pauli variant explicitly conjugates the intervening
            # target Z to X; otherwise both target Hadamards remain.
            if kind != "ccz" and self.h_policy != "grouped_pauli":
                self.emit("h", [t], module)
            finish = self.emit("barrier", q, module, after=z_ops)
            for wire in ports + (resources if self.reset_resources else []):
                self.emit("reset", [wire], module, after=[finish])
        self.gadgets.append({"kind": kind, "data": q, "resources": resources,
                             "start": start, "end": len(self.ops), "h_policy": self.h_policy})

    def segment(self, name, gates):
        module = len(self.modules)
        gates = list(gates)
        wires = sorted({q for gate in gates for q in gate[1:]})
        pairs = sorted({p for g in gates if g[0] in ("qand", "toffoli", "ccz")
                        for q in g[1:] for p in self.ports[q]})
        all_wires = sorted(set(wires + pairs))
        assert len(all_wires) < 255
        start = len(self.ops)
        self.emit("barrier", all_wires, module)
        for gate in gates:
            self.gate(gate, module)
        self.emit("barrier", all_wires, module)
        self.modules.append({"id": module, "name": name, "start": start, "end": len(self.ops),
                             "source_gates": gates, "wires": wires})


def specification(name, segments, patches, rows, cols, policy="module_homogeneous", h_policy="grouped", extra=None):
    # Reserve the same pair inventory for forward and inverse: QAND and its
    # matched erase exchange roles in the inverse direction.
    targets = {q for _, block in segments for g in block
               if g[0] in ("qand", "qand_dagger", "toffoli", "ccz") for q in g[1:]}
    patches, ports = add_resource_pairs(patches, targets, rows, cols)
    lowering = Lowerer(ports, h_policy)
    for label, block in segments:
        lowering.segment(label, block)
    spec = {"schema": "lsc-ccz-template-v1", "name": name, "router": "dijkstra",
            "scheduling_policy": policy, "feedback_latency": 0,
            "max_cycles": 1000000, "timeout_seconds": 1800,
            "layout": {"rows": rows, "cols": cols, "patches": patches},
            "ports": {str(k): v for k, v in ports.items()}, "modules": lowering.modules,
            "gadgets": lowering.gadgets, "operations": lowering.ops,
            "contract": {"CCZ_supply": "accepted correlated triple at specified sites before injection; no supply wait",
                "Bell": "vertical ZZ parity, X readouts, fixed padding; 3 cycles on adjacent pair",
                "corrections": "all routed conditional slots reserved; branches retain predicates",
                "data": "fixed source placement, original source SWAPs routed as three CNOTs",
                "area": "whole rectangle, including all idle and resource patches",
                "H": "3 cycles; full 2x2 deformation and restored orientation",
                "CNOT": "2 cycles; full directed path; vertex-disjoint concurrency",
                "reset": "1 cycle for both consumed resource and injection ports",
                "factories_delivery": "excluded under the same conditional supply convention",
                "optimization_imported": "none; no optimized paper paths, timings, matrices, or AES interfaces coordination"}}
    if extra:
        spec.update(extra)
    return spec
