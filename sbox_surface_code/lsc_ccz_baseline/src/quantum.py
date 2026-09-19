"""Independent branch maps and source recovery. No scheduler/router imports."""
from itertools import product
from math import sqrt
import json


def h(v, q):
    u = v.copy()
    mask = 1 << q
    for i in range(len(v)):
        if not i & mask:
            j = i | mask
            u[i], u[j] = (v[i] + v[j]) / sqrt(2), (v[i] - v[j]) / sqrt(2)
    return u


def cx(v, c, t):
    u = [0j] * len(v)
    for i, value in enumerate(v):
        u[i ^ ((1 << t) if (i >> c) & 1 else 0)] = value
    return u


def read_x(v, q, outcome):
    # Destructive X bra followed by a canonical |0> output placeholder.
    # Physical reset remains a separate charged operation in the trace.
    u = [0j] * len(v)
    for i in range(len(v)):
        if not (i >> q) & 1:
            u[i] = (v[i] + (-1) ** outcome * v[i | (1 << q)]) / sqrt(2)
    return u


def bell_bra(v, p, m, k, r):
    u = [0j] * len(v)
    for i in range(len(v)):
        if (i >> p) & 1 or (i >> m) & 1:
            continue
        for y in (0, 1):
            j = i | ((k ^ y) << p) | (y << m)
            u[i] += (-1) ** (r * y) * v[j] / sqrt(2)
    return u


def verify_bell_primitive():
    """ZZ parity + two X reads equals a Bell bra, including discarded bit."""
    error = 0.0
    for k, sa, sm in product((0, 1), repeat=3):
        r = sa ^ sm
        for a, m in product((0, 1), repeat=2):
            actual = (0 if a ^ m != k else (-1) ** (sa * a + sm * m) / 2)
            expected = (0 if a ^ m != k else (-1) ** (r * m + sa * k) / 2)
            error = max(error, abs(actual - expected))
    assert error < 1e-12
    return {"passed": True, "fine_grained_outcomes": 8, "basis_columns": 4,
            "coarse_outcomes": "k=ZZ outcome, r=X-readout parity",
            "unused_readout_bit_changes_only_global_phase": True, "max_error": error}


def simulate(block, gadget, basis, branch):
    data = gadget["data"]
    ids = {q: j for j, q in enumerate(data)}
    for op in block:
        for q in op["qubits"]:
            if q not in ids:
                ids[q] = len(ids)
    v = [0j] * (1 << len(ids))
    if gadget["kind"] == "qand_dagger":
        v[basis] = 1
    else:
        resource = [ids[q] for q in gadget["resources"]]
        for y in range(8):
            index = basis | sum(((y >> j) & 1) << q for j, q in enumerate(resource))
            v[index] = (-1 if y == 7 else 1) / sqrt(8)
    results = {}
    cursor = 0
    for op in block:
        kind = op["kind"]
        q = [ids[w] for w in op["qubits"]]
        if "guard" in op:
            g = op["guard"]
            assert g["bit"] in results, "correction before measurement"
            if results[g["bit"]][g["component"]] != g["equals"]:
                continue
        if kind == "h":
            v = h(v, q[0])
        elif kind == "cx":
            v = cx(v, *q)
        elif kind == "bell":
            k, r = branch[2 * cursor:2 * cursor + 2]
            cursor += 1
            results[op["id"]] = {"k": k, "r": r}
            v = bell_bra(v, q[0], q[1], k, r)
        elif kind == "mx":
            m = branch[0]
            results[op["id"]] = {"m": m}
            v = read_x(v, q[0], m)
        elif kind in ("z_cond", "x_cond"):
            bits = [results[b["event"]][b["component"]] for b in op["controls"]]
            if bits[0] ^ (bits[1] & bits[2]):
                if kind == "z_cond":
                    v = [-z if (i >> q[0]) & 1 else z for i, z in enumerate(v)]
                else:
                    v = [v[i ^ (1 << q[0])] for i in range(len(v))]
        elif kind == "x":
            v = [v[i ^ (1 << q[0])] for i in range(len(v))]
        elif kind == "z":
            v = [-z if (i >> q[0]) & 1 else z for i, z in enumerate(v)]
        else:
            assert kind in ("reset", "ccz_supply", "barrier"), kind
    return v


def verify_gadget(block, gadget):
    kind = gadget["kind"]
    if kind == "qand_dagger":
        bases, branches, magnitude = [0, 1, 2, 7], [(0,), (1,)], 1 / sqrt(2)
    else:
        bases = list(range(4)) if kind == "qand" else list(range(8))
        branches, magnitude = list(product((0, 1), repeat=6)), 1 / 8
    error = 0.0
    for branch in branches:
        scalar = None
        for basis in bases:
            actual = simulate(block, gadget, basis, branch)
            if kind == "ccz":
                target, sign = basis, (-1 if basis == 7 else 1)
            else:
                target, sign = basis ^ (4 if basis & 3 == 3 else 0), 1
            if scalar is None:
                scalar = actual[target] / sign
                assert abs(abs(scalar) - magnitude) < 1e-12, "wrong branch probability"
            error = max(error, max(abs(value - (sign * scalar if i == target else 0))
                                   for i, value in enumerate(actual)))
    assert error < 1e-12, f"wrong {kind} branch map: {error}"
    return {"passed": True, "kind": kind, "branches": len(branches), "basis_columns": len(bases),
            "max_Kraus_column_error": error, "arbitrary_reference_preserved_by_linear_map": True,
            "surviving_data_map_verified": True,
            "measured_systems_discarded_in_Kraus_representation": True,
            "physical_reset_lifecycle_checked_by_trace_audit": True}


def normalized(block, gadget):
    ids = {q: i for i, q in enumerate(gadget["data"])}
    for op in block:
        for q in op["qubits"]:
            if q not in ids:
                ids[q] = len(ids)
    rows = []
    for op in block:
        row = {k: v for k, v in op.items() if k not in ("id", "module", "result", "after", "qubits", "guard", "controls")}
        row["qubits"] = [ids[q] for q in op["qubits"]]
        if "after" in op:
            row["after"] = [x - gadget["start"] for x in op["after"]]
        if "guard" in op:
            row["guard"] = {**op["guard"], "bit": op["guard"]["bit"] - gadget["start"]}
        if "controls" in op:
            row["controls"] = [{**b, "event": b["event"] - gadget["start"]} for b in op["controls"]]
        rows.append(row)
    return json.dumps(rows, sort_keys=True)


def verify_templates(spec):
    templates, checked = {}, []
    for g in spec["gadgets"]:
        block = spec["operations"][g["start"]:g["end"]]
        key = (g["kind"], g["h_policy"])
        n = normalized(block, g)
        if key not in templates:
            checked.append(verify_gadget(block, g))
            templates[key] = n
        assert templates[key] == n, "inconsistent emitted gadget occurrence"
    starts = {g["start"]: g for g in spec["gadgets"]}
    recovered, index = [], 0
    while index < len(spec["operations"]):
        if index in starts:
            g = starts[index]
            recovered.append((g["kind"], *g["data"]))
            index = g["end"]
        else:
            op = spec["operations"][index]
            if op["kind"] != "barrier":
                assert op["kind"] in ("cx", "x"), "unaccounted operation outside templates"
                assert "guard" not in op
                recovered.append((op["kind"], *op["qubits"]))
            index += 1
    expected = []
    for module in spec["modules"]:
        for kind, *q in module["source_gates"]:
            if kind == "swap":
                a, b = q
                expected.extend((("cx", a, b), ("cx", b, a), ("cx", a, b)))
            else:
                expected.append((kind, *q))
    assert recovered == expected, "serialized program differs from source gates"
    return {"passed": True, "gadget_occurrences": len(spec["gadgets"]),
            "templates": checked, "serialized_program_recovered": True,
            "Bell_primitive": verify_bell_primitive()}
