"""Independent GF(2) and geometric verification; does not call the optimizer."""


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value, name, lower=0, upper=None):
    require(type(value) is int and value >= lower
            and (upper is None or value <= upper), f"invalid {name}")
    return value


def rows(instance):
    n = integer(instance["n"], "n", 4, 64)
    encoded = instance["target_rows_hex"]
    require(type(encoded) is list and len(encoded) == n, "matrix row count")
    require(all(type(value) is str and value.startswith("0x")
                and len(value) > 2 and all(c in "0123456789abcdefABCDEF" for c in value[2:])
                for value in encoded), "matrix rows must be hexadecimal strings")
    result = [int(value, 16) for value in encoded]
    require(all(value < (1 << n) for value in result), "matrix row exceeds n bits")
    work = result.copy()
    for column in range(n):
        pivot = next((r for r in range(column, n) if (work[r] >> column) & 1), None)
        require(pivot is not None, "target matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        for r in range(column + 1, n):
            if (work[r] >> column) & 1:
                work[r] ^= work[column]
    return result


def check_instance(instance):
    target = rows(instance)
    n, layout = instance["n"], instance["layout"]
    dr = integer(layout["data_rows"], "data rows", 2, 32)
    dc = integer(layout["data_cols"], "data columns", 2, 32)
    require(dr * dc == n, "layout must be a full data rectangle")
    require(type(layout["grid_rows"]) is int and type(layout["grid_cols"]) is int
            and layout["grid_rows"] == 2 * dr + 1 and layout["grid_cols"] == 2 * dc + 1,
            "grid shape does not match data shape")
    expected = [[2 * r + 1, 2 * c + 1] for r in range(dr) for c in range(dc)]
    require(layout["data_pos"] == expected, "data sites must be fixed row-major odd/odd sites")
    for position in layout["data_pos"]:
        require(all(type(x) is int for x in position), "non-integer data position")
    require(instance.get("output_policy", "free") == "free", "only free output is supported")
    return target


def verify(circuit, instance):
    target = check_instance(instance)
    n, layout = instance["n"], instance["layout"]
    require(type(circuit["n"]) is int and circuit["n"] == n, "wire count mismatch")
    require(circuit["mode"] == "vdp", "expected VDP circuit")
    require(circuit["layout"] == layout, "circuit changed the physical layout")
    require(check_instance(circuit) == target, "circuit changed the target matrix")
    permutation = circuit["output_permutation"]
    require(type(permutation) is list and len(permutation) == n
            and all(type(wire) is int for wire in permutation)
            and sorted(permutation) == list(range(n)), "invalid output permutation")
    terminals = {tuple(position) for position in layout["data_pos"]}
    state = [1 << wire for wire in range(n)]
    count = 0
    layers = circuit["layers"]
    require(type(layers) is list and len(layers) <= 4096, "invalid layer list")
    for index, layer in enumerate(layers):
        require(type(layer["cycles"]) is int and layer["cycles"] == 2,
                f"layer {index}: VDP layer must cost two cycles")
        gates = layer["operations"]
        require(type(gates) is list and 0 < len(gates) <= n // 2,
                f"layer {index}: invalid number of gates")
        used_wires, used_vertices = set(), set()
        for gate in gates:
            c = integer(gate["control"], "control", 0, n - 1)
            t = integer(gate["target"], "target", 0, n - 1)
            require(c != t and c not in used_wires and t not in used_wires,
                    f"layer {index}: repeated endpoint")
            used_wires.update((c, t))
            path = gate["path"]
            require(type(path) is list and 2 <= len(path) <= layout["grid_rows"] * layout["grid_cols"],
                    f"layer {index}: invalid path length")
            vertices = []
            for position in path:
                require(type(position) is list and len(position) == 2, "invalid path coordinate")
                r = integer(position[0], "path row", 0, layout["grid_rows"] - 1)
                col = integer(position[1], "path column", 0, layout["grid_cols"] - 1)
                vertices.append((r, col))
            require(path[0] == layout["data_pos"][c] and path[-1] == layout["data_pos"][t],
                    f"layer {index}: path endpoints do not match gate")
            require(len(set(vertices)) == len(vertices), f"layer {index}: self-intersecting path")
            require(not used_vertices.intersection(vertices), f"layer {index}: paths share a vertex")
            require(not terminals.intersection(vertices[1:-1]), f"layer {index}: path crosses a data site")
            require(all(abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1
                        for a, b in zip(vertices, vertices[1:])), f"layer {index}: nonadjacent path step")
            require(vertices[0][1] == vertices[1][1] and vertices[-1][0] == vertices[-2][0],
                    f"layer {index}: wrong control/target boundary orientation")
            used_vertices.update(vertices)
            state[t] ^= state[c]
            count += 1
    require(state == [target[wire] for wire in permutation], "CNOT replay does not implement P A")
    measured = {"cnots": count, "layers": len(layers), "cycles": 2 * len(layers)}
    require(circuit["stats"] == measured
            and all(type(x) is int for x in circuit["stats"].values()), "incorrect reported cost")
    return measured
