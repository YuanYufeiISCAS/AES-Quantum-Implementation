"""Fixed key-register geometry, route checks, and AES key recurrence."""
from dataclasses import dataclass
from typing import Any, Sequence
from .transport import verify as verify_teleportation_routes
from .common import require

ROWS = 149
KEY_COLS = 189
JOINT_COLS = 382
A = (0,3,5,4,2,1,7,6,10,9,15,14,8,11,13,12)
B = (0,5,3,4,10,15,9,14,2,7,1,6,8,13,11,12)
C = (0,7,1,4,10,13,11,14,2,5,3,6,8,15,9,12)
SOURCE_LABELS = (13,14,15,12)
ROUND_LOCAL_ORIGINS = ((122,104),(94,142),(94,104),(122,142))
FINAL_ORIGINS = ((3,15),(3,135),(3,55),(3,95))
Point = tuple[int,int]


def key_point(label: int, bit: int, placement: Sequence[int] = A) -> Point:
    slot = tuple(placement).index(label)
    return 30 + 28 * (slot // 4), 30 + 38 * (slot % 4) + 2 * bit


def cstar(module: int, wire: int, origins: Sequence[Point]) -> Point:
    row, col = origins[module]
    return row + 2 * (wire // 7), col + 2 * (wire % 7)


def footprint(origin: Point) -> set[Point]:
    row, col = origin
    return {
        (r, c)
        for r in range(row - 3, row + 12)
        for c in range(col - 3, col + 18)
    }


def global_key(point: Point) -> Point:
    return point[0], JOINT_COLS - 1 - point[1]


def adjacent(left: Point, right: Point) -> bool:
    return abs(left[0] - right[0]) + abs(left[1] - right[1]) == 1


@dataclass(frozen=True)
class Demand:
    source: Point
    target: Point
    move: bool = False


def as_paths(raw: Sequence[Sequence[Sequence[int]]]) -> tuple[tuple[Point, ...], ...]:
    return tuple(tuple((int(row), int(col)) for row, col in path) for path in raw)


def verify_batch(
    name: str,
    demands: Sequence[Demand],
    paths: Sequence[Sequence[Point]],
    live: set[Point],
    rows: int = ROWS,
    cols: int = KEY_COLS,
) -> None:
    require(len(demands) == len(paths), f"{name}: demand/path count")
    used: set[Point] = set()
    for index, (demand, raw_path) in enumerate(zip(demands, paths)):
        path = tuple(raw_path)
        require(len(path) >= 2, f"{name} path {index}: too short")
        require(path[0] == demand.source and path[-1] == demand.target,
                f"{name} path {index}: endpoints")
        require(len(path) == len(set(path)), f"{name} path {index}: not simple")
        require(all(0 <= row < rows and 0 <= col < cols for row, col in path),
                f"{name} path {index}: outside grid")
        require(all(adjacent(left, right) for left, right in zip(path, path[1:])),
                f"{name} path {index}: nonadjacent edge")
        require(path[1][1] == demand.source[1],
                f"{name} path {index}: source boundary")
        if demand.move:
            require(path[-2][1] == demand.target[1],
                    f"{name} path {index}: move target boundary")
            require((len(path) - 1) % 2 == 0,
                    f"{name} path {index}: move parity")
            first = tuple(zip(path[1::2], path[2::2]))
            second = tuple(zip(path[0::2], path[1::2]))
            for layer in (first, second):
                support = tuple(point for edge in layer for point in edge)
                require(len(support) == len(set(support)),
                        f"{name} path {index}: Bell layer is not a matching")
            covered = {point for edge in first + second for point in edge}
            require(covered == set(path), f"{name} path {index}: Bell support")
        else:
            require(path[-2][0] == demand.target[0],
                    f"{name} path {index}: CNOT target boundary")
        require(not ((set(path[1:-1]) & live) - {demand.source, demand.target}),
                f"{name} path {index}: crosses a live terminal")
        require(not (set(path) & used), f"{name}: paths are not vertex-disjoint")
        used.update(path)


def partition(count: int, batches: int, pattern: int) -> list[list[int]]:
    if pattern == 0:
        return [list(range(offset, count, batches)) for offset in range(batches)]
    if pattern == 1:
        return [
            [index for index in range(count) if (index // 8) % batches == offset]
            for offset in range(batches)
        ]
    raise ValueError(f"unknown partition pattern {pattern}")


def cstar_terminals(origins: Sequence[Point]) -> set[Point]:
    return {cstar(module, wire, origins) for module in range(4) for wire in range(28)}


def verify_origins(origins: Sequence[Point]) -> None:
    terminals = {key_point(label, bit) for label in range(16) for bit in range(8)}
    areas: list[set[Point]] = []
    for origin in origins:
        area = footprint(origin)
        require(all(0 <= row < ROWS and 0 <= col < KEY_COLS for row, col in area),
                "C* footprint outside key grid")
        require(not (area & terminals), "C* footprint intersects a key terminal")
        require(all(not (area & previous) for previous in areas),
                "C* footprints intersect")
        areas.append(area)


def round_local_interface_demands() -> list[Demand]:
    demands = [
        Demand(key_point(SOURCE_LABELS[module], bit),
               cstar(module, bit, ROUND_LOCAL_ORIGINS))
        for module in range(4)
        for bit in range(8)
    ]
    demands.extend(
        Demand(key_point(module, bit),
               cstar(module, 16 + bit, ROUND_LOCAL_ORIGINS), move=True)
        for module in range(4)
        for bit in range(8)
    )
    return demands


def word_demands(layer: int) -> list[Demand]:
    return [
        Demand(key_point(4 * layer + byte, bit),
               key_point(4 * (layer + 1) + byte, bit))
        for byte in range(4)
        for bit in range(8)
    ]


def addroundkey_demands() -> list[Demand]:
    return [
        Demand(global_key(key_point(label, bit)), key_point(label, bit))
        for label in range(16)
        for bit in range(8)
    ]


def verify_round_local_routes(
    interface: dict[str, Any], word_routes: dict[str, Any]
) -> None:
    verify_origins(ROUND_LOCAL_ORIGINS)
    key_terminals = {key_point(label, bit) for label in range(16) for bit in range(8)}
    module_terminals = cstar_terminals(ROUND_LOCAL_ORIGINS)
    live = key_terminals | module_terminals

    require(interface["kind"] == "canonical-cstar-interface", "round-local route kind")
    require(tuple(interface["grid"]) == (ROWS, KEY_COLS), "round-local route grid")
    require((interface["batches"], interface["pattern"]) == (2, 0),
            "round-local route grouping")
    demands = round_local_interface_demands()
    groups = partition(64, 2, 0)
    require(len(interface["groups"]) == 2, "round-local route batch count")
    for batch, (indices, record) in enumerate(zip(groups, interface["groups"])):
        selected = [demands[index] for index in indices]
        paths = as_paths(record["paths"])
        verify_batch(f"round-local-input-{batch}", selected, paths, live)

        reverse_demands: list[Demand] = []
        reverse_paths: list[tuple[Point, ...]] = []
        for index, demand, path in zip(indices, selected, paths):
            if index < 32:
                reverse_demands.append(demand)
                reverse_paths.append(path)
            else:
                reverse_demands.append(Demand(demand.target, demand.source, move=True))
                reverse_paths.append(tuple(reversed(path)))
        verify_batch(
            f"round-local-output-{batch}", reverse_demands, reverse_paths, live
        )

    require(word_routes["kind"] == "canonical-key-update", "key-word route kind")
    require(tuple(word_routes["grid"]) == (ROWS, KEY_COLS), "key-word route grid")
    require(tuple(tuple(origin) for origin in word_routes["module_origins"])
            == ROUND_LOCAL_ORIGINS, "key-word route origins")
    require(len(word_routes["layers"]) == 3, "key-word route layer count")
    groups = partition(32, 2, 1)
    for layer, record in enumerate(word_routes["layers"]):
        require((record["demands"], record["batches"]) == (32, 2),
                f"key-word layer {layer}: metadata")
        demands = word_demands(layer)
        for batch, (indices, group) in enumerate(zip(groups, record["groups"])):
            verify_batch(
                f"key-word-{layer}-{batch}",
                [demands[index] for index in indices],
                as_paths(group["paths"]),
                key_terminals,
            )


def verify_addroundkey_routes(payload: dict[str, Any]) -> None:
    require(payload["kind"] == "ring-key-ark", "AddRoundKey route kind")
    require(payload["placement"] == "A", "AddRoundKey placement")
    require(tuple(payload["grid"]) == (ROWS, JOINT_COLS), "AddRoundKey grid")
    require((payload["batches"], payload["pattern"]) == (2, 0),
            "AddRoundKey grouping")
    state_terminals = {key_point(label, bit) for label in range(16) for bit in range(8)}
    key_terminals = {global_key(point) for point in state_terminals}
    demands = addroundkey_demands()
    groups = partition(128, 2, 0)
    for batch, (indices, record) in enumerate(zip(groups, payload["groups"])):
        require(record["demands"] == len(indices), f"AddRoundKey batch {batch}: count")
        verify_batch(
            f"AddRoundKey-{batch}",
            [demands[index] for index in indices],
            as_paths(record["paths"]),
            state_terminals | key_terminals,
            rows=ROWS,
            cols=JOINT_COLS,
        )


def verify_round10_routes(
    lookahead: dict[str, Any],
    word_routes: dict[str, Any],
    addroundkey: dict[str, Any],
    transition: dict[str, Any],
) -> None:
    verify_origins(FINAL_ORIGINS)
    require(lookahead["kind"] == "round10-lookahead-chain", "Round 10 route kind")
    require(tuple(tuple(origin) for origin in lookahead["module_origins"])
            == FINAL_ORIGINS, "Round 10 C* origins")
    require(lookahead["canonical_tail"] == "key_word_routes.json",
            "Round 10 key-word route reference")
    require((lookahead["load_batches"], lookahead["g_to_w0_batches"],
             lookahead["canonical_tail_batches"]) == (4, 2, 6),
            "Round 10 route batch counts")

    key_terminals = {key_point(label, bit) for label in range(16) for bit in range(8)}
    module_terminals = cstar_terminals(FINAL_ORIGINS)
    live = key_terminals | module_terminals
    load = [
        Demand(key_point(SOURCE_LABELS[module], bit), cstar(module, bit, FINAL_ORIGINS))
        for module in range(4)
        for bit in range(8)
    ]
    for batch, (indices, raw) in enumerate(
        zip(lookahead["load_groups"], lookahead["load_paths"])
    ):
        verify_batch(
            f"Round10-input-copy-{batch}",
            [load[index] for index in indices],
            as_paths(raw),
            live,
        )

    g_to_w0 = [
        Demand(cstar(module, 16 + bit, FINAL_ORIGINS), key_point(module, bit))
        for module in range(4)
        for bit in range(8)
    ]
    for batch, (indices, raw) in enumerate(
        zip(lookahead["g_to_w0_groups"], lookahead["g_to_w0_paths"])
    ):
        verify_batch(
            f"Round10-G-to-w0-{batch}",
            [g_to_w0[index] for index in indices],
            as_paths(raw),
            live,
        )

    groups = partition(32, 2, 1)
    for layer, record in enumerate(word_routes["layers"]):
        demands = word_demands(layer)
        for batch, (indices, group) in enumerate(zip(groups, record["groups"])):
            verify_batch(
                f"Round10-key-word-{layer}-{batch}",
                [demands[index] for index in indices],
                as_paths(group["paths"]),
                live,
            )

    areas = [footprint(origin) for origin in FINAL_ORIGINS]
    for batch, group in enumerate(addroundkey["groups"]):
        key_side = {
            (row, JOINT_COLS - 1 - col)
            for path in group["paths"]
            for row, col in path
            if col >= KEY_COLS + 4
        }
        require(all(not (area & key_side) for area in areas),
                f"Round 9 AddRoundKey batch {batch} crosses an active C* region")

    verify_teleportation_routes(transition)
    require(transition["source_placement_name"] == "A", "final transition source")
    require(transition["target_placement_name"] == "C", "final transition target")
    require((transition["teleportation_batches"], transition["surface_latency"])
            == (2, 5), "final transition latency")
    for route in transition["routes"]:
        used = {tuple(point) for point in route["path"]}
        used.add(tuple(route["output_port"]))
        used.add(tuple(route["destination_data"]))
        require(not (used & module_terminals),
                "final key transition crosses a live Round 10 C* terminal")


def gf_mul(left: int, right: int) -> int:
    result = 0
    for _ in range(8):
        if right & 1:
            result ^= left
        left = ((left << 1) ^ (0x11B if left & 0x80 else 0)) & 0xFF
        right >>= 1
    return result


def gf_pow(value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = gf_mul(result, value)
        value = gf_mul(value, value)
        exponent >>= 1
    return result


def rotl8(value: int, amount: int) -> int:
    return ((value << amount) | (value >> (8 - amount))) & 0xFF


def aes_sbox(value: int) -> int:
    inverse = gf_pow(value, 254) if value else 0
    return (
        inverse
        ^ rotl8(inverse, 1)
        ^ rotl8(inverse, 2)
        ^ rotl8(inverse, 3)
        ^ rotl8(inverse, 4)
        ^ 0x63
    )


SBOX = tuple(aes_sbox(value) for value in range(256))


def xtime(value: int) -> int:
    return ((value << 1) ^ (0x11B if value & 0x80 else 0)) & 0xFF


def reference_next(key: Sequence[int], rcon: int) -> list[int]:
    result = list(key)
    g = [SBOX[key[13]] ^ rcon, SBOX[key[14]], SBOX[key[15]], SBOX[key[12]]]
    result[:4] = [key[index] ^ g[index] for index in range(4)]
    for word in range(1, 4):
        for byte in range(4):
            result[4 * word + byte] = (
                key[4 * word + byte] ^ result[4 * (word - 1) + byte]
            )
    return result
