"""Check fixed teleportation paths, retained outputs, and final placement."""
from typing import Iterable
from .common import require

A = (0,3,5,4,2,1,7,6,10,9,15,14,8,11,13,12)
B = (0,5,3,4,10,15,9,14,2,7,1,6,8,13,11,12)

SBOX_FINAL_PLACEMENT = (
    3, 20, 5, 22, 7, 10, 15, 21, 16, 8, 14, 6,
    19, 17, 18, 9, 13, 0, 1, 11, 2, 12, 4,
)
INTERFACE_SLOTS = {
    "canonical": tuple(range(8)),
    "sbox-scattered": tuple(SBOX_FINAL_PLACEMENT.index(bit) for bit in range(8)),
}


def shift_destination(index: int) -> int:
    row = index % 4
    column = index // 4
    return row + 4 * ((column - row) % 4)


C = tuple(shift_destination(index) for index in A)
PLACEMENTS = {"A": A, "B": B, "C": C}


def adjacent(left, right) -> bool:
    return abs(left[0] - right[0]) + abs(left[1] - right[1]) == 1


def boundary_measurement(left, right) -> tuple[str, str]:
    """Return the joint parity and destructive basis for one path edge."""
    assert adjacent(left, right)
    if left[1] == right[1]:
        return "ZZ", "X"
    return "XX", "Z"


def teleportation_matchings(vertices):
    """Expand q,a1,...,am into the two protected teleportation layers."""
    edge_count = len(vertices) - 1
    assert edge_count >= 3 and edge_count % 2 == 1
    prepare = [
        (vertices[index], vertices[index + 1], "bell_prepare")
        for index in range(1, edge_count - 1, 2)
    ]
    measure = [
        (vertices[index], vertices[index + 1], "bell_measure")
        for index in range(0, edge_count - 2, 2)
    ]
    measure.append((vertices[-2], vertices[-1], "move"))
    assert {vertex for left, right, _ in prepare for vertex in (left, right)} == set(vertices[1:-1])
    assert {vertex for left, right, _ in measure for vertex in (left, right)} == set(vertices)
    for layer in (prepare, measure):
        support = []
        for left, right, _ in layer:
            boundary_measurement(left, right)
            support.extend((left, right))
        assert len(support) == len(set(support))
    return prepare, measure


def verify(data: dict) -> dict:
    schema = data["schema"]
    assert schema in {
        "aes-shiftrows-teleport-restore-v1",
        "aes-shiftrows-teleport-restore-v2",
    }
    source_name = data["source_placement_name"]
    target_name = data["target_placement_name"]
    if schema.endswith("v1"):
        source_interface = target_interface = "canonical"
        source_local_slots = target_local_slots = INTERFACE_SLOTS["canonical"]
    else:
        source_interface = data["source_interface"]
        target_interface = data["target_interface"]
        assert source_interface in INTERFACE_SLOTS
        assert target_interface in INTERFACE_SLOTS
        source_local_slots = INTERFACE_SLOTS[source_interface]
        target_local_slots = INTERFACE_SLOTS[target_interface]
        assert tuple(data["source_logical_bit_local_slots"]) == source_local_slots
        assert tuple(data["target_logical_bit_local_slots"]) == target_local_slots
    source = PLACEMENTS[source_name]
    target = PLACEMENTS[target_name]
    assert tuple(data["source_placement"]) == source
    assert tuple(data["target_placement"]) == target
    batch_count = data["teleportation_batches"]
    assert batch_count in (1, 2, 3)
    assert data["primitive"] == (
        f"{batch_count}-batch alternating-Bell-chain permutation "
        "teleportation plus one parallel restore layer"
    )
    assert data["teleportation_schedule"] == "alternating-bell-chain-output-move-v1"
    assert data["boundary_convention"] == {"horizontal": "XX", "vertical": "ZZ"}
    assert data["teleportation_cycles"] == 2 * batch_count
    assert data["restore_cycles"] == 1
    assert data["surface_latency"] == 2 * batch_count + 1

    layout = data["layout"]
    rows, cols = layout["grid_rows"], layout["grid_cols"]
    row_pitch = layout["row_pitch"]
    col_pitch = layout["col_pitch"]
    margin = layout["margin"]
    assert rows == 2 * margin + 3 * row_pitch + 9
    assert cols == 2 * margin + 3 * col_pitch + 19
    assert row_pitch >= 12 and col_pitch >= 22 and margin >= 2
    modules = layout["modules"]
    assert len(modules) == 16
    data_coord = {}
    outer_vertices = []
    physical_slots = set()
    for module in modules:
        slot = module["physical_slot"]
        assert 0 <= slot < 16 and slot not in physical_slots
        physical_slots.add(slot)
        assert (module["macro_row"], module["macro_col"]) == divmod(slot, 4)
        origin = tuple(module["data_origin"])
        assert origin == (
            margin + 2 + module["macro_row"] * row_pitch,
            margin + 2 + module["macro_col"] * col_pitch,
        )
        assert module["transform"] == "translation_only"
        local_coordinates = module["local_data_coordinates"]
        assert len(local_coordinates) == 24
        assert len({tuple(local) for local in local_coordinates}) == 24
        for local_slot, local in enumerate(local_coordinates):
            assert tuple(local) == (2 * (local_slot // 8), 2 * (local_slot % 8))
            data_coord[slot, local_slot] = (origin[0] + local[0], origin[1] + local[1])
        if schema.endswith("v1"):
            assert module["canonical_output_local_slots"] == list(range(8))
            assert module["released_clean_local_slots"] == list(range(8, 24))
        else:
            assert module["logical_output_local_slots"] == list(source_local_slots)
            assert module["output_author_wires"] == list(range(8))
            assert tuple(module["sbox_final_placement"]) == SBOX_FINAL_PLACEMENT
            assert module["released_clean_local_slots"] == sorted(
                set(range(23)) - set(source_local_slots)
            )
            assert module["vacant_local_slots"] == [23]
        extent = module["outer_extent"]
        assert extent == [origin[0] - 2, origin[0] + 8, origin[1] - 2, origin[1] + 18]
        assert 0 <= extent[0] <= extent[1] < rows
        assert 0 <= extent[2] <= extent[3] < cols
        vertices = {
            (row, col)
            for row in range(extent[0], extent[1] + 1)
            for col in range(extent[2], extent[3] + 1)
        }
        assert all(vertices.isdisjoint(prior) for prior in outer_vertices)
        outer_vertices.append(vertices)
    assert physical_slots == set(range(16))
    live_data = {
        data_coord[slot, local_slot]
        for slot in range(16)
        for local_slot in source_local_slots
    }
    assert live_data == {tuple(vertex) for vertex in layout["live_data_vertices"]}

    expected = []
    outgoing_side = {}
    for destination_slot, label in enumerate(target):
        source_slot = source.index(label)
        expected.extend(
            (label, bit, source_slot, destination_slot)
            for bit in range(8)
            if source_slot != destination_slot
            or source_local_slots[bit] != target_local_slots[bit]
        )
    routes = data["routes"]
    assert len(routes) == data["route_count"] == len(expected)
    batch_sizes = [
        sum(route["batch"] == batch for route in routes)
        for batch in range(batch_count)
    ]
    assert sum(batch_sizes) == len(expected)
    assert all(size > 0 for size in batch_sizes)
    assert sorted((item["label"], item["bit"], item["source_slot"], item["destination_slot"]) for item in routes) == sorted(expected)

    used_by_batch = [set() for _ in range(batch_count)]
    source_seen = set()
    output_seen = set()
    restore_support = set()
    batch_area = [0] * batch_count
    batch_outputs = [set() for _ in range(batch_count)]
    layer_support = [[set(), set()] for _ in range(batch_count)]
    all_outputs_by_batch = [
        {
            tuple(route["output_port"])
            for route in routes
            if route["batch"] == batch
        }
        for batch in range(batch_count)
    ]
    for route in routes:
        bit = route["bit"]
        source_slot = route["source_slot"]
        destination_slot = route["destination_slot"]
        source_data = tuple(route["source_data"])
        sender_half = tuple(route["sender_half"])
        output_port = tuple(route["output_port"])
        destination_data = tuple(route["destination_data"])
        vertices = [tuple(vertex) for vertex in route["path"]]
        batch = route["batch"]
        assert 0 <= batch < batch_count
        assert source_data == data_coord[source_slot, source_local_slots[bit]]
        assert destination_data == data_coord[destination_slot, target_local_slots[bit]]
        side = sender_half[0] - source_data[0]
        assert side in (-1, 1) and sender_half[1] == source_data[1]
        output_side = output_port[0] - destination_data[0]
        assert output_side in (-1, 1) and output_port[1] == destination_data[1]
        prior = outgoing_side.setdefault((source_slot, bit), side)
        assert prior == side
        assert vertices[0] == source_data and vertices[1] == sender_half
        assert vertices[-1] == output_port
        assert route["source_interface"] == "vertical"
        assert route["restore_interface"] == "vertical"
        assert len(vertices) == len(set(vertices))
        assert route["path_edges"] == len(vertices) - 1
        prepare, measure = teleportation_matchings(vertices)
        assert boundary_measurement(*measure[0][:2]) == ("ZZ", "X")
        for layer_index, operations in enumerate((prepare, measure)):
            support = {
                vertex
                for left, right, _ in operations
                for vertex in (left, right)
            }
            assert layer_support[batch][layer_index].isdisjoint(support)
            layer_support[batch][layer_index].update(support)
        assert all(0 <= row < rows and 0 <= col < cols for row, col in vertices)
        assert all(adjacent(left, right) for left, right in zip(vertices, vertices[1:]))
        assert all(vertex not in live_data for vertex in vertices[1:])
        assert used_by_batch[batch].isdisjoint(vertices)
        used_by_batch[batch].update(vertices)
        earlier_outputs = set().union(*all_outputs_by_batch[:batch]) if batch else set()
        assert earlier_outputs.isdisjoint(vertices)
        assert source_data not in source_seen
        assert output_port not in output_seen
        source_seen.add(source_data)
        output_seen.add(output_port)
        batch_outputs[batch].add(output_port)
        assert adjacent(output_port, destination_data)
        assert boundary_measurement(output_port, destination_data) == ("ZZ", "X")
        move = {output_port, destination_data}
        assert restore_support.isdisjoint(move)
        restore_support.update(move)
        batch_area[batch] += len(vertices) - 1
    stored_before = [sum(len(outputs) for outputs in batch_outputs[:batch]) for batch in range(batch_count)]
    peak = max(batch_area[batch] + stored_before[batch] for batch in range(batch_count))
    assert batch_area == data["batch_reserved_path_patches"]
    assert stored_before == data["stored_outputs_before_batch"]
    assert peak == data["peak_reserved_nondata_patches"]

    initial = {
        data_coord[slot, source_local_slots[bit]]: (source[slot], bit)
        for slot in range(16)
        for bit in range(8)
    }
    replay = dict(initial)
    for route in routes:
        source_key = tuple(route["source_data"])
        destination_key = tuple(route["destination_data"])
        replay[destination_key] = initial[source_key]
    expected_replay = {
        data_coord[slot, target_local_slots[bit]]: (target[slot], bit)
        for slot in range(16)
        for bit in range(8)
    }
    assert all(replay[coord] == value for coord, value in expected_replay.items())
    if source_interface == target_interface:
        fixed = {slot for slot in range(16) if source[slot] == target[slot]}
        assert all(slot not in {route["source_slot"] for route in routes} for slot in fixed)
    for route in routes:
        key = (route["source_slot"], route["bit"])
        destination_key = (route["destination_slot"], route["bit"])
        output_side = route["output_port"][0] - route["destination_data"][0]
        assert output_side == -outgoing_side[destination_key]
        if route["source_slot"] != route["destination_slot"]:
            assert outgoing_side[destination_key] == -outgoing_side[key]
    if {source_name, target_name} == {"A", "B"}:
        destination_of_source = {
            route["source_slot"]: route["destination_slot"]
            for route in routes
            if route["bit"] == 0
        }
        assert all(destination_of_source[destination] == source_slot for source_slot, destination in destination_of_source.items())
    assert data["teleportation_cycles"] + data["restore_cycles"] == data["surface_latency"]
    return dict(passed=True, routes=len(routes), batches=batch_count, latency=data["surface_latency"])


SOURCE_SLOTS = (17,18,20,0,22,2,11,4)
BOUNDS = (-2,8,-2,18)

def slot_coord(slot: int) -> tuple[int, int]:
    return 2 * (slot // 8), 2 * (slot % 8)


def adjacent(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return abs(left[0] - right[0]) + abs(left[1] - right[1]) == 1


def check_matching(
    pairs: Iterable[tuple[tuple[int, int], tuple[int, int]]], name: str
) -> None:
    used: set[tuple[int, int]] = set()
    for left, right in pairs:
        require(adjacent(left, right), f"{name}: non-adjacent operation")
        require(left not in used and right not in used, f"{name}: not a matching")
        used.add(left)
        used.add(right)


def verify_adapter(payload: dict) -> dict[str, int]:
    require(payload["schema"] == "scatter-to-canonical-adapter-v1", "adapter schema")
    module = payload["module"]
    require(tuple(module["bounds_inclusive"]) == BOUNDS, "module bounds mismatch")
    require(tuple(module["source_slot_by_bit"]) == SOURCE_SLOTS, "source slots mismatch")
    require(module["canonical_target_slots"] == list(range(8)), "target slots mismatch")
    require(
        module["clean_source_slots"] == [s for s in range(23) if s not in SOURCE_SLOTS],
        "clean-slot set mismatch",
    )
    require(module["vacant_slot"] == 23, "vacant slot mismatch")

    schedule = payload["schedule"]
    require(schedule["teleportation_batches"] == 1, "expected one adapter batch")
    require(schedule["cycles_per_teleportation_batch"] == 2, "batch cost mismatch")
    require(schedule["parallel_final_moves"] == 1, "final move cost mismatch")
    require(schedule["adapter_latency_cycles"] == 3, "adapter latency mismatch")

    sources = tuple(slot_coord(slot) for slot in SOURCE_SLOTS)
    targets = tuple(slot_coord(bit) for bit in range(8))
    blocked = set(sources) | set(targets)
    records = schedule["paths"]
    require(len(records) == 8, "adapter does not contain eight routes")
    require(sorted(record["bit"] for record in records) == list(range(8)), "bit coverage")

    route_vertices: set[tuple[int, int]] = set()
    bell_prepare: list[tuple[tuple[int, int], tuple[int, int]]] = []
    bell_measure: list[tuple[tuple[int, int], tuple[int, int]]] = []
    final_moves: list[tuple[tuple[int, int], tuple[int, int]]] = []
    edge_total = 0
    for record in records:
        bit = int(record["bit"])
        source = tuple(record["source"])
        target = tuple(record["target"])
        sender = tuple(record["sender"])
        output = tuple(record["output_port"])
        route = [tuple(point) for point in record["path"]]
        require(record["source_slot"] == SOURCE_SLOTS[bit], f"bit {bit}: source slot")
        require(record["target_slot"] == bit, f"bit {bit}: target slot")
        require(source == sources[bit] and target == targets[bit], f"bit {bit}: endpoint")
        require(route[0] == source and route[1] == sender, f"bit {bit}: sender edge")
        require(route[-1] == output, f"bit {bit}: output port")
        require(sender[1] == source[1] and abs(sender[0] - source[0]) == 1,
                f"bit {bit}: source interface is not vertical")
        require(output[1] == target[1] and abs(output[0] - target[0]) == 1,
                f"bit {bit}: restore interface is not vertical")
        require(len(route) >= 2 and (len(route) - 1) % 2 == 1,
                f"bit {bit}: path must have odd edge count")
        require(len(set(route)) == len(route), f"bit {bit}: path is not simple")
        for point in route:
            require(BOUNDS[0] <= point[0] <= BOUNDS[1], f"bit {bit}: row out of bounds")
            require(BOUNDS[2] <= point[1] <= BOUNDS[3], f"bit {bit}: col out of bounds")
        for left, right in zip(route, route[1:]):
            require(adjacent(left, right), f"bit {bit}: non-adjacent path step")
        require(not any(point in blocked for point in route[1:]),
                f"bit {bit}: path crosses a protected data patch")
        require(not (route_vertices & set(route)), f"bit {bit}: paths are not VDP")
        route_vertices.update(route)
        bell_prepare.extend((route[i], route[i + 1]) for i in range(1, len(route) - 1, 2))
        bell_measure.extend((route[i], route[i + 1]) for i in range(0, len(route), 2))
        final_moves.append((output, target))
        edge_total += len(route) - 1

    check_matching(bell_prepare, "adapter Bell preparation layer")
    check_matching(bell_measure, "adapter Bell measurement layer")
    check_matching(final_moves, "adapter final-move layer")
    require({target for _, target in final_moves} == set(targets), "canonical output coverage")
    return {"paths": 8, "path_edges": edge_total, "latency": 3}
