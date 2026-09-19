"""Dependencies and patch lifetimes for concurrent native primitives."""

from collections import defaultdict

from .dependencies import derive_segment_constraints
from .native import NativeGeometry


def make_model(circuit, *, split_linear=False):
    parsed = derive_segment_constraints(circuit['operations'], circuit['logical_gadgets'])
    layout = circuit['layout']
    geom = NativeGeometry(layout['data_rows'], layout['data_cols'], circuit['logical_width'])
    nodes = []

    def node(identifier, kind, duration, wires, points, parent, start=0, **extra):
        wires = sorted(set(wires))
        points = set(map(tuple, points)) | {geom.coord(q) for q in wires}
        nodes.append(dict(id=identifier, kind=kind, duration=duration, wires=wires,
            footprint=[list(p) for p in sorted(points)], parent=parent,
            reference_start=start, **extra))

    # Each invocation prepares its local injection ports, including first use.
    node('boundary:port_reset', 'boundary_port_reset', 1, [], geom.ports,
         'boundary:port_reset')
    for op in circuit['operations']:
        sid, kind = op['id'], op['kind']
        if split_linear and kind == 'linear_cnot_layer':
            for k, child in enumerate(op['operations']):
                node(f'{sid}:child:{k}', 'linear_cnot', 2,
                     [child['control'], child['target']], child['path'], sid,
                     op['start'], cnot_count=1)
            continue
        wires = parsed['operation_supports'][sid]
        if kind == 'relabel':
            wires = list(range(circuit['logical_width']))
        node(sid, kind, op['duration'], wires, parsed['operation_footprints'][sid],
             sid, op['start'], cnot_count=
             len(op['operations']) if kind == 'linear_cnot_layer' else
             3 * len(op['gadgets']) if kind == 'native_ccz_data_to_port' else
             1 if kind == 'conditional_cnot' else 0)

    by_id = {n['id']: n for n in nodes}
    if len(by_id) != len(nodes):
        raise ValueError('duplicate event identifier')
    positions = {n['id']: i for i, n in enumerate(nodes)}
    edges, reasons = set(), defaultdict(set)

    def edge(a, b, reason):
        if a == b or positions[a] >= positions[b]:
            raise ValueError(f'reversed dependency: {a} -> {b}')
        edges.add((a, b))
        reasons[reason].add((a, b))

    last = {}
    for n in nodes:
        for q in n['wires']:
            if q in last:
                edge(last[q], n['id'], 'canonical_data_wire_order')
            last[q] = n['id']
    for consumer, outcomes in parsed['outcomes_by_consumer'].items():
        for outcome in outcomes:
            edge(parsed['producer_by_outcome'][outcome]['operation_id'], consumer,
                 'decoded_outcome')
    last_bell = {}
    for group in parsed['port_lifetimes']:
        acquire, bell, reset = group['acquire'], group['bell'], group['reset']
        edge(acquire, bell, 'injection_before_Bell')
        if reset:
            for point in by_id[reset]['footprint']:
                edge(last_bell[tuple(point)], reset, 'previous_Bell_before_reset')
            edge(reset, acquire, 'reset_before_injection')
        for point in group['ports'] + group['resources']:
            p = tuple(point)
            if p in last_bell:
                edge(last_bell[p], acquire, 'accepted_resource_reuse')
            last_bell[p] = bell
        edge('boundary:port_reset', acquire, 'boundary_ports_ready')
    return dict(nodes=nodes,
        dependencies=[list(e) for e in sorted(edges, key=lambda x: (positions[x[0]], positions[x[1]]))],
        dependency_counts={k: len(v) for k, v in reasons.items()},
        port_lifetimes=parsed['port_lifetimes'],
        producer_by_outcome=parsed['producer_by_outcome'],
        outcomes_by_consumer=parsed['outcomes_by_consumer'], split_linear=split_linear)


def reference_timeline(model):
    positions = {n['id']: i for i, n in enumerate(model['nodes'])}
    events = [dict(id=n['id'], start=n['reference_start'],
                   end=n['reference_start'] + n['duration']) for n in model['nodes']]
    events.sort(key=lambda e: (e['start'], positions[e['id']]))
    return dict(events=events, latency=max(e['end'] for e in events))
