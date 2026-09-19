"""Global reservations for fixed module calls and routed primitive batches."""
from collections import defaultdict
import heapq

from .common import require
from . import key_geometry as key
from .mixcolumns import BLOCK_SLOTS
from sbox_compile.native import NativeGeometry
from sbox_compile.parallel_model import make_model


def reservations(components, nodes):
    c = components
    records = []
    state_origins = [m['data_origin'] for m in c['shiftrows']['layout']['modules']]

    def emit(node, tag, offset, duration, points, **extra):
        points = sorted(set(map(tuple, points)))
        require(duration > 0 and offset >= 0 and offset+duration <= node['duration'],
                'primitive reservation exceeds its module interval')
        records.append(dict(id=node['id']+':'+tag, macro=node['id'],
            start=node['start']+offset, end=node['start']+offset+duration,
            footprint=[list(p) for p in points], **extra))

    def module_points(family):
        s = c[family]['source']
        g = NativeGeometry(s['layout']['data_rows'], s['layout']['data_cols'], s['logical_width'])
        points = g.vertices | g.resources
        model = make_model(s, split_linear=c[family]['split_linear'])
        require(all(set(map(tuple,n['footprint'])) <= points for n in model['nodes']),
                'S-box primitive outside module reservation')
        return points

    local_points = {name:module_points(name) for name in ('state','key')}

    def translated(points, origin, mirror=False):
        result = [(r+origin[0], col+origin[1]) for r,col in points]
        return [key.global_key(p) for p in result] if mirror else result

    def key_batches(node, groups, offset=0, tag='routes'):
        for i, group in enumerate(groups):
            paths = group['paths'] if isinstance(group,dict) else group
            emit(node, f'{tag}:{i}', offset+2*i, 2,
                 [key.global_key(p) for path in paths for p in path], kind='routed_key_batch')

    def teleport(node, payload, mirror=False):
        transform = key.global_key if mirror else tuple
        held = []
        batches = sorted({op['batch'] for op in payload['routes']})
        for batch in batches:
            selected = [op for op in payload['routes'] if op['batch'] == batch]
            points = [transform(p) for op in selected for p in op['path']]
            emit(node, f'teleport:{batch}', 2*batch, 2, points+held,
                 kind='teleportation_batch_with_retained_outputs')
            held += [transform(op['output_port']) for op in selected]
        emit(node, 'restore', 2*len(batches), 1,
             [transform(op[k]) for op in payload['routes'] for k in ('output_port','destination_data')],
             kind='parallel_restore')

    for node in nodes:
        kind = node['kind']
        if kind in ('state_sbox','key_sbox','final_sbox'):
            family = 'state' if kind == 'state_sbox' else 'key'
            origins = state_origins if family == 'state' else (
                key.FINAL_ORIGINS if kind == 'final_sbox' else key.ROUND_LOCAL_ORIGINS)
            for instance, origin in enumerate(origins):
                emit(node, f'instance:{instance}', 0, node['duration'],
                     translated(local_points[family], origin, family == 'key'),
                     kind='verified_sbox_module', component=family, instance=instance)
        elif kind == 'addroundkey':
            for i, group in enumerate(c['addroundkey']['groups']):
                emit(node, f'cnots:{i}', 2*i, 2,
                     [p for path in group['paths'] for p in path], kind='routed_cnot_batch')
        elif kind == 'shiftrows':
            teleport(node,c['shiftrows'])
        elif kind in ('initial_key','final_key'):
            teleport(node,c[kind],mirror=True)
        elif kind in ('key_input','key_output'):
            # Output moves reverse the same paths; CNOTs retain their direction.
            key_batches(node,c['key_interface']['groups'])
        elif kind == 'key_words':
            for i,layer in enumerate(c['key_words']['layers']):
                key_batches(node,layer['groups'],4*i,f'word:{i}')
        elif kind == 'final_copy':
            key_batches(node,c['lookahead']['load_paths'])
        elif kind == 'final_words':
            key_batches(node,c['lookahead']['g_to_w0_paths'],0,'G10')
            for i,layer in enumerate(c['key_words']['layers']):
                key_batches(node,layer['groups'],4+4*i,f'word:{i}')
        elif kind == 'mixcolumns':
            for block in range(4):
                r,col = state_origins[BLOCK_SLOTS[block][0]]
                for i,layer in enumerate(c['mixcolumns'][block%2]['layers']):
                    points = [(r-2+p[0],col-2+p[1]) for op in layer['operations'] for p in op['path']]
                    emit(node,f'block:{block}:layer:{i}',2*i,2,points,kind='routed_cnot_batch')
        elif kind == 'adapter':
            paths = c['adapter']['schedule']['paths']
            for instance,origin in enumerate(state_origins):
                emit(node,f'adapter:{instance}',0,2,
                     translated([p for op in paths for p in op['path']],origin),kind='teleportation_batch')
                emit(node,f'restore:{instance}',2,1,
                     translated([op[k] for op in paths for k in ('output_port','target')],origin),
                     kind='parallel_restore')
        else:
            raise ValueError(f'unknown forward module {kind}')
    return records


def verify_global(components, nodes, records):
    require(records == reservations(components,nodes), 'global reservations differ from fixed circuits')
    require(len({r['id'] for r in records}) == len(records), 'duplicate reservation')
    by_macro = defaultdict(list)
    for record in records:
        by_macro[record['macro']].append(record)
        require(all(type(r) is int and type(c) is int and 0 <= r < 149 and 0 <= c < 382
                    for r,c in record['footprint']), 'global footprint outside reserved grid')
    for node in nodes:
        require(min(r['start'] for r in by_macro[node['id']]) == node['start'], 'missing module start')
        require(max(r['end'] for r in by_macro[node['id']]) == node['end'], 'missing module end')

    active, occupied, peak = [], {}, 0
    for serial, record in sorted(enumerate(records), key=lambda p:(p[1]['start'],p[0])):
        while active and active[0][0] <= record['start']:
            _, old = heapq.heappop(active)
            for point in map(tuple, records[old]['footprint']):
                require(occupied.pop(point) == old, 'reservation release mismatch')
        for point in map(tuple,record['footprint']):
            require(point not in occupied,
                    f"global patch conflict at {point}: {record['id']} and "
                    f"{records[occupied[point]]['id'] if point in occupied else ''}")
            occupied[point] = serial
        heapq.heappush(active,(record['end'],serial))
        peak = max(peak,len(occupied))
    latency = max(r['end'] for r in records)
    return dict(passed=True, reservation_records=len(records), latency=latency,
                reserved_patches=149*382, reserved_patch_cycles=149*382*latency,
                peak_reserved_active_footprint=peak,
                sbox_footprints='whole module regions; internal event traces checked separately')
