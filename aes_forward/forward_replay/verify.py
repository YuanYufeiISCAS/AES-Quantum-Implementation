"""Replay local certificates, global reservations, interfaces, and AES action."""
from .common import SCHEMA, require
from . import key_geometry as key, mixcolumns, transport, schedule, semantics
from .circuit import CONTRACT,layout,ledger
from .geometry import verify_global
from sbox_compile.verify import require_verified

EXPECTED = dict(latency=5075,reserved_patches=56918,linear_cnot=92292,
                explicit_cnot=158052,ccz_states=10080,h=49200,
                qand=5280,qand_dagger=5280,ordinary_toffoli=4800)


def verify_components(c):
    require(set(c) == {'state','key','mixcolumns','shiftrows','initial_key','final_key',
                       'adapter','key_interface','key_words','addroundkey','lookahead'},
            'missing or unsupported circuit component')
    reports = {}
    for family,width in (('state',23),('key',28)):
        reports[family] = require_verified(c[family])
        source = c[family]['source']
        require(source['logical_width'] == width, 'S-box width')
        require(source['initial_placement'] == list(range(width)), 'S-box input placement')
    require(tuple(c['state']['source']['final_placement'].index(b) for b in range(8))
            == mixcolumns.SCATTERED, 'state S-box output interface')
    require(c['key']['source']['final_placement'] == list(range(28)), 'key S-box output interface')
    require(len(c['mixcolumns']) == 2,'expected two congruence classes of MixColumns circuits')
    reports['mixcolumns'] = [mixcolumns.verify(c['mixcolumns'][k%2],k) for k in range(4)]
    for name in ('shiftrows','initial_key','final_key'):
        reports[name] = transport.verify(c[name])
    require(c['shiftrows']['source_placement_name'] == 'A' and c['shiftrows']['target_placement_name'] == 'B',
            'state ShiftRows transition')
    require(c['shiftrows']['source_interface'] == c['shiftrows']['target_interface'] == 'sbox-scattered',
            'state ShiftRows uses the wrong local slots')
    require(c['initial_key']['source_placement_name'] == 'A' and c['initial_key']['target_placement_name'] == 'B',
            'initial key transition must be the A/B involution')
    reports['adapter'] = transport.verify_adapter(c['adapter'])
    key.verify_round_local_routes(c['key_interface'],c['key_words'])
    key.verify_addroundkey_routes(c['addroundkey'])
    key.verify_round10_routes(c['lookahead'],c['key_words'],c['addroundkey'],c['final_key'])
    require(len(c['addroundkey']['groups']) == 2, 'incomplete AddRoundKey groups')
    for layer in c['key_words']['layers']:
        require(len(layer['groups']) == 2, 'incomplete word-addition groups')
    for name,count in (('load',4),('g_to_w0',2)):
        groups,paths = c['lookahead'][name+'_groups'],c['lookahead'][name+'_paths']
        require(len(groups) == len(paths) == count, 'incomplete lookahead routing')
        require(sorted(i for group in groups for i in group) == list(range(32)),
                'lookahead does not cover all 32 bits exactly once')
    return reports


def verify(program,random_trials=256):
    require(__debug__, 'verification requires Python assertions; do not use -O')
    require(program['schema'] == SCHEMA, 'unsupported forward circuit schema')
    require(program['contract'] == CONTRACT, 'primitive or supply contract changed')
    require(program['layout'] == layout(), 'global layout or external interfaces changed')
    c,nodes = program['components'],program['macros']
    modules = verify_components(c)
    timing = schedule.check(nodes,c)
    geometry = verify_global(c,nodes,program['reservations'])
    resources = ledger(c,nodes)
    require(program['resources'] == resources, 'reported gate counts do not match the circuit')
    require(geometry['latency'] == timing['latency'], 'global trace and module latency differ')
    for name,value in EXPECTED.items():
        actual = timing['latency'] if name == 'latency' else (
            geometry['reserved_patches'] if name == 'reserved_patches' else resources[name])
        require(actual == value, f'paper result mismatch for {name}: {actual} != {value}')
    functional = semantics.verify(c,nodes,random_trials)
    for recorded,field in (('linear_cnot','linear_cnot'),('qand','qand'),
                           ('qand_dagger','qand_dagger'),('toffoli','ordinary_toffoli')):
        require(functional['counts'][recorded] == resources[field], 'executed gate counts differ')
    return dict(passed=True,timing=timing,geometry=geometry,resources=resources,functional=functional,
        sboxes={name:dict(latency=modules[name]['cost']['latency'],
                         checked_inputs=modules[name]['checked_inputs']) for name in ('state','key')},
        mixcolumns=modules['mixcolumns'],
        scope='forward AES-128 under the stated protected primitives and ready CCZ inputs')
