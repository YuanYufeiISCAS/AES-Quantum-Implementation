"""Assemble the supplied circuit witnesses into one forward implementation."""
from collections import Counter
from copy import deepcopy

from .common import SCHEMA
from . import key_geometry as key
from .geometry import reservations
from .schedule import graph


CONTRACT = {
    'CNOT_cycles': 2,
    'restoring_H_cycles': 3,
    'resource_Bell_cycles': 3,
    'measurement_reset_cycles': 1,
    'physical_Pauli_and_feedback': 'explicit, available before dependent operations; zero leading d-scaled cycles',
    'CCZ_supply': 'accepted CCZ states ready at the specified resource sites before consumption',
    'excluded': ['factories','buffering','delivery','ciphertext comparison','diffusion','uncomputation'],
}


def layout():
    return dict(rows=149,cols=382,state_columns=189,separating_columns=4,key_columns=189,
                initial_byte_placement=list(key.B),ciphertext_byte_placement=list(key.C),
                key_reflection='(row,col) -> (row,381-col)',
                round_local_key_origins=[list(p) for p in key.ROUND_LOCAL_ORIGINS],
                final_key_origins=[list(p) for p in key.FINAL_ORIGINS])


def ledger(c,nodes):
    calls = Counter(n['kind'] for n in nodes)
    state_calls = 16*calls['state_sbox']
    key_calls = 4*(calls['key_sbox']+calls['final_sbox'])
    module_calls = ((state_calls,c['state']),(key_calls,c['key']))
    def total(field):
        return sum(n*m['cost'][field] for n,m in module_calls)
    roles = Counter()
    for n,module in module_calls:
        for gadget in module['source']['logical_gadgets']:
            roles[gadget['role']] += n
    route_count = lambda groups: sum(len(group['paths']) for group in groups)
    word_cnots = sum(route_count(layer['groups']) for layer in c['key_words']['layers'])
    other = dict(
        mixcolumns=calls['mixcolumns']*2*sum(sum(len(layer['operations']) for layer in m['layers'])
                                           for m in c['mixcolumns']),
        key_interfaces=(calls['key_input']+calls['key_output'])*len(key.round_local_interface_demands()[:32]),
        final_key_copy=calls['final_copy']*sum(map(len,c['lookahead']['load_paths'])),
        key_word_additions=(calls['key_words']+calls['final_words'])*word_cnots
                          +calls['final_words']*sum(map(len,c['lookahead']['g_to_w0_paths'])),
        addroundkey=calls['addroundkey']*route_count(c['addroundkey']['groups']))
    linear = total('N_CNOT_linear')+sum(other.values())
    ccz = total('N_CCZ')
    corrections = total('N_conditional_CNOT')
    return dict(state_sbox_calls=state_calls,key_sbox_calls=key_calls,
        ccz_states=ccz,linear_cnot=linear,linear_cnot_outside_sboxes=other,
        injection_cnot=3*ccz,conditional_cnot=corrections,
        explicit_cnot=linear+3*ccz+corrections,h=total('N_H'),
        qand=roles['qand'],qand_dagger=roles['qand_dagger'],ordinary_toffoli=roles['toffoli'],
        qand_x_measurements=total('N_QAND_dagger'),qand_resets=total('N_QAND_dagger'))


def assemble(components):
    nodes = graph(components)
    return dict(schema=SCHEMA,contract=deepcopy(CONTRACT),layout=layout(),components=components,
                macros=nodes,reservations=reservations(components,nodes),
                resources=ledger(components,nodes))
