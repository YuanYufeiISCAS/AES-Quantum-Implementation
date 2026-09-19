"""Execute the compiled gate action on physical data sites, without synthesis."""
from collections import Counter, defaultdict
import random
import subprocess

from .common import require
from . import key_geometry as key
from .mixcolumns import BLOCK_SLOTS, contract


def physical_nct(source):
    gadgets = {g['source_id']:g for g in source['logical_gadgets']}
    result = []
    for op in source['operations']:
        if op['kind'] == 'linear_cnot_layer':
            result.extend(('cx',o['control'],o['target']) for o in op['operations'])
        elif op['kind'] == 'x':
            result.append(('x',op['wire']))
        elif op['kind'] == 'native_ccz_data_to_port':
            result.extend((gadgets[g['source_id']]['role'],*g['physical_indices']) for g in op['gadgets'])
        elif op['kind'] == 'qand_measure_x':
            declared = gadgets[op['id'].removesuffix(':measureX')]
            result.append(('qand_dagger',*declared['physical_indices']))
    return result


def verify(components, nodes, random_trials=256):
    require(type(random_trials) is int and random_trials >= 0, 'invalid random-trial count')
    c = components
    native = {name:physical_nct(c[name]['source']) for name in ('state','key')}
    rng = random.Random(2026091728)
    cases = [(bytes.fromhex('00112233445566778899aabbccddeeff'),bytes(range(16))),
             (bytes(16),bytes(16)),(bytes([255])*16,bytes([255])*16)]
    cases += [(rng.randbytes(16),rng.randbytes(16)) for _ in range(random_trials)]
    ones = (1 << len(cases))-1
    memory, counts = defaultdict(int), Counter()
    origins = [m['data_origin'] for m in c['shiftrows']['layout']['modules']]

    def state_point(module,slot):
        r,col = origins[module]
        return r+2*(slot//8),col+2*(slot%8)

    def key_point(label,bit,placement=key.A):
        return key.global_key(key.key_point(label,bit,placement))

    def cstar_points(module,final=False):
        places = key.FINAL_ORIGINS if final else key.ROUND_LOCAL_ORIGINS
        return [key.global_key(key.cstar(module,q,places)) for q in range(28)]

    def cx(a,b):
        memory[b] ^= memory[a]
        counts['linear_cnot'] += 1

    def move_many(pairs):
        sources,targets = [a for a,_ in pairs],[b for _,b in pairs]
        require(len(set(sources)) == len(sources) and len(set(targets)) == len(targets),
                'duplicate transfer endpoint')
        require(all(memory[b] == 0 for b in targets), 'nonclean transport target')
        values = [memory[a] for a in sources]
        for a in sources:
            memory[a] = 0
        for b,value in zip(targets,values):
            memory[b] = value
        counts['state_transfer_steps'] += len(pairs)

    def transport(payload,mirror=False):
        transform = key.global_key if mirror else tuple
        for batch in sorted({op['batch'] for op in payload['routes']}):
            move_many([(transform(op['source_data']),transform(op['output_port']))
                       for op in payload['routes'] if op['batch'] == batch])
        move_many([(transform(op['output_port']),transform(op['destination_data']))
                   for op in payload['routes']])

    def sbox_call(points,family):
        for kind,*q in native[family]:
            if kind == 'cx':
                cx(points[q[0]],points[q[1]])
            elif kind == 'x':
                memory[points[q[0]]] ^= ones
                counts['affine_x'] += 1
            else:
                a,b,t = (points[i] for i in q)
                product = memory[a] & memory[b]
                if kind == 'qand':
                    require(memory[t] == 0, 'dirty QAND target')
                elif kind == 'qand_dagger':
                    require(memory[t] == product, 'invalid measured QAND erasure')
                else:
                    require(kind == 'toffoli', 'unknown nonlinear source action')
                memory[t] ^= product
                counts[kind] += 1

    def word_chain():
        for layer in range(3):
            for demand in key.word_demands(layer):
                cx(key.global_key(demand.source),key.global_key(demand.target))

    constants = [1]
    for _ in range(9):
        constants.append(key.xtime(constants[-1]))

    def rcon(r,final=False):
        points = cstar_points(0,True)[16:24] if final else [key_point(0,b) for b in range(8)]
        for bit,p in enumerate(points):
            if constants[r-1] >> bit & 1:
                memory[p] ^= ones
                counts['affine_x'] += 1

    for module,label in enumerate(key.B):
        for bit in range(8):
            memory[state_point(module,bit)] = sum(((p[label] >> bit)&1) << i for i,(p,_) in enumerate(cases))
            memory[key_point(label,bit,key.B)] = sum(((k[label] >> bit)&1) << i for i,(_,k) in enumerate(cases))

    round_key_checks = 0
    for node in nodes:
        kind = node['kind']
        if kind == 'addroundkey':
            for demand in key.addroundkey_demands():
                cx(demand.source,demand.target)
        elif kind == 'state_sbox':
            for m in range(16):
                sbox_call([state_point(m,q) for q in range(23)],'state')
        elif kind == 'shiftrows':
            transport(c['shiftrows'])
        elif kind in ('initial_key','final_key'):
            transport(c[kind],mirror=True)
        elif kind == 'mixcolumns':
            for block in range(4):
                r,col = origins[BLOCK_SLOTS[block][0]]
                local,_,live,_ = contract(block)
                points = [(r-2+a,col-2+b) for a,b in local]
                for layer in c['mixcolumns'][block%2]['layers']:
                    for op in layer['operations']:
                        cx(points[op['control']],points[op['target']])
                require(all(memory[p] == 0 for q,p in enumerate(points) if q not in live),
                        'unclean MixColumns workspace')
        elif kind in ('key_input','key_output'):
            for demand in key.round_local_interface_demands():
                a,b = key.global_key(demand.source),key.global_key(demand.target)
                if demand.move:
                    move_many([(b,a) if kind == 'key_output' else (a,b)])
                else:
                    cx(a,b)
        elif kind in ('key_sbox','final_sbox'):
            for m in range(4):
                sbox_call(cstar_points(m,kind == 'final_sbox'),'key')
        elif kind == 'key_words':
            rcon(node['round'])
            word_chain()
            # Check every round key, not just the final ciphertext.
            for trial,(_,raw_key) in enumerate(cases):
                expected = list(raw_key)
                for constant in constants[:node['round']]:
                    expected = key.reference_next(expected,constant)
                actual = bytes(sum(((memory[key_point(label,b)] >> trial)&1) << b for b in range(8))
                               for label in range(16))
                require(actual == bytes(expected), 'incorrect intermediate round key')
                round_key_checks += 1
        elif kind == 'final_copy':
            for m,label in enumerate(key.SOURCE_LABELS):
                for b in range(8):
                    cx(key_point(label,b),cstar_points(m,True)[b])
        elif kind == 'final_words':
            rcon(10,final=True)
            for m in range(4):
                for b in range(8):
                    cx(cstar_points(m,True)[16+b],key_point(m,b))
            word_chain()
        elif kind == 'adapter':
            for m,(r,col) in enumerate(origins):
                paths = c['adapter']['schedule']['paths']
                triples = [(state_point(m,op['source_slot']),
                            (r+op['output_port'][0],col+op['output_port'][1]),
                            state_point(m,op['target_slot'])) for op in paths]
                move_many([(a,p) for a,p,b in triples])
                move_many([(p,b) for a,p,b in triples])
        else:
            raise ValueError(f'unknown forward module {kind}')

    def read_bytes(points,trial):
        return bytes(sum(((memory[p] >> trial)&1) << b for b,p in enumerate(byte)) for byte in points)

    outputs = [[state_point(key.C.index(label),b) for b in range(8)] for label in range(16)]
    final_key = [[key_point(label,b,key.C) for b in range(8)] for label in range(16)]
    xs = [cstar_points(m,True)[:8] for m in range(4)]
    ys = [cstar_points(m,True)[16:24] for m in range(4)]
    known = []
    for trial,(plaintext,raw_key) in enumerate(cases):
        reference = subprocess.run(['openssl','enc','-aes-128-ecb','-K',raw_key.hex(),'-nopad','-nosalt'],
            input=plaintext,capture_output=True,check=True,timeout=10).stdout
        actual = read_bytes(outputs,trial)
        require(actual == reference,f'OpenSSL ciphertext mismatch: case {trial}')
        current = list(raw_key)
        for constant in constants[:9]:
            current = key.reference_next(current,constant)
        require(read_bytes(final_key,trial) == bytes(key.reference_next(current,constants[9])), 'final round key')
        copied = bytes(current[label] for label in key.SOURCE_LABELS)
        g = bytes(key.SBOX[x] ^ (constants[9] if m == 0 else 0) for m,x in enumerate(copied))
        require(read_bytes(xs,trial) == copied and read_bytes(ys,trial) == g, 'retained final lookahead registers')
        if trial < 3:
            known.append(actual.hex())
    require(known[0] == '69c4e0d86a7b0430d8cdb78070b4c55a','AES-128 known-answer test')
    live = {p for byte in outputs+final_key+xs+ys for p in byte}
    require(all(value == 0 for p,value in memory.items() if p not in live), 'nonzero final scratch')
    return dict(passed=True,cases=len(cases),random_cases=random_trials,known_ciphertexts=known,
                intermediate_round_key_checks=round_key_checks,OpenSSL_crosscheck=True,
                scratch_zero=True,retained_lookahead_checked=True,counts=dict(counts))
