"""Replay the two fixed 52-wire CNOT circuits against AES-derived targets."""
from .common import require
from .key_geometry import gf_mul

SLOTS = (0,1,2,3,4,5,6,7,11,17,18,20,22)
SCATTERED = (17,18,20,0,22,2,11,4)
ORIGINS = ((2,2),(2,40),(30,2),(30,40))
BLOCK_SLOTS = ((0,1,4,5),(2,3,6,7),(8,9,12,13),(10,11,14,15))
RHO = ((0,1,2,3),(3,0,1,2),(2,3,0,1),(1,2,3,0))
PI = ((0,3,2,1),(1,0,3,2),(2,1,0,3),(3,2,1,0))


def contract(block):
    positions = [(r+2*(s//8), c+2*(s%8)) for r,c in ORIGINS for s in SLOTS]
    inputs = [13*m+SLOTS.index(s) for m in range(4) for s in SCATTERED]
    outputs = [13*m+SLOTS.index(s) for m in range(4) for s in range(8)]
    coefficients = ((2,3,1,1),(1,2,3,1),(1,1,2,3),(3,1,1,2))
    matrix = [0]*32
    for byte in range(4):
        for bit in range(8):
            for out in range(4):
                value = gf_mul(1 << bit, coefficients[PI[block][out]][RHO[block][byte]])
                for k in range(8):
                    if value >> k & 1:
                        matrix[8*out+k] |= 1 << (8*byte+bit)
    return positions, inputs, outputs, matrix


def verify(circuit, block):
    require(set(circuit) == {'n','mode','layout','output_permutation','layers','stats'},
            'unsupported MixColumns field')
    require(circuit['n'] == 52 and circuit['mode'] == 'vdp', 'MixColumns width/mode')
    require(circuit['output_permutation'] == list(range(52)), 'unpaid output permutation')
    positions, inputs, outputs, matrix = contract(block)
    require(circuit['layout'] == dict(grid_rows=39, grid_cols=59,
            data_positions=[r*59+c for r,c in positions]), 'MixColumns terminals changed')
    rows = [0]*52
    for bit, q in enumerate(inputs):
        rows[q] = 1 << bit
    count = 0
    terminals = set(positions)
    for layer in circuit['layers']:
        require(set(layer) == {'mode','surface_depth','operations'}, 'unsupported layer field')
        require(layer['mode'] == 'vdp' and layer['surface_depth'] == 2, 'CNOT layer duration')
        used = set()
        for op in layer['operations']:
            require(set(op) == {'control','target','path'}, 'unsupported CNOT field')
            a, b = op['control'], op['target']
            require(type(a) is int and type(b) is int and 0 <= a < 52 and 0 <= b < 52 and a != b,
                    'invalid CNOT endpoints')
            path = tuple(map(tuple, op['path']))
            require(len(path) >= 3 and path[0] == positions[a] and path[-1] == positions[b],
                    'CNOT route endpoints')
            require(len(set(path)) == len(path), 'repeated path vertex')
            require(all(type(r) is int and type(c) is int and 0 <= r < 39 and 0 <= c < 59
                        for r,c in path), 'route outside MixColumns grid')
            require(all(abs(r-u)+abs(c-v) == 1 for (r,c),(u,v) in zip(path,path[1:])),
                    'nonadjacent route step')
            require(path[1][1] == path[0][1] and path[-2][0] == path[-1][0],
                    'fixed boundary directions violated')
            require(not set(path[1:-1]) & terminals, 'route crosses a static data terminal')
            require(not used.intersection(path), 'non-VDP CNOT layer')
            used.update(path)
            rows[b] ^= rows[a]
            count += 1
    target = [0]*52
    for bit, q in enumerate(outputs):
        target[q] = matrix[bit]
    require(rows == target, 'AES linear map or output cleanup is incorrect')
    stats = dict(cnots=count, layers=len(circuit['layers']), surface_depth=2*len(circuit['layers']))
    require(circuit['stats'] == stats, 'MixColumns cost mismatch')
    return dict(passed=True, block=block, checked_basis_inputs=32, clean_outputs=20, **stats)
