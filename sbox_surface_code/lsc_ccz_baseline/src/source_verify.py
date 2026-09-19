"""Bit-sliced AES tests and QAND preconditions, separate from physical timing."""
from functools import lru_cache
import random
import subprocess


def require(condition, message):
    if not condition:
        raise ValueError(message)


def multiply(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        a = ((a << 1) ^ (0x11b if a & 128 else 0)) & 255
        b >>= 1
    return result


@lru_cache(maxsize=256)
def sbox(value):
    inverse = 0 if value == 0 else 1
    for _ in range(254):
        inverse = multiply(inverse, value)
    result = inverse ^ 0x63
    for shift in range(1, 5):
        result ^= ((inverse << shift) | (inverse >> (8 - shift))) & 255
    return result


def mix_column(values):
    a, b, c, d = values
    return [multiply(a, 2) ^ multiply(b, 3) ^ c ^ d,
            a ^ multiply(b, 2) ^ multiply(c, 3) ^ d,
            a ^ b ^ multiply(c, 2) ^ multiply(d, 3),
            multiply(a, 3) ^ b ^ c ^ multiply(d, 2)]


def reference_encrypt(plaintext, key):
    state, current_key = list(plaintext), list(key)
    state = [a ^ b for a, b in zip(state, current_key)]
    rcon = 1
    for round_number in range(1, 11):
        state = [sbox(x) for x in state]
        state = [state[4 * ((c + r) % 4) + r] for c in range(4) for r in range(4)]
        if round_number < 10:
            state = sum((mix_column(state[c:c+4]) for c in range(0, 16, 4)), [])
        word = [current_key[i] for i in (13, 14, 15, 12)]
        sub = [sbox(x) for x in word]
        sub[0] ^= rcon
        for i in range(4):
            current_key[i] ^= sub[i]
        for i in range(4, 16):
            current_key[i] ^= current_key[i-4]
        state = [a ^ b for a, b in zip(state, current_key)]
        rcon = multiply(rcon, 2)
    return bytes(state), bytes(current_key), word, sub


def execute(gates, state, mask):
    for position, (kind, *q) in enumerate(gates):
        if kind == 'cx':
            state[q[1]] ^= state[q[0]]
        elif kind in ('toffoli', 'qand', 'qand_dagger'):
            product = state[q[0]] & state[q[1]]
            if kind == 'qand':
                require(state[q[2]] == 0, f'dirty QAND target at gate {position}')
            if kind == 'qand_dagger':
                require(state[q[2]] == product, f'invalid QAND erasure at gate {position}')
            state[q[2]] ^= product
        elif kind == 'swap':
            state[q[0]], state[q[1]] = state[q[1]], state[q[0]]
        elif kind == 'x':
            state[q[0]] ^= mask
        else:
            raise ValueError(f'unsupported source gate: {kind}')


def verify_sources(forward, inverse, trials=256):
    """Check 259 inputs, including the standard AES known-answer example.

    These are functional tests, not exhaustive verification over 256 input bits.
    Branch identities of the emitted quantum gadgets are checked separately.
    """
    interface = forward['interface']
    require(interface == inverse['interface'], 'forward/inverse interface mismatch')
    require(forward['layout'] == inverse['layout'], 'forward/inverse layout mismatch')
    gates = [gate for _, block in forward['segments'] for gate in block]
    inverse_gates = [gate for _, block in inverse['segments'] for gate in block]
    expected_inverse = [
        [{'qand': 'qand_dagger', 'qand_dagger': 'qand'}.get(g[0], g[0]), *g[1:]]
        for g in reversed(gates)]
    require(inverse_gates == expected_inverse, 'inverse source is not the matched inverse')
    rng = random.Random(20260916)
    cases = [(bytes.fromhex('00112233445566778899aabbccddeeff'), bytes(range(16))),
             (bytes(16), bytes(16)), (bytes([255])*16, bytes([255])*16)]
    cases += [(rng.randbytes(16), rng.randbytes(16)) for _ in range(trials)]
    state = [0] * interface['width']
    for trial, (plain, key) in enumerate(cases):
        for wires, value in ((interface['initial_state'], plain), (interface['initial_key'], key)):
            for bit, wire in enumerate(wires):
                state[wire] |= ((value[bit // 8] >> (bit % 8)) & 1) << trial
    original = state.copy()
    mask = (1 << len(cases)) - 1
    execute(gates, state, mask)

    def read_bytes(wires, trial):
        bits = [(state[q] >> trial) & 1 for q in wires]
        return bytes(sum(bits[start+b] << b for b in range(8)) for start in range(0, len(bits), 8))

    require(all(state[q] == 0 for q in interface['clean_scratch']), 'unclean AES scratch')
    for trial, (plain, key) in enumerate(cases):
        ciphertext, key10, lx, ly = reference_encrypt(plain, key)
        if trial == 0:
            require(ciphertext.hex() == '69c4e0d86a7b0430d8cdb78070b4c55a', 'known-answer ciphertext')
            require(key10.hex() == '13111d7fe3944a17f307a78b4d2b30c5', 'known-answer final key')
        external = subprocess.run(['openssl', 'enc', '-aes-128-ecb', '-K', key.hex(), '-nopad', '-nosalt'],
                                  input=plain, capture_output=True, check=True, timeout=10).stdout
        require(ciphertext == external, f'reference/OpenSSL mismatch {trial}')
        require(read_bytes(interface['final_state'], trial) == external, f'ciphertext {trial}')
        require(read_bytes(interface['final_key'], trial) == key10, f'final key {trial}')
        for byte, wires in enumerate(interface['final_lookahead_modules']):
            require(read_bytes(wires[:8], trial) == bytes([lx[byte]]), 'lookahead input')
            require(read_bytes(wires[16:24], trial) == bytes([ly[byte]]), 'lookahead output')
    execute(inverse_gates, state, mask)
    require(state == original, 'inverse failed to restore all input/work wires')
    return dict(passed=True, cases=len(cases), random_cases=trials,
                OpenSSL_crosscheck=True, QAND_preconditions=True,
                clean_scratch=True, inverse_restores_all_wires=True)
