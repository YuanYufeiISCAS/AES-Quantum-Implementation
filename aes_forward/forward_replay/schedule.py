"""Instantiate the paper's fixed module dependency graph; no search."""
from .common import require


def durations(components):
    c = components
    return dict(state=max(e['end'] for e in c['state']['schedule']['events']),
        key=max(e['end'] for e in c['key']['schedule']['events']),
        mixcolumns=max(2*len(m['layers']) for m in c['mixcolumns']),
        shiftrows=2*len({r['batch'] for r in c['shiftrows']['routes']})+1,
        initial_key=2*len({r['batch'] for r in c['initial_key']['routes']})+1,
        final_key=2*len({r['batch'] for r in c['final_key']['routes']})+1,
        adapter=2*c['adapter']['schedule']['teleportation_batches']+1,
        key_interface=2*len(c['key_interface']['groups']),
        key_words=2*sum(len(layer['groups']) for layer in c['key_words']['layers']),
        final_copy=2*len(c['lookahead']['load_groups']),
        final_add=2*len(c['lookahead']['g_to_w0_groups']),
        addroundkey=2*len(c['addroundkey']['groups']))


def graph(components):
    d = durations(components)
    nodes, by_name = [], {}

    def add(name, kind, duration, predecessors, round_number):
        start = max((by_name[p]['end'] for p in predecessors), default=0)
        node = dict(id=name, kind=kind, round=round_number, duration=duration,
                    predecessors=predecessors, start=start, end=start+duration)
        nodes.append(node)
        by_name[name] = node
        return name

    prior = add('ARK0', 'addroundkey', d['addroundkey'], [], 0)
    for r in range(1,10):
        state = add(f'SB{r}', 'state_sbox', d['state'], [prior], r)
        if r != 1:
            state = add(f'SR{r}', 'shiftrows', d['shiftrows'], [state], r)
        state = add(f'MC{r}', 'mixcolumns', d['mixcolumns'], [state], r)
        key = prior
        if r == 1:
            key = add('key_B_to_A', 'initial_key', d['initial_key'], [key], r)
        key = add(f'key_input{r}', 'key_input', d['key_interface'], [key], r)
        key = add(f'Cstar{r}', 'key_sbox', d['key'], [key], r)
        key = add(f'key_output{r}', 'key_output', d['key_interface'], [key], r)
        key = add(f'key_words{r}', 'key_words', d['key_words'], [key], r)
        prior = add(f'ARK{r}', 'addroundkey', d['addroundkey'], [state,key], r)
    final = add('key10_copy', 'final_copy', d['final_copy'], [key], 10)
    final = add('Cstar10', 'final_sbox', d['key'], [final], 10)
    final = add('key10_words', 'final_words', d['final_add']+d['key_words'], [final,prior], 10)
    final = add('key_A_to_C', 'final_key', d['final_key'], [final], 10)
    state = add('SB10', 'state_sbox', d['state'], [prior], 10)
    state = add('final_state_adapter', 'adapter', d['adapter'], [state], 10)
    add('ARK10', 'addroundkey', d['addroundkey'], [state,final], 10)
    return nodes


def check(nodes, components):
    # Regenerate dependencies from the specified AES architecture, not from
    # success flags or the expected final latency.
    expected = graph(components)
    require(nodes == expected, 'forward module dependencies or timestamps differ')
    by_name = {n['id']:n for n in nodes}
    for node in nodes:
        require(node['end'] == node['start']+node['duration'], 'macro duration mismatch')
        for predecessor in node['predecessors']:
            require(by_name[predecessor]['end'] <= node['start'], 'unmet module dependency')
    ends = [by_name[f'ARK{r}']['end'] for r in range(11)]
    return dict(passed=True, macro_count=len(nodes), latency=ends[-1],
                round_ends=ends, round_durations=[ends[0]]+[b-a for a,b in zip(ends,ends[1:])],
                primitive_durations=durations(components))
