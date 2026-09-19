"""Independent timeline replay. Does not import the scheduler.

Every reordered quantum instrument retains the order on each physical data
wire and the producer-before-consumer order for measurement outcomes. Thus
reordering is by disjoint-support exchanges, including on entangled inputs.
Complete physical footprints and continuous port lifetimes are checked anew.
"""
from collections import Counter, defaultdict
from .parallel_model import make_model


def audit(artifact, timeline, *, split_linear=False):
    if not __debug__:
        raise ValueError('timeline verification requires Python assertions; do not use -O')
    model=make_model(artifact,split_linear=split_linear)
    nodes={n['id']:n for n in model['nodes']}
    events=timeline['events']
    by_id={e['id']:e for e in events}
    assert len(events)==len(nodes)==len(by_id) and set(by_id)==set(nodes), 'event cover'
    order={e['id']:i for i,e in enumerate(events)}
    occupancy={}; wire_occupancy={}; kinds_by_tick=defaultdict(set)
    for e in events:
        sid=e['id']; n=nodes[sid]; start,end=e['start'],e['end']
        assert type(start) is int and type(end) is int and start>=0 and end==start+n['duration'], 'duration'
        for tick in range(start,end):
            kinds_by_tick[tick].add(n['kind'])
            for point in map(tuple,n['footprint']):
                assert (tick,point) not in occupancy, 'overlapping physical footprint'
                occupancy[tick,point]=sid
            for q in n['wires']:
                assert (tick,q) not in wire_occupancy, 'overlapping data support'
                wire_occupancy[tick,q]=sid
        if start==end:
            for point in map(tuple,n['footprint']):
                for other in events:
                    if other['start']<start<other['end']:
                        assert point not in map(tuple,nodes[other['id']]['footprint']), 'Pauli inside occupied interval'
    # Independently check the canonical projection onto each data wire.
    previous={}
    for n in model['nodes']:
        e=by_id[n['id']]
        for q in n['wires']:
            if q in previous:
                p=previous[q]
                assert p['end']<=e['start'], 'data wire order'
                if p['start']==p['end']==e['start']:
                    assert order[p['id']]<order[e['id']], 'zero-cycle data order'
            previous[q]=e
    for a,b in model['dependencies']:
        assert by_id[a]['end']<=by_id[b]['start'], 'dependency'
        if by_id[a]['start']==by_id[a]['end']==by_id[b]['start']:
            assert order[a]<order[b], 'zero-cycle dependency'
    for consumer,outcomes in model['outcomes_by_consumer'].items():
        for outcome in outcomes:
            p=model['producer_by_outcome'][outcome]
            assert by_id[p['operation_id']]['end']<=by_id[consumer]['start'], 'feedback readiness'
    live={}
    for g in model['port_lifetimes']:
        a,b=by_id[g['acquire']],by_id[g['bell']]
        assert a['end']<=b['start'], 'Bell before injection'
        own={g['acquire'],g['bell']}
        for tick in range(a['start'],b['end']):
            for point in map(tuple,g['ports']+g['resources']):
                assert (tick,point) not in live, 'overlapping CCZ lifetime'
                live[tick,point]=g['group']
                assert occupancy.get((tick,point)) in own|{None}, 'live CCZ/port overwritten'
    latency=max(e['end'] for e in events)
    assert type(timeline['latency']) is int and latency==timeline['latency'], 'makespan'
    mixed=[t for t,ks in kinds_by_tick.items() if len(ks)>1]
    h_cnot=[t for t,ks in kinds_by_tick.items() if 'h' in ks and
            ks&{'linear_cnot','linear_cnot_layer','conditional_cnot','native_ccz_data_to_port'}]
    return dict(passed=True,latency=latency,operations=len(nodes),dependencies=len(model['dependencies']),
                conditional_CNOTs=sum(n['kind']=='conditional_cnot' for n in nodes.values()),
                H=sum(n['kind']=='h' for n in nodes.values()),
                CNOT_slots=sum(n.get('cnot_count',0) for n in nodes.values()),
                mixed_primitive_cycles=len(mixed),H_CNOT_overlap_cycles=len(h_cnot),
                checked_CCZ_lifetimes=len(model['port_lifetimes']),
                correctness='original branch/geometry proofs plus data-order/outcome-preserving disjoint-instrument exchanges',
                operation_counts=dict(Counter(n['kind'] for n in nodes.values())))
