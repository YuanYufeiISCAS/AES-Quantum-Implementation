"""Final scheduling of verified circuits without changing gates or routes."""

from copy import deepcopy

from .model import ROOT, read
from .parallel_model import make_model, reference_timeline
from .parallel_scheduler import compact, search
from .parallel_verify import audit

SCHEMA = 'sbox-surface-code-parallel-v1'


def _cost(source, timeline):
    # Stage-duration sums are not a latency formula for concurrent execution.
    result = {k: v for k, v in source['cost'].items()
              if k.startswith('N_') or k == 'R_CCZ'}
    result['latency'] = max(e['end'] for e in timeline['events'])
    return result


def materialize(source, timeline, *, split_linear=False):
    from .verify import require_verified

    require_verified(source)
    audit(source, timeline, split_linear=split_linear)
    space = deepcopy(source['space_time'])
    cost = _cost(source, timeline)
    space['reserved_patch_cycles'] = space['reserved_module_patches'] * cost['latency']
    return dict(schema=SCHEMA, case=source['case'], source=deepcopy(source),
                split_linear=split_linear, schedule=deepcopy(timeline),
                cost=cost, space_time=space)


def verify(circuit):
    from .verify import require_verified

    if circuit['schema'] != SCHEMA or type(circuit['split_linear']) is not bool:
        raise ValueError('invalid concurrent circuit schema')
    source = circuit['source']
    if source.get('schema') == SCHEMA:
        raise ValueError('nested concurrent circuit')
    original = require_verified(source)
    if circuit['case'] != source['case']:
        raise ValueError('concurrent circuit case mismatch')
    timing = audit(source, circuit['schedule'], split_linear=circuit['split_linear'])
    cost = _cost(source, circuit['schedule'])
    if circuit['cost'] != cost:
        raise ValueError('concurrent cost differs from the complete event trace')
    space = deepcopy(source['space_time'])
    space['reserved_patch_cycles'] = space['reserved_module_patches'] * cost['latency']
    if circuit['space_time'] != space:
        raise ValueError('concurrent space-time cost mismatch')
    if timing['H'] != cost['N_H'] or timing['conditional_CNOTs'] != cost['N_conditional_CNOT']:
        raise ValueError('missing nonlinear operation')
    if timing['CNOT_slots'] != cost['N_CNOT_linear'] + 3*cost['N_CCZ'] + cost['N_conditional_CNOT']:
        raise ValueError('missing CNOT operation')
    return dict(passed=True, case=circuit['case'], cost=cost,
                checked_inputs=original['checked_inputs'], scratch_zero=original['scratch_zero'],
                checked_operations=original['checked_operations'],
                checked_events=timing['operations'], timing=timing, errors=[])


def schedule(source, *, starts=2048):
    """Run the same two finite scheduling portfolios used for the paper."""
    from .verify import require_verified

    if type(starts) is not int or starts < 1:
        raise ValueError('starts must be a positive integer')
    if source.get('schema') == SCHEMA:
        source = source['source']
    require_verified(source)
    variants = []
    for split in (False, True):
        model = make_model(source, split_linear=split)
        initial = reference_timeline(model)
        audit(source, initial, split_linear=split)
        result = search(model, starts=starts, initial=initial)
        checked = audit(source, result['best'], split_linear=split)
        variants.append(dict(split_linear=split, result=result, audit=checked))
    selected = min(variants, key=lambda v: v['audit']['latency'])
    circuit = materialize(source, selected['result']['best'], split_linear=selected['split_linear'])
    report = dict(starts_per_variant=starts,
                  source_latency=source['cost']['latency'],
                  latency=circuit['cost']['latency'],
                  variants=[dict(split_linear=v['split_linear'],
                    latency=v['audit']['latency'], cpu_seconds=v['result']['cpu_seconds'])
                    for v in variants])
    return circuit, report


def paper_cases(case=None):
    entries = read(ROOT / 'paper_results/index.json')
    if case:
        entries = [e for e in entries if e['case'] == case and e['role'] == 'reference']
    for entry in entries:
        circuit = read(ROOT / 'paper_results' / entry['file'])
        report = verify(circuit)
        if report['cost'] != entry['cost'] or circuit['case'] != entry['case']:
            raise ValueError('paper index differs from the verified circuit')
        yield entry, circuit, report


def reproduce(circuit):
    """Recompute start times from a checked published resource-use order."""
    verify(circuit)
    split = circuit['split_linear']
    rebuilt = compact(make_model(circuit['source'], split_linear=split), circuit['schedule'])
    result = materialize(circuit['source'], rebuilt, split_linear=split)
    if result['cost'] != circuit['cost']:
        raise ValueError('published schedule is not fully compacted')
    return result
