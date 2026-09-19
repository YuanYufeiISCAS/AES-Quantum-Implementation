"""Reproduce or independently replay the adapted LSC-CCZ AES baseline."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))
from audit import audit
from inputs import DIRECTIONS, build_spec, load_source
from quantum import verify_templates
from source_verify import verify_sources


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, sort_keys=True, separators=(',', ':'))
        stream.write('\n')


def main(argv=None):
    if not __debug__:
        raise ValueError('verification requires Python assertions; do not use -O')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('reproduce', 'verify'))
    parser.add_argument('--direction', choices=(*DIRECTIONS, 'both'), default='both')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--runner', type=Path, default=ROOT / 'build/lsc_ccz_runner')
    args = parser.parse_args(argv)
    directions = DIRECTIONS if args.direction == 'both' else (args.direction,)
    sources = {d: load_source(d) for d in DIRECTIONS}
    functional = verify_sources(sources['forward'], sources['inverse'])
    output = args.output.resolve()
    expected = json.loads((ROOT / 'expected.json').read_text())
    if args.command == 'reproduce':
        if not args.runner.is_file():
            raise ValueError('build lsc_ccz_runner before reproducing the baseline')
        for protected in ('src', 'vendor', 'inputs', 'tests', 'build'):
            if output.is_relative_to(ROOT / protected):
                raise ValueError('output must be separate from source, inputs, and build files')
        output.mkdir(parents=True, exist_ok=False)
        save(output / 'source_verification.json', functional)
    summaries = {}
    for direction in directions:
        spec = build_spec(sources[direction])
        quantum = verify_templates(spec)
        path = output / direction
        if args.command == 'reproduce':
            path.mkdir()
            save(path / 'input.json', spec)
            with (path / 'compiler.log').open('x') as log:
                process = subprocess.run([str(args.runner.resolve()), str(path / 'input.json'),
                    str(path / 'trace.jsonl')], stdout=subprocess.PIPE, stderr=log,
                    text=True, check=True, timeout=spec['timeout_seconds'] + 60)
            result = json.loads(process.stdout)
            save(path / 'result.json', result)
        else:
            saved = json.loads((path / 'input.json').read_text())
            if saved != json.loads(json.dumps(spec)):
                raise ValueError('saved input differs from source-derived lowering')
            result = json.loads((path / 'result.json').read_text())
        checked = audit(spec, path / 'trace.jsonl', result)
        summary = {k: checked[k] for k in expected[direction]}
        if summary != expected[direction]:
            raise ValueError(f'{direction}: independently checked costs differ from expected.json: {summary}')
        if args.command == 'reproduce':
            save(path / 'audit.json', checked)
            save(path / 'quantum.json', quantum)
        summaries[direction] = summary
        print(json.dumps(dict(direction=direction, passed=True, **summary)), flush=True)
    if len(summaries) == 2:
        total = dict(latency_logical_cycles=sum(s['latency_logical_cycles'] for s in summaries.values()),
                     CCZ_states=sum(s['CCZ_states'] for s in summaries.values()),
                     reserved_patches=max(s['reserved_patches'] for s in summaries.values()))
        if args.command == 'reproduce':
            save(output / 'summary.json', dict(passed=True, directions=summaries, total=total))
        print(json.dumps(dict(total=total, passed=True)), flush=True)
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (AssertionError, ValueError, OSError, KeyError, subprocess.SubprocessError) as error:
        print(f'lsc_ccz_baseline: {error}', file=sys.stderr)
        raise SystemExit(1)
