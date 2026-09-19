"""Reproduce the fixed forward AES-128 implementation; no search or synthesis."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from forward_replay.common import ROOT,load_components,read,write
from forward_replay.circuit import assemble
from forward_replay.verify import verify


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='new directory for the complete circuit and verification result')
    parser.add_argument('--verify',type=Path,help='independently replay an existing circuit.json.gz')
    parser.add_argument('--random-trials',type=int,default=256)
    args = parser.parse_args(argv)
    if args.output and args.verify:
        parser.error('choose reproduction or verification, not both')
    if args.output and args.output.exists():
        parser.error('output directory already exists')
    if args.output:
        output = args.output.resolve()
        for name in ('data','forward_replay','tests'):
            if output.is_relative_to(ROOT/name):
                parser.error('output must be separate from inputs and code')
    program = read(args.verify) if args.verify else assemble(load_components())
    report = verify(program,args.random_trials)
    if args.output:
        output.mkdir(parents=True,exist_ok=False)
        write(output/'circuit.json.gz',program)
        write(output/'verification.json',report)
    print(json.dumps(dict(passed=True,latency=report['timing']['latency'],
        reserved_patches=report['geometry']['reserved_patches'],
        resources=report['resources'],checked_AES_inputs=report['functional']['cases'],
        known_ciphertext=report['functional']['known_ciphertexts'][0]),indent=2))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (AssertionError,ValueError,KeyError,IndexError,TypeError,OSError,subprocess.SubprocessError) as error:
        print(f'aes_forward: {error}',file=sys.stderr)
        raise SystemExit(1)
