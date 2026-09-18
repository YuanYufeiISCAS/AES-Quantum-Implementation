"""Public command line: verify, reproduce, search and inspect."""

import argparse
from pathlib import Path
import sys
import subprocess

from .io import ROOT, load
from .reproduce import reproduce
from .search import search
from .verify import verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("verify", help="independently verify the saved circuit")
    check.add_argument("circuit", type=Path, nargs="?", default=ROOT / "circuits/aes.json")
    check.add_argument("--instance", type=Path, default=ROOT / "instances/aes.json")
    replay = commands.add_parser("reproduce", help="recompute AES from the matrix and fixed recipe")
    replay.add_argument("--backend", type=Path, default=ROOT / "build/linear_surface_backend")
    replay.add_argument("--instance", type=Path, default=ROOT / "instances/aes.json")
    replay.add_argument("--recipe", type=Path, default=ROOT / "recipes/aes_sota.json")
    replay.add_argument("--output", type=Path, default=ROOT / "runs/reproduce")
    fresh = commands.add_parser("search", help="fresh bounded search; does not read saved circuits or recipes")
    fresh.add_argument("--instance", type=Path, default=ROOT / "instances/aes.json")
    fresh.add_argument("--seed", type=int, default=1)
    fresh.add_argument("--restarts", type=int, default=8)
    fresh.add_argument("--backend", type=Path, default=ROOT / "build/linear_surface_backend")
    fresh.add_argument("--output", type=Path, default=ROOT / "runs/search")
    show = commands.add_parser("inspect", help="print the gate layers and output permutation")
    show.add_argument("circuit", type=Path, nargs="?", default=ROOT / "circuits/aes.json")
    show.add_argument("--instance", type=Path, default=ROOT / "instances/aes.json")
    show.add_argument("--paths", action="store_true", help="also print the routed grid vertices")
    args = parser.parse_args()
    try:
        if args.command in ("reproduce", "search") and not args.backend.is_file():
            raise ValueError(f"backend not found: {args.backend}; build with CMake first")
        if args.command == "reproduce":
            reproduce(args.backend, args.instance, args.recipe, args.output)
            print(f"saved: {args.output / 'aes.json'}")
        elif args.command == "search":
            search(args.backend, args.instance, args.seed, args.restarts, args.output)
            print(f"saved: {args.output / 'best.json'}")
        else:
            circuit = load(args.circuit)
            stats = verify(circuit, load(args.instance))
            print(f'PASS: {stats["cnots"]} CNOT / {stats["layers"]} VDP layers / {stats["cycles"]} cycles')
            if args.command == "inspect":
                print("output permutation (zero-based):", circuit["output_permutation"])
                for index, layer in enumerate(circuit["layers"]):
                    print(f'layer {index + 1}: ' + " ".join(
                        f'{op["control"]}->{op["target"]}' for op in layer["operations"]))
                    if args.paths:
                        for op in layer["operations"]:
                            print(f'  {op["control"]}->{op["target"]}: {op["path"]}')
    except (ValueError, OSError, KeyError, TypeError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"error: {error}\n")
    except KeyboardInterrupt:
        parser.exit(130, "interrupted; completed stage files were preserved\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
