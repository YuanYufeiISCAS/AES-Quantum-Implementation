"""Submission entry points for synthesis, concurrent scheduling, and replay."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .model import CASES, ROOT, fixture, objective, read, write
from .replay import replay
from .verify import require_verified


def publish(directory, circuit, *, recipe=None, report=None):
    directory = Path(directory).resolve()
    for protected in ("results", "parents", "recipes", "paper_results", "lsc_ccz_baseline", "sbox_compile", "tests"):
        if directory.is_relative_to(ROOT / protected):
            raise ValueError("output must be separate from submission inputs and source")
    require_verified(circuit)
    directory.mkdir(parents=True, exist_ok=False)
    write(directory / "circuit.json.gz", circuit)
    if recipe is not None:
        write(directory / "recipe.json.gz", recipe)
    if report is not None:
        write(directory / "report.json", report)


def main(argv=None):
    if not __debug__:
        raise ValueError('verification requires Python assertions; do not use -O')
    parser = argparse.ArgumentParser(
        description="Joint AES S-box compilation under the strict native cost model."
    )
    commands = parser.add_subparsers(dest="command", required=True)
    verify_parser = commands.add_parser(
        "verify", help="independently verify saved complete circuits"
    )
    verify_parser.add_argument("circuit", nargs="?", type=Path)
    verify_parser.add_argument("--case", choices=CASES)
    verify_parser.add_argument("--paper", action="store_true", help="verify the paper's concurrent schedules")
    reproduce_parser = commands.add_parser("reproduce", help="replay the paper's concurrent schedules")
    reproduce_parser.add_argument("--case", choices=CASES)
    reproduce_parser.add_argument("--output", type=Path)
    schedule_parser = commands.add_parser("schedule", help="schedule independent native primitives concurrently")
    schedule_parser.add_argument("circuit", nargs="?", type=Path)
    schedule_parser.add_argument("--case", choices=CASES)
    schedule_parser.add_argument("--paper", action="store_true", help="use the integrated paper reference circuit")
    schedule_parser.add_argument("--starts", type=int, default=2048)
    schedule_parser.add_argument("--output", type=Path, required=True)
    replay_parser = commands.add_parser("replay", help="rebuild saved circuit transformations")
    replay_parser.add_argument("--case", choices=CASES, required=True)
    replay_parser.add_argument("--output", type=Path)
    search_parser = commands.add_parser(
        "search", help="search from a parent without loading results or recipes"
    )
    search_parser.add_argument("--case", choices=CASES, required=True)
    search_parser.add_argument(
        "--parent", type=Path, help="optional independently verified parent circuit"
    )
    search_parser.add_argument("--seed", type=int, default=2026914746)
    search_parser.add_argument("--iterations", type=int, default=32)
    search_parser.add_argument(
        "--budget", type=int, default=100, help="deterministic work multiplier"
    )
    search_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "reproduce" or (args.command == "verify" and args.paper):
        from .parallel import paper_cases, reproduce

        if getattr(args, "circuit", None):
            parser.error("--paper cannot be combined with a circuit path")
        for entry, circuit, report in paper_cases(args.case):
            if args.command == "reproduce":
                circuit = reproduce(circuit)
                if args.output:
                    name = entry['file'].removesuffix('.json.gz')
                    publish(args.output / name, circuit)
            print(json.dumps(dict(case=entry['case'], role=entry['role'],
                                  objective=objective(circuit), passed=True)))
        return 0
    if args.command == "schedule":
        from .parallel import paper_cases, schedule

        if bool(args.circuit) == bool(args.case):
            parser.error("provide exactly one circuit path or --case")
        if args.paper and not args.case:
            parser.error("--paper requires --case")
        if args.output.exists():
            parser.error("output directory already exists")
        if args.paper:
            _, stored, _ = next(paper_cases(args.case))
            source = stored['source']
        else:
            source = read(args.circuit) if args.circuit else fixture(args.case)
        circuit, report = schedule(source, starts=args.starts)
        publish(args.output, circuit, report=report)
        print(json.dumps(dict(case=circuit['case'], objective=objective(circuit), passed=True)))
        return 0
    if args.command == "verify":
        if args.circuit and args.case:
            parser.error("choose a circuit path or a case")
        circuits = (
            [read(args.circuit)]
            if args.circuit
            else [fixture(c) for c in ([args.case] if args.case else CASES)]
        )
        if not args.circuit and not args.case:
            circuits.append(read(ROOT / "results/inplace_3row_cover.json.gz"))
        for circuit in circuits:
            report = require_verified(circuit)
            print(
                json.dumps(
                    {
                        "case": circuit["case"],
                        "objective": objective(circuit),
                        "checked_inputs": report["checked_inputs"],
                        "passed": True,
                    }
                )
            )
        return 0
    if args.command == "replay":
        parent = fixture(args.case, parent=True)
        recipe = read(ROOT / "recipes" / f"{args.case}.json.gz")
        circuit = replay(parent, recipe, progress=lambda x: print(json.dumps(x), flush=True))
        if args.output:
            publish(args.output, circuit)
        print(json.dumps({"case": args.case, "objective": objective(circuit), "passed": True}))
        return 0
    from .search import optimize

    parent = read(args.parent) if args.parent else fixture(args.case, parent=True)
    if parent["case"] != args.case:
        parser.error("parent belongs to another case")
    if args.output.exists():
        parser.error("output directory already exists")
    result = optimize(
        parent,
        seed=args.seed,
        iterations=args.iterations,
        budget=args.budget,
        progress=lambda x: print(json.dumps(x), flush=True),
    )
    publish(args.output, result["circuit"], recipe=result["recipe"], report=result["report"])
    print(
        json.dumps({"case": args.case, "objective": objective(result["circuit"]), "passed": True})
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, KeyError) as error:
        print(f"sbox_compile: {error}", file=sys.stderr)
        raise SystemExit(1)
