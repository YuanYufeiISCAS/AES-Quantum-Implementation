"""Matrix-to-circuit runs; saved witnesses are never synthesis inputs."""

import json
from pathlib import Path
import subprocess
import time

from .io import load, normalize, request, save
from .verify import check_instance, integer, require, verify


def score(circuit):
    stats = circuit["stats"]
    return stats["cycles"], stats["layers"], stats["cnots"]


def backend_run(backend, arguments, instance, circuit=None, free_output=True, pool=None):
    started = time.monotonic()
    completed = subprocess.run([str(Path(backend).resolve()), *map(str, arguments)],
                               input=request(instance, circuit=circuit, free_output=free_output),
                               text=True, stdout=subprocess.PIPE, check=True)
    # The private stream consists of complete JSON values; whitespace may
    # include newlines inside a witness. Incomplete output is an error.
    decoder, offset, candidates, work = json.JSONDecoder(), 0, [], {}
    while offset < len(completed.stdout):
        while offset < len(completed.stdout) and completed.stdout[offset].isspace():
            offset += 1
        if offset == len(completed.stdout):
            break
        event, offset = decoder.raw_decode(completed.stdout, offset)
        if event["kind"] in ("candidate", "exploration"):
            candidate = normalize(event["result"])
            verify(candidate, instance)
            if circuit is not None and not free_output:
                if candidate["output_permutation"] != circuit["output_permutation"]:
                    raise ValueError("fixed-output stage changed the permutation")
            if event["kind"] == "candidate":
                candidates.append(candidate)
            if pool is not None:
                pool.admit(candidate)
        elif event["kind"] == "complete":
            work = {key: event[key] for key in ("proposals", "routes", "rounds")}
        else:
            raise ValueError("unexpected backend event")
    if not candidates:
        raise ValueError("backend produced no verified complete circuit")
    return min(candidates, key=score), {**work, "seconds": round(time.monotonic() - started, 3)}


def refine(backend, instance, circuit, stage, pool=None):
    arguments = ["refine", *(stage[key] for key in ("seed", "proposals", "routes", "rounds", "growth", "slack")),
                 int(stage["free_output"])]
    return backend_run(backend, arguments, instance, circuit, stage["free_output"], pool)


def check_recipe(recipe):
    require(set(recipe) == {"synthesis_seed", "initial_expected", "stages"}, "unexpected recipe fields")
    integer(recipe["synthesis_seed"], "synthesis seed", 0, (1 << 64) - 1)
    require(type(recipe["stages"]) is list and 1 <= len(recipe["stages"]) <= 64, "invalid recipe stages")
    limits = {"seed": (0, (1 << 64) - 1), "proposals": (1, 262144), "routes": (1, 4096),
              "rounds": (1, 4096), "growth": (0, 8), "slack": (0, 8)}
    for stage in recipe["stages"]:
        require(set(stage) == set(limits) | {"free_output", "expected"}, "unexpected stage fields")
        for key, (lower, upper) in limits.items():
            integer(stage[key], key, lower, upper)
        require(type(stage["free_output"]) is bool, "invalid output policy")
    for expected in [recipe["initial_expected"], *(stage["expected"] for stage in recipe["stages"])]:
        require(set(expected) == {"cnots", "layers", "cycles"}, "invalid expected cost")
        for key in expected:
            integer(expected[key], key)
        require(expected["cycles"] == 2 * expected["layers"], "inconsistent expected cycles")


def reproduce(backend, instance_path, recipe_path, output):
    instance, recipe = load(instance_path), load(recipe_path)
    check_instance(instance)
    check_recipe(recipe)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    records = []
    best, work = backend_run(backend, ["synth", recipe["synthesis_seed"]], instance)
    stages = [None, *recipe["stages"]]
    for index, stage in enumerate(stages):
        if stage is not None:
            best, work = refine(backend, instance, best, stage)
        expected = recipe["initial_expected"] if stage is None else stage["expected"]
        save(output / f"stage_{index:02d}.json", best)
        records.append({"stage": index, **work, **best["stats"]})
        save(output / "run.json", {"recipe": recipe, "stages": records})
        print(f'stage {index}: {best["stats"]}; {work["seconds"]} s', flush=True)
        if best["stats"] != expected:
            raise ValueError(f"stage {index}: expected {expected}, obtained {best['stats']}")
    save(output / "aes.json", best)
    return best
