"""Bounded multi-start search, independent of saved circuits and recipes."""

from collections import Counter
from copy import deepcopy
from pathlib import Path
import random

from .io import load, save
from .reproduce import backend_run, refine, score
from .verify import check_instance, integer, verify


def diversity_key(circuit):
    gates = sorted((op["control"], op["target"])
                   for layer in circuit["layers"] for op in layer["operations"])
    return tuple(gates), tuple(circuit["output_permutation"])


class Population:
    """At most eight verified circuits; always retain the lexicographic best."""

    def __init__(self, instance):
        self.instance = instance
        self.entries = {}
        self.selections = 0

    @property
    def best(self):
        return min(self.entries.values(), key=lambda entry: (score(entry["circuit"]), entry["key"]))["circuit"]

    def admit(self, circuit):
        verify(circuit, self.instance)
        key = diversity_key(circuit)
        previous = self.entries.get(key)
        if previous is not None and score(previous["circuit"]) <= score(circuit):
            return
        self.entries[key] = {"key": key, "circuit": deepcopy(circuit),
                             "uses": previous["uses"] if previous else 0}
        best = self.best
        selected = [diversity_key(best)]
        available = [k for k, entry in self.entries.items()
                     if k != selected[0]
                     and entry["circuit"]["stats"]["cycles"] <= best["stats"]["cycles"] + 2
                     and entry["circuit"]["stats"]["cnots"] <= best["stats"]["cnots"] + 8]

        def distance(a, b):
            left, right = Counter(a[0]), Counter(b[0])
            return (sum(abs(left[g] - right[g]) for g in left.keys() | right.keys())
                    + sum(x != y for x, y in zip(a[1], b[1])))

        while available and len(selected) < 8:
            chosen = min(available, key=lambda k: (-min(distance(k, other) for other in selected),
                                                  score(self.entries[k]["circuit"]), k))
            selected.append(chosen)
            available.remove(chosen)
        self.entries = {k: self.entries[k] for k in selected}

    def select(self):
        best_key = diversity_key(self.best)
        others = [entry for key, entry in self.entries.items() if key != best_key]
        chosen = (min(others, key=lambda entry: (entry["uses"], entry["key"]))
                  if self.selections % 2 and others else self.entries[best_key])
        chosen["uses"] += 1
        self.selections += 1
        return deepcopy(chosen["circuit"])


def search(backend, instance_path, seed, restarts, output):
    integer(seed, "seed", 0, (1 << 64) - 1)
    integer(restarts, "restarts", 0, 4096)
    instance = load(instance_path)
    check_instance(instance)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    rng = random.Random(seed)
    pool = Population(instance)
    initial, work = backend_run(backend, ["synth", seed], instance, pool=pool)
    records = [{"kind": "synthesis", "seed": seed, **work, **initial["stats"]}]
    save(output / "best.json", pool.best)
    # Three fixed local search policies: greedy, free-output walk, and
    # fixed-output walk. Only budgets/seeds vary; no experimental switches.
    profiles = ((0, 0, False), (6, 4, True), (4, 2, False))
    for index in range(restarts):
        if index > 0 and index % 4 == 0:
            synthesis_seed = rng.getrandbits(64)
            candidate, work = backend_run(backend, ["synth", synthesis_seed], instance, pool=pool)
            records.append({"kind": "synthesis", "seed": synthesis_seed, **work, **candidate["stats"]})
        growth, slack, free_output = profiles[index % len(profiles)]
        stage = {"seed": rng.getrandbits(64), "proposals": 65536, "routes": 512,
                 "rounds": 128, "growth": growth, "slack": slack, "free_output": free_output}
        candidate, work = refine(backend, instance, pool.select(), stage, pool=pool)
        records.append({"kind": "refinement", "parameters": stage, **work, **candidate["stats"]})
        save(output / "best.json", pool.best)
        save(output / "run.json", {"seed": seed, "restarts": restarts, "attempts": records})
        print(f'restart {index + 1}: best {pool.best["stats"]}', flush=True)
    save(output / "run.json", {"seed": seed, "restarts": restarts, "attempts": records})
    return pool.best
