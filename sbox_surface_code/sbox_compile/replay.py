"""Rebuild a saved transformation chain without reading a final circuit."""

from __future__ import annotations

import copy

from .compiler import compile_circuit
from .geometry import _require
from .lower import routes_from_artifact
from .model import objective
from .verify import require_verified
from .windows import move_placement


def transition(old, new):
    fields = ("linear_candidates", "logical_schedule", "logical_gadgets", "search", "schedule")
    inputs = {k: copy.deepcopy(new[k]) for k in fields if old[k] != new[k]}
    routes = routes_from_artifact(new)
    if routes != routes_from_artifact(old):
        inputs["routes"] = routes
    inputs["expected_cost"] = copy.deepcopy(new["cost"])
    return inputs


def apply_inputs(source, inputs):
    allowed = {
        "linear_candidates",
        "logical_schedule",
        "logical_gadgets",
        "search",
        "schedule",
        "routes",
        "expected_cost",
    }
    _require(set(inputs) <= allowed, "unexpected replay input")
    front = copy.deepcopy(source)
    for key in ("linear_candidates", "logical_schedule", "logical_gadgets", "search"):
        if key in inputs:
            front[key] = copy.deepcopy(inputs[key])
    stages = inputs.get("schedule", source["schedule"])["stages"]
    result = compile_circuit(
        front, routes=inputs.get("routes", routes_from_artifact(source)), stages=stages
    )
    _require(result["cost"] == inputs["expected_cost"], "replay cost mismatch")
    return result


def replay(parent, recipe, *, progress=None):
    require_verified(parent)
    _require(recipe["case"] == parent["case"], "recipe case mismatch")
    current = copy.deepcopy(parent)
    for index, step in enumerate(recipe["steps"]):
        before = current
        start, stop = step["start"], step["stop"]
        boundaries = step["boundaries"]
        n = parent["logical_width"]
        _require(0 <= start <= stop < len(current["linear_candidates"]), "invalid recipe window")
        _require(len(boundaries) == stop - start + 2, "incomplete boundary sequence")
        _require(boundaries[0] == boundaries[-1] == list(range(n)), "unpaid external boundary")
        window = apply_inputs(current, step["window"])
        for j, (old, new) in enumerate(
            zip(current["linear_candidates"], window["linear_candidates"], strict=True)
        ):
            if start <= j <= stop:
                for field, boundary in (
                    ("input_placement", boundaries[j - start]),
                    ("output_placement", boundaries[j - start + 1]),
                ):
                    _require(
                        new["choice"][field] == move_placement(old["choice"][field], boundary),
                        "recipe boundary does not match materialized matrix",
                    )
            else:
                _require(old == new, "recipe changes a matrix outside its window")
        current = window
        if step["nonlinear"] is not None:
            current = apply_inputs(window, step["nonlinear"])
            _require(
                current["linear_candidates"] == window["linear_candidates"]
                and current["logical_schedule"] == window["logical_schedule"]
                and current["logical_gadgets"] == window["logical_gadgets"],
                "nonlinear replay changes the linear/logical problem",
            )
        if progress:
            progress(
                {
                    "step": index,
                    "before": objective(before),
                    "window": objective(window),
                    "after": objective(current),
                }
            )
    _require(current["cost"] == recipe["expected_cost"], "final replay cost mismatch")
    _require(objective(current) <= objective(parent), "replayed incumbent regressed")
    require_verified(current)
    return current
