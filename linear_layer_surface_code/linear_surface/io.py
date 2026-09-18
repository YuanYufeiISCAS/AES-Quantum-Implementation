"""JSON files and the private C++ process protocol (standard library only)."""

import json
import os
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def load(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON field: {key}")
            result[key] = value
        return result

    with Path(path).open(encoding="utf-8") as stream:
        return json.load(stream, object_pairs_hook=unique)


def save(path, value):
    """Atomically replace a generated JSON file; never leave a partial witness."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix=f".{path.name}.", delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write("\n")
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def normalize(result):
    """Export only the mathematical instance and complete routed witness."""
    if result["mode"] != "vdp" or any(layer.get("mode", "vdp") != "vdp" for layer in result["layers"]):
        raise ValueError("backend returned a non-VDP circuit")
    return {
        "n": result["n"], "mode": result["mode"], "layout": result["layout"],
        "target_rows_hex": result["target_rows_hex"],
        "output_permutation": result["output_permutation"],
        "stats": {"cnots": result["stats"]["cnots"],
                  "layers": result["stats"]["layers"],
                  "cycles": result["stats"]["cycles"]},
        "layers": [
            {"cycles": layer["cycles"],
             "operations": [{key: gate[key] for key in ("control", "target", "path")}
                            for gate in layer["operations"]]}
            for layer in result["layers"]
        ],
    }


def request(instance, *, circuit=None, free_output=True):
    layout = instance["layout"]
    columns = layout["grid_cols"]
    lines = [
        f'{instance["n"]} {layout["data_rows"]} {layout["data_cols"]} '
        f'{"free" if free_output else "fixed"}',
        " ".join(str(r * columns + c) for r, c in layout["data_pos"]),
        " ".join(instance["target_rows_hex"]),
        "-" if free_output else " ".join(map(str, circuit["output_permutation"])),
    ]
    if circuit is not None:
        lines.append(str(len(circuit["layers"])))
        for layer in circuit["layers"]:
            fields = [len(layer["operations"])]
            for gate in layer["operations"]:
                path = gate["path"]
                fields.extend([gate["control"], gate["target"], len(path)])
                fields.extend(r * columns + c for r, c in path)
            lines.append(" ".join(map(str, fields)))
    return "\n".join(lines) + "\n"
