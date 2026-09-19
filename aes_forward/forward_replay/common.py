"""Portable input/output; reuse the repository's native S-box verifier."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parent
sys.path.insert(0, str(REPOSITORY / 'sbox_surface_code'))
from sbox_compile.model import read, write

SCHEMA = 'aes128-fixed-forward-v1'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load_components():
    components = read(ROOT / 'data/components.json.gz')
    for name, case in (('state', 'inplace_3row'), ('key', 'cstar_4row')):
        components[name] = read(REPOSITORY / 'sbox_surface_code/paper_results' / (case + '.json.gz'))
    return components
