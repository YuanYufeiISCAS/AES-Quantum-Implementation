"""Load source gates and lower them without external repositories or circuits."""
import gzip
import json
from pathlib import Path

from lower import specification

ROOT = Path(__file__).resolve().parents[1]
DIRECTIONS = ('forward', 'inverse')


def load_source(direction):
    if direction not in DIRECTIONS:
        raise ValueError('unknown AES direction')
    path = ROOT / 'inputs' / ('aes_' + direction + '.json.gz')
    source = json.loads(gzip.decompress(path.read_bytes()))
    if source['schema'] != 'lsc-ccz-source-v1':
        raise ValueError('unsupported source format')
    return source


def build_spec(source):
    layout = source['layout']
    return specification(source['name'], source['segments'], layout['patches'],
        layout['rows'], layout['cols'], policy='dependency_overlap',
        h_policy='grouped_pauli', extra={'interface': source['interface']})
