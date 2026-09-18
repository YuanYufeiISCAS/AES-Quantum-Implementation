"""layout."""

# Adapted from the authors' sbox/surface/ccz_scheduler.py; campaign infrastructure removed.
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple, Iterable

Coord = tuple[int, int]


@dataclass(frozen=True)
class FoundryLayout:
    """Surface-code foundry with a data rectangle and one outer port ring."""

    data_rows: int
    data_cols: int

    @property
    def n_data(self) -> int:
        return self.data_rows * self.data_cols

    @property
    def grid_rows(self) -> int:
        return 2 * self.data_rows + 1

    @property
    def grid_cols(self) -> int:
        return 2 * self.data_cols + 1

    @property
    def min_row(self) -> int:
        return -2

    @property
    def max_row(self) -> int:
        return 2 * self.data_rows + 2

    @property
    def min_col(self) -> int:
        return -2

    @property
    def max_col(self) -> int:
        return 2 * self.data_cols + 2

    def data_coord(self, physical_index: int) -> Coord:
        if physical_index < 0 or physical_index >= self.n_data:
            raise ValueError(f"physical data index out of range: {physical_index}")
        return (2 * (physical_index // self.data_cols), 2 * (physical_index % self.data_cols))

    def data_coords(self) -> List[Coord]:
        return [self.data_coord(i) for i in range(self.n_data)]

    def all_vertices(self) -> Set[Coord]:
        return {
            (r, c)
            for r in range(self.min_row, self.max_row + 1)
            for c in range(self.min_col, self.max_col + 1)
        }

    def neighbors(self, coord: Coord) -> Iterable[Coord]:
        r, c = coord
        for nr, nc in ((r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)):
            if self.min_row <= nr <= self.max_row and self.min_col <= nc <= self.max_col:
                yield (nr, nc)

    def port_coords(self) -> List[Coord]:
        """Return the Section 5.2 ring order, starting at the upper-left port."""
        ports: List[Coord] = []
        top = self.min_row
        bottom = self.max_row
        left = self.min_col
        right = self.max_col
        inner_max_row = 2 * self.data_rows
        inner_max_col = 2 * self.data_cols
        ports.extend(((top, c) for c in range(0, inner_max_col + 1, 2)))
        ports.extend(((r, right) for r in range(0, inner_max_row, 2)))
        ports.extend(((bottom, c) for c in range(inner_max_col, -1, -2)))
        ports.extend(((r, left) for r in range(inner_max_row - 2, -1, -2)))
        return ports

    def port_triples(self) -> List[Tuple[Coord, Coord, Coord]]:
        ports = self.port_coords()
        return [(ports[i], ports[i + 1], ports[i + 2]) for i in range(0, len(ports) - 2, 3)]

    def to_dict(self) -> Dict[str, object]:
        return {
            "data_rows": self.data_rows,
            "data_cols": self.data_cols,
            "grid_rows": self.grid_rows,
            "grid_cols": self.grid_cols,
            "extended_rows": [self.min_row, self.max_row],
            "extended_cols": [self.min_col, self.max_col],
            "port_count": len(self.port_coords()),
            "port_triple_count": len(self.port_triples()),
            "ring_order": "upper-left clockwise, every three consecutive ports form one |CCZ> state",
        }


def resource_coord_for_port(layout: FoundryLayout, port: Coord) -> Coord:
    r, c = port
    if r == layout.min_row:
        return (r - 1, c)
    if r == layout.max_row:
        return (r + 1, c)
    if c == layout.min_col:
        return (r, c - 1)
    if c == layout.max_col:
        return (r, c + 1)
    return (r - 1, c)
