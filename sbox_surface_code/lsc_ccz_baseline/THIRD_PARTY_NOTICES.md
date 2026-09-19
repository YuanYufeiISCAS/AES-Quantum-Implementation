# Third-party sources

## Lattice Surgery Compiler

The bundled subset of `liblsqecc` comes from
`https://github.com/latticesurgery-com/liblsqecc`, revision
`fddaecf0d929b0afa0ae72a1adc1df865fab4e18`.
Its GNU GPL version 3 license is preserved in `vendor/liblsqecc/LICENSE`.

The following files include our dispatch adaptation, distributed in this
submission on September 18, 2026:

- `include/lsqecc/ls_instructions/ls_instructions.hpp`: common-model instruction variant and dependencies.
- `src/ls_instructions/ls_instructions.cpp`: printing support for that variant.
- `src/patches/dense_patch_computation.cpp`: dispatch to the primitive adapter.
- `src/scheduler/wave_scheduler.cpp`: recognition of the adapted instruction.

The corresponding adapter source is supplied under `src/`. The graph-search
and wave-selection algorithms are retained. The backend therefore must be
identified as **adapted LSC-CCZ**, not unmodified LSC. Complete source for this
build is included; no prebuilt library from another experiment is required.

## Header dependencies

- nlohmann/json: `vendor/liblsqecc/external/json/LICENSE.MIT`.
- Tessil ordered-map: `vendor/liblsqecc/external/ordered-map/LICENSE`.
- cppitertools: `vendor/liblsqecc/external/include/cppitertools/LICENSE.md`.
- InfInt: Mozilla Public License 2.0; the original notice and license location
  are retained in `vendor/liblsqecc/external/include/infint/InfInt.h`.

Original copyright and license notices in the vendored headers are retained.
