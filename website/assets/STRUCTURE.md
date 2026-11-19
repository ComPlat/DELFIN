# fac-Ir(ppy)3 — structure provenance

The downloadable `ir-ppy3.xyz` is an actual 61-atom MANTA construction of
tris(2-phenylpyridinato)iridium(III), C33H24IrN3. It is a starting geometry,
not an experimental crystal structure or a DFT-optimised geometry.

Source SMILES: the `Ir(ppy)3` entry in `tests/test_isomer_benchmark.py`.
DELFIN source commit: `a7e22dcda398f3fcf55c9d4d5f38e91e84f73871`.
Generated on 2026-11-16 with the local DELFIN environment (RDKit 2025.09.6).

```sh
python -m delfin.cli_manta \
  '[N+]12=CC=CC=C1C(C=CC=C3)=C3[Ir-3]24([N+]5=CC=CC=C5C6=C4C=CC=C6)([N+]7=CC=CC=C78)C9=C8C=CC=C9' \
  -o preview --max-isomers 1 --quality fast --seeds 1 --num-confs 1 \
  --construction builder --no-binding-modes
```

The command emitted three structures and reported its isomer-cap warning
(exit 2). `000__fac.xyz` was selected unchanged for the illustration. This
is deliberately one displayed geometry, not a complete conformer ensemble.
The capped preview run is not an energetic ranking or a production calculation.

`molecule.js` embeds the exact XYZ coordinates to allow offline/file previews.
Display bonds are inferred when distance < 1.15 × the sum of covalent radii
(H 0.31, C 0.76, N 0.71, Ir 1.41 Å). They indicate connectivity, not bond order.
The resulting Ir coordination is three C and three N donors. The viewer centres,
rotates and scales the coordinates only; it does not deform the structure.

Colours: Ir gold, N blue, C slate, H white. The viewer supports pointer drag,
arrow-key rotation, +/- zoom, reset, hiding hydrogens and global motion pause.
Canvas rendering is used so a WebGL-capable GPU is not required.
