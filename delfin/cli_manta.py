"""Command-line interface for MANTA — the deterministic, complete
coordination-isomer × conformer manifold generator.

    delfin-manta "SMILES" -o out/ [--rank] [--method gfn2|gfnff] ...

From a (metal) SMILES it constructs the coordination-isomer manifold and writes
one XYZ per emitted isomer/conformer plus a ``manifest.json``.

Energy ranking (GFN2/GFN-FF, requires ``xtb`` on PATH) is **OFF by default**, so
the output is byte-identical to the library default; pass ``--rank`` to enable it.
Charge is taken from the SMILES (e.g. ``[Co+3]``).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

# "All": a guard high above any real coordination-isomer set (octahedral MABCDEF
# = 30; with chelate conformers still only low hundreds). The enumeration itself
# is provably complete (Burnside-Pólya); generation terminates naturally once the
# finite isomer set is exhausted, so this only bounds pathological blow-ups. If it
# is ever actually hit we WARN — completeness is never silently violated.
_ALL_ISOMERS = 100_000


def _safe_name(label: str, idx: int) -> str:
    """Filesystem-safe ``NNN__label.xyz`` from an isomer label."""
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", (label or "").strip()).strip("_")
    return f"{idx:03d}__{s or 'isomer'}.xyz"


def _atom_lines(block: str) -> list:
    """Non-empty coordinate lines, skipping any existing count/comment header."""
    lines = [ln for ln in block.splitlines() if ln.strip()]
    if len(lines) >= 2 and lines[0].strip().isdigit():
        return lines[2:]  # already standard XYZ -> drop count + comment
    return lines


def _to_xyz(block: str, comment: str) -> str:
    """Wrap a bare coordinate block in a valid standard XYZ file (count + comment)."""
    atoms = _atom_lines(block)
    comment = (comment or "MANTA").replace("\n", " ").strip() or "MANTA"
    return f"{len(atoms)}\n{comment}\n" + "\n".join(atoms) + "\n"


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="delfin-manta",
        description=(
            "MANTA: construct the complete coordination-isomer × conformer "
            "manifold from a (metal) SMILES. Deterministic, force-field-free at "
            "the metal, license-clean. Output is a UFF-quality STARTING geometry "
            "— relax with xTB/DFT before production calculations."
        ),
    )
    p.add_argument("smiles", help="input SMILES (charge encoded in the SMILES, e.g. '[Co+3]')")
    p.add_argument("-o", "--out", default=Path("manta_out"), type=Path,
                   help="output directory (default: ./manta_out)")
    p.add_argument("--rank", action="store_true",
                   help="energy-rank the ensemble with xtb (default: off / byte-identical)")
    p.add_argument("--method", choices=["gfn2", "gfnff", "gfn1", "gfn0"], default="gfn2",
                   help="ranking Hamiltonian when --rank is set (default: gfn2)")
    p.add_argument("--max-isomers", type=int, default=None, dest="max_isomers",
                   help="optionally cap the number of isomers (default: ALL — "
                        "the complete enumerated set)")
    p.add_argument("--quality", choices=["fast", "normal", "max", "extreme"], default=None,
                   help="conformer-depth preset (seeds x templates x cap): "
                        "fast(12 seeds) | normal(20) | max(40) | extreme(60). More seeds = "
                        "more conformers = higher chance the GLOBAL MINIMUM is in the manifold. "
                        "(default: library default ~normal)")
    p.add_argument("--seeds", type=int, default=None, dest="seeds",
                   help="override the ETKDG conformer-seed count directly (the key "
                        "completeness/speed switch; overrides --quality's seed count). "
                        "Higher = more conformers/rotamers, slower.")
    p.add_argument("--num-confs", type=int, default=None, dest="num_confs",
                   help="conformers embedded per isomer (default: 200)")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="suppress the per-structure listing")
    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)

    # Ranking is env-gated inside the converter; flip the switches BEFORE import.
    if args.rank:
        os.environ["DELFIN_FFFREE_GFNFF_RANK"] = "1"
        os.environ["DELFIN_CONF_RANK_METHOD"] = args.method

    from delfin.smiles_converter import smiles_to_xyz_isomers

    cap = args.max_isomers if args.max_isomers is not None else _ALL_ISOMERS
    kwargs = {"max_isomers": cap}
    if args.quality is not None:
        kwargs["quality_mode"] = args.quality
    if args.seeds is not None:
        kwargs["seeds_override"] = args.seeds
    if args.num_confs is not None:
        kwargs["num_confs"] = args.num_confs

    result = smiles_to_xyz_isomers(args.smiles, **kwargs)
    if isinstance(result, tuple) and len(result) == 2:
        isomers, error = result
    else:
        isomers, error = result, None

    if error:
        print(f"delfin-manta: error: {error}", file=sys.stderr)
        return 1
    if not isomers:
        print("delfin-manta: error: no structures generated", file=sys.stderr)
        return 1

    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    manifest = []
    for i, (xyz, label) in enumerate(isomers):
        fname = _safe_name(label, i)
        comment = f"{label}  |  {args.smiles}" if label else args.smiles
        (out / fname).write_text(_to_xyz(xyz, comment))
        natoms = len(_atom_lines(xyz))
        manifest.append({"index": i, "label": label, "file": fname, "natoms": natoms})
        if not args.quiet:
            print(f"  [{i:03d}] {label or '(single)'}  ({natoms} atoms)  -> {fname}")

    (out / "manifest.json").write_text(json.dumps(
        {
            "smiles": args.smiles,
            "count": len(manifest),
            "ranked": bool(args.rank),
            "method": args.method if args.rank else None,
            "isomers": manifest,
        },
        indent=2,
    ))

    print(f"\ndelfin-manta: {len(manifest)} structure(s) written to {out}/  (manifest.json)")
    if len(isomers) >= cap:
        print(
            f"delfin-manta: WARNING: output reached the cap of {cap} isomers and may "
            f"be INCOMPLETE — raise --max-isomers for the full enumerated set.",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
