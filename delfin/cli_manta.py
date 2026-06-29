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

# Construction configs (env-flag sets). 'champion' = the full SHIP-31 rich
# construction (maximum coverage + conformer/motif handling); 'builder' = lean
# core + reach; 'default' = library/legacy. champion is the MANTA default.
_CHAMPION_FLAGS = (
    "AROM_PLANARIZE", "ARYL_RING_SIZE", "BACKBONE_REEMBED", "CHELATE_BACKBONE",
    "CN_EXTEND", "CONF_COMPLETE", "CONFORMER_COVERAGE", "CONFORMER_SEATING",
    "DIATOMIC_ORIENT", "DONOR_BEND", "HAPTO_AXIS_ROT", "HAPTO_HALFSANDWICH_GATE",
    "INTERLIG_RANK", "INTERLIG_VDW_GATE", "JOINT_DECLASH", "LIGAND_RIGID",
    "MD_CONTEXT", "METALLOID_DONOR", "MULTIBOND_EXEMPT", "NHC_CARBENE",
    "PI_COPLANAR_M", "PLANAR_MER", "RIGID_HAPTO", "SIGMA_ENSEMBLE",
    "TOPOLOGY_GATE", "TORSION_RELAX", "XH_COLLAPSE", "KAPPA4", "CONF_ENERGY_RANK",
)
_BUILDER_FLAGS = ("KAPPA4", "SIGMA_ENSEMBLE", "CONF_ENERGY_RANK")


def _apply_construction_env(config: str) -> None:
    """Set the DELFIN_FFFREE_* construction env for the chosen config (before import)."""
    if config == "default":
        return
    os.environ["DELFIN_FFFREE_BUILDER"] = "1"
    os.environ["DELFIN_FRAME_RANK_FIX"] = "1"
    flags = _CHAMPION_FLAGS if config == "champion" else _BUILDER_FLAGS
    for f in flags:
        os.environ["DELFIN_FFFREE_" + f] = "1"


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
    p.add_argument("--quality", choices=["fast", "normal", "max", "extreme"], default="extreme",
                   help="conformer-depth preset (seeds x templates x cap): "
                        "fast(12 seeds) | normal(20) | max(40) | extreme(60, DEFAULT). "
                        "Convergence study: fast/normal MISS the GFN2 global minimum by "
                        "~2.5 kcal/mol on multi-isomer systems; extreme captures it for all "
                        "tested -> extreme is the default for guaranteed global-min coverage. "
                        "Lower it only for quick previews.")
    p.add_argument("--seeds", type=int, default=None, dest="seeds",
                   help="override the ETKDG conformer-seed count directly (the key "
                        "completeness/speed switch; overrides --quality's seed count). "
                        "Higher = more conformers/rotamers, slower.")
    p.add_argument("--num-confs", type=int, default=None, dest="num_confs",
                   help="conformers embedded per isomer (default: 200)")
    p.add_argument("--construction", choices=["champion", "builder", "default"],
                   default="champion",
                   help="construction config: champion (full SHIP-31 rich, DEFAULT) | "
                        "builder (lean core+reach) | default (library/legacy)")
    p.add_argument("--collapse-variants", dest="collapse", action="store_true", default=False,
                   help="merge label-variant isomers (default: OFF — keep every variant = "
                        "maximum richness)")
    p.add_argument("--no-binding-modes", dest="binding_modes", action="store_false", default=True,
                   help="disable alternative binding-mode isomers (default: ON)")
    p.add_argument("--hapto-approx", choices=["auto", "on", "off"], default="auto",
                   help="hapto (eta) approximation: auto (default) | on | off")
    p.add_argument("--no-uff", dest="apply_uff", action="store_false", default=True,
                   help="skip the UFF cleanup of emitted geometries (default: UFF on)")
    p.add_argument("--no-deterministic", dest="deterministic", action="store_false", default=True,
                   help="allow non-deterministic embedding (default: deterministic)")
    p.add_argument("--charge", type=int, default=None,
                   help="override the metal/complex charge (default: read from the SMILES)")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="suppress the per-structure listing")
    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)

    # Construction config + ranking are env-gated; set ALL switches BEFORE import.
    _apply_construction_env(args.construction)
    if args.rank:
        os.environ["DELFIN_FFFREE_GFNFF_RANK"] = "1"
        os.environ["DELFIN_CONF_RANK_METHOD"] = args.method
    if args.charge is not None:
        os.environ["DELFIN_GFNFF_CHARGE"] = str(int(args.charge))

    from delfin.smiles_converter import smiles_to_xyz_isomers

    cap = args.max_isomers if args.max_isomers is not None else _ALL_ISOMERS
    kwargs = {
        "max_isomers": cap,
        "collapse_label_variants": bool(args.collapse),
        "include_binding_mode_isomers": bool(args.binding_modes),
        "apply_uff": bool(args.apply_uff),
        "deterministic": bool(args.deterministic),
    }
    if args.quality is not None:
        kwargs["quality_mode"] = args.quality
    if args.seeds is not None:
        kwargs["seeds_override"] = args.seeds
    if args.num_confs is not None:
        kwargs["num_confs"] = args.num_confs
    if args.hapto_approx != "auto":
        kwargs["hapto_approx"] = (args.hapto_approx == "on")

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
