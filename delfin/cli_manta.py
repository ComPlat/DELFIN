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

# Construction configs (env-flag sets). 'champion' = the DE-BLOATED best config; 'builder' = lean
# core + reach; 'default' = library/legacy. champion is the MANTA default.
#
# 2026-07-04 DE-BLOAT: the old ~29-flag champion stack (dense conformer generators + geometry
# correctors + chelate cap-raisers) was NET-NEGATIVE. Measured on ~450 region-balanced Batch.txt
# systems via the canonical whole-manifold eye: the full old stack scored 13.7% topology-correct
# vs 33.6% for zero-flags legacy, WORSE in every donor class (P/S 4x, NO/O 3x, ...). Removing the
# 20 proven-negative flags -> 2.1x, never-worse in every class, while KEEPING the edge-class
# builders (hapto/carbonyl/NHC/metalloid/kappa4/arom) so under-sampled edge chemistry cannot
# regress, plus main's own never-worse CAGE_MD_GUARD + CN4_BOTH additions. Removed flags stay in
# git history + remain individually env-gated (DELFIN_FFFREE_<flag>). Spot-verified: KITNEJ/QIDGEP/
# VULMOE (50-junk topo-wrong manifolds under the old stack) now build the crystal topology.
_CHAMPION_FLAGS = (
    "DET_CLASSIFY",       # DETERMINISM: classify conformers without sharing one RDKit mol across
                          # threads (its ring/property caches init lazily and are not thread-safe, so
                          # the SCORES raced -> another conformer won its fingerprint -> same label,
                          # different geometry; and the FILTERS raced -> a different frame COUNT).
                          # landed 2026-07-10, whole-pool byte-determinism, full-fidelity, all pillars:
                          #   main            108 systems, byte_identical 101, nondeterministic 7
                          #   main+DET_CLASSIFY                        108,                    0
                          # capability never-worse: valid 85->85, cap_LOST=0, build_lost=0, rt 0 lost
    "AROM_PLANARIZE", "ARYL_RING_SIZE", "DIATOMIC_ORIENT", "HAPTO_AXIS_ROT",
    "HAPTO_HALFSANDWICH_GATE", "METALLOID_DONOR", "NHC_CARBENE", "RIGID_HAPTO", "KAPPA4",
    "METALLOID_MD_LEN",   # correct M-metalloid bond length (no row-offset overshoot); landed 2026-07-09:
                          # heavy_donor cap-never-worse (+8 valid frames, 0 lost), byte-id off metalloids
    "JOINT_DECLASH", "DECLASH_METALLOID_LIGAND",  # general inter-ligand azimuthal declash (hard M-D
                          # invariant + rollback) + metalloid-ligand awareness; landed 2026-07-09:
                          # full:120 whole-space cap-never-worse, +3 valid frames, mean -0.3
    "SP3C_TET_SEAT",      # tetrahedral seating of PENDANT sp3-alkyl C donors (M-C-X ~114 deg, not 180);
                          # landed 2026-07-09: class:sigma_coord:100 cap-never-worse (n=94, 0 lost,
                          # mean -2.009) + full:120 whole-space (n=105, 0 lost, +1 perfect manifold).
                          # Scoped by element+graph only: not aromatic, not in-ring (sp2 carbene/aryl
                          # have no vacant tetrahedral slot; ring C orientation is fixed by the scaffold)
    "CAGE_MD_GUARD",   # main's user-approved net-positive cage/over-coord M-D guard (kept)
    "CN4_BOTH",        # main's native-additive CN4 dual-geometry completeness, never-worse (kept)
    "METALLOID_MD_CLAMP", # post-UFF re-snap of M-metalloid bonds (Sb/As/Bi/Te/Se/Ge/Sn/Pb) to ideal;
                          # OB-UFF has no params for these soft donors and collapses them (QIGFOF Ag-Sb
                          # 2.01 vs 2.65).  Completes METALLOID_MD_LEN (which sets the target; this holds
                          # it through UFF).  landed 2026-07-11 under the NOISE-BAND gate: champion+CLAMP
                          # vs champion never-worse BOTH directions, mean -0.204, 0 regressions, whole-
                          # pool byte-determinism 109/109 (TOCSAI needed per_timeout 1800, then identical)
    "STEREOCENTER_ENUM",  # #18: coordination-created STEREOCENTRE completeness -- additively builds every
                          # buildable +/- fold of a pnictogen (N/P/As/Sb/Bi) X-H or X-R stereocentre (the
                          # metal locks the pyramidal inversion), so the crystal's N-H fold is no longer
                          # missing.  Eye: ccdc_isomer_realized false->true (USEMOW [..N+N-N+N-], OFAJOV).
                          # landed 2026-07-14, complete eye+battery gate (NO mogul, per user), full:120
                          # never-worse: cap/build/poly/isomers/roundtrip/realism/battery/good_regr ALL 0,
                          # 0 ccdc-breaks, mean_delta ~0, whole-pool byte-determinism 108/108 (TOCSAI a
                          # no-op timeout).  Additive+deterministic; scoped to group-15 (O/C are NOT stable
                          # centres and broke the ccdc floor when included -> restricted).  Impl:
                          # delfin/manta/_stereocenter_enum.py; env DELFIN_(FFFREE_)STEREOCENTER_ENUM.
    "D8_SQ_ADD",          # #19: ADDITIVE d8 CN4 square-planar completeness.  A UFF-SP-4 octahedron... no:
                          # a UFF SQUARE frame is added as a PURELY ADDITIVE sibling (force_d8_sq) next to
                          # the normal build; the PRIMARY frame is untouched, so a bulky ligand keeps its
                          # valid tetrahedral primary (NO broken_regressed -- the earlier REPLACE-the-frame
                          # SP-4 cost HAKQES a frame) while a d8-no-valid system GAINS its valid SP-4 frame.
                          # landed 2026-07-14, full:120 never-worse: cap_LOST=0, good_regressions=0,
                          # broken_regressed=0, no ccdc_isomer_lost / poly_cshm_regressed, mean_delta -0.116
                          # (slightly BETTER); rescues the d8-no-valid tail (+2 on pool_d8novalid,
                          # ITANUN/XETKAI 0->topo).  Byte-deterministic (29/29).  Only blockers were the
                          # TOCSAI build-TIMEOUT (additive is a bit slower -> raise per-timeout) and an
                          # axis-record artefact.  Env DELFIN_FFFREE_D8_SQ_ADD; impl smiles_converter batch
                          # loop.  CN4 SQ has ONE polyhedron so no OC/TPR isomer ambiguity (unlike CN6).
    "CN6_OH_ADD",         # #20: ADDITIVE CN6 octahedron completeness -- adds the OC-6 octahedron isomer as a
                          # PURELY ADDITIVE sibling for CN6 systems missing it (AVEDOW iso_theory=2: coverage
                          # 50%->100%, builds the theory-expected 2nd isomer).  landed 2026-07-16, smart-1000
                          # (full:1000 --ab): 243 affected + 728 byte-identical, cap_LOST=0, cap_gained=4
                          # (ALUDAO/APIKER/QUHWAT/TADYIJ), good_regressions=0, poly_lost=0, ccdc_isomer_lost=0,
                          # mean_delta -0.402 (BETTER); never_worse_ok=False ONLY because REALISM+RT not
                          # measured (--no-mogul), topology_floor=True.  ENABLED BY two two-sided-verified
                          # eye/gate sharpenings the user's per-frame push surfaced (agent_workspace): (1) the
                          # polyhedron floor now reads the BEST (min CShM to crystal shape) over ANY realistic
                          # frame, not the single min-RMSD best_valid frame -- RMSD-selection had hijacked
                          # AVEDOW/MOYDOV onto the added-isomer sibling -> false poly_lost on a completeness
                          # WIN; firing-side kept (16 systems still poly_cshm>2, max 33).  (2) ccdc_isomer_lost
                          # valid-baseline-scoped like ccdc_pucker_lost (SUPTON: NO valid frame either way ->
                          # its 'realisation' was on a broken frame; firing-side kept -- FONKOL still fires).
                          # Env DELFIN_FFFREE_CN6_OH_ADD; impl smiles_converter energy-outlier-cut root (keep
                          # inf-energy topo-valid isomer primaries; VOYWUD 6->5->6).
    "AROM_SEAT",          # #21: construction-time FREE-aromatic bond-length seat toward the CCDC delocalised
                          # targets (C-C 1.387 / C-N 1.349 / C-O 1.369 / C-S 1.732 / N-N 1.364 = the eye's
                          # org_prior "ar" mu).  refine() aromatic-aware (coordinated ring SYSTEMS EXCLUDED ->
                          # coordination byte-identical, never-worse by construction) + UFF AddDistanceConstraint
                          # on the metal-free path + bond_decollapse target.  landed 2026-07-17 (night):
                          # class:sigma_coord:100 never-worse (n=3 affected, cap_LOST=0, mean -13.5) + full:120
                          # never-worse (n=3, within band).  Reach = FREE aromatics only; coordinated aromatics
                          # deferred to embed-time seating (AROM_EMBED_SEAT, separate).  Env DELFIN_FFFREE_AROM_SEAT.
)
_BUILDER_FLAGS = ("KAPPA4", "SIGMA_ENSEMBLE", "CONF_ENERGY_RANK")


def _apply_construction_env(config: str) -> None:
    """Set the DELFIN_FFFREE_* construction env for the chosen config (before import)."""
    if config == "default":
        return
    os.environ["DELFIN_FFFREE_BUILDER"] = "1"
    os.environ["DELFIN_FRAME_RANK_FIX"] = "1"
    os.environ["DELFIN_CHIRAL_ENUM"] = "1"   # Lambda/Delta enantiomer enumeration (>=2 chelate pairs)
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
                   help="ranking/opt Hamiltonian when --rank/--opt is set (default: gfn2)")
    p.add_argument("--opt", type=int, default=None, metavar="N",
                   help="POST-PROCESSING (opt-in): xtb geometry-optimise the top-N emitted "
                        "structures (0 = the whole manifold). Default: off — the manifold is "
                        "emitted with its construction geometry, unchanged. Independent of --rank.")
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


def _geometry_opt_topn(isomers, topn, method, charge):
    """POST-PROCESSING (opt-in, --opt): xtb geometry-optimise the top-N emitted structures
    (0 = all).  ``isomers`` = ``[(xyz, label), ...]``; returns the same shape with optimised
    geometries substituted for the head.  Best-effort: any structure whose optimisation fails
    keeps its construction geometry.  The manifold is UNCHANGED when --opt is omitted."""
    if not isomers or topn is None or int(topn) < 0:
        return isomers
    try:
        from delfin.manta import _gfnff_rank as _gff
    except Exception:
        return isomers
    if not _gff.available():
        print("delfin-manta: --opt requested but xtb was not found on PATH; "
              "keeping construction geometry.", file=sys.stderr)
        return isomers
    import concurrent.futures as _cf
    n = len(isomers) if int(topn) == 0 else min(int(topn), len(isomers))
    head, tail = list(isomers[:n]), list(isomers[n:])

    def _opt_one(item):
        xyz, label = item
        try:
            r = _gff.gfnff_optimize_autospin(xyz, charge=int(charge), method=method)
            new_xyz = r[0] if isinstance(r, tuple) else r
            return (new_xyz or xyz, label)
        except Exception:
            return (xyz, label)
    with _cf.ThreadPoolExecutor(max_workers=max(1, min(n, (os.cpu_count() or 4)))) as ex:
        head = list(ex.map(_opt_one, head))
    return head + tail


def _charge_for_opt(smiles, cli_charge):
    """Charge for --opt: the explicit --charge if given, else RDKit formal charge of the SMILES
    (same derivation the Submit tab uses), else 0."""
    if cli_charge is not None:
        return int(cli_charge)
    try:
        from rdkit import Chem as _Chem
        _m = _Chem.MolFromSmiles(smiles, sanitize=False)
        if _m is not None:
            return int(_Chem.GetFormalCharge(_m))
    except Exception:
        pass
    return 0


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

    # POST-PROCESSING (opt-in): geometry-optimise the top-N.  Omitted -> pure manifold, unchanged.
    if args.opt is not None and args.opt >= 0:
        _chg = _charge_for_opt(args.smiles, args.charge)
        isomers = _geometry_opt_topn(isomers, args.opt, args.method, _chg)

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
