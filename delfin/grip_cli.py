#!/usr/bin/env python3
"""delfin-grip — Geometric Refinement by Internal Priors (CLI).

The standalone command-line frontend for the GRIP method.  Takes an XYZ
structure (organic, organometallic, transition-metal complex, ...) and
produces a CCDC-grounded, force-field-free, geometrically-polished XYZ.

Usage examples::

    # Basic post-xTB refinement
    delfin-grip input.xyz -o output.xyz

    # Explicit metal + donor specification (skips the heuristic)
    delfin-grip input.xyz -o output.xyz --metal 0 --donors 1,2,3,4

    # Use a SMILES hint to validate connectivity and ensure correct topology
    delfin-grip input.xyz -o output.xyz --smiles "[Pd](Cl)(Cl)(N)N"

    # Custom library + tighter convergence
    delfin-grip input.xyz -o output.xyz --library grip_lib_v2.npz --max-iter 500

    # Ensemble mode: enumerate Pólya isomers, return top-k
    delfin-grip input.xyz -o output.xyz --ensemble --k 10 --smiles "[Pd](Cl)(Cl)(N)N"

    # Refine every frame of a trajectory
    delfin-grip trajectory.xyz -o refined.xyz --all-frames

    # Skip the optional post-corrector
    delfin-grip input.xyz -o output.xyz --no-corrector

    # Show library / version info
    delfin-grip --info

    # Validate-only (no writing)
    delfin-grip --validate input.xyz

Design contract
---------------
* Default-OFF byte-identical to HEAD — this is a NEW file; existing GRIP
  callers are not affected.
* Deterministic — ``PYTHONHASHSEED=0`` is set at startup; no RNG inside
  the CLI itself.
* Universal — graph-only, no SMILES-specific gating.  Works on organics,
  TMCs, hapto-π systems, multi-metal clusters with explicit ``--metal``.
* Returns ``exit 0`` on a clean polish, ``2`` on validation failure
  (severity-non-improvement / constraint violation), ``1`` on hard errors
  (bad XYZ, missing library, RDKit failure, ...).
"""
from __future__ import annotations

import argparse
import logging
import os
# Deterministic seed BEFORE any numpy / RDKit import.
os.environ.setdefault("PYTHONHASHSEED", "0")

import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

from delfin import __version__ as _DELFIN_VERSION
from delfin.grip_io import (
    Frame,
    detect_metal_and_donors,
    match_smiles_to_xyz,
    perceive_bonds_from_xyz,
    read_xyz,
    write_xyz,
)

_LOG = logging.getLogger("delfin-grip")

__all__ = [
    "main",
    "build_parser",
    "refine_frame",
    "GripCliResult",
]


# ---------------------------------------------------------------------------
# Exit codes
# ---------------------------------------------------------------------------
EXIT_OK = 0
EXIT_HARD_ERROR = 1
EXIT_VALIDATION_FAIL = 2


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="delfin-grip",
        description=(
            "GRIP — CCDC-grounded, force-field-free geometric refinement. "
            "Refines XYZ structures via Mahalanobis projection onto the "
            "Mogul fragment manifold under hard coordination constraints."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  delfin-grip input.xyz -o output.xyz\n"
            "  delfin-grip in.xyz -o out.xyz --metal 0 --donors 1,2,3,4\n"
            "  delfin-grip in.xyz -o out.xyz --smiles '[Pd](Cl)(Cl)(N)N'\n"
            "  delfin-grip --info\n"
            "  delfin-grip --validate input.xyz\n"
        ),
    )
    p.add_argument(
        "input_xyz",
        nargs="?",
        type=str,
        help="Input XYZ file (single frame or multi-frame trajectory).",
    )
    p.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Output XYZ path.  Default: input.grip.xyz next to the input.",
    )
    p.add_argument(
        "--metal",
        type=int,
        default=None,
        help=(
            "Atom index (0-based) of the metal centre.  If omitted, the "
            "CLI auto-detects the first transition-metal atom."
        ),
    )
    p.add_argument(
        "--donors",
        type=str,
        default=None,
        help=(
            "Comma-separated list of donor-atom indices (0-based).  If "
            "omitted, donors are auto-detected from covalent-radius proximity."
        ),
    )
    p.add_argument(
        "--smiles",
        type=str,
        default=None,
        help=(
            "Optional SMILES hint.  When provided, the CLI cross-checks "
            "the XYZ connectivity and uses the SMILES atom ordering as a "
            "tiebreaker for ambiguous mappings (see --strict-mapping)."
        ),
    )
    p.add_argument(
        "--charge",
        type=int,
        default=0,
        help="Net molecular charge for bond perception (default: 0).",
    )
    p.add_argument(
        "--library",
        type=str,
        default=None,
        help=(
            "Path to the GRIP Mogul library (.npz).  Defaults to the "
            "release-pinned grip_lib_v1.npz (or $DELFIN_GRIP_LIB_PATH "
            "if set)."
        ),
    )
    p.add_argument(
        "--max-iter",
        type=int,
        default=200,
        help="L-BFGS max iterations (default: 200).",
    )
    p.add_argument(
        "--md-tol",
        type=float,
        default=0.05,
        help="M-D invariant half-width in Ångström (default: 0.05).",
    )
    p.add_argument(
        "--clash-weight",
        type=float,
        default=None,
        help=(
            "Pauli-floor clash penalty weight.  Default: env "
            "DELFIN_FFFREE_GRIP_CLASH_WEIGHT or 5.0."
        ),
    )
    p.add_argument(
        "--ensemble",
        action="store_true",
        help=(
            "Enumerate Pólya × Cremer-Pople candidates and rank.  Requires "
            "--smiles for the SMILES-based decomposition."
        ),
    )
    p.add_argument(
        "--k", "--top-k",
        dest="top_k",
        type=int,
        default=3,
        help="Top-K candidates to emit in --ensemble mode (default: 3).",
    )
    p.add_argument(
        "--ensemble-full",
        action="store_true",
        help=(
            "Emit the full ranked ensemble (overrides --k).  Use with care "
            "for combinatorial-heavy SMILES."
        ),
    )
    p.add_argument(
        "--all-frames",
        action="store_true",
        help=(
            "Refine every frame of a multi-frame XYZ.  Default: refine only "
            "the FIRST frame."
        ),
    )
    p.add_argument(
        "--no-corrector",
        action="store_true",
        help=(
            "Skip the post-GRIP corrector (sp2-flatten, h-axis rotation, "
            "M-D-short repulsion).  Default: corrector ON."
        ),
    )
    p.add_argument(
        "--skip-hapto",
        action="store_true",
        help=(
            "Skip the hapto-π protection (debug/forensic only).  Default: "
            "hapto skip ON (matches HEAD behaviour)."
        ),
    )
    p.add_argument(
        "--strict-mapping",
        action="store_true",
        help=(
            "Reject mappings with any element mismatch.  Default: warn-and-"
            "continue."
        ),
    )
    p.add_argument(
        "--info",
        action="store_true",
        help="Print version + library coverage stats and exit.",
    )
    p.add_argument(
        "--validate",
        action="store_true",
        help=(
            "Run GRIP but only REPORT severity-before/-after and constraint "
            "status; do not write an output file."
        ),
    )
    p.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress non-error output.",
    )
    p.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    p.add_argument(
        "--version",
        action="version",
        version=f"delfin-grip {_DELFIN_VERSION}",
    )
    return p


# ---------------------------------------------------------------------------
# Refinement result container
# ---------------------------------------------------------------------------
class GripCliResult:
    """Holds the per-frame refinement diagnostics for the CLI summary."""

    def __init__(
        self,
        frame_index: int,
        accepted: bool,
        severity_before: float,
        severity_after: float,
        n_iter: int,
        n_terms: int,
        rollback_reason: str = "",
        corrector_diag: Optional[dict] = None,
    ):
        self.frame_index = frame_index
        self.accepted = accepted
        self.severity_before = severity_before
        self.severity_after = severity_after
        self.n_iter = n_iter
        self.n_terms = n_terms
        self.rollback_reason = rollback_reason
        self.corrector_diag = corrector_diag or {}

    def summary(self) -> str:
        verdict = "ACCEPTED" if self.accepted else "ROLLBACK"
        delta = self.severity_before - self.severity_after
        rel = ""
        if (
            np.isfinite(self.severity_before)
            and self.severity_before > 1e-12
            and np.isfinite(delta)
        ):
            rel = f" ({100.0 * delta / self.severity_before:+.1f}%)"
        lines = [
            f"frame {self.frame_index}: {verdict}",
            f"  severity before: {self.severity_before:.4f}",
            f"  severity after:  {self.severity_after:.4f}{rel}",
            f"  L-BFGS iter: {self.n_iter} / fragments: {self.n_terms}",
        ]
        if self.rollback_reason:
            lines.append(f"  rollback reason: {self.rollback_reason}")
        if self.corrector_diag and self.corrector_diag.get("any_enabled"):
            cd = self.corrector_diag
            lines.append(
                f"  post-corrector: flatten={cd.get('flatten_n', 0)} "
                f"haxis={cd.get('haxis_n', 0)} mdshort={cd.get('mdshort_n', 0)}"
            )
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Single-frame refinement
# ---------------------------------------------------------------------------
def refine_frame(
    frame: Frame,
    *,
    frame_index: int = 0,
    metal_arg: Optional[int] = None,
    donors_arg: Optional[Sequence[int]] = None,
    smiles_hint: Optional[str] = None,
    charge: int = 0,
    library_path: Optional[str] = None,
    max_iter: int = 200,
    md_tol: float = 0.05,
    clash_weight: Optional[float] = None,
    run_corrector: bool = True,
    skip_hapto: bool = True,
    strict_mapping: bool = False,
) -> Tuple[Frame, GripCliResult]:
    """Refine a single :class:`Frame` and return the polished frame + diagnostics.

    The function:

    1. Perceives bonds from the XYZ via :func:`perceive_bonds_from_xyz`.
    2. (Optional) validates the perceived graph against ``smiles_hint``.
    3. Resolves metal + donor indices (explicit args > auto-detect).
    4. Calls :func:`grip_polish` with ``return_diagnostics=True``.
    5. (Optional) runs :func:`post_grip_corrections`.

    On any HARD failure (missing metal, RDKit error) raises ``RuntimeError``.
    A clean accept-if-better rejection is NOT an error — the input frame is
    returned with ``accepted=False`` in the result.
    """
    # Honour the library override env var so child imports see it.
    if library_path:
        os.environ["DELFIN_GRIP_LIB_PATH"] = str(library_path)

    # Honour the hapto skip env var for the polish.
    if skip_hapto:
        os.environ.setdefault("DELFIN_FFFREE_GRIP_SKIP_HAPTO", "1")
    else:
        os.environ["DELFIN_FFFREE_GRIP_SKIP_HAPTO"] = "0"

    # --- 1. Bond perception ---------------------------------------------
    try:
        xyz_mol = perceive_bonds_from_xyz(
            frame.elements, frame.coordinates, charge=charge,
        )
    except Exception as exc:
        raise RuntimeError(f"bond perception failed: {exc!r}") from exc

    # --- 2. Optional SMILES validation -----------------------------------
    if smiles_hint:
        mapping = match_smiles_to_xyz(
            smiles_hint, xyz_mol, frame.coordinates,
            strict_element=strict_mapping,
        )
        if mapping is None and strict_mapping:
            raise RuntimeError(
                "SMILES->XYZ mapping failed under --strict-mapping; "
                "check that --smiles atom count + elements match the XYZ "
                "(remember to write metal-donor bonds explicitly, e.g. "
                "'[Pd](Cl)(Cl)(N)N' not 'Pd Cl Cl N N')."
            )
        if mapping is None:
            _LOG.warning(
                "SMILES->XYZ mapping failed; proceeding with perceived "
                "connectivity only."
            )

    # --- 3. Resolve metal + donors --------------------------------------
    if metal_arg is not None:
        metal_idx = int(metal_arg)
        if not (0 <= metal_idx < xyz_mol.GetNumAtoms()):
            raise RuntimeError(
                f"--metal {metal_idx} out of range [0, {xyz_mol.GetNumAtoms()})"
            )
        if donors_arg is not None:
            donors = [int(d) for d in donors_arg]
            if not donors:
                # Empty list explicit -> fall back to auto-detect
                _, donors = detect_metal_and_donors(xyz_mol, frame.coordinates)
        else:
            _, donors = detect_metal_and_donors(xyz_mol, frame.coordinates)
    else:
        metal_idx, donors = detect_metal_and_donors(xyz_mol, frame.coordinates)
        if donors_arg is not None:
            donors = [int(d) for d in donors_arg]

    if metal_idx is None:
        # Organic-only input: GRIP still works, but with metal=donors=empty
        # (no constraint stack on coordination, no polyhedron).  Use sentinel
        # values that grip_polish accepts gracefully.
        _LOG.info(
            "no transition-metal detected; running GRIP in organic-only "
            "mode (no M-D constraints)."
        )
        metal_idx = 0
        donors = []

    # --- 4. GRIP polish ---------------------------------------------------
    try:
        from delfin.fffree.grip_polish import grip_polish, GripPolishResult
    except Exception as exc:
        raise RuntimeError(f"cannot import grip_polish: {exc!r}") from exc

    try:
        result = grip_polish(
            frame.coordinates,
            xyz_mol,
            metal=int(metal_idx),
            donors=list(donors),
            geom="",   # geometry tag unused in standalone CLI (no decompose)
            max_iter=int(max_iter),
            md_tol=float(md_tol),
            clash_weight=clash_weight,
            return_diagnostics=True,
        )
    except Exception as exc:
        raise RuntimeError(f"grip_polish raised: {exc!r}") from exc

    if isinstance(result, GripPolishResult):
        P_refined = result.P
        accepted = bool(result.accepted)
        sev_before = float(result.severity_before)
        sev_after = float(result.severity_after)
        n_iter = int(result.n_iter)
        n_terms = int(result.n_terms)
        rollback_reason = str(result.rollback_reason)
    else:  # defensive — older API returned an ndarray
        P_refined = np.asarray(result, dtype=np.float64)
        accepted = not np.allclose(P_refined, frame.coordinates, atol=1e-9)
        sev_before = sev_after = float("nan")
        n_iter = 0
        n_terms = 0
        rollback_reason = ""

    if P_refined is None:
        P_refined = frame.coordinates.copy()

    # --- 5. Post-GRIP corrector (optional) -------------------------------
    corrector_diag: dict = {}
    if run_corrector and donors:
        try:
            from delfin.fffree.post_grip_corrector import post_grip_corrections
            # The post-corrector reads its own env-flags to decide which
            # sub-correctors fire.  We do not unconditionally enable them —
            # only the user's process env governs them.
            P_corr, corrector_diag = post_grip_corrections(
                frame.elements, P_refined,
                metal_idx=int(metal_idx),
                donor_idxs=list(donors),
            )
            if P_corr is not None and np.all(np.isfinite(P_corr)):
                P_refined = np.asarray(P_corr, dtype=np.float64)
        except Exception as exc:
            _LOG.debug("post_grip_corrector failed: %r (ignored)", exc)

    refined_frame = Frame(
        elements=list(frame.elements),
        coordinates=np.asarray(P_refined, dtype=np.float64),
        comment=(frame.comment + " [GRIP refined]").strip(),
    )
    cli_result = GripCliResult(
        frame_index=frame_index,
        accepted=accepted,
        severity_before=sev_before,
        severity_after=sev_after,
        n_iter=n_iter,
        n_terms=n_terms,
        rollback_reason=rollback_reason,
        corrector_diag=corrector_diag,
    )
    return refined_frame, cli_result


# ---------------------------------------------------------------------------
# Ensemble mode
# ---------------------------------------------------------------------------
def _refine_ensemble(
    smiles: str,
    *,
    top_k: int = 3,
    emit_full: bool = False,
    library_path: Optional[str] = None,
    skip_hapto: bool = True,
) -> List[Tuple[Frame, dict]]:
    """Run the GRIP ensemble enumeration on ``smiles`` and return
    refined frames + diagnostics for each top-ranked candidate.
    """
    if library_path:
        os.environ["DELFIN_GRIP_LIB_PATH"] = str(library_path)
    if skip_hapto:
        os.environ.setdefault("DELFIN_FFFREE_GRIP_SKIP_HAPTO", "1")
    # Ensure the ensemble module's env flag is on for grip_ensemble_enumerate's
    # internal hooks (it self-gates on this).
    os.environ["DELFIN_FFFREE_GRIP_ENSEMBLE"] = "1"

    try:
        from delfin.fffree.grip_ensemble import grip_ensemble_enumerate
    except Exception as exc:
        raise RuntimeError(f"cannot import grip_ensemble: {exc!r}") from exc

    result = grip_ensemble_enumerate(
        smiles,
        top_k=top_k,
    )
    if result.skip_reason:
        raise RuntimeError(
            f"ensemble skipped: {result.skip_reason}"
        )

    out: List[Tuple[Frame, dict]] = []
    candidates = result.candidates if emit_full else result.top_k
    for c in candidates:
        frame = Frame(
            elements=list(c.syms),
            coordinates=np.asarray(c.P, dtype=np.float64),
            comment=(
                f"GRIP-Ensemble {c.label} isomer={c.isomer_id} "
                f"conformer={c.conformer_id} severity={c.severity:.3f} "
                f"clash={c.clash_count} cshm={c.cshm:.3f} score={c.score:.3f}"
            ),
        )
        out.append((
            frame,
            {
                "label": c.label,
                "isomer_id": c.isomer_id,
                "conformer_id": c.conformer_id,
                "severity": c.severity,
                "clash_count": c.clash_count,
                "cshm": c.cshm,
                "score": c.score,
                "accepted": c.accepted,
            },
        ))
    return out


# ---------------------------------------------------------------------------
# --info handler
# ---------------------------------------------------------------------------
def _print_info(library_path: Optional[str] = None) -> int:
    """Print version + library statistics.  Returns shell exit code."""
    print(f"delfin-grip — DELFIN v{_DELFIN_VERSION}")
    try:
        from delfin.fffree.grip_mogul_lookup import (
            DEFAULT_LIB_PATH,
            DELFIN_GRIP_LIB_PATH_ENV,
            GripLibrary,
        )
        if library_path:
            os.environ[DELFIN_GRIP_LIB_PATH_ENV] = str(library_path)
        path_env = os.environ.get(DELFIN_GRIP_LIB_PATH_ENV, "").strip()
        resolved = path_env or str(DEFAULT_LIB_PATH)
        print(f"library path: {resolved}")
        try:
            lib = GripLibrary.get(Path(resolved))
            print(f"library version: v{lib.version}")
            print(f"library size:")
            print(f"  master entries: {lib.n_master:,}")
            print(f"  original (no-fallback) entries: {lib.n_orig:,}")
        except FileNotFoundError:
            print(f"library NOT FOUND at {resolved}")
            print(
                "  (delfin-grip will still run, but with an empty fragment "
                "library — i.e. a no-op polish).  Build via "
                "scripts/grip_build_mogul_lib.py."
            )
            return EXIT_HARD_ERROR
    except Exception as exc:
        print(f"library inspection failed: {exc!r}")
        return EXIT_HARD_ERROR
    return EXIT_OK


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point — returns shell exit code (0/1/2)."""
    parser = build_parser()
    args = parser.parse_args(argv)

    # --- Logging level --------------------------------------------------
    if args.quiet:
        level = logging.ERROR
    elif args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(
        level=level, format="%(levelname)s [%(name)s] %(message)s"
    )

    # --- --info short-circuit ------------------------------------------
    if args.info:
        return _print_info(library_path=args.library)

    if not args.input_xyz:
        parser.error("input_xyz is required (or use --info / --validate).")

    in_path = Path(args.input_xyz)
    if not in_path.exists():
        _LOG.error("input file not found: %s", in_path)
        return EXIT_HARD_ERROR

    # --- Output path -------------------------------------------------------
    if args.validate:
        out_path: Optional[Path] = None
    elif args.output:
        out_path = Path(args.output)
    else:
        out_path = in_path.with_suffix("")
        out_path = out_path.with_name(out_path.name + ".grip.xyz")

    # --- Parse donors argument ----------------------------------------
    donors_arg: Optional[List[int]] = None
    if args.donors is not None:
        try:
            donors_arg = [
                int(x.strip()) for x in args.donors.split(",") if x.strip()
            ]
        except ValueError as exc:
            _LOG.error("--donors parse error: %s", exc)
            return EXIT_HARD_ERROR

    # --- ENSEMBLE mode --------------------------------------------------
    if args.ensemble:
        if not args.smiles:
            _LOG.error(
                "--ensemble requires --smiles for SMILES-based Pólya "
                "decomposition.  Aborting."
            )
            return EXIT_HARD_ERROR
        try:
            ensemble_frames = _refine_ensemble(
                args.smiles,
                top_k=int(args.top_k),
                emit_full=bool(args.ensemble_full),
                library_path=args.library,
                skip_hapto=not args.skip_hapto,
            )
        except RuntimeError as exc:
            _LOG.error("ensemble failed: %s", exc)
            return EXIT_HARD_ERROR
        if not ensemble_frames:
            _LOG.error("ensemble produced no valid candidates.")
            return EXIT_VALIDATION_FAIL
        if args.validate:
            _LOG.info("--validate: %d candidates generated:", len(ensemble_frames))
            for i, (_, diag) in enumerate(ensemble_frames):
                _LOG.info(
                    "  rank %d: severity=%.3f clash=%d cshm=%.3f score=%.3f",
                    i + 1, diag["severity"], diag["clash_count"],
                    diag["cshm"], diag["score"],
                )
            return EXIT_OK
        try:
            write_xyz(
                out_path, [fr for fr, _ in ensemble_frames],
                comment_prefix="",  # comments already carry per-candidate detail
            )
        except Exception as exc:
            _LOG.error("write failed: %s", exc)
            return EXIT_HARD_ERROR
        _LOG.info(
            "wrote %d ensemble candidates -> %s",
            len(ensemble_frames), out_path,
        )
        return EXIT_OK

    # --- Single / trajectory mode --------------------------------------
    try:
        frames = read_xyz(in_path)
    except Exception as exc:
        _LOG.error("read_xyz failed: %s", exc)
        return EXIT_HARD_ERROR

    if not args.all_frames:
        frames = frames[:1]

    refined: List[Frame] = []
    results: List[GripCliResult] = []
    any_accepted = False
    for i, frame in enumerate(frames):
        try:
            rf, res = refine_frame(
                frame,
                frame_index=i,
                metal_arg=args.metal,
                donors_arg=donors_arg,
                smiles_hint=args.smiles,
                charge=int(args.charge),
                library_path=args.library,
                max_iter=int(args.max_iter),
                md_tol=float(args.md_tol),
                clash_weight=args.clash_weight,
                run_corrector=not args.no_corrector,
                skip_hapto=not args.skip_hapto,
                strict_mapping=bool(args.strict_mapping),
            )
        except RuntimeError as exc:
            _LOG.error("frame %d: %s", i, exc)
            return EXIT_HARD_ERROR
        results.append(res)
        refined.append(rf)
        if res.accepted:
            any_accepted = True
        if not args.quiet:
            print(res.summary())

    if args.validate:
        # --validate: no output file, exit code reflects the verdict.
        return EXIT_OK if any_accepted else EXIT_VALIDATION_FAIL

    try:
        write_xyz(out_path, refined)
    except Exception as exc:
        _LOG.error("write failed: %s", exc)
        return EXIT_HARD_ERROR
    if not args.quiet:
        print(f"wrote {len(refined)} frame(s) -> {out_path}")

    return EXIT_OK if any_accepted else EXIT_VALIDATION_FAIL


if __name__ == "__main__":
    sys.exit(main())
