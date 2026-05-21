"""Batch GUPPY sampling across many SMILES.

Reads a CSV or plain-text file with one SMILES per row and runs GUPPY
(`delfin.guppy_sampling.run_sampling`) in isolated subdirectories per entry.
Charge is re-derived from each SMILES — no charge column is honored.

Can be invoked in three ways:

- Sequentially on the full batch (default) — all rows processed in one call.
- Single-row mode via ``--row N`` (1-based) for SLURM job-array dispatch.
- As an imported module via :func:`run_sampling_batch`.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from delfin.common.logging import get_logger
from delfin.guppy_sampling import (
    _ALLOWED_START_STRATEGIES,
    _derive_charge_from_smiles,
    _read_smiles_lines,
    run_sampling,
)

logger = get_logger(__name__)

_SLUG_SANITIZE_RE = re.compile(r"[^A-Za-z0-9._+-]+")


def _sanitize_slug(name: str) -> str:
    cleaned = _SLUG_SANITIZE_RE.sub("_", name).strip("_")
    return cleaned[:64] if cleaned else ""


def _derive_slug(smiles: str, name: Optional[str], index: int) -> str:
    if name:
        slug = _sanitize_slug(name)
        if slug:
            return f"{index:04d}_{slug}"
    digest = hashlib.sha1(smiles.encode("utf-8")).hexdigest()[:10]
    return f"{index:04d}_{digest}"


@dataclass
class BatchEntryResult:
    slug: str
    smiles: str
    resolved_charge: int
    status: str  # "ok" | "failed" | "error"
    best_xtb_eh: Optional[float] = None
    best_goat_eh: Optional[float] = None
    winner_source: Optional[str] = None
    walltime_s: float = 0.0
    workdir: Optional[str] = None
    error: Optional[str] = None


@dataclass
class BatchSharedArgs:
    runs: int = 20
    pal: int = 16
    maxcore: int = 6000
    parallel_jobs: int = 4
    method: str = "XTB2"
    seed: int = 31
    allow_partial: bool = True
    goat_topk: int = 0
    goat_parallel_jobs: Optional[int] = None
    start_strategy: str = "isomers"
    max_isomers: int = 100
    rmsd_cutoff: float = 0.3
    energy_window_kcal: float = 25.0
    output_name: str = "GUPPY_try.xyz"
    workdir_name: str = "GUPPY"
    extra: Dict[str, str] = field(default_factory=dict)


def _write_input_txt(entry_dir: Path, smiles: str) -> Path:
    path = entry_dir / "input.txt"
    path.write_text(f"{smiles}\n", encoding="utf-8")
    return path


def _read_summary_json(entry_dir: Path, shared: BatchSharedArgs) -> Optional[dict]:
    summary_path = entry_dir / shared.workdir_name / "guppy_goat_summary.json"
    if not summary_path.exists():
        return None
    try:
        return json.loads(summary_path.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        logger.warning("Could not parse GUPPY summary %s: %s", summary_path, exc)
        return None


def _best_energy_from_xyz(xyz_path: Path) -> Optional[float]:
    if not xyz_path.exists():
        return None
    try:
        with xyz_path.open("r", encoding="utf-8") as handle:
            # Skip atom count line
            handle.readline()
            comment = handle.readline().strip()
    except Exception:
        return None
    tokens = comment.split()
    for tok in tokens[1:]:
        try:
            return float(tok)
        except ValueError:
            continue
    return None


def _run_one_entry(
    index: int,
    smiles: str,
    name: Optional[str],
    *,
    base_workdir: Path,
    shared: BatchSharedArgs,
) -> BatchEntryResult:
    slug = _derive_slug(smiles, name, index)
    entry_dir = base_workdir / slug
    entry_dir.mkdir(parents=True, exist_ok=True)

    resolved_charge = _derive_charge_from_smiles(smiles)
    logger.info(
        "Batch entry %d: slug=%s SMILES=%s charge=%+d", index, slug, smiles, resolved_charge
    )

    input_file = _write_input_txt(entry_dir, smiles)
    guppy_workdir = entry_dir / shared.workdir_name
    output_file = entry_dir / shared.output_name

    t_start = time.time()
    try:
        rc = run_sampling(
            input_file=input_file,
            runs=shared.runs,
            charge=None,  # always derive from SMILES
            pal=shared.pal,
            maxcore=shared.maxcore,
            parallel_jobs=shared.parallel_jobs,
            method=shared.method,
            output_file=output_file,
            workdir=guppy_workdir,
            seed=shared.seed,
            allow_partial=shared.allow_partial,
            goat_topk=shared.goat_topk,
            goat_parallel_jobs=shared.goat_parallel_jobs,
            start_strategy=shared.start_strategy,
            max_isomers=shared.max_isomers,
            rmsd_cutoff=shared.rmsd_cutoff,
            energy_window_kcal=shared.energy_window_kcal,
        )
    except Exception as exc:  # noqa: BLE001
        elapsed = time.time() - t_start
        logger.exception("Batch entry %d (%s) crashed: %s", index, slug, exc)
        return BatchEntryResult(
            slug=slug, smiles=smiles, resolved_charge=resolved_charge,
            status="error", walltime_s=elapsed, workdir=str(entry_dir), error=str(exc),
        )

    elapsed = time.time() - t_start

    summary = _read_summary_json(entry_dir, shared)
    best_xtb = _best_energy_from_xyz(output_file.with_name("best_coordniation_pre_goat.xyz"))
    best_final = _best_energy_from_xyz(output_file.with_name("best_coordniation.xyz"))
    winner_source = summary.get("winner_source") if summary else ("xtb" if best_final is not None else None)
    best_goat = None
    if summary and summary.get("winner", {}).get("source") == "goat":
        best_goat = best_final

    status = "ok" if rc == 0 else "failed"
    return BatchEntryResult(
        slug=slug,
        smiles=smiles,
        resolved_charge=resolved_charge,
        status=status,
        best_xtb_eh=best_xtb,
        best_goat_eh=best_goat,
        winner_source=winner_source,
        walltime_s=elapsed,
        workdir=str(entry_dir),
    )


def _write_batch_results(base_workdir: Path, results: List[BatchEntryResult]) -> Path:
    path = base_workdir / "batch_results.json"
    payload = {
        "total_entries": len(results),
        "ok": sum(1 for r in results if r.status == "ok"),
        "failed": sum(1 for r in results if r.status != "ok"),
        "entries": [
            {
                "slug": r.slug,
                "smiles": r.smiles,
                "resolved_charge": r.resolved_charge,
                "status": r.status,
                "best_xtb_eh": r.best_xtb_eh,
                "best_goat_eh": r.best_goat_eh,
                "winner_source": r.winner_source,
                "walltime_s": round(r.walltime_s, 3),
                "workdir": r.workdir,
                "error": r.error,
            }
            for r in results
        ],
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    _write_batch_summary_md(base_workdir, results)
    return path


def _write_batch_summary_md(base_workdir: Path, results: List[BatchEntryResult]) -> Path:
    """Publication-grade Markdown summary for the batch."""
    md_path = base_workdir / "batch_summary.md"
    total = len(results)
    ok = sum(1 for r in results if r.status == "ok")
    runtimes = sorted(r.walltime_s for r in results if r.walltime_s > 0)
    median_runtime = runtimes[len(runtimes) // 2] if runtimes else 0.0

    lines: List[str] = []
    lines.append("# GUPPY Batch Summary")
    lines.append("")
    lines.append(f"- Total entries: **{total}**")
    lines.append(f"- Successful: **{ok}** ({(ok / total * 100):.1f}%)" if total else "- Successful: 0")
    lines.append(f"- Failed: **{total - ok}**")
    lines.append(f"- Median walltime per entry: **{median_runtime:.1f} s**")
    lines.append("")
    lines.append("| # | slug | SMILES | charge | status | best (Eh) | winner | walltime (s) |")
    lines.append("|---|------|--------|--------|--------|-----------|--------|-------------:|")
    for idx, r in enumerate(results, start=1):
        best_energy = r.best_goat_eh if r.best_goat_eh is not None else r.best_xtb_eh
        best_str = f"{best_energy:.6f}" if best_energy is not None else "—"
        winner_str = r.winner_source or "—"
        smiles_disp = (r.smiles[:40] + "…") if len(r.smiles) > 41 else r.smiles
        lines.append(
            f"| {idx} | `{r.slug}` | `{smiles_disp}` | {r.resolved_charge:+d} | "
            f"{r.status} | {best_str} | {winner_str} | {r.walltime_s:.1f} |"
        )
    lines.append("")
    lines.append(
        "Reproduce: same input CSV + `--seed` + `--start-strategy` + `--max-isomers` "
        "yields byte-identical outputs (deterministic isomer enumeration)."
    )
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return md_path


def run_sampling_batch(
    smiles_entries: Sequence[Tuple[str, Optional[str]]],
    *,
    base_workdir: Path,
    shared_args: BatchSharedArgs,
    row: Optional[int] = None,
) -> int:
    """Iterate batch entries. Returns 0 if all succeed, 1 otherwise.

    If ``row`` is set (1-based), only that entry is processed — used by
    SLURM job arrays where each array task handles one SMILES.
    """
    base_workdir.mkdir(parents=True, exist_ok=True)

    if row is not None:
        if row < 1 or row > len(smiles_entries):
            logger.error(
                "Row %d out of range for %d entries.", row, len(smiles_entries)
            )
            return 1
        smiles, name = smiles_entries[row - 1]
        result = _run_one_entry(
            row, smiles, name, base_workdir=base_workdir, shared=shared_args
        )
        single_results = [result]
        _write_batch_results(base_workdir / result.slug, single_results)
        logger.info(
            "Row %d done: status=%s, best_final_eh=%s",
            row,
            result.status,
            result.best_goat_eh if result.best_goat_eh is not None else result.best_xtb_eh,
        )
        return 0 if result.status == "ok" else 1

    results: List[BatchEntryResult] = []
    for idx, (smi, name) in enumerate(smiles_entries, start=1):
        result = _run_one_entry(
            idx, smi, name, base_workdir=base_workdir, shared=shared_args
        )
        results.append(result)

    path = _write_batch_results(base_workdir, results)
    ok = sum(1 for r in results if r.status == "ok")
    logger.info(
        "Batch complete: %d/%d entries OK. Summary: %s",
        ok,
        len(results),
        path,
    )
    return 0 if ok == len(results) else 1


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="delfin-guppy-batch",
        description=(
            "Run GUPPY sampling for every SMILES in a batch file. Charge is "
            "always derived from each SMILES and flows into XTB and GOAT."
        ),
    )
    parser.add_argument(
        "input_file",
        help="Batch file: .csv or plain-text with one SMILES per line.",
    )
    parser.add_argument(
        "--workdir",
        default=os.environ.get("GUPPY_BATCH_WORKDIR", "GUPPY_BATCH"),
        help="Base directory for per-SMILES subdirectories (default: GUPPY_BATCH).",
    )
    parser.add_argument(
        "--row",
        type=int,
        default=None,
        help="Process only this 1-based row (for SLURM --array task index).",
    )
    parser.add_argument("--runs", type=int, default=int(os.environ.get("GUPPY_RUNS", "20")))
    parser.add_argument(
        "--pal",
        type=int,
        default=int(os.environ.get("GUPPY_PAL", os.environ.get("SLURM_CPUS_PER_TASK", "16"))),
    )
    parser.add_argument(
        "--maxcore",
        type=int,
        default=int(os.environ.get("GUPPY_MAXCORE", os.environ.get("DELFIN_MAXCORE", "6000"))),
    )
    parser.add_argument(
        "--parallel-jobs",
        type=int,
        default=int(os.environ.get("GUPPY_PARALLEL_JOBS", "4")),
    )
    parser.add_argument("--method", default=os.environ.get("GUPPY_XTB_METHOD", "XTB2"))
    parser.add_argument("--seed", type=int, default=int(os.environ.get("GUPPY_SEED", "31")))
    parser.add_argument(
        "--goat-topk",
        type=int,
        default=int(os.environ.get("GUPPY_GOAT_TOPK", "0")),
    )
    parser.add_argument(
        "--goat-parallel-jobs",
        type=int,
        default=int(
            os.environ.get("GUPPY_GOAT_PARALLEL_JOBS", os.environ.get("GUPPY_PARALLEL_JOBS", "4"))
        ),
    )
    parser.add_argument(
        "--start-strategy",
        choices=list(_ALLOWED_START_STRATEGIES),
        default=os.environ.get("GUPPY_START_STRATEGY", "isomers"),
    )
    parser.add_argument(
        "--max-isomers",
        type=int,
        default=int(os.environ.get("GUPPY_MAX_ISOMERS", "100")),
    )
    parser.add_argument(
        "--rmsd-cutoff",
        type=float,
        default=float(os.environ.get("GUPPY_RMSD_CUTOFF", "0.3")),
    )
    parser.add_argument(
        "--energy-window-kcal",
        type=float,
        default=float(os.environ.get("GUPPY_ENERGY_WINDOW_KCAL", "25.0")),
    )
    parser.add_argument("--allow-partial", action="store_true", default=True)
    parser.add_argument("--no-allow-partial", dest="allow_partial", action="store_false")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    input_path = Path(args.input_file)
    smiles_entries = _read_smiles_lines(input_path)
    logger.info("Loaded %d SMILES entries from %s", len(smiles_entries), input_path)

    shared = BatchSharedArgs(
        runs=args.runs,
        pal=args.pal,
        maxcore=args.maxcore,
        parallel_jobs=args.parallel_jobs,
        method=str(args.method).strip() or "XTB2",
        seed=args.seed,
        allow_partial=args.allow_partial,
        goat_topk=args.goat_topk,
        goat_parallel_jobs=args.goat_parallel_jobs,
        start_strategy=args.start_strategy,
        max_isomers=args.max_isomers,
        rmsd_cutoff=args.rmsd_cutoff,
        energy_window_kcal=args.energy_window_kcal,
    )

    row = args.row
    if row is None:
        env_task = os.environ.get("SLURM_ARRAY_TASK_ID")
        if env_task is not None:
            try:
                row = int(env_task)
            except ValueError:
                logger.warning("Ignoring non-integer SLURM_ARRAY_TASK_ID=%r", env_task)

    return run_sampling_batch(
        smiles_entries,
        base_workdir=Path(args.workdir),
        shared_args=shared,
        row=row,
    )


if __name__ == "__main__":
    raise SystemExit(main())
