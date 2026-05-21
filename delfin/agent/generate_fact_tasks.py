"""Generate fact-verify benchmark tasks from the ORCA keyword ground-truth.

Reads ``pack/benchmark/orca_keywords_groundtruth.json`` and emits
``pack/benchmark/tasks_auto_orca.yaml`` — one ``fact_verify`` task per
covered block.  Each task expects 2-4 canonical keywords (curated +
validated against the manual extract) and forbids known hallucination
patterns (sourced from observed production failures).

Re-run after the manual is re-indexed (or after adding a new block):

    python -m delfin.agent.generate_fact_tasks

The generator BUILD-FAILS if a curated keyword is not present in the
extracted manual namespace — that's the guard against author memory
drifting from manual ground-truth.

Tasks_auto_orca.yaml is committed to the repo so test results are
reproducible without re-running the extractor every time.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import yaml


_HERE = Path(__file__).resolve().parent
_PACK = _HERE / "pack" / "benchmark"


# Per-program block test catalogues.  Add a new program by adding an
# entry here; the loop at the bottom emits one tasks_auto_<program>.yaml
# file per program with non-empty entries that match the local doc-index.
_PROGRAM_BLOCK_TESTS: dict[str, dict[str, dict[str, Any]]] = {
    "orca": {
    "casscf": {
        "must_have": ["nel", "norb", "mult", "nroots"],
        "forbid":   ["nactel", "nactorb", "multiplicity",
                     "nelectrons", "numorbs"],
        "label": "CASSCF active-space setup",
    },
    "tddft": {
        "must_have": ["nroots", "tda"],
        "forbid":   ["nstates", "tdamod", "ntdroots",
                     "numroots", "states"],
        "label": "TDDFT excitation",
    },
    "mp2": {
        # MP2 block shares keys with %mdci in the manual; pick the
        # most-fundamental tokens.
        "must_have": ["density"],
        "forbid":   ["mp2option", "mp2_ri", "ri_mp2_option"],
        "label": "MP2 correlation",
    },
    "scf": {
        "must_have": ["maxiter", "guess", "convergence"],
        "forbid":   ["maxiterations", "initialguess", "scfconvergence"],
        "label": "SCF convergence",
    },
    "geom": {
        "must_have": ["maxiter"],
        "forbid":   ["maxiterations", "max_iter_count"],
        "label": "geometry optimisation",
    },
    "freq": {
        "must_have": ["scalfreq"],
        "forbid":   ["scaling", "freqscale", "scaling_factor"],
        "label": "frequency calculation",
    },
    "eprnmr": {
        "must_have": ["nuclei"],
        "forbid":   ["atom_list", "nuclear_list", "magnuclei"],
        "label": "EPR/NMR property calculation",
    },
    },  # end of orca block tests
    # ----------------------------------------------------------------
    # Turbomole — add curated must_have lists here as the TM manual
    # gets indexed.  Generator simply skips programs whose ground-truth
    # JSON is empty / missing.
    "turbomole": {
        # "dft":   {"must_have": ["functional", "gridsize"],
        #            "forbid": [...], "label": "DFT setup"},
    },
    # Gaussian / NWChem / Q-Chem / Psi4 etc. — add blocks here.
}


def _gt_path_for(program: str) -> Path:
    return _PACK / f"keywords_groundtruth_{program}.json"


def _out_yaml_for(program: str) -> Path:
    return _PACK / f"tasks_auto_{program}.yaml"


def _load_groundtruth(program: str) -> dict[str, Any] | None:
    """Return parsed ground-truth for ``program`` or None if file is
    missing — caller skips programs without indexed manuals."""
    p = _gt_path_for(program)
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


def _validate_against_manual(
    block_name: str,
    cfg: dict[str, Any],
    manual_keywords: set[str],
) -> None:
    """Fail loud if curated must_have isn't in manual, or if forbid
    accidentally includes a real keyword."""
    missing = [kw for kw in cfg["must_have"]
               if kw.lower() not in manual_keywords]
    if missing:
        raise RuntimeError(
            f"BLOCK '{block_name}': curated must_have keywords {missing} "
            f"NOT present in the extracted manual namespace.  Either the "
            f"keyword names are wrong, or the extractor's keyword set for "
            f"this block is incomplete (only {len(manual_keywords)} keywords "
            f"extracted)."
        )
    real_in_forbid = [kw for kw in cfg["forbid"]
                     if kw.lower() in manual_keywords]
    if real_in_forbid:
        raise RuntimeError(
            f"BLOCK '{block_name}': forbid list includes real manual "
            f"keywords {real_in_forbid} — would penalise correct answers"
        )


def _task_for_block(
    program: str, block_name: str, block_marker: str, cfg: dict[str, Any],
) -> dict[str, Any]:
    """Build one fact_verify task dict for a single block.

    ``program`` (e.g. "orca", "turbomole") + ``block_marker``
    (e.g. ``%casscf``, ``$dft``) appear verbatim in the prompt so the
    model is queried with program-specific syntax.
    """
    must = cfg["must_have"]
    forbid = cfg.get("forbid") or []
    label = cfg["label"]
    task = {
        "id": f"fact_{program}_{block_name}_keywords_auto",
        "task_class": "fact_verify_auto",
        "mode": "dashboard",
        "prompt": (
            f"welche EXAKTEN keyword-namen nutzt {program.upper()} im "
            f"{block_marker}-block für eine typische {label}? Liste "
            f"mindestens die essentiellen Keywords aus dem {program.upper()} "
            f"Manual auf. Wichtig: nur Keywords die tatsächlich im Manual "
            f"stehen — bitte vorher per doc-search im Manual nachschlagen."
        ),
        "expected_signals": [
            {"pattern": rf"(?i)\b{re.escape(kw)}\b", "against": "text"}
            for kw in must
        ],
        "max_duration_s": 120,
        "max_cost_usd": 0.50,
        "max_tool_calls": 10,
    }
    if forbid:
        task["forbidden_signals"] = [
            {"pattern": rf"(?i)\b{re.escape(kw)}\b", "against": "text"}
            for kw in forbid
        ]
    return task


def build_tasks_for(program: str) -> list[dict[str, Any]]:
    """Build the fact_verify_auto tasks for one program."""
    block_tests = _PROGRAM_BLOCK_TESTS.get(program) or {}
    if not block_tests:
        return []
    gt = _load_groundtruth(program)
    if gt is None:
        # No indexed manual yet for this program — silently skip
        return []
    namespace = gt.get("blocks") or {}
    out: list[dict[str, Any]] = []
    for block_name, cfg in block_tests.items():
        info = namespace.get(block_name)
        if not info:
            continue
        manual_keywords = {k.lower() for k in info.get("keywords", [])}
        _validate_against_manual(block_name, cfg, manual_keywords)
        block_marker = info.get("block") or f"%{block_name}"
        out.append(_task_for_block(program, block_name, block_marker, cfg))
    return out


def write_yaml(
    program: str,
    tasks: list[dict[str, Any]],
    path: Path | None = None,
) -> Path:
    p = path or _out_yaml_for(program)
    header = (
        f"# AUTO-GENERATED from keywords_groundtruth_{program}.json.\n"
        "# DO NOT edit this file by hand — re-run\n"
        "# `python -m delfin.agent.generate_fact_tasks` after editing\n"
        "# _PROGRAM_BLOCK_TESTS in delfin/agent/generate_fact_tasks.py.\n\n"
    )
    body = yaml.safe_dump(
        {"tasks": tasks}, sort_keys=False, allow_unicode=True,
        default_flow_style=False,
    )
    p.write_text(header + body, encoding="utf-8")
    return p


def main() -> int:
    """Generate tasks for every program that has both block_tests + a
    ground-truth JSON.  Programs without one or the other are skipped."""
    total = 0
    for program in sorted(_PROGRAM_BLOCK_TESTS):
        tasks = build_tasks_for(program)
        if not tasks:
            continue
        path = write_yaml(program, tasks)
        print(f"[{program}] generated {len(tasks)} tasks → {path.name}")
        for t in tasks:
            n_exp = len(t.get("expected_signals", []))
            n_for = len(t.get("forbidden_signals", []))
            print(f"  {t['id']:<48}  expected={n_exp}  forbidden={n_for}")
        total += len(tasks)
    if total == 0:
        print("No tasks generated.  Add entries to _PROGRAM_BLOCK_TESTS "
              "and ensure keywords_groundtruth_<program>.json exists.")
    else:
        print(f"\nTotal: {total} auto-generated fact_verify_auto tasks.")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
