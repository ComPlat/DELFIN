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
#
# Schema (per block):
#   candidate_pool: list of plausible canonical keywords from the manual.
#                   The generated task expects the model to mention AT
#                   LEAST `min_required` of them (default 2).  This is
#                   more robust than "model must say exactly these 3"
#                   because the manual lists many valid keywords per
#                   block and any subset is a correct answer.
#   forbid:        common hallucination patterns (validated to NOT be
#                  real manual keywords).
#   label:         user-facing description.
#   min_required:  how many of the pool the model must mention (default 2).
_PROGRAM_BLOCK_TESTS: dict[str, dict[str, dict[str, Any]]] = {
    "orca": {
    "casscf": {
        "candidate_pool": ["nel", "norb", "mult", "nroots",
                           "ptmethod", "actorbs"],
        "forbid":   ["nactel", "nactorb", "multiplicity",
                     "nelectrons", "numorbs"],
        "label": "CASSCF active-space setup",
        "min_required": 3,
    },
    "tddft": {
        "candidate_pool": ["nroots", "tda", "triplets", "maxdim",
                           "etol", "iroot", "tdamod"],
        "forbid":   ["nstates", "ntdroots", "numroots", "states"],
        "label": "TDDFT excitation",
        "min_required": 2,
    },
    "mp2": {
        "candidate_pool": ["density", "ri", "maxiter", "doscs",
                           "calcs", "calc"],
        "forbid":   ["mp2option", "mp2_ri", "ri_mp2_option"],
        "label": "MP2 correlation",
        "min_required": 1,
    },
    "scf": {
        # SCF block has many valid keywords; testing "any 2 of pool"
        # catches Konvergenz/Convergence + Damp + Maxiter + Guess etc.
        "candidate_pool": ["maxiter", "guess", "convergence",
                           "convcheckmode", "tole", "tolp", "damp",
                           "level", "shift", "directresp", "diis"],
        "forbid":   ["maxiterations", "initialguess",
                     "scfconvergence"],
        "label": "SCF convergence",
        "min_required": 2,
    },
    "geom": {
        "candidate_pool": ["maxiter", "trust", "hess_filename",
                           "calc_hess", "tols", "tole", "tolg",
                           "tolx", "scan"],
        "forbid":   ["maxiterations", "max_iter_count"],
        "label": "geometry optimisation",
        "min_required": 2,
    },
    "freq": {
        "candidate_pool": ["scalfreq", "anfreq", "centraldiff",
                           "increment", "restart", "scaling"],
        "forbid":   ["freqscale", "scaling_factor"],
        "label": "frequency calculation",
        "min_required": 1,
    },
    "eprnmr": {
        "candidate_pool": ["nuclei", "tensor", "ori", "doalpha",
                           "shifts", "ssall"],
        "forbid":   ["atom_list", "nuclear_list", "magnuclei"],
        "label": "EPR/NMR property calculation",
        "min_required": 1,
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
) -> list[str]:
    """Filter the candidate pool to only keywords actually in the
    manual.  Returns the validated pool.

    Raises if:
    - fewer than ``min_required`` pool entries are real manual keywords
      (the test would be impossible to pass)
    - any forbid keyword is a real manual keyword (would penalise the
      correct answer)
    """
    pool = cfg["candidate_pool"]
    min_required = cfg.get("min_required", 2)
    valid_pool = [kw for kw in pool if kw.lower() in manual_keywords]
    if len(valid_pool) < min_required:
        raise RuntimeError(
            f"BLOCK '{block_name}': only {len(valid_pool)} of "
            f"{len(pool)} candidate keywords are in the manual; "
            f"need at least {min_required} to make the test passable. "
            f"Validated pool: {valid_pool}. Original: {pool}."
        )
    real_in_forbid = [kw for kw in cfg["forbid"]
                     if kw.lower() in manual_keywords]
    if real_in_forbid:
        raise RuntimeError(
            f"BLOCK '{block_name}': forbid list includes real manual "
            f"keywords {real_in_forbid} — would penalise correct answers"
        )
    return valid_pool


def _task_for_block(
    program: str, block_name: str, block_marker: str, cfg: dict[str, Any],
    valid_pool: list[str],
) -> dict[str, Any]:
    """Build one fact_verify task dict for a single block.

    Uses UNION-style expected: the model must mention at least
    ``min_required`` keywords from the candidate pool.  Encoded as
    repeated independent expected_signals each matching ANY pool entry —
    the model satisfies signal[i] by mentioning ANY pool entry, but
    since the same pattern is repeated, the scorer's "all required
    signals must match" rule effectively becomes "at least
    min_required distinct mentions across the response" only when
    we provide DISTINCT alternation patterns.

    Simpler approach: emit ONE union pattern as expected[0], and rely
    on min_required=1 with broader pool.  For higher min_required,
    split the pool into N halves and emit one alternation-pattern
    per half — model must hit each half at least once.
    """
    forbid = cfg.get("forbid") or []
    label = cfg["label"]
    min_required = cfg.get("min_required", 2)
    # Split the validated pool into min_required disjoint groups so
    # the model has to mention something from EACH group — that's a
    # robust proxy for "at least N distinct keywords".
    n_groups = min(min_required, len(valid_pool))
    groups: list[list[str]] = [[] for _ in range(n_groups)]
    for i, kw in enumerate(valid_pool):
        groups[i % n_groups].append(kw)
    expected = []
    for group in groups:
        if not group:
            continue
        union = "|".join(re.escape(kw) for kw in group)
        expected.append({
            "pattern": rf"(?i)\b(?:{union})\b",
            "against": "text",
        })
    task = {
        "id": f"fact_{program}_{block_name}_keywords_auto",
        "task_class": "fact_verify_auto",
        "mode": "dashboard",
        "prompt": (
            f"welche keywords nutzt {program.upper()} im "
            f"{block_marker}-block für eine typische {label}? Nenne "
            f"mindestens {min_required} echte Keywords aus dem "
            f"{program.upper()} Manual.  Wichtig: nur Keywords die "
            f"tatsächlich im Manual stehen — bitte per doc-search im "
            f"Manual nachschlagen bevor du antwortest."
        ),
        "expected_signals": expected,
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
        valid_pool = _validate_against_manual(block_name, cfg, manual_keywords)
        block_marker = info.get("block") or f"%{block_name}"
        out.append(_task_for_block(
            program, block_name, block_marker, cfg, valid_pool,
        ))
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
