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
_GROUNDTRUTH = _HERE / "pack" / "benchmark" / "orca_keywords_groundtruth.json"
_OUT_YAML = _HERE / "pack" / "benchmark" / "tasks_auto_orca.yaml"


# Per-block test definitions.
#
# must_have:   keywords the model must mention (validated against the
#              extracted manual namespace — generator aborts if any of
#              these is NOT actually in the manual)
# forbid:      hallucination patterns observed in real production OR
#              plausible-looking synonyms that ARE NOT in the manual
#              (validator confirms forbidden list does not accidentally
#              include real keywords)
# label:       short user-facing description of what the block is for
_BLOCK_TESTS: dict[str, dict[str, Any]] = {
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
}


def _load_groundtruth() -> dict[str, Any]:
    if not _GROUNDTRUTH.exists():
        raise FileNotFoundError(
            f"ORCA keyword ground-truth missing: {_GROUNDTRUTH}.  "
            "Run extract_keyword_namespace + write JSON first."
        )
    return json.loads(_GROUNDTRUTH.read_text(encoding="utf-8"))


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


def _task_for_block(block_name: str, cfg: dict[str, Any]) -> dict[str, Any]:
    """Build one fact_verify task dict for a single block."""
    must = cfg["must_have"]
    forbid = cfg["forbid"]
    label = cfg["label"]
    task = {
        "id": f"fact_orca_{block_name}_keywords_auto",
        "task_class": "fact_verify_auto",
        "mode": "dashboard",
        "prompt": (
            f"welche EXAKTEN keyword-namen nutzt ORCA im %{block_name} "
            f"block für eine typische {label}? Liste mindestens die "
            f"essentiellen Keywords aus dem ORCA Manual auf. Wichtig: "
            f"nur Keywords die tatsächlich im Manual stehen — bitte "
            f"vorher per doc-search im Manual nachschlagen."
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


def build_tasks() -> list[dict[str, Any]]:
    gt = _load_groundtruth()
    namespace = gt.get("blocks") or {}
    out: list[dict[str, Any]] = []
    for block_name, cfg in _BLOCK_TESTS.items():
        info = namespace.get(block_name)
        if not info:
            continue
        manual_keywords = {k.lower() for k in info.get("keywords", [])}
        _validate_against_manual(block_name, cfg, manual_keywords)
        out.append(_task_for_block(block_name, cfg))
    return out


def write_yaml(tasks: list[dict[str, Any]], path: Path | None = None) -> Path:
    p = path or _OUT_YAML
    header = (
        "# AUTO-GENERATED from orca_keywords_groundtruth.json.\n"
        "# DO NOT edit this file by hand — re-run\n"
        "# `python -m delfin.agent.generate_fact_tasks` after editing\n"
        "# _BLOCK_TESTS in delfin/agent/generate_fact_tasks.py.\n\n"
    )
    body = yaml.safe_dump(
        {"tasks": tasks}, sort_keys=False, allow_unicode=True,
        default_flow_style=False,
    )
    p.write_text(header + body, encoding="utf-8")
    return p


def main() -> int:
    tasks = build_tasks()
    path = write_yaml(tasks)
    print(f"Generated {len(tasks)} fact_verify_auto tasks → {path}")
    for t in tasks:
        n_exp = len(t.get("expected_signals", []))
        n_for = len(t.get("forbidden_signals", []))
        print(f"  {t['id']:<48}  expected={n_exp}  forbidden={n_for}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
