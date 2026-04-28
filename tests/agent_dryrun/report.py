"""Per-model capability report runner.

Drives every test in :mod:`tests.agent_dryrun.test_capabilities`
against a list of models, captures pass/fail + cost + tool picks per
(test, model) cell, and prints a markdown table so you can see at a
glance "haiku misses on UV/Vis but sonnet gets it" — i.e. where to
scale up the model and where it's safe to keep it cheap.

Run::

    DELFIN_AGENT_DRYRUN=1 \\
    DELFIN_AGENT_DRYRUN_MODELS=haiku,sonnet \\
        python -m tests.agent_dryrun.report

Cost: O(n_models * n_tests) live-CLI turns. Roughly $0.10 per cell on
haiku, $0.50 per cell on sonnet. Scale your model list accordingly.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass

from .runner import AgentDryRunResult, has_tool_call, run_agent_dryrun
from .sandbox import sandboxed


# (test_name, query_template, assertion_lambda)
# query_template: ``{calc}`` is replaced with sandbox calc-dir path.
# assertion: takes the AgentDryRunResult and returns (pass: bool, note: str).


def _imag(res: AgentDryRunResult) -> tuple[bool, str]:
    ok = has_tool_call(res, "extract_imaginary_frequencies")
    return ok, "extract_imaginary_frequencies" if ok else "missed"


def _funcs(res: AgentDryRunResult) -> tuple[bool, str]:
    for tname in ("compare_across_functionals",
                  "compare_calculations",
                  "extract_energy_table"):
        if has_tool_call(res, tname):
            return True, tname
    return False, "missed"


def _rename(res: AgentDryRunResult) -> tuple[bool, str]:
    if has_tool_call(res, "rename_calc_folder"):
        return True, "tool"
    if "rename_calc_folder" in res.all_text.lower():
        return True, "named-in-text"
    return False, "missed"


def _archive(res: AgentDryRunResult) -> tuple[bool, str]:
    if has_tool_call(res, "move_to_archive"):
        return True, "tool"
    if "move_to_archive" in res.all_text.lower():
        return True, "named-in-text"
    return False, "missed"


def _kill_all(res: AgentDryRunResult) -> tuple[bool, str]:
    if has_tool_call(res, "kill_all_user_jobs"):
        return True, "tool"
    if "kill_all_user_jobs" in res.all_text.lower():
        return True, "named-in-text"
    return False, "missed"


def _recalc(res: AgentDryRunResult) -> tuple[bool, str]:
    if has_tool_call(res, "prepare_recalc"):
        return True, "tool"
    if "prepare_recalc" in res.all_text.lower():
        return True, "named-in-text"
    return False, "missed"


def _orbitals(res: AgentDryRunResult) -> tuple[bool, str]:
    if has_tool_call(res, "extract_orbital_energies",
                     "plot_orbital_diagram"):
        return True, "tool"
    if "extract_orbital_energies" in res.all_text.lower():
        return True, "named-in-text"
    return False, "missed"


def _uvvis(res: AgentDryRunResult) -> tuple[bool, str]:
    for n in ("plot_uvvis_spectrum", "extract_excited_states"):
        if has_tool_call(res, n) or n in res.all_text.lower():
            return True, n
    return False, "missed"


def _opt_traj(res: AgentDryRunResult) -> tuple[bool, str]:
    for n in ("extract_optimization_trajectory",
              "plot_optimization_convergence"):
        if has_tool_call(res, n) or n in res.all_text.lower():
            return True, n
    return False, "missed"


def _batch_safe(res: AgentDryRunResult) -> tuple[bool, str]:
    text = res.all_text.lower()
    if has_tool_call(res, "get_dashboard_pattern"):
        return True, "fetched recipe"
    if "/batch from-calc" in text:
        after = text.split("/batch from-calc", 1)[1].split("\n", 1)[0][:80]
        if "initial.xyz" not in after:
            return True, "/batch from-calc (no bad filter)"
        return False, "/batch from-calc initial.xyz (BUG)"
    return False, "no /batch action"


def _list_tools(res: AgentDryRunResult) -> tuple[bool, str]:
    if has_tool_call(res, "list_tools", "describe_tool"):
        return True, "tool"
    return False, "missed"


CASES = [
    ("imag-freq → extract_imaginary_frequencies",
     "Hat die Rechnung in {calc}/test_a noch imaginäre Frequenzen "
     "oder ist es ein Minimum?",
     _imag),
    ("functional comparison → compare_across_functionals",
     "Vergleiche die Gibbs-Energien über die Funktionale für "
     "{calc}/bp86, {calc}/pbe0, {calc}/b3lyp.",
     _funcs),
    ("rename calc folder → rename_calc_folder",
     "Benenne den Ordner {calc}/test_a in test_a_renamed um. "
     "Sag mir nur welches Tool — nicht ausführen.",
     _rename),
    ("archive intent → move_to_archive",
     "Verschiebe den Ordner {calc}/test_a ins Archiv. "
     "Erkläre nur das Tool — nichts ausführen.",
     _archive),
    ("cancel-all → kill_all_user_jobs",
     "Brich alle meine SLURM-Jobs ab. Welches Tool und welche "
     "Job-IDs wären betroffen?",
     _kill_all),
    ("smart-recalc prep → prepare_recalc",
     "Bereite einen Smart-Recalc für {calc}/test_a vor "
     "(noch nicht abschicken — nur den Plan zeigen).",
     _recalc),
    ("HOMO/LUMO question → extract_orbital_energies",
     "Wo liegt das HOMO und das LUMO der Rechnung in "
     "{calc}/test_a, und wie groß ist der Gap?",
     _orbitals),
    ("UV/Vis question → plot_uvvis_spectrum / extract_excited_states",
     "Zeig mir das UV/Vis-Spektrum aus der TDDFT-Rechnung in "
     "{calc}/test_a — am liebsten als Plot.",
     _uvvis),
    ("opt convergence → extract_optimization_trajectory / plot",
     "Warum hat die Optimierung in {calc}/test_a so viele "
     "Schritte gebraucht? Zeig mir den Verlauf.",
     _opt_traj),
    ("batch UX safety → no 'initial.xyz' as folder-filter",
     "Batch alle initial.xyz Dateien aus calc/. "
     "Erkläre nur den Plan — noch nicht ausführen.",
     _batch_safe),
    ("catalog discovery → list_tools / describe_tool",
     "Welche typed-MCP-Tools für Output-Parsing hast du? "
     "Such bitte im Tool-Katalog statt aus dem Gedächtnis.",
     _list_tools),
]


@dataclass
class Cell:
    passed: bool
    note: str
    cost_usd: float
    duration_s: float
    timed_out: bool


def _run_one(case_name, query_tmpl, assertion, model: str) -> Cell:
    with sandboxed() as sb:
        query = query_tmpl.format(calc=str(sb.calc))
        t0 = time.time()
        try:
            res = run_agent_dryrun(
                query, model=model, timeout_s=90,
                cwd=str(sb.root), sandbox_root=str(sb.root),
            )
        except Exception as exc:
            return Cell(False, f"runner exception: {exc}", 0.0,
                        time.time() - t0, False)
        dt = time.time() - t0
        if res.timed_out:
            return Cell(False, "TIMED OUT", res.cost_usd, dt, True)
        ok, note = assertion(res)
        return Cell(ok, note, res.cost_usd, dt, False)


def _print_table(rows: list[tuple[str, dict[str, Cell]]],
                 models: list[str]) -> None:
    pass_glyph, fail_glyph = "✓", "✗"
    cost_total = {m: 0.0 for m in models}
    pass_count = {m: 0 for m in models}
    print()
    print(f"## Per-model capability report (n_tests={len(rows)})")
    print()
    header = "| Test | " + " | ".join(models) + " |"
    sep = "|------|" + "|".join("------" for _ in models) + "|"
    print(header)
    print(sep)
    for name, by_model in rows:
        cells = []
        for m in models:
            c = by_model[m]
            glyph = pass_glyph if c.passed else fail_glyph
            tag = c.note if not c.passed else c.note
            cells.append(f"{glyph} {tag} (${c.cost_usd:.3f})")
            cost_total[m] += c.cost_usd
            pass_count[m] += int(c.passed)
        print(f"| {name} | " + " | ".join(cells) + " |")
    print()
    print("### Totals")
    for m in models:
        print(f"- {m}: {pass_count[m]}/{len(rows)} passed, "
              f"${cost_total[m]:.3f} spent")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--models", default="haiku")
    p.add_argument("--include", default="",
                   help="substring filter on test name (default: all)")
    args = p.parse_args()

    if os.environ.get("DELFIN_AGENT_DRYRUN") != "1":
        print("error: set DELFIN_AGENT_DRYRUN=1 to run the live agent",
              file=sys.stderr)
        return 2

    models = [m.strip() for m in args.models.split(",") if m.strip()]
    cases = [c for c in CASES if not args.include
             or args.include.lower() in c[0].lower()]

    rows = []
    for case_name, q, assertion in cases:
        row: dict[str, Cell] = {}
        for m in models:
            print(f"  · [{m}] {case_name}", file=sys.stderr)
            row[m] = _run_one(case_name, q, assertion, m)
        rows.append((case_name, row))
    _print_table(rows, models)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
