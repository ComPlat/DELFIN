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


def _make_picks(*names: str):
    """Factory: builds an assertion that passes if ANY of the named tools
    were attempted (substring match, case-insensitive). Falls back to
    text-mention so the assertion isn't punished by plan-mode blocks."""
    def _check(res: AgentDryRunResult) -> tuple[bool, str]:
        if has_tool_call(res, *names):
            return True, names[0]
        text_low = res.all_text.lower()
        for n in names:
            if n in text_low:
                return True, f"named:{n}"
        return False, "missed"
    return _check


_thermo = _make_picks("extract_thermochem")
_dipole = _make_picks("extract_dipole")
_energy_table = _make_picks("extract_energy_table", "find_calculation_extreme")
_find_errors = _make_picks("find_orca_errors")
_lowest_gibbs = _make_picks("find_calculation_extreme", "extract_energy_table")
_diff_two = _make_picks("compare_calculations", "compare_across_functionals")
_validate_inp = _make_picks("validate_orca_input")
_submit = _make_picks("submit_calculation")
_list_jobs = _make_picks("list_active_calculations")
_get_pattern = _make_picks("get_dashboard_pattern", "list_dashboard_patterns")
_widgets = _make_picks("list_dashboard_widgets", "get_widget_options")
_explain_feature = _make_picks("explain_delfin_feature", "list_delfin_features")
_orca_manual = _make_picks(
    "check_orca_manual_indexed",
    "search_docs",
    "read_pdf",
)
_pdf_read = _make_picks("read_pdf", "search_pdf_local", "extract_pdf_section")
_calc_options = _make_picks("list_calc_options")
_create_folder = _make_picks("create_calc_folder")
_move_calc = _make_picks("move_calc_folder")
_delete_folder = _make_picks("delete_calc_folder")
_plot_distribution = _make_picks("plot_energy_distribution")
_plot_correlation = _make_picks("plot_energy_correlation")
_plot_orbital_diag = _make_picks("plot_orbital_diagram")
_pipeline_run = _make_picks("pipeline_run", "pipeline_prepare")
_cleanup = _make_picks("cleanup")
_describe_tool = _make_picks("describe_tool")
_ssh_jobs = _make_picks("list_ssh_transfer_jobs")


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
    # ---------- Phase A extras ----------
    ("thermochem question → extract_thermochem",
     "Gib mir die Thermochemie (T, P, ZPE, Enthalpie, Entropie, G) "
     "von {calc}/test_a.",
     _thermo),
    ("dipole question → extract_dipole",
     "Wie groß ist das Dipolmoment der Rechnung in {calc}/test_a "
     "(Vektor + Magnitude in Debye)?",
     _dipole),
    ("multi-folder energy table → extract_energy_table",
     "Gib mir eine Tabelle mit Gibbs- und Single-Point-Energien für "
     "{calc}/bp86, {calc}/pbe0, {calc}/b3lyp, {calc}/test_a.",
     _energy_table),
    ("ORCA error scan → find_orca_errors",
     "Such in {calc}/test_a nach bekannten ORCA-Fehlern (SCF-Diverge, "
     "OOM, Multiplicity, …).",
     _find_errors),
    ("lowest-Gibbs question → find_calculation_extreme",
     "Welcher Ordner unter {calc}/bp86, {calc}/pbe0, {calc}/b3lyp "
     "hat die niedrigste Gibbs-Energie?",
     _lowest_gibbs),
    ("two-folder diff → compare_calculations",
     "Vergleiche {calc}/bp86 und {calc}/pbe0 head-to-head: gleicher "
     "Functional? gleiche Basis? Δ Gibbs in kcal/mol?",
     _diff_two),
    # ---------- ORCA-Builder + lifecycle ----------
    ("INP validation → validate_orca_input",
     "Prüf bitte einen ORCA-Input-Block:\n"
     "! BP86 def2-SVP Opt Freq\n%pal nprocs 8 end\n%maxcore 2000\n"
     "*xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n"
     "Gibt es Probleme?",
     _validate_inp),
    ("submit one folder → submit_calculation",
     "Submitte den Ordner {calc}/bp86 als ORCA-Job (nur erklären "
     "welches Tool — keine Mutation).",
     _submit),
    ("list active jobs → list_active_calculations",
     "Welche meiner Jobs laufen gerade auf SLURM?",
     _list_jobs),
    ("ssh transfer status → list_ssh_transfer_jobs",
     "Welche SSH-Transfers laufen gerade ins Remote-Archiv?",
     _ssh_jobs),
    # ---------- Discovery + recipes ----------
    ("ask for batch recipe → get_dashboard_pattern",
     "Wie genau funktioniert /batch from-calc — zeig mir das Rezept.",
     _get_pattern),
    ("widget catalog question → list_dashboard_widgets",
     "Welche Widgets gibt's auf dem ORCA-Builder-Tab und welche "
     "Werte akzeptieren sie?",
     _widgets),
    ("DELFIN concept question → explain_delfin_feature",
     "Was macht der OCCUPIER in DELFIN? Erklär bitte mit Quellen "
     "aus dem DELFIN-Feature-Explainer.",
     _explain_feature),
    ("describe one tool → describe_tool",
     "Beschreibe das Tool extract_imaginary_frequencies — "
     "welche Argumente, was returnst es?",
     _describe_tool),
    # ---------- Literature ----------
    ("ORCA-syntax question → check_orca_manual_indexed",
     "Wie funktioniert in ORCA das %tddft-Block für CIS-Singulett-"
     "Berechnung? Schau im ORCA-Manual nach (nicht raten).",
     _orca_manual),
    ("ad-hoc PDF read → read_pdf",
     "Lies bitte aus der PDF /tmp/papers/example.pdf die Seiten 1-3 "
     "und sag mir, was drin steht.",
     _pdf_read),
    # ---------- Calc-options dropdown + folder mgmt ----------
    ("options-menu lookup → list_calc_options",
     "Welche Options stehen für eine ausgewählte CONTROL.txt zur "
     "Verfügung? Schau bitte im typed Calc-Options-Map nach.",
     _calc_options),
    ("create new folder → create_calc_folder",
     "Leg bitte einen neuen Ordner {calc}/new_run an. Sag nur welches "
     "Tool du verwenden würdest.",
     _create_folder),
    ("move folder → move_calc_folder",
     "Verschiebe {calc}/test_a nach {calc}/bp86 (calc → calc). "
     "Welches Tool?",
     _move_calc),
    ("delete folder (3-lock) → delete_calc_folder",
     "Lösch bitte den Ordner {calc}/test_a permanent. Welches Tool, "
     "und welche Sicherheits-Locks?",
     _delete_folder),
    # ---------- Plots beyond Phase D ----------
    ("energy distribution histogram → plot_energy_distribution",
     "Plotte ein Histogramm der Gibbs-Energien über alle Ordner "
     "in calc/ ({calc}/bp86, {calc}/pbe0, {calc}/b3lyp).",
     _plot_distribution),
    ("energy correlation scatter → plot_energy_correlation",
     "Plotte ein Streudiagramm Gibbs vs. Single-Point-Energy für "
     "{calc}/bp86, {calc}/pbe0, {calc}/b3lyp, {calc}/test_a.",
     _plot_correlation),
    ("orbital diagram plot → plot_orbital_diagram",
     "Zeig mir das Orbitaldiagramm um HOMO/LUMO für {calc}/test_a "
     "als Plot.",
     _plot_orbital_diag),
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
    p.add_argument("--out", default="",
                   help="write the markdown table here (default: stdout)")
    args = p.parse_args()

    if os.environ.get("DELFIN_AGENT_DRYRUN") != "1":
        print("error: set DELFIN_AGENT_DRYRUN=1 to run the live agent",
              file=sys.stderr)
        return 2

    models = [m.strip() for m in args.models.split(",") if m.strip()]
    cases = [c for c in CASES if not args.include
             or args.include.lower() in c[0].lower()]

    rows = []
    for i, (case_name, q, assertion) in enumerate(cases, 1):
        row: dict[str, Cell] = {}
        for m in models:
            print(f"  ({i}/{len(cases)}) [{m}] {case_name}",
                  file=sys.stderr, flush=True)
            row[m] = _run_one(case_name, q, assertion, m)
            # Per-cell flush so output survives a crash mid-sweep.
            if args.out:
                with open(args.out + ".partial", "a") as fh:
                    print(
                        f"  ({i}) {m} | {case_name} | "
                        f"{'PASS' if row[m].passed else 'FAIL'} | "
                        f"{row[m].note} | ${row[m].cost_usd:.3f}",
                        file=fh, flush=True,
                    )
        rows.append((case_name, row))

    out_path = args.out or None
    if out_path:
        import io
        buf = io.StringIO()
        _print_table_to(rows, models, fh=buf)
        with open(out_path, "w") as fh:
            fh.write(buf.getvalue())
        print(f"\nwrote table to {out_path}", file=sys.stderr)
        # Also stream to stdout so callers without --out still see it
        sys.stdout.write(buf.getvalue())
    else:
        _print_table_to(rows, models, fh=sys.stdout)
    return 0


def _print_table_to(rows, models, *, fh) -> None:
    """Variant of _print_table that takes an explicit file handle."""
    pass_glyph, fail_glyph = "✓", "✗"
    cost_total = {m: 0.0 for m in models}
    pass_count = {m: 0 for m in models}
    print(file=fh)
    print(f"## Per-model capability report (n_tests={len(rows)})", file=fh)
    print(file=fh)
    header = "| Test | " + " | ".join(models) + " |"
    sep = "|------|" + "|".join("------" for _ in models) + "|"
    print(header, file=fh)
    print(sep, file=fh)
    for name, by_model in rows:
        cells = []
        for m in models:
            c = by_model[m]
            glyph = pass_glyph if c.passed else fail_glyph
            cells.append(f"{glyph} {c.note} (${c.cost_usd:.3f})")
            cost_total[m] += c.cost_usd
            pass_count[m] += int(c.passed)
        print(f"| {name} | " + " | ".join(cells) + " |", file=fh)
    print(file=fh)
    print("### Totals", file=fh)
    for m in models:
        print(f"- {m}: {pass_count[m]}/{len(rows)} passed, "
              f"${cost_total[m]:.3f} spent", file=fh)


if __name__ == "__main__":
    raise SystemExit(main())
