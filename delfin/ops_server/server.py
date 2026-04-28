"""MCP server for DELFIN runtime actions.

Exposes ``delfin.api`` functions as MCP tools so an agent can invoke DELFIN
workflows in a structured, typed way (no fragile ACTION-string parsing).

Tools provided
--------------
Read-only (safe):
    - ``qm_check``         — inspect QM tool resolution
    - ``csp_check``        — CSP tool availability
    - ``mlp_check``        — MLP backend availability
    - ``analysis_check``   — analysis tool availability
    - ``stop_dry_run``     — list DELFIN processes that *would* be signaled

Mutating (require explicit ``allow_mutate=True``):
    - ``cleanup``          — remove scratch artifacts (supports dry_run)
    - ``stop``             — signal DELFIN processes
    - ``pipeline_run``     — full DELFIN pipeline
    - ``pipeline_prepare`` — generate CONTROL template
    - ``run_orca_input``   — run ORCA on a .inp file
    - ``co2``              — CO2 Coordinator workflow
    - ``tadf_xtb``         — TADF xTB workflow
    - ``hyperpol``         — hyperpolarisability workflow

Mutating tools refuse to execute unless the ``allow_mutate`` parameter is
passed explicitly.  This matches the dashboard's rule that destructive
agent actions require user confirmation.

The tool functions are defined at module level so they can be imported and
tested without the optional FastMCP dependency.  ``run_server`` only loads
FastMCP when actually starting the stdio server.
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Any

from delfin import api as delfin_api


# ---------------------------------------------------------------------------
# Result formatting
# ---------------------------------------------------------------------------

_STDOUT_TRUNC = 16_000  # chars; large enough for normal CLI output


def _format_result(rc: delfin_api.CommandResult, *, action: str, dry_run: bool = False) -> str:
    """Format a ``CommandResult`` as a JSON string for the MCP tool reply."""
    stdout = rc.stdout or ""
    stderr = rc.stderr or ""
    if len(stdout) > _STDOUT_TRUNC:
        stdout = stdout[:_STDOUT_TRUNC] + (
            f"\n... [truncated {len(rc.stdout) - _STDOUT_TRUNC} chars]"
        )
    if len(stderr) > _STDOUT_TRUNC:
        stderr = stderr[:_STDOUT_TRUNC] + (
            f"\n... [truncated {len(rc.stderr) - _STDOUT_TRUNC} chars]"
        )
    payload: dict[str, Any] = {
        "action": action,
        "dry_run": dry_run,
        "returncode": rc.returncode,
        "ok": rc.ok,
        "stdout": stdout,
        "stderr": stderr,
        "argv": rc.argv,
    }
    return json.dumps(payload, indent=2, ensure_ascii=False)


def _refuse_mutation(action: str) -> str:
    """Standard refusal payload when ``allow_mutate`` is not set."""
    return json.dumps({
        "action": action,
        "ok": False,
        "error": "mutation_blocked",
        "message": (
            f"The tool '{action}' modifies state. Pass allow_mutate=True "
            "to execute. The agent should ask the user before doing so."
        ),
    }, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# Read-only tool implementations
# ---------------------------------------------------------------------------

def tool_qm_check(tools: str = "") -> str:
    """Check QM tool resolution (xtb, crest, xtb4stda, std2, stda, dftb+).

    Args:
        tools: optional comma-separated subset (e.g. "xtb,crest"). Empty = all.
    """
    names = [t.strip() for t in tools.split(",") if t.strip()] or None
    rc = delfin_api.qm_check(tools=names)
    return _format_result(rc, action="qm_check")


def tool_csp_check() -> str:
    """Check CSP (genarris) tool availability."""
    return _format_result(delfin_api.csp_check(), action="csp_check")


def tool_mlp_check() -> str:
    """Check MLP backend availability (torchani, AIMNet2, MACE)."""
    return _format_result(delfin_api.mlp_check(), action="mlp_check")


def tool_analysis_check() -> str:
    """Check analysis tools (Multiwfn, CENSO, ANMR, morfeus)."""
    return _format_result(delfin_api.analysis_check(), action="analysis_check")


def tool_list_dashboard_patterns() -> str:
    """List the names of operational-pattern recipes available on demand.

    Each name maps to a concrete slash-chain recipe for one dashboard
    workflow (batch jobs, smart recalc, ORCA submit, …). Use the
    returned list to pick a name, then call get_dashboard_pattern(name)
    to fetch only the recipe you need — none are pre-loaded into the
    system prompt, so this is the only way to see them.
    """
    names = delfin_api.list_dashboard_patterns()
    return "Available dashboard pattern names:\n- " + "\n- ".join(names)


def tool_get_dashboard_pattern(name: str) -> str:
    """Return the slash-chain recipe for a named dashboard workflow.

    Available names (call list_dashboard_patterns to see them all):
    batch, control_edit, smart_recalc, submit_orca, analyze, recalc,
    cancel.

    Use this when the user asks for one of those workflows and you
    aren't 100% sure of the exact ACTION: chain — the recipe contains
    the verbatim slash commands you should emit, plus the don't-
    reinvent rules (e.g. never hand-roll batch text, always use
    /batch from-calc).

    Args:
        name: the pattern to fetch (case-insensitive; "-" / " " also OK).
    """
    return delfin_api.get_dashboard_pattern(name)


# ---------------------------------------------------------------------------
# P1 — Output parsing tools (read-only, structured returns)
# ---------------------------------------------------------------------------


def _orca_parse_to_dict(parsed) -> dict:
    """Render an OrcaParseResult as a stable JSON-friendly dict."""
    return {
        "path": parsed.path,
        "final_single_point": parsed.final_single_point,
        "gibbs_free_energy": parsed.gibbs_free_energy,
        "zpe": parsed.zpe,
        "scf_converged": parsed.scf_converged,
        "opt_converged": parsed.opt_converged,
        "imag_freq_count": parsed.imag_freq_count,
        "walltime_s": parsed.walltime_s,
        "n_atoms": parsed.n_atoms,
        "functional": parsed.functional,
        "basis": parsed.basis,
        "error_summary": parsed.error_summary,
    }


def tool_parse_orca_output(path: str) -> str:
    """Parse one ORCA .out file and return a structured snapshot.

    Returns a JSON object with: final_single_point (Hartree),
    gibbs_free_energy (Hartree), zpe (Hartree), scf_converged (bool),
    opt_converged (bool), imag_freq_count (int), walltime_s (float),
    n_atoms (int), functional (str), basis (str), error_summary (str).
    Missing values are null. Use this BEFORE writing a Python script
    to grep the file — one tool call replaces dozens of regexes.

    Args:
        path: absolute path to the ORCA .out file.
    """
    import json as _json
    parsed = delfin_api.parse_orca_output(path)
    return _json.dumps(_orca_parse_to_dict(parsed), indent=2)


def tool_find_orca_errors(folder: str) -> str:
    """Scan all *.out files in ``folder`` for known ORCA error patterns.

    Returns a JSON list of {type, message, line_number, suggestion}
    entries. Empty list = no patterns matched (NOT proof of success —
    use parse_orca_output for that).

    Detected error types: scf_diverge, oom, basis, multiplicity, mpi,
    timeout, other.

    Args:
        folder: absolute path to the calc folder containing .out files.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    errors = delfin_api.find_orca_errors(folder)
    return _json.dumps([_asdict(e) for e in errors], indent=2)


def tool_extract_thermochem(folder: str) -> str:
    """Extract the full thermochemistry block from an ORCA Freq output.

    Picks the first .out in ``folder`` containing thermochemistry data.
    Returns JSON with: temperature_k, pressure_atm, zpe, thermal_corr,
    enthalpy_corr, entropy_total, gibbs_corr, final_gibbs (all in
    Hartree except T and P).

    Args:
        folder: absolute path to a calc folder.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    result = delfin_api.extract_thermochem(folder)
    return _json.dumps(_asdict(result), indent=2)


def tool_extract_energy_table(
    folders: str,
    properties: str = "",
) -> str:
    """Walk a list of folders and collect energies into rows.

    Returns a JSON list of rows. Each row has ``folder``, ``status``
    ("ok" / "missing" / "no_output"), and one entry per requested
    property. Rows with status != "ok" carry None for properties.

    Recognised properties: gibbs, zpe, single_point, scf_converged,
    opt_converged, imag_freqs, walltime_s.

    Args:
        folders: comma-separated absolute paths (or a single path).
        properties: comma-separated property names. Empty → defaults
            to "gibbs,zpe,single_point".
    """
    import json as _json
    folder_list = [f.strip() for f in folders.split(",") if f.strip()]
    prop_list = [p.strip() for p in properties.split(",") if p.strip()]
    rows = delfin_api.extract_energy_table(
        folder_list, properties=prop_list or None,
    )
    return _json.dumps(rows, indent=2)


def tool_plot_energy_distribution(
    folders: str,
    properties: str = "gibbs,single_point",
    plot_type: str = "histogram",
    title: str = "",
    bins: int = 30,
) -> str:
    """Plot energy distributions across calculations and write a PNG.

    Reads ``properties`` from each folder's largest .out, then renders:
    - histogram (default) — one panel per property with mean + median lines.
    - bar — one bar per folder per property (good for small N).
    - boxplot — distribution summary side-by-side.

    Output PNG lands in agent_workspace/ where the dashboard's inline-
    artifact hook displays it in the chat automatically. Returns JSON
    with: path, n_points, title, properties, statistics (per-property
    {n, min, max, mean, range}), error.

    Args:
        folders: comma-separated absolute paths.
        properties: comma-separated subset of gibbs/zpe/single_point.
        plot_type: histogram | bar | boxplot.
        title: figure title (auto-generated if empty).
        bins: histogram bin count (only used for plot_type=histogram).
    """
    import json as _json
    from dataclasses import asdict as _asdict
    folder_list = [f.strip() for f in folders.split(",") if f.strip()]
    prop_list = [p.strip() for p in properties.split(",") if p.strip()]
    result = delfin_api.plot_energy_distribution(
        folder_list,
        properties=prop_list or None,
        plot_type=plot_type,
        title=title,
        bins=int(bins),
    )
    return _json.dumps(_asdict(result), indent=2)


def tool_list_tools(category: str = "", query: str = "") -> str:
    """Browse the typed-tool catalog without paying for full schemas.

    Returns JSON list of {name, category, summary} entries. Use this
    BEFORE making up a tool name — the catalog is the source of truth.
    Filter by category (checks/workflow/parsing/plotting/jobs/...) or
    query (case-insensitive substring match).

    Args:
        category: optional category to filter by.
        query: optional substring to match against name/summary/category.
    """
    import json as _json
    return _json.dumps(
        delfin_api.list_tools(category=category, query=query),
        indent=2,
    )


def tool_describe_tool(name: str) -> str:
    """Return the full description + signature of one typed tool.

    Use AFTER list_tools to read the docstring for one specific tool
    before deciding to call it. Unknown name → JSON with hint.

    Args:
        name: exact tool name (case-insensitive).
    """
    import json as _json
    return _json.dumps(delfin_api.describe_tool(name), indent=2)


def tool_list_dashboard_widgets(tab: str = "") -> str:
    """Catalog of widgets the agent can drive via /ui ACTION commands.

    Returns JSON list of {name, tab, type, purpose} entries plus, for
    dropdowns, an explicit ``values`` list — no /ui options round-trip
    needed. Filter by ``tab`` (submit/orca/calc/agent) for fewer rows.

    Args:
        tab: optional tab name to filter by.
    """
    import json as _json
    return _json.dumps(
        delfin_api.list_dashboard_widgets(tab=tab),
        indent=2,
    )


def tool_get_widget_options(name: str) -> str:
    """Return the allowed values for a dropdown widget.

    Use BEFORE setting a value so /ui doesn't reject it. Empty list
    when the widget isn't a dropdown (or doesn't exist).

    Args:
        name: widget name (e.g. "orca-method").
    """
    import json as _json
    return _json.dumps(delfin_api.get_widget_options(name), indent=2)


def tool_validate_orca_input(inp_text: str) -> str:
    """Sanity-check the text of an ORCA .inp and report issues.

    Use this when the user asks "passt alles im ORCA Builder?" — read
    the orca-preview widget value first (via /ui orca-preview show),
    then pass it here. Returns JSON list of
    {severity, code, message, suggestion} entries.

    Severities: error (definitely broken), warning (probably wrong),
    info (style hint). Empty list = no obvious problems detected (NOT
    proof the calculation is correct).

    Args:
        inp_text: raw ORCA .inp file content.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    issues = delfin_api.validate_orca_input(inp_text)
    return _json.dumps([_asdict(i) for i in issues], indent=2)


def tool_submit_calculation(
    folder: str,
    job_name: str = "",
    mode: str = "delfin",
    time_limit: str = "12:00:00",
    pal: int = 12,
    maxcore: int = 6000,
    allow_mutate: bool = False,
) -> str:
    """Submit a folder via the live backend (DESTRUCTIVE — needs allow_mutate).

    Default ``allow_mutate=False`` returns a "would submit" preview
    with the full args so you can confirm with the user. Set
    allow_mutate=True ONLY after the user explicitly says yes.

    Returns JSON with: job_id, submitted (bool), folder, backend,
    message, error.

    Args:
        folder: absolute path to the job folder.
        job_name: defaults to folder basename.
        mode: "delfin" (full pipeline) or "orca" (single-step).
        time_limit: SLURM HH:MM:SS.
        pal: cores. maxcore: per-core memory in MB.
        allow_mutate: must be True for the submit to actually run.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    result = delfin_api.submit_calculation(
        folder, job_name=job_name, mode=mode, time_limit=time_limit,
        pal=int(pal), maxcore=int(maxcore),
        allow_mutate=bool(allow_mutate),
    )
    return _json.dumps(_asdict(result), indent=2)


def tool_cancel_calculation(
    job_id: str,
    allow_mutate: bool = False,
) -> str:
    """Cancel a running job by id (DESTRUCTIVE — needs allow_mutate).

    Default returns a "would cancel" hint without scancel. Confirm
    with the user, then re-call with allow_mutate=True.

    Args:
        job_id: SLURM job id (or local PID).
        allow_mutate: must be True for actual cancel.
    """
    import json as _json
    return _json.dumps(
        delfin_api.cancel_calculation(
            job_id, allow_mutate=bool(allow_mutate),
        ),
        indent=2,
    )


# ---------------------------------------------------------------------------
# Calculations-tab folder management (mutating, allow_mutate-gated)
# ---------------------------------------------------------------------------


def _split_csv(value: str) -> list[str]:
    return [v.strip() for v in (value or "").split(",") if v.strip()]


def tool_rename_calc_folder(
    src: str,
    new_name: str,
    allowed_roots: str = "",
    allow_mutate: bool = False,
) -> str:
    """Rename one calc folder in place (DESTRUCTIVE).

    new_name must be a basename (no slashes). Default-blocks until
    allow_mutate=True; first call returns a "would rename" hint with
    src/dst so you can confirm with the user.

    Args:
        src: absolute path to the folder.
        new_name: new basename (no slashes, no spaces).
        allowed_roots: comma-separated list of writable roots. Empty =
            inferred from cwd (calculations/ + calc/, never archive/).
        allow_mutate: must be True for the actual rename.
    """
    import json as _json
    return _json.dumps(
        delfin_api.rename_calc_folder(
            src, new_name,
            allowed_roots=_split_csv(allowed_roots) or None,
            allow_mutate=bool(allow_mutate),
        ),
        indent=2,
    )


def tool_create_calc_folder(
    parent: str,
    name: str,
    allowed_roots: str = "",
    allow_mutate: bool = False,
) -> str:
    """mkdir one new sub-folder under ``parent`` (DESTRUCTIVE).

    Args:
        parent: absolute path to the parent directory.
        name: new folder basename.
        allowed_roots: comma-separated writable roots. Empty = inferred.
        allow_mutate: required for actual mkdir.
    """
    import json as _json
    return _json.dumps(
        delfin_api.create_calc_folder(
            parent, name,
            allowed_roots=_split_csv(allowed_roots) or None,
            allow_mutate=bool(allow_mutate),
        ),
        indent=2,
    )


def tool_move_calc_folder(
    src: str,
    dst_parent: str,
    allowed_roots: str = "",
    allow_mutate: bool = False,
) -> str:
    """Move one calc folder into another calc/ location (DESTRUCTIVE).

    Use for re-organizing the calculations tree. For sending things to
    archive use ``move_to_archive`` instead — it enforces direction.

    Args:
        src: absolute path to the folder being moved.
        dst_parent: absolute path to the new parent directory.
        allowed_roots: comma-separated writable roots. Empty = inferred.
        allow_mutate: must be True for the actual move.
    """
    import json as _json
    return _json.dumps(
        delfin_api.move_calc_folder(
            src, dst_parent,
            allowed_roots=_split_csv(allowed_roots) or None,
            allow_mutate=bool(allow_mutate),
        ),
        indent=2,
    )


def tool_move_to_archive(
    src: str,
    archive_root: str = "",
    allowed_roots: str = "",
    allow_mutate: bool = False,
) -> str:
    """Move a folder calc/ -> archive/ (DESTRUCTIVE, direction enforced).

    Args:
        src: absolute path to the folder being archived.
        archive_root: absolute path to the archive directory. Empty =
            uses ``./archive`` next to the cwd.
        allowed_roots: comma-separated source-side roots. Empty = inferred.
        allow_mutate: required for the actual move.
    """
    import json as _json
    return _json.dumps(
        delfin_api.move_to_archive(
            src, archive_root=archive_root,
            allowed_roots=_split_csv(allowed_roots) or None,
            allow_mutate=bool(allow_mutate),
        ),
        indent=2,
    )


def tool_delete_calc_folder(
    folder: str,
    confirm_token: str = "",
    allowed_roots: str = "",
    allow_mutate: bool = False,
) -> str:
    """Permanently delete a calc folder (DESTRUCTIVE — 3-lock gate).

    Three locks: ``allow_mutate=True``, target inside calc roots
    (NEVER archive/), and ``confirm_token`` MUST equal the folder's
    basename verbatim. Missing any lock → safe refusal.

    Args:
        folder: absolute path to the folder.
        confirm_token: must equal the folder's basename.
        allowed_roots: comma-separated writable roots.
        allow_mutate: must be True for actual rmtree.
    """
    import json as _json
    return _json.dumps(
        delfin_api.delete_calc_folder(
            folder, confirm_token=confirm_token,
            allowed_roots=_split_csv(allowed_roots) or None,
            allow_mutate=bool(allow_mutate),
        ),
        indent=2,
    )


def tool_read_pdf(
    path: str,
    pages: str = "",
    max_chars: int = 50000,
) -> str:
    """Read plain text from a PDF without indexing it first.

    Use for ad-hoc PDFs the user dropped in the Literature tab when
    you don't need persistent search. For repeated lookups, prefer
    index_new_pdf + search_docs.

    Args:
        path: absolute path to the PDF.
        pages: optional 1-based page selector — "5" / "3-7" /
            "1,5-10,15". Empty → entire document.
        max_chars: cap on returned text size (default 50000).

    Returns JSON: {path, n_pages_total, n_pages_read, text,
    truncated, error}.
    """
    import json as _json
    return _json.dumps(
        delfin_api.read_pdf(path, pages=pages, max_chars=int(max_chars)),
        indent=2,
    )


def tool_search_pdf_local(
    path: str,
    query: str,
    context_lines: int = 3,
    max_hits: int = 20,
    case_sensitive: bool = False,
) -> str:
    """Substring-search inside ONE PDF; return matching paragraphs.

    Cheaper than indexing for one-off lookups. Each hit carries
    surrounding context so the agent can quote the relevant passage.

    Args:
        path: absolute path to the PDF.
        query: substring to look for.
        context_lines: lines above/below each hit (default 3).
        max_hits: cap on returned hits (default 20).
        case_sensitive: default False.

    Returns JSON: {path, query, hits: [{page, line, snippet}], error}.
    """
    import json as _json
    return _json.dumps(
        delfin_api.search_pdf_local(
            path, query,
            context_lines=int(context_lines),
            max_hits=int(max_hits),
            case_sensitive=bool(case_sensitive),
        ),
        indent=2,
    )


def tool_extract_pdf_section(
    path: str,
    heading: str,
    max_chars: int = 8000,
) -> str:
    """Pull a single section from a PDF by heading.

    Searches for the heading line, then returns everything from there
    until the next heading-like line. Useful when search_docs returned
    a hit and you want the full section without round-tripping through
    read_section.

    Args:
        path: absolute path to the PDF.
        heading: heading text to match (case-insensitive substring).
        max_chars: cap on returned text (default 8000).

    Returns JSON: {path, heading_found, text, next_heading, error}.
    """
    import json as _json
    return _json.dumps(
        delfin_api.extract_pdf_section(
            path, heading, max_chars=int(max_chars),
        ),
        indent=2,
    )


def tool_list_literature_files(folder: str = "") -> str:
    """List PDFs / MDs / TXTs in the Literature folder.

    Empty folder → uses the Literature directory the indexer
    auto-detects (next to the DELFIN repo or ~/literature). Recursive.

    Returns JSON list of {path, name, size_bytes, mtime_iso, ext}.

    Args:
        folder: optional explicit folder to scan.
    """
    import json as _json
    return _json.dumps(
        delfin_api.list_literature_files(folder=folder),
        indent=2,
    )


def tool_check_orca_manual_indexed() -> str:
    """Quick-check whether the ORCA manual is in the doc-search index.

    ALWAYS call this FIRST when the user asks an ORCA-specific question
    (keyword syntax, %block, methodology, basis pairing). If the
    ``indexed`` field comes back False, surface the ``hint`` to the
    user verbatim — it asks them to drop the manual into the
    Literature tab and call index_new_pdf.

    Returns JSON: {indexed, doc_ids, hint}.
    """
    import json as _json
    return _json.dumps(delfin_api.check_orca_manual_indexed(), indent=2)


def tool_index_new_pdf(
    path: str,
    doc_id: str = "",
    title: str = "",
) -> str:
    """Add a PDF to the doc-search index.

    Use after the user drops a new manual / paper into the Literature
    tab so subsequent search_docs calls find it. Re-builds the
    existing index plus the new PDF.

    Returns JSON: {ok, doc_id, title, sections_indexed, index_path,
    error}.

    Args:
        path: absolute path to the PDF on disk.
        doc_id: explicit identifier (defaults to filename stem).
        title: human-readable title (defaults to filename stem).
    """
    import json as _json
    return _json.dumps(
        delfin_api.index_new_pdf(
            path, doc_id=doc_id, title=title,
        ),
        indent=2,
    )


def tool_list_delfin_features(category: str = "") -> str:
    """Browse the curated catalog of DELFIN concepts/features.

    Use before explain_delfin_feature when you're not sure of the
    exact concept name. Each entry: {name, category, summary
    (truncated to 140 chars)}.

    Args:
        category: optional filter (config / workflow / module /
            agent / methodology).
    """
    import json as _json
    return _json.dumps(
        delfin_api.list_delfin_features(category=category),
        indent=2,
    )


def tool_explain_delfin_feature(name: str) -> str:
    """Explain a DELFIN concept (CONTROL keys, Smart Recalc, OCCUPIER, …).

    Use when the user asks "wie funktioniert X in DELFIN?" or "was
    macht <feature>?". Returns curated prose + source-file pointers
    so you can read deeper if needed. Unknown name → JSON with
    candidates and available list.

    Available concepts: control, control_keys, relativistic_methods,
    pipeline, smart_recalc, occupier, guppy, csp, mlp, modes,
    permissions, co2, tadf_xtb, hyperpol, outcomes, session_boot,
    live_state.

    Args:
        name: concept name (case-insensitive; "-" / " " also ok).
    """
    import json as _json
    return _json.dumps(
        delfin_api.explain_delfin_feature(name),
        indent=2,
    )


def tool_list_active_calculations() -> str:
    """Live list of running/pending jobs (read-only).

    Returns JSON list of {job_id, name, status, runtime_s, directory}
    entries. Empty when no jobs are active.
    """
    import json as _json
    return _json.dumps(delfin_api.list_active_calculations(), indent=2)


def tool_plot_energy_correlation(
    folders: str,
    x: str = "single_point",
    y: str = "gibbs",
    title: str = "",
) -> str:
    """Scatter plot one energy property against another across folders.

    Pearson correlation + linear best-fit line are drawn on top so the
    user can see at a glance whether ``y`` tracks ``x`` linearly.
    Returns JSON with path, n_points, title, properties=[x, y], and
    statistics={n, pearson_r, slope, intercept}.

    Args:
        folders: comma-separated absolute paths.
        x: gibbs | zpe | single_point.
        y: gibbs | zpe | single_point.
        title: figure title (auto-generated if empty).
    """
    import json as _json
    from dataclasses import asdict as _asdict
    folder_list = [f.strip() for f in folders.split(",") if f.strip()]
    result = delfin_api.plot_energy_correlation(
        folder_list, x=x, y=y, title=title,
    )
    return _json.dumps(_asdict(result), indent=2)


def tool_find_calculation_extreme(
    folders: str,
    property: str = "gibbs",
    extreme: str = "min",
    n: int = 5,
) -> str:
    """Return the N folders with the lowest/highest value of a property.

    Direct answer to "find the .out with the lowest Gibbs energy"
    type questions. Folders that fail to parse the property are
    excluded from the ranking, so a clean list is returned.

    Args:
        folders: comma-separated absolute paths.
        property: gibbs (default) | zpe | single_point | imag_freqs |
            walltime_s.
        extreme: "min" (lowest, default) or "max" (highest).
        n: how many top entries to return (default 5).
    """
    import json as _json
    folder_list = [f.strip() for f in folders.split(",") if f.strip()]
    rows = delfin_api.find_calculation_extreme(
        folder_list, property=property, extreme=extreme, n=int(n),
    )
    return _json.dumps(rows, indent=2)


def tool_extract_imaginary_frequencies(folder: str) -> str:
    """List imaginary modes (n_imag, mode_index, frequency_cm) for one folder.

    Reads the largest ``*.out`` containing a VIBRATIONAL FREQUENCIES
    block and returns:

    - ``n_imag``: count of imaginary modes
    - ``modes``: list of {mode_index, frequency_cm}
    - ``most_negative``: most negative cm**-1 value (None if no imag)
    - ``is_minimum`` (n_imag == 0), ``is_ts`` (n_imag == 1)
    - ``error``: filled when no Freq output is present

    Use this when the user asks "ist X ein Minimum?" / "TS suchen" /
    "imaginäre Frequenzen vergleichen". Combine with
    ``compare_across_functionals`` to scan many folders.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.extract_imaginary_frequencies(folder)),
        indent=2,
    )


def tool_compare_calculations(folder_a: str, folder_b: str) -> str:
    """Side-by-side diff of two calculation folders.

    Returns: method/basis match flags, both functionals + bases,
    Gibbs/SPE for each, ``delta_*_kcal`` (B - A in kcal/mol), imag-freq
    counts, and human notes (e.g. "different functional", "B has 1 imag
    freq — not a minimum").

    Use this when the user wants to compare two single calculations
    head-to-head. For more than 2 folders use
    ``compare_across_functionals``.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.compare_calculations(folder_a, folder_b)),
        indent=2,
    )


def tool_compare_across_functionals(
    folders: str,
    include_imag: bool = True,
    sort_by: str = "gibbs",
) -> str:
    """Multi-folder comparison table grouped by functional/basis.

    Returns one row per folder with: functional, basis, gibbs,
    single_point, zpe, n_imag, is_minimum, status. Sortable by gibbs
    (default), single_point, zpe, functional, or folder.

    Direct answer to "vergleiche imaginäre Frequenzen über Funktionale"
    or "Welcher Functional gibt das tiefste Minimum?".

    Args:
        folders: comma-separated absolute paths.
        include_imag: if True (default), also extract imaginary-freq
            counts (slightly slower but usually wanted).
        sort_by: gibbs | single_point | zpe | functional | folder.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    folder_list = [f.strip() for f in folders.split(",") if f.strip()]
    rows = delfin_api.compare_across_functionals(
        folder_list, include_imag=include_imag, sort_by=sort_by,
    )
    return _json.dumps([_asdict(r) for r in rows], indent=2)


def tool_stop_dry_run(workspace: str) -> str:
    """List DELFIN processes that would be signaled (no actual signal sent)."""
    rc = delfin_api.stop(workspace=workspace, dry_run=True)
    return _format_result(rc, action="stop_dry_run", dry_run=True)


# ---------------------------------------------------------------------------
# Mutating tool implementations (require allow_mutate=True)
# ---------------------------------------------------------------------------

def tool_cleanup(
    orca: bool = False,
    dry_run: bool = True,
    workspace: str = "",
    scratch: str = "",
    allow_mutate: bool = False,
) -> str:
    """Remove DELFIN scratch artifacts.

    Defaults to dry_run=True for safety.  To actually delete, pass
    dry_run=False AND allow_mutate=True.
    """
    if not dry_run and not allow_mutate:
        return _refuse_mutation("cleanup")
    rc = delfin_api.cleanup(
        orca=orca, dry_run=dry_run, workspace=workspace,
        scratch=scratch or None,
    )
    return _format_result(rc, action="cleanup", dry_run=dry_run)


def tool_stop(
    signal_name: str = "INT",
    workspace: str = "",
    dry_run: bool = True,
    cleanup_after: bool = False,
    wait_seconds: float = 3.0,
    allow_mutate: bool = False,
) -> str:
    """Signal running DELFIN processes."""
    if not dry_run and not allow_mutate:
        return _refuse_mutation("stop")
    rc = delfin_api.stop(
        workspace=workspace, signal_name=signal_name, dry_run=dry_run,
        cleanup_after=cleanup_after, wait_seconds=wait_seconds,
    )
    return _format_result(rc, action="stop", dry_run=dry_run)


def tool_pipeline_prepare(
    control_file: str = "CONTROL.txt",
    overwrite: bool = False,
    allow_mutate: bool = False,
) -> str:
    """Generate a CONTROL.txt template (``delfin --define``)."""
    if not allow_mutate:
        return _refuse_mutation("pipeline_prepare")
    rc = delfin_api.pipeline_prepare(control_file=control_file, overwrite=overwrite)
    return _format_result(rc, action="pipeline_prepare")


def tool_pipeline_run(
    control_file: str = "CONTROL.txt",
    cleanup: bool = True,
    recalc: bool = False,
    overwrite: bool = False,
    define: str = "",
    extra_args: str = "",
    allow_mutate: bool = False,
) -> str:
    """Run the full DELFIN pipeline."""
    if not allow_mutate:
        return _refuse_mutation("pipeline_run")
    extras = extra_args.split() if extra_args else None
    rc = delfin_api.pipeline_run(
        control_file=control_file, cleanup=cleanup, recalc=recalc,
        overwrite=overwrite, define=define or None, extra_args=extras,
    )
    return _format_result(rc, action="pipeline_run")


def tool_run_orca_input(
    input_file: str = "",
    output: str = "",
    allow_mutate: bool = False,
) -> str:
    """Run ORCA on a .inp file via DELFIN's ORCA resolver."""
    if not allow_mutate:
        return _refuse_mutation("run_orca_input")
    rc = delfin_api.run_orca_input(
        input_file=input_file or None, output=output or None,
    )
    return _format_result(rc, action="run_orca_input")


def tool_co2(
    define: bool = False,
    force: bool = False,
    recalc: bool = False,
    charge: int = 0,
    multiplicity: int = 0,
    solvent: str = "",
    metal: str = "",
    broken_sym: str = "",
    allow_mutate: bool = False,
) -> str:
    """CO2 Coordinator workflow."""
    if not allow_mutate:
        return _refuse_mutation("co2")
    rc = delfin_api.co2(
        define=define, force=force, recalc=recalc,
        charge=charge if charge else None,
        multiplicity=multiplicity if multiplicity else None,
        solvent=solvent or None, metal=metal or None,
        broken_sym=broken_sym or None,
    )
    return _format_result(rc, action="co2")


def tool_tadf_xtb(extra_args: str = "", allow_mutate: bool = False) -> str:
    """Run the TADF xTB workflow."""
    if not allow_mutate:
        return _refuse_mutation("tadf_xtb")
    extras = extra_args.split() if extra_args else None
    rc = delfin_api.tadf_xtb(extra_args=extras)
    return _format_result(rc, action="tadf_xtb")


def tool_hyperpol(extra_args: str = "", allow_mutate: bool = False) -> str:
    """Run the hyperpolarisability workflow."""
    if not allow_mutate:
        return _refuse_mutation("hyperpol")
    extras = extra_args.split() if extra_args else None
    rc = delfin_api.hyperpol(extra_args=extras)
    return _format_result(rc, action="hyperpol")


# ---------------------------------------------------------------------------
# Server bootstrap
# ---------------------------------------------------------------------------

def run_server(argv: list[str] | None = None) -> None:
    """Start the MCP server on stdio."""
    parser = argparse.ArgumentParser(prog="delfin-ops-server")
    parser.add_argument(
        "--workspace",
        default=os.getcwd(),
        help="Default workspace directory used when tools take a workspace arg.",
    )
    args = parser.parse_args(argv)
    default_workspace = args.workspace

    from mcp.server.fastmcp import FastMCP

    mcp = FastMCP(
        "delfin-ops",
        instructions=(
            "DELFIN operations server. Run DELFIN workflows (pipeline, "
            "ORCA, TADF, hyperpol, CO2), inspect tool availability, and "
            "control running processes. Read-only tools are always safe; "
            "mutating tools require allow_mutate=True and should be "
            "confirmed with the user first."
        ),
    )

    # Read-only — register module functions directly
    mcp.tool()(tool_qm_check)
    mcp.tool()(tool_csp_check)
    mcp.tool()(tool_mlp_check)
    mcp.tool()(tool_analysis_check)
    mcp.tool()(tool_list_dashboard_patterns)
    mcp.tool()(tool_get_dashboard_pattern)
    # P1 — output parsing (read-only, structured returns)
    mcp.tool()(tool_parse_orca_output)
    mcp.tool()(tool_find_orca_errors)
    mcp.tool()(tool_extract_thermochem)
    mcp.tool()(tool_extract_energy_table)
    mcp.tool()(tool_find_calculation_extreme)
    # Imaginary-frequency + functional-comparison helpers
    mcp.tool()(tool_extract_imaginary_frequencies)
    mcp.tool()(tool_compare_calculations)
    mcp.tool()(tool_compare_across_functionals)
    # P1 — statistical plots (PNG → agent_workspace, auto-displayed)
    mcp.tool()(tool_plot_energy_distribution)
    mcp.tool()(tool_plot_energy_correlation)
    # Tool / widget catalogs (cheap on-demand discovery)
    mcp.tool()(tool_list_tools)
    mcp.tool()(tool_describe_tool)
    mcp.tool()(tool_list_dashboard_widgets)
    mcp.tool()(tool_get_widget_options)
    # ORCA Builder validation
    mcp.tool()(tool_validate_orca_input)
    # Job lifecycle (read-only list + mutating submit/cancel)
    mcp.tool()(tool_list_active_calculations)
    mcp.tool()(tool_submit_calculation)
    mcp.tool()(tool_cancel_calculation)
    # Calc folder management (mutating, allow_mutate-gated)
    mcp.tool()(tool_rename_calc_folder)
    mcp.tool()(tool_create_calc_folder)
    mcp.tool()(tool_move_calc_folder)
    mcp.tool()(tool_move_to_archive)
    mcp.tool()(tool_delete_calc_folder)
    # ORCA-manual lookup + literature indexing
    mcp.tool()(tool_check_orca_manual_indexed)
    mcp.tool()(tool_index_new_pdf)
    # PDF on-demand reading (no pre-indexing)
    mcp.tool()(tool_read_pdf)
    mcp.tool()(tool_search_pdf_local)
    mcp.tool()(tool_extract_pdf_section)
    mcp.tool()(tool_list_literature_files)
    # DELFIN-feature explainer
    mcp.tool()(tool_list_delfin_features)
    mcp.tool()(tool_explain_delfin_feature)

    # stop_dry_run needs the default workspace closed over
    @mcp.tool(name="stop_dry_run", description=tool_stop_dry_run.__doc__)
    def _stop_dry_run(workspace: str = "") -> str:
        return tool_stop_dry_run(workspace or default_workspace)

    # Mutating wrappers default workspace = server cwd
    @mcp.tool(name="cleanup", description=tool_cleanup.__doc__)
    def _cleanup(
        orca: bool = False,
        dry_run: bool = True,
        workspace: str = "",
        scratch: str = "",
        allow_mutate: bool = False,
    ) -> str:
        return tool_cleanup(
            orca=orca, dry_run=dry_run,
            workspace=workspace or default_workspace,
            scratch=scratch, allow_mutate=allow_mutate,
        )

    @mcp.tool(name="stop", description=tool_stop.__doc__)
    def _stop(
        signal_name: str = "INT",
        workspace: str = "",
        dry_run: bool = True,
        cleanup_after: bool = False,
        wait_seconds: float = 3.0,
        allow_mutate: bool = False,
    ) -> str:
        return tool_stop(
            signal_name=signal_name,
            workspace=workspace or default_workspace,
            dry_run=dry_run, cleanup_after=cleanup_after,
            wait_seconds=wait_seconds, allow_mutate=allow_mutate,
        )

    mcp.tool()(tool_pipeline_prepare)
    mcp.tool()(tool_pipeline_run)
    mcp.tool()(tool_run_orca_input)
    mcp.tool()(tool_co2)
    mcp.tool()(tool_tadf_xtb)
    mcp.tool()(tool_hyperpol)

    mcp.run(transport="stdio")
