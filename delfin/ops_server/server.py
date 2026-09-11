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

Mutating (require the host's grant):
    - ``cleanup``          — remove scratch artifacts (supports dry_run)
    - ``stop``             — signal DELFIN processes
    - ``pipeline_run``     — full DELFIN pipeline
    - ``pipeline_prepare`` — generate CONTROL template
    - ``run_orca_input``   — run ORCA on a .inp file
    - ``co2``              — CO2 Coordinator workflow
    - ``tadf_xtb``         — TADF xTB workflow
    - ``hyperpol``         — hyperpolarisability workflow

Mutating tools refuse to execute unless the process that started this server
set ``DELFIN_OPS_ALLOW_MUTATE``.  The consent used to be an ``allow_mutate``
parameter ON the mutating tool, which meant the caller being gated supplied
its own permission — self-attested consent, invisible to every gate
upstream.  The MCP-facing wrappers therefore drop the parameter from the
schema entirely: there is nothing for a caller to pass.  The Python
functions keep the keyword so in-process callers (dashboard, tests) can
still pass a decision they actually made.

The tool functions are defined at module level so they can be imported and
tested without the optional MCP SDK installed.  ``run_server`` only loads
the SDK's server class when actually starting the stdio server, and reaches
it through ``delfin.mcp_compat`` because the class changed module and name
between the SDK's 1.x and 2.x lines.
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
    return json.dumps(payload, separators=(",", ":"), ensure_ascii=False)


def _refuse_mutation(action: str) -> str:
    """Standard refusal payload when the host has not granted mutation.

    The message deliberately does NOT tell the caller to pass a flag. It
    used to, and that was the whole gate: ``allow_mutate`` was a parameter
    on the tool the model was calling, so the consent for
    delete-calc-folder, kill-all-user-jobs, pipeline-run and move-to-archive
    was supplied by the same party the consent was protecting against, and
    the permission gate upstream never saw it. Consent now comes from the
    process that STARTED this server.
    """
    return json.dumps({
        "action": action,
        "ok": False,
        "error": "mutation_blocked",
        "message": (
            f"The tool '{action}' modifies state and this server was not "
            "started with mutation granted. Ask the USER to approve the "
            f"action; there is no argument that turns it on."
        ),
    }, separators=(",", ":"), ensure_ascii=False)


# Environment variable the host sets when the user has granted this server
# permission to change things. Read at call time, never taken from the tool
# arguments, so nothing the model emits can set it.
_MUTATION_ENV = "DELFIN_OPS_ALLOW_MUTATE"


def host_grants_mutation() -> bool:
    """Whether the process that started this server granted mutation.

    Out of band by construction: an MCP tool argument travels with the
    model's message, so a boolean parameter named ``allow_mutate`` is the
    model attesting its own consent. This is read from the environment the
    host controls.
    """
    raw = (os.environ.get(_MUTATION_ENV) or "").strip().lower()
    return raw in ("1", "true", "yes", "on")


def _host_gated(fn):
    """Expose *fn* without its ``allow_mutate`` parameter.

    The wrapper keeps the same call, drops the flag from the schema the
    model sees, and supplies the value from :func:`host_grants_mutation`.
    Anything the model sends under that name is discarded.
    """
    import functools
    import inspect

    sig = inspect.signature(fn)
    params = [p for n, p in sig.parameters.items() if n != "allow_mutate"]

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        kwargs.pop("allow_mutate", None)
        return fn(*args, allow_mutate=host_grants_mutation(), **kwargs)

    wrapper.__signature__ = sig.replace(parameters=params)
    wrapper.__doc__ = _strip_allow_mutate_doc(fn.__doc__)
    return wrapper


def _strip_allow_mutate_doc(doc: str | None) -> str:
    """Drop any ``allow_mutate`` line from a docstring shown to the model.

    A description that tells the model to pass a flag it no longer has
    teaches it to retry with an argument that is thrown away.
    """
    if not doc:
        return ""
    kept = [ln for ln in doc.splitlines()
            if "allow_mutate" not in ln]
    return "\n".join(kept).rstrip()


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


def _dumps(obj) -> str:
    """Every tool result, as compact JSON.

    The engine caps a tool result at 5000 chars for a model whose profile
    sets no larger cap, and cuts the middle out of anything longer. Over
    the nine-run archive fixture the pretty-printed ranking was 6244
    chars and the comparison 5863: the model read a table with its
    middle missing, marked as truncated (2026-09-11). Compact, the same
    tables are 4338 and 3871. Indentation is for people; the reader here
    is a model, and every char of it is a token it pays for.
    """
    import json as _json
    return _json.dumps(obj, ensure_ascii=False, separators=(",", ":"))


def _safe(name: str, fn):
    """A tool that raises answers with the reason, as JSON.

    The MCP layer turns an uncaught exception into 'Error executing tool
    <name>' and nothing else; a model that got that for a figure it was
    told to draw could not tell a missing directory from a missing
    library, and drew the figure by hand (2026-09-11). The signature and
    docstring are the wrapped function's, so the schema the model sees
    is unchanged.
    """
    import functools

    @functools.wraps(fn)
    def _call(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except Exception as exc:               # noqa: BLE001 - reported, not hidden
            return _dumps({"error": f"{type(exc).__name__}: {exc}", "tool": name})
    return _call


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
        "status": getattr(parsed, "status", "ok"),
        "outcome": getattr(parsed, "outcome", ""),
    }


def tool_parse_orca_output(path: str) -> str:
    """Parse one ORCA .out file and return a structured snapshot.

    Returns a JSON object with: final_single_point (Hartree),
    gibbs_free_energy (Hartree), zpe (Hartree), scf_converged (bool),
    opt_converged (bool), imag_freq_count (int), walltime_s (float),
    n_atoms (int), functional (str), basis (str), error_summary (str),
    status ("ok" | "no_output" | "missing" | "read_error") and outcome
    (how the run ended, e.g. "failed (exit code 1025)" -- the same
    phrase calc_status gives). status says whether there was an output
    to parse at all; "no_output" is a run that has not written yet, not
    a parse failure. functional/basis come from the output when it
    states them, else from the folder (DELFIN_Data.json, CONTROL.txt,
    .inp) -- a DELFIN run's output does not name its method. Missing
    values are null. Use this BEFORE writing a Python script to grep
    the file — one tool call replaces dozens of regexes.

    Args:
        path: absolute path to the ORCA .out file, or to the calculation
            folder (its largest .out is parsed).
    """
    parsed = delfin_api.parse_orca_output(path)
    return _dumps(_orca_parse_to_dict(parsed))


def tool_calc_status(folder: str) -> str:
    """Did this calculation succeed, fail, or is it still running -- with the evidence.

    Returns {"folder", "state", "outcome", "method", "evidence",
    "last_activity", "last_activity_age_s"}.
    state: succeeded | failed | finished | running | stalled | pending |
    unknown | missing. "stalled" is a run with no exit code whose files
    have not been written for over six hours (last_activity_age_s says
    how long); "running" is one written recently. outcome: the phrase
    with its source ("failed (exit code
    1025)", "running or crashed (run log present, no exit code)", "no
    output yet (input present; not started or still running)",
    "finished per ORCA output (no exit code file)"). evidence: every
    file that had a say -- the exit-code marker, the run log's last
    line, the state file's status, the output's termination line -- and
    what it said, so the answer can be cited. One call per folder; for
    a table over many folders use extract_energy_table, whose outcome
    column is the same phrase. list_active_calculations is NOT this: it
    lists the scheduler's jobs, not the state of a folder on disk.

    Args:
        folder: absolute path to the calculation folder.
    """
    from dataclasses import asdict as _asdict
    return _dumps(_asdict(delfin_api.calculation_status(folder)))


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
    from dataclasses import asdict as _asdict
    errors = delfin_api.find_orca_errors(folder)
    return _dumps([_asdict(e) for e in errors])


def tool_extract_thermochem(folder: str) -> str:
    """Extract the full thermochemistry block from an ORCA Freq output.

    Picks the first .out in ``folder`` containing thermochemistry data.
    Returns JSON with: temperature_k, pressure_atm, zpe, thermal_corr,
    enthalpy_corr, entropy_total, gibbs_corr, final_gibbs (all in
    Hartree except T and P).

    Args:
        folder: absolute path to a calc folder.
    """
    from dataclasses import asdict as _asdict
    result = delfin_api.extract_thermochem(folder)
    return _dumps(_asdict(result))


def _common_root(paths: list) -> str:
    """The directory every path shares, or "" when there is none."""
    import os as _os
    clean = [str(p) for p in paths if p]
    if not clean:
        return ""
    try:
        root = _os.path.commonpath(clean)
    except ValueError:
        return ""
    if root and any(_os.path.normpath(p) == root for p in clean):
        root = _os.path.dirname(root)
    return root


def _relative_to_root(rows: list, root: str) -> list:
    """Rows with their folder spelled relative to *root*.

    An absolute path repeated once per row and once per skipped folder
    was what pushed a nine-run table past the tool-result cap -- and
    by how much depended on where the archive happened to live. The
    root is said once; a row names its folder the way a person does.
    """
    import os as _os
    if not root:
        return rows
    out = []
    for row in rows:
        if isinstance(row, dict) and row.get("folder"):
            row = dict(row)
            try:
                row["folder"] = _os.path.relpath(str(row["folder"]), root)
            except ValueError:
                pass
        out.append(row)
    return out


def _grouped_by_method(rows: list) -> dict:
    """Rows -> {"note", "root", "groups": [{"method", "rows"}]}, in row
    order, folders relative to the shared root.

    The rule travels with the data: a caller that reads the JSON reads
    that energies compare only within a method before it reads a number.
    """
    root = _common_root([r.get("folder") for r in rows if isinstance(r, dict)])
    groups: list = []
    index: dict = {}
    for row in _relative_to_root(rows, root):
        method = row.get("method") if isinstance(row, dict) else None
        if method not in index:
            index[method] = len(groups)
            groups.append({"method": method, "rows": []})
        groups[index[method]]["rows"].append(row)
    return {"note": delfin_api.METHOD_NOTE, "root": root, "groups": groups}


def tool_extract_energy_table(
    folders: str,
    properties: str = "",
) -> str:
    """Walk a list of folders and collect energies into rows.

    Returns a JSON list of rows. Each row has ``folder``, ``status``
    ("ok" / "missing" / "no_output"), ``method`` ("PBE0/def2-SVP/DMF":
    functional, dispersion, basis and solvent -- a total energy compares
    only within one method, and the gas phase is not a solvent), ``outcome``
    (succeeded / failed (exit code N) / running or crashed / unknown --
    "no_output" alone does not say which), ``state`` (one word by the
    same rule calc_status uses: succeeded / failed / finished / running /
    stalled / pending / unknown -- "running" when the scheduler has a
    job on the folder or it was written recently, "stalled" after six
    hours without a write), ``last_activity`` and ``last_activity_age_s``
    (when the newest file in the folder was written), and one entry per
    requested property. Rows with status != "ok" carry None for properties.

    Recognised properties: gibbs, zpe, single_point, scf_converged,
    opt_converged, imag_freqs, walltime_s.

    Args:
        folders: comma-separated absolute paths (or a single path).
        properties: comma-separated property names. Empty → defaults
            to "gibbs,zpe,single_point".
    
    """
    folder_list = [f.strip() for f in folders.split(",") if f.strip()]
    prop_list = [p.strip() for p in properties.split(",") if p.strip()]
    rows = delfin_api.extract_energy_table(
        folder_list, properties=prop_list or None,
    )
    return _dumps(rows)


def tool_plot_energy_distribution(
    folders: str,
    properties: str = "gibbs,single_point",
    plot_type: str = "histogram",
    title: str = "",
    bins: int = 30,
    output_path: str = "",
) -> str:
    """Plot energy distributions across calculations and write a PNG.

    Reads ``properties`` from each folder's largest .out, then renders:
    - histogram (default) — one panel per property with mean + median lines.
    - bar — one bar per folder per property (good for small N).
    - boxplot — distribution summary side-by-side.
    - bar_by_method — one bar per folder, grouped and coloured by method
      (functional/dispersion/basis/solvent): the figure for "which run is
      lowest", since a total energy compares only within a method.
      statistics carries the lowest per group and the excluded folders.

    Output PNG lands in agent_workspace/ where the dashboard's inline-
    artifact hook displays it in the chat automatically. Returns JSON
    with: path, n_points, title, properties, statistics (per-property
    {n, min, max, mean, range}), error.

    Args:
        folders: comma-separated absolute paths.
        properties: comma-separated subset of gibbs/zpe/single_point.
        plot_type: histogram | bar | boxplot | bar_by_method.
        title: figure title (auto-generated if empty).
        bins: histogram bin count (only used for plot_type=histogram).
        output_path: where to write the PNG (default: the agent workspace).
    """
    from dataclasses import asdict as _asdict
    folder_list = [f.strip() for f in folders.split(",") if f.strip()]
    prop_list = [p.strip() for p in properties.split(",") if p.strip()]
    result = delfin_api.plot_energy_distribution(
        folder_list,
        properties=prop_list or None,
        plot_type=plot_type,
        title=title,
        bins=int(bins), output_path=output_path)
    return _dumps(_asdict(result))


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
        separators=(",", ":"),
    )


def tool_describe_tool(name: str) -> str:
    """Return the full description + signature of one typed tool.

    Use AFTER list_tools to read the docstring for one specific tool
    before deciding to call it. Unknown name → JSON with hint.

    Args:
        name: exact tool name (case-insensitive).
    """
    return _dumps(delfin_api.describe_tool(name))


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
        separators=(",", ":"),
    )


def tool_get_widget_options(name: str) -> str:
    """Return the allowed values for a dropdown widget.

    Use BEFORE setting a value so /ui doesn't reject it. Empty list
    when the widget isn't a dropdown (or doesn't exist).

    Args:
        name: widget name (e.g. "orca-method").
    """
    return _dumps(delfin_api.get_widget_options(name))


def tool_validate_orca_input(inp_text: str) -> str:
    """Sanity-check the text of an ORCA .inp and report issues.

    Use this when the user asks "is everything OK in the ORCA Builder?" — read
    the orca-preview widget value first (via /ui orca-preview show),
    then pass it here. Returns JSON list of
    {severity, code, message, suggestion} entries.

    Severities: error (definitely broken), warning (probably wrong),
    info (style hint). Empty list = no obvious problems detected (NOT
    proof the calculation is correct).

    Args:
        inp_text: raw ORCA .inp file content.
    """
    from dataclasses import asdict as _asdict
    issues = delfin_api.validate_orca_input(inp_text)
    return _dumps([_asdict(i) for i in issues])


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
    from dataclasses import asdict as _asdict
    result = delfin_api.submit_calculation(
        folder, job_name=job_name, mode=mode, time_limit=time_limit,
        pal=int(pal), maxcore=int(maxcore),
        allow_mutate=bool(allow_mutate),
    )
    return _dumps(_asdict(result))


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
        separators=(",", ":"),
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
        separators=(",", ":"),
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
        separators=(",", ":"),
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
        separators=(",", ":"),
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
        separators=(",", ":"),
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
        separators=(",", ":"),
    )


# ---------------------------------------------------------------------------
# Bulk job control + recalc preparation + (Options) dispatcher
# ---------------------------------------------------------------------------


def tool_list_ssh_transfer_jobs(limit: int = 8) -> str:
    """List queued/running/finished SSH transfer jobs (read-only).

    Returns up to ``limit`` most-recently-updated entries. The actual
    ``run_transfer_job`` step is dashboard-only — drive it via
    ACTION: /calc options ssh-transfer.

    Args:
        limit: max entries (default 8, most recent first).
    """
    import json as _json
    return _json.dumps(
        delfin_api.list_ssh_transfer_jobs(limit=int(limit)),
        separators=(",", ":"),
    )


def tool_kill_all_user_jobs(
    only_running: bool = False,
    allow_mutate: bool = False,
) -> str:
    """Cancel every active job for the current user (DESTRUCTIVE).

    Default-blocks: returns a list of job_ids that WOULD be cancelled
    so you can confirm with the user. Pass ``only_running=True`` to
    skip PENDING jobs.

    Args:
        only_running: skip PENDING jobs.
        allow_mutate: required for actual scancel.
    """
    import json as _json
    return _json.dumps(
        delfin_api.kill_all_user_jobs(
            only_running=bool(only_running),
            allow_mutate=bool(allow_mutate),
        ),
        separators=(",", ":"),
    )


def tool_prepare_recalc(
    folder: str,
    mode: str = "smart",
    time_limit: str = "24:00:00",
    override: str = "",
    allow_mutate: bool = False,
) -> str:
    """Pre-flight + submit a Smart-Recalc / classic Recalc / Override.

    Reads PAL/maxcore from CONTROL.txt or the largest .inp, validates
    the time-limit format, and either returns a dry-run plan
    (allow_mutate=False) or submits via the live backend.

    Args:
        folder: absolute path to the calc folder.
        mode: 'smart' (delfin-recalc), 'classic' (delfin-recalc-classic),
            or 'override' (delfin-recalc-override; needs override=).
        time_limit: SLURM time spec HH:MM:SS.
        override: only used when mode='override' (STAGE=INDEX[,...]).
        allow_mutate: required for actual submission.
    """
    from dataclasses import asdict as _asdict
    plan = delfin_api.prepare_recalc(
        folder, mode=mode, time_limit=time_limit,
        override=override, allow_mutate=bool(allow_mutate),
    )
    return _dumps(_asdict(plan))


def tool_list_calc_options(filename: str) -> str:
    """Return the (Options)-dropdown items for one selected basename.

    Mirrors the dashboard's context-aware menu; the file_type field
    shows which classification was applied (control_txt / occupier_txt
    / inp / out / xyz / csv_complete / final_interp / other).
    """
    import json as _json
    return _json.dumps(
        delfin_api.list_calc_options(filename),
        separators=(",", ":"),
    )


def tool_run_calc_option(
    folder: str,
    option: str,
    target_file: str = "",
    time_limit: str = "24:00:00",
    override: str = "",
    allow_mutate: bool = False,
) -> str:
    """Dispatch one (Options) action.

    Recalc / Smart Recalc / Override route to ``prepare_recalc``;
    everything else (Visualize, MO Plot, Print NMR, Plot Trajectory,
    Preselection, RMSD, Build Batch from XYZ, Calc NMR, Calc CENSO/ANMR,
    hyperpol_xtb, tadf_xtb) returns a hint to drive it via ACTION:
    slash-commands (those run inside the dashboard UI).

    Args:
        folder: absolute path to the calc folder.
        option: dropdown text — see list_calc_options output.
        target_file: optional basename when the option needs one.
        time_limit: SLURM time for recalcs.
        override: STAGE=INDEX when option='Override'.
        allow_mutate: required for any submission.
    """
    import json as _json
    return _json.dumps(
        delfin_api.run_calc_option(
            folder, option,
            target_file=target_file,
            time_limit=time_limit,
            override=override,
            allow_mutate=bool(allow_mutate),
        ),
        separators=(",", ":"),
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
        separators=(",", ":"),
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
        separators=(",", ":"),
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
        separators=(",", ":"),
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
        separators=(",", ":"),
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
    return _dumps(delfin_api.check_orca_manual_indexed())


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
        separators=(",", ":"),
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
        separators=(",", ":"),
    )


def tool_explain_delfin_feature(name: str) -> str:
    """Explain a DELFIN concept (CONTROL keys, Smart Recalc, OCCUPIER, …).

    Use when the user asks "how does X work in DELFIN?" or "what
    does <feature> do?". Returns curated prose + source-file pointers
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
        separators=(",", ":"),
    )


def tool_list_active_calculations() -> str:
    """Live list of the scheduler's running/pending jobs (read-only).

    Returns JSON list of {job_id, name, status, runtime_s, directory}
    entries. Empty when no jobs are active -- which says nothing about
    folders on disk: a calculation submitted elsewhere, or one that
    crashed, is not in this list. For "which of these folders is still
    running" read the ``outcome`` column of extract_energy_table
    (running or crashed / no output yet / failed (exit code N)).
    """
    return _dumps(delfin_api.list_active_calculations())


def tool_plot_energy_correlation(
    folders: str,
    x: str = "single_point",
    y: str = "gibbs",
    title: str = "",
    output_path: str = "",
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
        output_path: where to write the PNG (default: the agent workspace).
    """
    from dataclasses import asdict as _asdict
    folder_list = [f.strip() for f in folders.split(",") if f.strip()]
    result = delfin_api.plot_energy_correlation(
        folder_list, x=x, y=y, title=title, output_path=output_path)
    return _dumps(_asdict(result))


def tool_plot_orbital_diagram(
    folder: str,
    n_below: int = 5,
    n_above: int = 5,
    title: str = "",
    output_path: str = "",
) -> str:
    """Render an MO level diagram around HOMO/LUMO; PNG → workspace.

    Reads the LAST ORBITAL ENERGIES block from the folder's largest
    .out, then draws horizontal lines at each orbital's eV value.
    HOMO highlighted blue, LUMO red. JSON return carries the PNG path
    so the dashboard's inline-render hook displays it in chat.

    Args:
        folder: absolute path to one calculation folder.
        n_below: occupied orbitals to show below HOMO (default 5).
        n_above: virtuals above LUMO (default 5).
        title: optional figure title.
        output_path: where to write the PNG (default: the agent workspace).
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.plot_orbital_diagram(
            folder, n_below=n_below, n_above=n_above, title=title, output_path=output_path)),
        separators=(",", ":"),
    )


def tool_plot_optimization_convergence(
    folder: str,
    title: str = "",
    output_path: str = "",
) -> str:
    """Render energy + ΔE per cycle for an Opt run; PNG → workspace.

    Two panels: absolute energy (Eh) and ΔE (kcal/mol on symlog).
    Title carries the converged/not-converged status. Direct answer
    to 'why did the optimization take 50 steps?'.
        output_path: where to write the PNG (default: the agent workspace).
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.plot_optimization_convergence(
            folder, title=title, output_path=output_path)),
        separators=(",", ":"),
    )


def tool_plot_uvvis_spectrum(
    folder: str,
    fwhm_nm: float = 20.0,
    wavelength_min: float = 200.0,
    wavelength_max: float = 800.0,
    n_points: int = 1000,
    title: str = "",
    output_path: str = "",
) -> str:
    """Gaussian-broadened UV/Vis from TDDFT; PNG → workspace.

    Uses the LAST ABSORPTION SPECTRUM block from the folder's .out.
    Stick spectrum drawn underneath the broadened curve.

    Args:
        folder: absolute path to a TDDFT calc folder.
        fwhm_nm: broadening FWHM in nm (default 20).
        wavelength_min, wavelength_max: plot window in nm.
        n_points: convolution grid resolution.
        title: optional figure title.
        output_path: where to write the PNG (default: the agent workspace).
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.plot_uvvis_spectrum(
            folder,
            fwhm_nm=fwhm_nm,
            wavelength_min=wavelength_min,
            wavelength_max=wavelength_max,
            n_points=n_points,
            title=title, output_path=output_path)),
        separators=(",", ":"),
    )


def tool_plot_scf_convergence(
    folder: str,
    cycle_index: int = -1,
    title: str = "",
    output_path: str = "",
) -> str:
    """Plot SCF iteration curves (energy vs iteration); PNG → workspace.

    For multi-cycle output (geom optimization), draws one curve per
    geom step. Set cycle_index >= 0 to plot a single cycle.
    Default = -1 (treat as None → all cycles overlaid).

    Direct answer to 'why does the SCF not converge?'.
        output_path: where to write the PNG (default: the agent workspace).
    """
    import json as _json
    from dataclasses import asdict as _asdict
    ci = None if int(cycle_index) < 0 else int(cycle_index)
    return _json.dumps(
        _asdict(delfin_api.plot_scf_convergence(
            folder, cycle_index=ci, title=title, output_path=output_path)),
        separators=(",", ":"),
    )


def tool_plot_population_charges(
    folder: str,
    method: str = "mulliken",
    title: str = "",
    output_path: str = "",
) -> str:
    """Bar chart of atomic charges (Mulliken or Loewdin); PNG → workspace.

    Args:
        folder: absolute path to a calc folder.
        method: 'mulliken' (default) or 'loewdin'.
        title: optional figure title.
        output_path: where to write the PNG (default: the agent workspace).
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.plot_population_charges(
            folder, method=method, title=title, output_path=output_path)),
        separators=(",", ":"),
    )


def tool_plot_vibrational_spectrum(
    folder: str,
    fwhm_cm: float = 12.0,
    freq_min: float = 0.0,
    freq_max: float = 4000.0,
    n_points: int = 1500,
    title: str = "",
    output_path: str = "",
) -> str:
    """Lorentzian-broadened IR spectrum from full mode list; PNG.

    Direct answer to 'what does the IR spectrum look like?' / 'which
    vibrational modes are IR-active?'.

    Args:
        folder: absolute path to a Freq calc folder.
        fwhm_cm: Lorentzian FWHM in cm-1 (default 12).
        freq_min, freq_max: window in cm-1 (default 0-4000).
        n_points: convolution grid resolution.
        title: optional figure title.
        output_path: where to write the PNG (default: the agent workspace).
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.plot_vibrational_spectrum(
            folder,
            fwhm_cm=fwhm_cm,
            freq_min=freq_min,
            freq_max=freq_max,
            n_points=n_points,
            title=title, output_path=output_path)),
        separators=(",", ":"),
    )


def tool_find_calculation_extreme(
    folders: str,
    property: str = "gibbs",
    extreme: str = "min",
    n: int = 5,
) -> str:
    """The N lowest/highest folders by a property, PER METHOD.

    Returns {"note", "root", "property_requested", "property_used",
    "groups": [{"method", "rows"}], "skipped": [{"folder", "reason",
    "outcome"}]}. Folders are given relative to "root", the directory
    they all share, so a long path is said once:
    within each method (functional/basis) the top n rows, ranked; groups
    are not ranked against each other, because a total energy compares
    only within one method. "Find the .out with the lowest Gibbs energy"
    is answered per method.

    Every folder that is not in a group is in "skipped" with the reason
    -- missing, no output (the outcome says whether it is still running,
    crashed or never started), or the property is not in its output.
    When NO folder carries the requested property but single point
    energies are there, the ranking uses single_point and says so in
    property_used and in the note; an archive of single points has no
    Gibbs energy to rank by.

    Args:
        folders: comma-separated absolute paths.
        property: gibbs (default) | zpe | single_point | imag_freqs |
            walltime_s.
        extreme: "min" (lowest, default) or "max" (highest).
        n: how many top entries per method to return (default 5).
    
    """
    folder_list = [f.strip() for f in folders.split(",") if f.strip()]
    res = delfin_api.find_calculation_extreme_explained(
        folder_list, property=property, extreme=extreme, n=int(n),
    )
    root = _common_root([r.get("folder") for r in res["rows"] + res["skipped"]])
    out = _grouped_by_method(res["rows"])
    out["root"] = root
    for g in out["groups"]:
        g["rows"] = _relative_to_root([dict(r, folder=os.path.join(root, r["folder"]) if root and r.get("folder") and not os.path.isabs(str(r["folder"])) else r.get("folder")) for r in g["rows"]], root)
    out["property_requested"] = res["property_requested"]
    out["property_used"] = res["property_used"]
    out["skipped"] = _relative_to_root(res["skipped"], root)
    if res["property_used"] != res["property_requested"]:
        out["note"] += (f" No folder carries {res['property_requested']}; "
                        f"ranked by {res['property_used']} instead.")
    if not res["rows"]:
        out["note"] += (" Nothing to rank: 'skipped' says why each "
                        "folder is left out.")
    return _dumps(out)


def tool_extract_imaginary_frequencies(folder: str) -> str:
    """List imaginary modes (n_imag, mode_index, frequency_cm) for one folder.

    Reads the largest ``*.out`` containing a VIBRATIONAL FREQUENCIES
    block and returns:

    - ``n_imag``: count of imaginary modes
    - ``modes``: list of {mode_index, frequency_cm}
    - ``most_negative``: most negative cm**-1 value (None if no imag)
    - ``is_minimum`` (n_imag == 0), ``is_ts`` (n_imag == 1)
    - ``error``: filled when no Freq output is present

    Use this when the user asks "is X a minimum?" / "find the TS" /
    "compare imaginary frequencies". Combine with
    ``compare_across_functionals`` to scan many folders.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.extract_imaginary_frequencies(folder)),
        separators=(",", ":"),
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
        separators=(",", ":"),
    )


def tool_compare_across_functionals(
    folders: str,
    include_imag: bool = True,
    sort_by: str = "gibbs",
) -> str:
    """Multi-folder comparison table grouped by method (functional/basis).

    Returns {"note", "groups": [{"method", "rows"}]}: one row per folder
    with functional, basis, gibbs, single_point, zpe, n_imag,
    is_minimum, status. Sorting by gibbs (default), single_point or zpe
    orders rows WITHIN a method; groups are never ranked against each
    other, because a total energy compares only within one method.
    Answers "which run is lowest within each method" and "compare
    imaginary frequencies across functionals" -- not "which functional
    gives the lowest minimum", which has no answer.

    Args:
        folders: comma-separated absolute paths.
        include_imag: if True (default), also extract imaginary-freq
            counts (slightly slower but usually wanted).
        sort_by: gibbs | single_point | zpe | functional | folder.
    
    """
    from dataclasses import asdict as _asdict
    folder_list = [f.strip() for f in folders.split(",") if f.strip()]
    rows = delfin_api.compare_across_functionals(
        folder_list, include_imag=include_imag, sort_by=sort_by,
    )
    return _dumps(_grouped_by_method([_asdict(r) for r in rows]))


def tool_extract_orbital_energies(folder: str) -> str:
    """Pull the LAST ORBITAL ENERGIES block + HOMO/LUMO/gap.

    Returns the full orbital list (index, occupation, energy_eh,
    energy_ev) plus homo_index/lumo_index/homo_ev/lumo_ev/gap_ev.
    Direct answer to "where is the HOMO?" / "how large is the gap?".
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.extract_orbital_energies(folder)),
        separators=(",", ":"),
    )


def tool_extract_excited_states(folder: str) -> str:
    """Pull the TDDFT/CIS transition table from a folder's .out.

    Selects the LAST ABSORPTION SPECTRUM block (post-optimization).
    Each row has state_from, state_to, energy_ev, energy_cm,
    wavelength_nm, fosc. Use for UV/Vis spectrum analysis. Also answers
    the questions the table alone did not: first_bright (the lowest
    bright transition), brightest, brightest_visible (the strongest
    inside visible_range_nm), each as {index, state, wavelength_nm,
    energy_ev, fosc, in_visible} -- with the definitions stated
    (bright_threshold_fosc, visible_range_nm) so the answer can be
    quoted with its rule. null when no transition qualifies.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.extract_excited_states(folder)),
        separators=(",", ":"),
    )


def tool_extract_dipole(folder: str) -> str:
    """Pull the dipole-moment vector + magnitude from a folder's .out.

    Returns dx/dy/dz (a.u.), magnitude_au, and magnitude_debye
    (computed if ORCA didn't print Debye explicitly).
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.extract_dipole(folder)),
        separators=(",", ":"),
    )


def tool_extract_optimization_trajectory(folder: str) -> str:
    """Per-cycle energies + convergence flag for an Opt run.

    Each row carries cycle, energy_eh, delta_e (Δ to previous cycle).
    ``converged`` follows ORCA's OPTIMIZATION RUN DONE / HAS
    CONVERGED markers. Useful for "why did the optimization take 50
    steps?" and trajectory plotting.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.extract_optimization_trajectory(folder)),
        separators=(",", ":"),
    )


def tool_extract_scf_convergence(folder: str) -> str:
    """Per-iteration SCF history for every geom-step in a folder.

    Returns one entry per SCF cycle with the full iteration table
    (iteration #, energy, ΔE, max density change). Useful for
    'why does the SCF not converge?' diagnostics + plotting
    energy(iter) curves to spot oscillation / divergence.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.extract_scf_convergence(folder)),
        separators=(",", ":"),
    )


def tool_extract_mulliken_charges(folder: str) -> str:
    """Mulliken atomic charges + spin populations from .out.

    Returns the LAST Mulliken block (post-optimization). Open-shell
    calculations include spin_population per atom; closed-shell omit
    it. Direct answer to 'what is the charge distribution?'.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.extract_mulliken_charges(folder)),
        separators=(",", ":"),
    )


def tool_extract_loewdin_charges(folder: str) -> str:
    """Loewdin atomic charges (basis-set-stable alternative to Mulliken).

    Loewdin charges are the preferred metric for cross-method
    comparison (Mulliken charges depend on basis-set choice). Same
    return shape as extract_mulliken_charges.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.extract_loewdin_charges(folder)),
        separators=(",", ":"),
    )


def tool_extract_vibrational_modes(folder: str) -> str:
    """All vibrational modes (real + imaginary) + IR intensities.

    Each row: mode_index, frequency_cm, ir_intensity (km/mol if
    present), is_imaginary. Use for full vibrational spectrum
    analysis or to pick a mode for animation.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.extract_vibrational_modes(folder)),
        separators=(",", ":"),
    )


def tool_extract_delfin_json(folder: str) -> str:
    """Read the structured DELFIN_Data.json pipeline state.

    Returns workflow_stages, per-stage energies + timings, total
    cost (USD), and raw_keys (for schema-version discovery).
    Direct answer to 'what has the pipeline done so far?'.
    """
    import json as _json
    from dataclasses import asdict as _asdict
    return _json.dumps(
        _asdict(delfin_api.extract_delfin_json(folder)),
        separators=(",", ":"),
    )


def tool_extract_calc_summary_table(folders: str) -> str:
    """Multi-property comparison row per folder.

    For each folder: functional, basis, gibbs, single_point, zpe,
    HOMO/LUMO/gap, n_imag, dipole_debye, walltime. Useful for
    benchmark tables across methods/molecules.

    Args:
        folders: comma-separated absolute paths.
    """
    from dataclasses import asdict as _asdict
    folder_list = [f.strip() for f in folders.split(",") if f.strip()]
    rows = delfin_api.extract_calc_summary_table(folder_list)
    return _dumps([_asdict(r) for r in rows])


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

    from ..mcp_compat import load_server_class

    mcp = load_server_class()(
        "delfin-ops",
        instructions=(
            "DELFIN operations server. Run DELFIN workflows (pipeline, "
            "ORCA, TADF, hyperpol, CO2), inspect tool availability, and "
            "control running processes. Read-only tools are always safe; "
            "mutating tools need the host's grant, which no tool argument "
            "can set — ask the user, then let them grant it."
        ),
    )

    # Read-only — register module functions directly
    mcp.tool(name="qm_check")(_safe("qm_check", tool_qm_check))
    mcp.tool(name="csp_check")(_safe("csp_check", tool_csp_check))
    mcp.tool(name="mlp_check")(_safe("mlp_check", tool_mlp_check))
    mcp.tool(name="analysis_check")(_safe("analysis_check", tool_analysis_check))
    mcp.tool(name="list_dashboard_patterns")(_safe("list_dashboard_patterns", tool_list_dashboard_patterns))
    mcp.tool(name="get_dashboard_pattern")(_safe("get_dashboard_pattern", tool_get_dashboard_pattern))
    # P1 — output parsing (read-only, structured returns)
    mcp.tool(name="parse_orca_output")(_safe("parse_orca_output", tool_parse_orca_output))
    mcp.tool(name="calc_status")(_safe("calc_status", tool_calc_status))
    mcp.tool(name="find_orca_errors")(_safe("find_orca_errors", tool_find_orca_errors))
    mcp.tool(name="extract_thermochem")(_safe("extract_thermochem", tool_extract_thermochem))
    mcp.tool(name="extract_energy_table")(_safe("extract_energy_table", tool_extract_energy_table))
    mcp.tool(name="find_calculation_extreme")(_safe("find_calculation_extreme", tool_find_calculation_extreme))
    # Imaginary-frequency + functional-comparison helpers
    mcp.tool(name="extract_imaginary_frequencies")(_safe("extract_imaginary_frequencies", tool_extract_imaginary_frequencies))
    mcp.tool(name="compare_calculations")(_safe("compare_calculations", tool_compare_calculations))
    mcp.tool(name="compare_across_functionals")(_safe("compare_across_functionals", tool_compare_across_functionals))
    # Output-analysis depth (orbitals / TDDFT / dipole / opt trajectory)
    mcp.tool(name="extract_orbital_energies")(_safe("extract_orbital_energies", tool_extract_orbital_energies))
    mcp.tool(name="extract_excited_states")(_safe("extract_excited_states", tool_extract_excited_states))
    mcp.tool(name="extract_dipole")(_safe("extract_dipole", tool_extract_dipole))
    mcp.tool(name="extract_optimization_trajectory")(_safe("extract_optimization_trajectory", tool_extract_optimization_trajectory))
    # Phase E parsers (SCF / population / vib modes / DELFIN json / summary table)
    mcp.tool(name="extract_scf_convergence")(_safe("extract_scf_convergence", tool_extract_scf_convergence))
    mcp.tool(name="extract_mulliken_charges")(_safe("extract_mulliken_charges", tool_extract_mulliken_charges))
    mcp.tool(name="extract_loewdin_charges")(_safe("extract_loewdin_charges", tool_extract_loewdin_charges))
    mcp.tool(name="extract_vibrational_modes")(_safe("extract_vibrational_modes", tool_extract_vibrational_modes))
    mcp.tool(name="extract_delfin_json")(_safe("extract_delfin_json", tool_extract_delfin_json))
    mcp.tool(name="extract_calc_summary_table")(_safe("extract_calc_summary_table", tool_extract_calc_summary_table))
    # P1 — statistical plots (PNG → agent_workspace, auto-displayed)
    mcp.tool(name="plot_energy_distribution")(_safe("plot_energy_distribution", tool_plot_energy_distribution))
    mcp.tool(name="plot_energy_correlation")(_safe("plot_energy_correlation", tool_plot_energy_correlation))
    # Phase D plots (orbitals / opt convergence / UV/Vis spectrum)
    mcp.tool(name="plot_orbital_diagram")(_safe("plot_orbital_diagram", tool_plot_orbital_diagram))
    mcp.tool(name="plot_optimization_convergence")(_safe("plot_optimization_convergence", tool_plot_optimization_convergence))
    mcp.tool(name="plot_uvvis_spectrum")(_safe("plot_uvvis_spectrum", tool_plot_uvvis_spectrum))
    # Phase E plots (SCF / charges / vibrational IR spectrum)
    mcp.tool(name="plot_scf_convergence")(_safe("plot_scf_convergence", tool_plot_scf_convergence))
    mcp.tool(name="plot_population_charges")(_safe("plot_population_charges", tool_plot_population_charges))
    mcp.tool(name="plot_vibrational_spectrum")(_safe("plot_vibrational_spectrum", tool_plot_vibrational_spectrum))
    # Tool / widget catalogs (cheap on-demand discovery)
    mcp.tool(name="list_tools")(_safe("list_tools", tool_list_tools))
    mcp.tool(name="describe_tool")(_safe("describe_tool", tool_describe_tool))
    mcp.tool(name="list_dashboard_widgets")(_safe("list_dashboard_widgets", tool_list_dashboard_widgets))
    mcp.tool(name="get_widget_options")(_safe("get_widget_options", tool_get_widget_options))
    # ORCA Builder validation
    mcp.tool(name="validate_orca_input")(_safe("validate_orca_input", tool_validate_orca_input))
    # Job lifecycle (read-only list + mutating submit/cancel)
    mcp.tool(name="list_active_calculations")(_safe("list_active_calculations", tool_list_active_calculations))
    mcp.tool(name="submit_calculation")(_host_gated(tool_submit_calculation))
    mcp.tool(name="cancel_calculation")(_host_gated(tool_cancel_calculation))
    # Calc folder management (mutating, allow_mutate-gated)
    mcp.tool(name="rename_calc_folder")(_host_gated(tool_rename_calc_folder))
    mcp.tool(name="create_calc_folder")(_host_gated(tool_create_calc_folder))
    mcp.tool(name="move_calc_folder")(_host_gated(tool_move_calc_folder))
    mcp.tool(name="move_to_archive")(_host_gated(tool_move_to_archive))
    mcp.tool(name="delete_calc_folder")(_host_gated(tool_delete_calc_folder))
    # Bulk job control + recalc preparation + Options dispatcher
    mcp.tool(name="kill_all_user_jobs")(_host_gated(tool_kill_all_user_jobs))
    mcp.tool(name="prepare_recalc")(_host_gated(tool_prepare_recalc))
    mcp.tool(name="list_calc_options")(_safe("list_calc_options", tool_list_calc_options))
    mcp.tool(name="run_calc_option")(_host_gated(tool_run_calc_option))
    mcp.tool(name="list_ssh_transfer_jobs")(_safe("list_ssh_transfer_jobs", tool_list_ssh_transfer_jobs))
    # ORCA-manual lookup + literature indexing
    mcp.tool(name="check_orca_manual_indexed")(_safe("check_orca_manual_indexed", tool_check_orca_manual_indexed))
    mcp.tool(name="index_new_pdf")(_safe("index_new_pdf", tool_index_new_pdf))
    # PDF on-demand reading (no pre-indexing)
    mcp.tool(name="read_pdf")(_safe("read_pdf", tool_read_pdf))
    mcp.tool(name="search_pdf_local")(_safe("search_pdf_local", tool_search_pdf_local))
    mcp.tool(name="extract_pdf_section")(_safe("extract_pdf_section", tool_extract_pdf_section))
    mcp.tool(name="list_literature_files")(_safe("list_literature_files", tool_list_literature_files))
    # DELFIN-feature explainer
    mcp.tool(name="list_delfin_features")(_safe("list_delfin_features", tool_list_delfin_features))
    mcp.tool(name="explain_delfin_feature")(_safe("explain_delfin_feature", tool_explain_delfin_feature))

    # stop_dry_run needs the default workspace closed over
    @mcp.tool(name="stop_dry_run", description=tool_stop_dry_run.__doc__)
    def _stop_dry_run(workspace: str = "") -> str:
        return tool_stop_dry_run(workspace or default_workspace)

    # Mutating wrappers default workspace = server cwd. No allow_mutate
    # parameter: the value comes from the host, never from the caller.
    @mcp.tool(name="cleanup",
              description=_strip_allow_mutate_doc(tool_cleanup.__doc__))
    def _cleanup(
        orca: bool = False,
        dry_run: bool = True,
        workspace: str = "",
        scratch: str = "",
    ) -> str:
        return tool_cleanup(
            orca=orca, dry_run=dry_run,
            workspace=workspace or default_workspace,
            scratch=scratch, allow_mutate=host_grants_mutation(),
        )

    @mcp.tool(name="stop",
              description=_strip_allow_mutate_doc(tool_stop.__doc__))
    def _stop(
        signal_name: str = "INT",
        workspace: str = "",
        dry_run: bool = True,
        cleanup_after: bool = False,
        wait_seconds: float = 3.0,
    ) -> str:
        return tool_stop(
            signal_name=signal_name,
            workspace=workspace or default_workspace,
            dry_run=dry_run, cleanup_after=cleanup_after,
            wait_seconds=wait_seconds, allow_mutate=host_grants_mutation(),
        )

    mcp.tool(name="pipeline_prepare")(_host_gated(tool_pipeline_prepare))
    mcp.tool(name="pipeline_run")(_host_gated(tool_pipeline_run))
    mcp.tool(name="run_orca_input")(_host_gated(tool_run_orca_input))
    mcp.tool(name="co2")(_host_gated(tool_co2))
    mcp.tool(name="tadf_xtb")(_host_gated(tool_tadf_xtb))
    mcp.tool(name="hyperpol")(_host_gated(tool_hyperpol))

    mcp.run(transport="stdio")
