"""Programmatic Python API for DELFIN.

Wraps the CLI subcommands (``delfin/cli.py``) with typed Python functions
returning structured ``CommandResult`` objects.  Intended for embedding in
notebooks, agent tooling (MCP server), tests, and downstream automation.

The internal pattern is:

1. Build the argv list from typed kwargs.
2. Run ``cli_main(argv)`` with captured stdout/stderr.
3. Return a ``CommandResult`` containing exit code + captured streams.

This keeps the CLI as the single source of truth and avoids logic duplication.
"""
from __future__ import annotations

import contextlib
import io
import sys
from dataclasses import dataclass, field
from typing import Sequence


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

@dataclass
class CommandResult:
    """Structured result of a DELFIN CLI subcommand.

    Attributes
    ----------
    returncode : int
        Exit code (0 = success).  ``cli_main`` returns this directly; on
        ``SystemExit`` we capture ``code`` from the exception.
    stdout : str
        Captured standard output.
    stderr : str
        Captured standard error.
    argv : list[str]
        The argv passed to ``cli_main`` (useful for debugging / logging).
    """

    returncode: int = 0
    stdout: str = ""
    stderr: str = ""
    argv: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.returncode == 0


# ---------------------------------------------------------------------------
# Internal CLI runner
# ---------------------------------------------------------------------------

def _has_flag(argv: Sequence[str], *flags: str) -> bool:
    for flag in flags:
        if flag in argv:
            return True
        prefix = f"{flag}="
        if any(arg.startswith(prefix) for arg in argv):
            return True
    return False


def _run_cli(argv: list[str], *, capture: bool = True) -> CommandResult:
    """Invoke ``delfin.cli.main`` with a built argv and capture the streams.

    Parameters
    ----------
    argv : list[str]
        The argv passed to ``cli_main`` (without the leading ``delfin``).
    capture : bool
        If True (default), capture stdout/stderr into the result.  If False,
        streams pass through to the caller's terminal.
    """
    from delfin.cli import main as cli_main

    if not capture:
        try:
            rc = cli_main(argv)
        except SystemExit as exc:
            rc = int(exc.code) if exc.code is not None else 0
        return CommandResult(returncode=rc, argv=list(argv))

    out_buf = io.StringIO()
    err_buf = io.StringIO()
    rc: int = 0
    with contextlib.redirect_stdout(out_buf), contextlib.redirect_stderr(err_buf):
        try:
            rc = cli_main(argv)
        except SystemExit as exc:
            rc = int(exc.code) if exc.code is not None else 0
    return CommandResult(
        returncode=int(rc or 0),
        stdout=out_buf.getvalue(),
        stderr=err_buf.getvalue(),
        argv=list(argv),
    )


# ---------------------------------------------------------------------------
# Pipeline (run / prepare) — kept rückwärtskompatibel
# ---------------------------------------------------------------------------

def run(
    control_file: str = "CONTROL.txt",
    *,
    cleanup: bool = True,
    recalc: bool = False,
    overwrite: bool = False,
    define: str | None = None,
    extra_args: Sequence[str] | None = None,
) -> int:
    """Execute the DELFIN pipeline programmatically using CLI semantics.

    Returns the raw exit code (int) for backwards compatibility.
    Use :func:`pipeline_run` for a structured ``CommandResult``.
    """
    return pipeline_run(
        control_file=control_file,
        cleanup=cleanup,
        recalc=recalc,
        overwrite=overwrite,
        define=define,
        extra_args=extra_args,
        capture=False,
    ).returncode


def prepare(control_file: str = "CONTROL.txt", overwrite: bool = False) -> int:
    """Convenience wrapper for ``delfin --define`` to create CONTROL templates."""
    return pipeline_prepare(
        control_file=control_file, overwrite=overwrite, capture=False
    ).returncode


def pipeline_run(
    control_file: str = "CONTROL.txt",
    *,
    cleanup: bool = True,
    recalc: bool = False,
    overwrite: bool = False,
    define: str | None = None,
    extra_args: Sequence[str] | None = None,
    capture: bool = True,
) -> CommandResult:
    """Execute the DELFIN pipeline and return a ``CommandResult``."""
    argv: list[str] = list(extra_args or [])

    if recalc and not _has_flag(argv, "--recalc"):
        argv.append("--recalc")
    if not cleanup and not _has_flag(argv, "--no-cleanup"):
        argv.append("--no-cleanup")
    if overwrite and not _has_flag(argv, "--overwrite"):
        argv.append("--overwrite")

    if not _has_flag(argv, "--control", "-F") and control_file != "CONTROL.txt":
        argv.extend(["--control", control_file])

    if define is not None and not _has_flag(argv, "--define", "-D"):
        argv.extend(["--define", define])

    return _run_cli(argv, capture=capture)


def pipeline_prepare(
    control_file: str = "CONTROL.txt",
    overwrite: bool = False,
    *,
    capture: bool = True,
) -> CommandResult:
    """Generate a CONTROL template via ``delfin --define <file>``."""
    argv = ["--define", control_file]
    if overwrite:
        argv.append("--overwrite")
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# Tool / runtime checks (read-only, safe)
# ---------------------------------------------------------------------------

def qm_check(
    tools: Sequence[str] | None = None,
    *,
    capture: bool = True,
) -> CommandResult:
    """Inspect QM tool resolution (xtb, crest, etc.).  Read-only."""
    argv = ["qm_check"]
    if tools:
        argv.extend(list(tools))
    return _run_cli(argv, capture=capture)


def csp_check(*, capture: bool = True) -> CommandResult:
    """Check CSP (genarris) tool availability.  Read-only."""
    return _run_cli(["csp_check"], capture=capture)


def mlp_check(*, capture: bool = True) -> CommandResult:
    """Check MLP (torchani / AIMNet2 / MACE) backend availability.  Read-only."""
    return _run_cli(["mlp_check"], capture=capture)


def analysis_check(*, capture: bool = True) -> CommandResult:
    """Check analysis tools (Multiwfn, CENSO, ANMR, morfeus).  Read-only."""
    return _run_cli(["analysis_check"], capture=capture)


# ---------------------------------------------------------------------------
# QM tool runner
# ---------------------------------------------------------------------------

def qm_run(
    tool: str,
    tool_args: Sequence[str] | None = None,
    *,
    cwd: str = ".",
    capture_tool: bool = True,
    capture: bool = True,
) -> CommandResult:
    """Run a single QM tool via DELFIN's resolver.

    Parameters
    ----------
    tool : str
        Tool name (e.g. ``xtb``, ``crest``, ``xtb4stda``, ``stda``, ``dftb+``).
    tool_args : sequence of str
        Arguments forwarded to the tool.
    cwd : str
        Working directory for the tool process.
    capture_tool : bool
        Pass ``--capture`` so ORCA prints stdout/stderr after completion.
    capture : bool
        Capture DELFIN's own output into the ``CommandResult``.
    """
    argv = ["qm_run", tool, "--cwd", cwd]
    if capture_tool:
        argv.append("--capture")
    if tool_args:
        argv.append("--")
        argv.extend(list(tool_args))
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# Process control (cleanup / stop)
# ---------------------------------------------------------------------------

def cleanup(
    *,
    orca: bool = False,
    dry_run: bool = False,
    workspace: str = ".",
    scratch: str | None = None,
    capture: bool = True,
) -> CommandResult:
    """Remove DELFIN scratch artifacts and optionally stop ORCA."""
    argv = ["cleanup", "--workspace", workspace]
    if orca:
        argv.append("--orca")
    if dry_run:
        argv.append("--dry-run")
    if scratch:
        argv.extend(["--scratch", scratch])
    return _run_cli(argv, capture=capture)


def stop(
    *,
    workspace: str = ".",
    signal_name: str = "INT",
    dry_run: bool = False,
    cleanup_after: bool = False,
    wait_seconds: float = 3.0,
    capture: bool = True,
) -> CommandResult:
    """Send a signal (INT/TERM/KILL) to running DELFIN processes.

    Parameters
    ----------
    signal_name : {"INT", "TERM", "KILL"}
        Signal to send.  Defaults to graceful SIGINT.
    cleanup_after : bool
        Run cleanup after signaling (only when not dry-run).
    """
    if signal_name not in {"INT", "TERM", "KILL"}:
        raise ValueError(
            f"signal_name must be one of INT/TERM/KILL, got {signal_name!r}"
        )
    argv = ["stop", "--workspace", workspace, "--signal", signal_name,
            "--wait-seconds", str(wait_seconds)]
    if dry_run:
        argv.append("--dry-run")
    if cleanup_after:
        argv.append("--cleanup")
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# ORCA invocation
# ---------------------------------------------------------------------------

def run_orca_input(
    input_file: str | None = None,
    output: str | None = None,
    *,
    capture: bool = True,
) -> CommandResult:
    """Run ORCA on a ``.inp`` file (uses DELFIN's ORCA resolver)."""
    argv: list[str] = ["run_orca"]
    if input_file:
        argv.append(input_file)
    if output:
        argv.extend(["--output", output])
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# CO2 Coordinator
# ---------------------------------------------------------------------------

def co2(
    *,
    define: bool = False,
    force: bool = False,
    recalc: bool = False,
    charge: int | None = None,
    multiplicity: int | None = None,
    solvent: str | None = None,
    metal: str | None = None,
    broken_sym: str | None = None,
    capture: bool = True,
) -> CommandResult:
    """CO2 Coordinator workflow."""
    argv: list[str] = ["co2"]
    if define:
        argv.append("--define")
    if force:
        argv.append("--force")
    if recalc:
        argv.append("--recalc")
    if charge is not None:
        argv.extend(["--charge", str(charge)])
    if multiplicity is not None:
        argv.extend(["--multiplicity", str(multiplicity)])
    if solvent is not None:
        argv.extend(["--solvent", solvent])
    if metal is not None:
        argv.extend(["--metal", metal])
    if broken_sym is not None:
        argv.extend(["--broken_sym", broken_sym])
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# Specialised workflows: TADF (xTB), hyperpolarisability
# ---------------------------------------------------------------------------

def tadf_xtb(
    extra_args: Sequence[str] | None = None,
    *,
    capture: bool = True,
) -> CommandResult:
    """Run the TADF xTB workflow (``delfin tadf_xtb``)."""
    argv: list[str] = ["tadf_xtb"]
    if extra_args:
        argv.extend(list(extra_args))
    return _run_cli(argv, capture=capture)


def hyperpol(
    extra_args: Sequence[str] | None = None,
    *,
    capture: bool = True,
) -> CommandResult:
    """Run the hyperpolarisability workflow (``delfin hyperpol``)."""
    argv: list[str] = ["hyperpol"]
    if extra_args:
        argv.extend(list(extra_args))
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# P1 — Output parsing (read-only structured returns for the agent)
# ---------------------------------------------------------------------------


@dataclass
class OrcaParseResult:
    """Structured ORCA-output snapshot.

    Numeric fields are ``None`` when the corresponding pattern wasn't
    found in the output. ``scf_converged`` / ``opt_converged`` are
    tri-state via ``None`` when the marker is absent — don't mistake
    "no marker" for "did not converge".
    """
    path: str
    final_single_point: float | None = None
    gibbs_free_energy: float | None = None
    zpe: float | None = None
    scf_converged: bool | None = None
    opt_converged: bool | None = None
    imag_freq_count: int | None = None
    walltime_s: float | None = None
    n_atoms: int | None = None
    functional: str = ""
    basis: str = ""
    error_summary: str = ""


def parse_orca_output(path: str) -> OrcaParseResult:
    """Parse one ORCA output file and return a structured snapshot.

    Wraps the per-property helpers in :mod:`delfin.energies` so the
    agent can ask "what happened in this run?" with a single tool
    call instead of grepping the .out file line by line.

    Missing fields → ``None`` (not zero — distinguishes "absent"
    from "found but zero-valued"). Read errors don't raise; the
    result's ``error_summary`` carries the diagnosis.
    """
    from pathlib import Path as _P
    from delfin import energies as _e
    p = _P(path)
    if not p.exists():
        return OrcaParseResult(path=str(p), error_summary="file not found")
    if not p.is_file():
        return OrcaParseResult(path=str(p), error_summary="not a file")
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except Exception as exc:
        return OrcaParseResult(path=str(p), error_summary=f"read error: {exc}")

    out = OrcaParseResult(path=str(p))
    try:
        out.final_single_point = _e.find_electronic_energy(str(p))
    except Exception:
        pass
    try:
        out.gibbs_free_energy = _e.find_gibbs_energy(str(p))
    except Exception:
        pass
    try:
        out.zpe = _e.find_ZPE(str(p))
    except Exception:
        pass

    # SCF convergence: ORCA writes a positive marker on success.
    if "SCF CONVERGED AFTER" in text or "SCF CONVERGED" in text:
        out.scf_converged = True
    elif "SCF NOT CONVERGED" in text or "SCF iterations did not converge" in text:
        out.scf_converged = False

    # Optimization convergence (geom_opt jobs only).
    if ("OPTIMIZATION RUN DONE" in text
            or "THE OPTIMIZATION HAS CONVERGED" in text):
        out.opt_converged = True
    elif ("OPTIMIZATION DID NOT CONVERGE" in text
            or "FAILED TO CONVERGE THE GEOMETRY OPTIMIZATION" in text):
        out.opt_converged = False

    # Imaginary frequency count.
    import re as _re
    m_imag = _re.search(
        r"Number of imaginary frequencies\s*\.\.\.\s*(\d+)", text,
    )
    if m_imag:
        out.imag_freq_count = int(m_imag.group(1))
    else:
        # Count negative-frequency lines: "  4:    -123.45 cm**-1"
        imag_lines = _re.findall(
            r"^\s*\d+:\s+(-\d+\.\d+)\s+cm\*\*-1", text, _re.MULTILINE,
        )
        if imag_lines:
            out.imag_freq_count = len(imag_lines)

    # Number of atoms.
    m_atoms = _re.search(r"Number of atoms\s+\.\.\.\s+(\d+)", text)
    if m_atoms:
        out.n_atoms = int(m_atoms.group(1))

    # Walltime: "TOTAL RUN TIME: 0 days 1 hours 23 minutes 45 seconds 678 msec"
    m_wall = _re.search(
        r"TOTAL RUN TIME:\s+(\d+)\s+days\s+(\d+)\s+hours\s+(\d+)\s+minutes\s+"
        r"(\d+)\s+seconds\s+(\d+)\s+msec",
        text,
    )
    if m_wall:
        d, h, mi, s, ms = (int(g) for g in m_wall.groups())
        out.walltime_s = d * 86400 + h * 3600 + mi * 60 + s + ms / 1000.0

    # Functional / basis from the input keyword line.
    m_kw = _re.search(r"^\s*\|\s*\d+>\s*!\s*(.+)$", text, _re.MULTILINE)
    if m_kw:
        for token in m_kw.group(1).split():
            T = token.upper()
            if T in ("BP86", "PBE0", "B3LYP", "TPSS", "WB97X", "WB97X-D3",
                    "CAM-B3LYP", "R2SCAN", "PBE", "M06", "M062X", "BLYP",
                    "REVTPSS"):
                out.functional = token
            if ("DEF2-" in T or "PCSSEG" in T or "MA-DEF2-" in T
                    or "SARC-" in T or "X2C-" in T):
                out.basis = token

    if (out.final_single_point is None
            and out.scf_converged is not True
            and "ABORTING THE RUN" in text):
        out.error_summary = "ORCA aborted the run"

    return out


@dataclass
class OrcaError:
    """One detected error pattern in an ORCA output."""
    type: str
    message: str
    line_number: int | None = None
    suggestion: str = ""


_ORCA_ERROR_PATTERNS: list[tuple[str, str, str]] = [
    (r"SCF NOT CONVERGED|SCF iterations did not converge",
     "scf_diverge",
     "Add SlowConv or tighten SCFconv; try DIIS off or KDIIS; check spin guess."),
    (r"oom_kill|Killed by signal 9|out of memory",
     "oom",
     "Increase --mem (sbatch) or lower %maxcore in the .inp."),
    (r"Number of basis functions.*exceeds",
     "basis",
     "Reduce basis (def2-TZVP → def2-SVP) or split atom regions."),
    (r"Wrong multiplicity|Multiplicity\s+\d+\s+impossible",
     "multiplicity",
     "Recheck charge/multiplicity — total electrons must match."),
    (r"MPI ERROR|MPI_ABORT|orterun.*detected.*aborted",
     "mpi",
     "MPI process died — check sbatch --ntasks vs nprocs in %pal."),
    (r"DUE TO TIME LIMIT|Job exceeded|TIMEOUT",
     "timeout",
     "Bump SLURM time limit; or restart from .gbw + .opt with maxiter=N."),
    (r"FATAL ERROR|ABORTING THE RUN",
     "other",
     "ORCA aborted — read the output around this line for the cause."),
]


def find_orca_errors(folder: str) -> list[OrcaError]:
    """Scan ``*.out`` files in ``folder`` for known ORCA error patterns.

    Non-recursive (one folder at a time keeps results focused). Each
    match yields an OrcaError with a short suggestion. An empty list
    is NOT proof of success — use ``parse_orca_output`` for that.
    """
    from pathlib import Path as _P
    import re as _re
    errors: list[OrcaError] = []
    p = _P(folder)
    if not p.exists() or not p.is_dir():
        return errors
    for out_file in sorted(p.glob("*.out")):
        try:
            lines = out_file.read_text(
                encoding="utf-8", errors="replace",
            ).splitlines()
        except Exception:
            continue
        for ln_no, line in enumerate(lines, start=1):
            for pat, etype, suggestion in _ORCA_ERROR_PATTERNS:
                if _re.search(pat, line, _re.IGNORECASE):
                    errors.append(OrcaError(
                        type=etype,
                        message=line.strip()[:200],
                        line_number=ln_no,
                        suggestion=suggestion,
                    ))
                    break
    return errors


@dataclass
class ThermochemResult:
    """Thermochemistry block extracted from an ORCA Freq output."""
    path: str
    temperature_k: float | None = None
    pressure_atm: float | None = None
    zpe: float | None = None
    thermal_corr: float | None = None
    enthalpy_corr: float | None = None
    entropy_total: float | None = None
    gibbs_corr: float | None = None
    final_gibbs: float | None = None


def extract_thermochem(folder: str) -> ThermochemResult:
    """Extract the full thermochemistry block from an ORCA Freq output.

    Picks the first ``*.out`` containing thermochemistry data in the
    folder. All energies in Hartree, T in K, P in atm.
    """
    from pathlib import Path as _P
    import re as _re
    p = _P(folder)
    if not p.exists():
        return ThermochemResult(path=str(p))
    out_files = sorted(p.glob("*.out")) if p.is_dir() else [p]
    for out_file in out_files:
        try:
            text = out_file.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue
        if "THERMOCHEMISTRY" not in text.upper():
            continue
        result = ThermochemResult(path=str(out_file))
        FLOAT = r"(-?\d+\.\d+)"
        for attr, pattern in [
            ("temperature_k",  rf"Temperature\s+\.\.\.\s+{FLOAT}\s*K"),
            ("pressure_atm",   rf"Pressure\s+\.\.\.\s+{FLOAT}\s*atm"),
            ("zpe",            rf"Zero point energy\s+\.\.\.\s+{FLOAT}\s*Eh"),
            ("thermal_corr",   rf"Total thermal correction\s+{FLOAT}\s*Eh"),
            ("enthalpy_corr",  rf"Total Enthalpy\s+\.\.\.\s+{FLOAT}\s*Eh"),
            ("entropy_total",  rf"Final entropy term\s+\.\.\.\s+{FLOAT}\s*Eh"),
            ("gibbs_corr",     rf"G-E\(el\)\s+\.\.\.\s+{FLOAT}\s*Eh"),
            ("final_gibbs",
             rf"Final Gibbs free energy\s+\.\.\.\s+{FLOAT}\s*Eh"),
        ]:
            m = _re.search(pattern, text)
            if m:
                setattr(result, attr, float(m.group(1)))
        return result
    return ThermochemResult(path=str(p))


def extract_energy_table(
    folders: list[str] | str,
    properties: list[str] | None = None,
) -> list[dict]:
    """Walk a list of folders and collect energies into rows.

    Each row is a dict with at least ``folder`` and ``status`` plus
    one entry per requested property. Defaults to
    ``["gibbs", "zpe", "single_point"]``.

    Recognised properties: ``gibbs``, ``zpe``, ``single_point``,
    ``scf_converged``, ``opt_converged``, ``imag_freqs``,
    ``walltime_s``.

    Folders without ORCA output get ``status="no_output"`` and
    ``None`` for every property — easier to filter than missing rows.
    """
    if properties is None:
        properties = ["gibbs", "zpe", "single_point"]
    if isinstance(folders, str):
        folders = [folders]

    from pathlib import Path as _P
    rows: list[dict] = []
    for folder in folders:
        p = _P(folder)
        if not p.exists() or not p.is_dir():
            row: dict = {"folder": str(folder), "status": "missing"}
            for prop in properties:
                row[prop] = None
            rows.append(row)
            continue
        out_files = sorted(p.glob("*.out"))
        if not out_files:
            row = {"folder": str(folder), "status": "no_output"}
            for prop in properties:
                row[prop] = None
            rows.append(row)
            continue
        target = max(out_files, key=lambda f: f.stat().st_size)
        parsed = parse_orca_output(str(target))
        row = {
            "folder": str(folder),
            "status": "ok",
            "output_file": target.name,
        }
        for prop in properties:
            if prop == "gibbs":
                row[prop] = parsed.gibbs_free_energy
            elif prop == "zpe":
                row[prop] = parsed.zpe
            elif prop == "single_point":
                row[prop] = parsed.final_single_point
            elif prop == "scf_converged":
                row[prop] = parsed.scf_converged
            elif prop == "opt_converged":
                row[prop] = parsed.opt_converged
            elif prop == "imag_freqs":
                row[prop] = parsed.imag_freq_count
            elif prop == "walltime_s":
                row[prop] = parsed.walltime_s
            else:
                row[prop] = None
        rows.append(row)
    return rows


def find_calculation_extreme(
    folders: list[str] | str,
    *,
    property: str = "gibbs",
    extreme: str = "min",
    n: int = 5,
) -> list[dict]:
    """Return the N folders with the lowest/highest value of a property.

    Direct answer to "open the .out with the lowest Gibbs energy"
    type questions: pass the candidate folders, get back the top
    ``n`` rows sorted ascending (``extreme="min"``) or descending
    (``extreme="max"``).

    Folders that fail to parse the property are excluded from the
    ranking.

    Args:
        folders: list of paths (or single path) to scan.
        property: one of ``gibbs``, ``zpe``, ``single_point``,
            ``imag_freqs``, ``walltime_s``.
        extreme: ``"min"`` (lowest) or ``"max"`` (highest).
        n: how many top rows to return (clipped to len(rows)).
    """
    rows = extract_energy_table(folders, properties=[property])
    valid = [r for r in rows if r.get(property) is not None]
    if not valid:
        return []
    reverse = (str(extreme).lower() == "max")
    valid.sort(key=lambda r: float(r[property]), reverse=reverse)
    return valid[: max(1, int(n))]


# ---------------------------------------------------------------------------
# On-demand operational-pattern lookup
# ---------------------------------------------------------------------------
#
# Concrete slash-chain recipes for common dashboard workflows.  These
# used to be baked into ``dashboard_agent.md`` but that paid token cost
# every turn — even for conversations that have nothing to do with
# batch jobs or recalc.  Now they live behind a typed lookup the agent
# only calls when it actually needs a workflow recipe.
#
# Adding a new pattern: drop an entry into ``_DASHBOARD_PATTERNS`` and
# add the name to the docstring for discoverability.

_DASHBOARD_PATTERNS: dict[str, str] = {
    "batch": (
        "## Batch jobs (Submit-Tab)\n\n"
        "- Build a batch from EVERY calc-folder's initial geometry:\n"
        "    `ACTION: /batch from-calc`\n"
        "- Build from a glob subset (e.g. all Casagrande folders):\n"
        "    `ACTION: /batch from-calc Casagrande*`\n"
        "- Add one SMILES line manually:\n"
        "    `ACTION: /batch add Name;SMILES;charge=…`\n"
        "- Show / clear the current batch text:\n"
        "    `ACTION: /batch show`  ·  `ACTION: /batch clear`\n\n"
        "**Never** assemble batch text by reading XYZ files yourself —\n"
        "`/batch from-calc` already collects every initial.xyz (or\n"
        "fallback input.txt / coords.xyz) across `calculations/` and\n"
        "writes a properly formatted block into the Submit-Tab textarea.\n\n"
        "Example:\n"
        "  > User: \"Bau einen Batch aus allen XYZ in calc/.\"\n"
        "  > Agent:\n"
        "  >\n"
        "  >     ACTION: /batch from-calc\n"
    ),
    "control_edit": (
        "## CONTROL.txt edits (Submit-Tab)\n\n"
        "- One key change → `ACTION: /control key <key> <value>`\n"
        "- Replace whole content (rare) → `ACTION: /control set <multi-line>`\n"
        "- Validate before submit → `ACTION: /control validate`\n\n"
        "Never paste the full CONTROL into chat — use `/control key` per change.\n"
    ),
    "smart_recalc": (
        "## Smart Recalc (Calc browser)\n\n"
        "Selecting Smart Recalc auto-loads CONTROL into the editor, so\n"
        "the right idiom is `/ui calc-editor replace <old> <new>`:\n\n"
        "  > User: \"smart recalc von Foo/bar mit PAL=1, 1 h\"\n"
        "  > Agent:\n"
        "  >\n"
        "  >     ACTION: /calc cd Foo/bar\n"
        "  >     ACTION: /calc select CONTROL.txt\n"
        "  >     ACTION: /ui calc-options value Smart Recalc\n"
        "  >     ACTION: /ui calc-editor replace PAL=40 PAL=1\n"
        "  >     ACTION: /ui calc-override-time value 01:00:00\n"
        "  >\n"
        "  > Then ASK before submitting:\n"
        "  >\n"
        "  >     ACTION: /ui calc-override-btn click\n\n"
        "`calc-override-btn` is the submit; `calc-override-time` is the\n"
        "time field. Don't confuse with `calc-recalc-btn` /\n"
        "`calc-submit-recalc-btn` — those are NOT the Smart-Recalc panel.\n"
    ),
    "submit_orca": (
        "## Submit a single ORCA job (ORCA Builder)\n\n"
        "- Set fields → `ACTION: /ui orca-method value PBE0`,\n"
        "  `/ui orca-basis value def2-TZVP`, `/ui orca-charge value 0`, …\n"
        "- Switch to the tab → `ACTION: /tab orca`\n"
        "- Submit (after the user explicitly OKs!) → `ACTION: /orca submit`\n"
    ),
    "analyze": (
        "## Analyze existing calculations\n\n"
        "- One folder, full → `ACTION: /analyze <dir>`\n"
        "- Just energies → `ACTION: /analyze energy <dir>`\n"
        "- SCF convergence → `ACTION: /analyze convergence <dir>`\n"
        "- Error scan → `ACTION: /analyze errors <dir>`\n"
        "- All folders overview → `ACTION: /analyze status`\n\n"
        "For multi-folder energy tables prefer\n"
        "`mcp__delfin-ops__extract_energy_table(folders=[...])` — it\n"
        "returns structured data the agent can format directly.\n"
    ),
    "recalc": (
        "## Recalc check / submit\n\n"
        "- Check one folder (safe, read-only) → `ACTION: /recalc check <dir>`\n"
        "- Scan everything (safe) → `ACTION: /recalc check-all`\n"
        "- Submit recalc (DESTRUCTIVE — needs explicit user OK):\n"
        "    `ACTION: /recalc <dir>` and confirm\n"
        "- Bulk auto-recalc only on explicit \"alle neuberechnen\":\n"
        "    `ACTION: /recalc auto`\n"
    ),
    "cancel": (
        "## Cancel jobs\n\n"
        "- One job → `ACTION: /cancel <job_id>` (after user OK)\n"
        "- All — only on explicit \"cancel all\":\n"
        "    `ACTION: /cancel all`\n"
    ),
}


def list_dashboard_patterns() -> list[str]:
    """Return the names of every operational pattern available for lookup.

    Use this when the agent isn't sure which pattern matches the user's
    request — a names-only list is cheap and lets the agent pick.
    """
    return sorted(_DASHBOARD_PATTERNS)


def get_dashboard_pattern(name: str) -> str:
    """Return the slash-chain recipe for a named dashboard workflow.

    Available names: ``batch``, ``control_edit``, ``smart_recalc``,
    ``submit_orca``, ``analyze``, ``recalc``, ``cancel``.

    Names are matched case-insensitively. An unknown name returns a
    short error string with the available choices so the agent can
    self-correct.
    """
    if not name:
        return (
            "No pattern requested. Call with one of: "
            + ", ".join(list_dashboard_patterns())
        )
    key = str(name).strip().lower().replace("-", "_").replace(" ", "_")
    text = _DASHBOARD_PATTERNS.get(key)
    if text is None:
        return (
            f"Unknown pattern: {name!r}. Available: "
            + ", ".join(list_dashboard_patterns())
        )
    return text


# ---------------------------------------------------------------------------
# Module surface
# ---------------------------------------------------------------------------

__all__ = [
    "CommandResult",
    "OrcaParseResult",
    "OrcaError",
    "ThermochemResult",
    "run",
    "prepare",
    "pipeline_run",
    "pipeline_prepare",
    "qm_check",
    "csp_check",
    "mlp_check",
    "analysis_check",
    "qm_run",
    "cleanup",
    "stop",
    "run_orca_input",
    "co2",
    "tadf_xtb",
    "hyperpol",
    "list_dashboard_patterns",
    "get_dashboard_pattern",
    # P1 — output parsing
    "parse_orca_output",
    "find_orca_errors",
    "extract_thermochem",
    "extract_energy_table",
    "find_calculation_extreme",
]
