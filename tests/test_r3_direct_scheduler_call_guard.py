"""R3 phase 3: no DIRECT scheduler call may exist outside the allow-list.

Finding 1's throttle is only worth anything if a future module cannot simply
add a new ``subprocess.run(["squeue", ...])`` somewhere else and bypass the
scheduler client. This guard is an AST scan over every ``.py`` under
``delfin/``: a DIRECT call is one whose first argument is a list-literal (or a
bare string) whose first element equals a scheduler command name, i.e. the
command is spelled out at the call site and executed directly rather than
routed through :func:`delfin.scheduler_client.query`.

Every file the scanner finds must be in :data:`ALLOWED`. Those files are the
sanctioned scheduler callers:

* ``job_monitor`` — routed through the throttled scheduler_client (owned);
* the others — pre-existing callers outside this package's ownership that
  migrate separately (listed for the operator) or use the command once
  (``--wait`` farm helper, ``sinfo`` probe).

Anything not in :data:`ALLOWED` fails the scan. A bypass that hides the
command in a variable (``[bin, ...]`` bound from ``shutil.which("sacct")``)
or concatenates it is deliberately out of scope here — those are the
operator's migration files — and the scanner documents that gap.
"""
from __future__ import annotations

import ast
from pathlib import Path

# The scheduler commands whose throttling finding 1 is about.
_COMMANDS = {"squeue", "sacct", "scontrol", "sbatch", "sinfo"}

# Files that are allowed to run a scheduler command directly. Every detected
# file MUST be in here; a file without a detected call must NOT be here.
ALLOWED = {
    "delfin/agent/job_monitor.py",
    "delfin/agent/pack/benchmark/accept/chem_opt_is_a_minimum.py",
    "delfin/cluster_utils.py",
    "delfin/dashboard/backend_slurm.py",
    "delfin/dashboard/tab_calculations_browser.py",
    "delfin/mlp_tools/__init__.py",
    "delfin/slurm_submit.py",
}


def _is_command_token(node) -> bool:
    """A bare string (or f-string constant) equal to a scheduler command."""
    if isinstance(node, ast.Constant):
        return isinstance(node.value, str) and node.value in _COMMANDS
    if isinstance(node, ast.JoinedStr) and node.values:
        first = node.values[0]
        return isinstance(first, ast.Constant) and first.value in _COMMANDS
    return False


def _call_func_name(node: ast.Call) -> str:
    """The plain name of the function being called, e.g. 'which' or 'run'."""
    if isinstance(node.func, ast.Name):
        return node.func.id
    if isinstance(node.func, ast.Attribute):
        return node.func.attr
    return ""


def _directly_runs(node: ast.Call) -> bool:
    """True if ``node`` runs a scheduler command spelled out at the call.

    ``shutil.which("squeue")`` and ``_sh.which(...)`` are *presence probes*,
    not runs — they resolve a path, they never execute the command, so they
    are not the per-widget query load finding 1 objects to. ``which`` calls
    are therefore not classified as direct runs.
    """
    if _call_func_name(node) == "which":
        return False
    if not node.args:
        return False
    first = node.args[0]
    if isinstance(first, ast.List) and first.elts and _is_command_token(first.elts[0]):
        return True
    return _is_command_token(first)  # bare subprocess.run("squeue")


def _scan_source(text: str) -> list[int]:
    """Line numbers of direct scheduler calls in ``text``."""
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return []
    hits = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and _directly_runs(node):
            hits.append(node.lineno)
    return sorted(hits)


def _scan_file(path: Path) -> list[int]:
    try:
        return _scan_source(path.read_text())
    except OSError:
        return []


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def test_scanner_flags_a_literal_squeue_call():
    src = "subprocess.run(['squeue', '-u', 'me'], text=True)"
    assert _scan_source(src) == [1]


def test_scanner_flags_bare_string_and_slayout_forms():
    assert _scan_source("run_fn(['scontrol', 'show', 'job', x])") == [1]
    assert _scan_source("subprocess.run('squeue', shell=True)") == [1]


def test_scanner_does_not_flag_docstrings_or_regex_mentions():
    src = (
        "def f():\n"
        "    \"\"\"the squeue state codes\"\"\"\n"
        "    pat = re.compile(r'^\\s*(?:squeue|sacct)\\b')\n"
        "    return pat\n"
    )
    assert _scan_source(src) == []


def test_scanner_does_not_flag_which_presence_probe():
    # shutil.which(...) is a presence probe, not a run — not throttled.
    assert _scan_source("_sh.which('squeue')") == []
    assert _scan_source("shutil.which('sbatch') and shutil.which('squeue')") == []


def test_scanner_does_not_flag_bash_strings_not_executed():
    # A here-doc / template built for another process is not a direct call.
    src = "script = '\\n'.join(['squeue -u %s' % u, '#SBATCH --time=1'])\n"
    assert _scan_source(src) == []


def test_whole_tree_has_no_direct_call_outside_allow_list():
    """The guard: every file with a direct scheduler call is allow-listed,
    and every allow-listed file actually has one (no bloat)."""
    root = _repo_root() / "delfin"
    detected: dict[str, list[int]] = {}
    for p in sorted(root.rglob("*.py")):
        hits = _scan_file(p)
        if hits:
            detected[str(p.relative_to(_repo_root()))] = hits

    allowed_with_calls = {f for f in ALLOWED if f in detected}
    assert detected.keys() <= ALLOWED, (
        "direct scheduler call outside the allow-list — route it through "
        "delfin.scheduler_client instead: "
        f"{sorted(detected.keys() - ALLOWED)}"
    )
    # No bloat: an allow-listed file with no detected call means the guard
    # would silently accept a future direct call in it.
    assert allowed_with_calls == ALLOWED, (
        "allow-list bloat: no direct call detected in "
        f"{sorted(ALLOWED - allowed_with_calls)}"
    )
