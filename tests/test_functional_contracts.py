"""Functional contract tests — guard that user-visible DELFIN behavior
doesn't silently change under refactoring.

Where `test_no_regression_undefined_names.py` says "the code imports
and doesn't reference undefined names", THIS file says "the functions
users rely on still have the expected signatures and still behave the
way users have learned to expect". Each contract is one pytest case,
so the failure message points at exactly which guarantee was broken.

Motivated by the refactor bug wave: every large commit on main
(`cf78545`, `3ed6566`, `e6761e4`, `cceea3e`, `39b1919`, …) silently
broke a code path users depend on. The affected users (Tilman, Jan)
only noticed after queued SLURM jobs failed hours later. These tests
shift that discovery from "user-hits-it-in-production" to
"pre-commit".

Contracts pinned here (add more as surprises come up):

  1. CLI entry points still import and expose a callable `main`.
  2. Critical functions have the signature callers assume.
  3. Core workflow primitives accept the kwargs used at call sites.
  4. Load-settings survives legacy / minimal settings.json shapes.
  5. Default workflow shapes match what the CLI consumes.
"""

from __future__ import annotations

import importlib
import inspect
from pathlib import Path

import pytest


# ===========================================================================
# 1. CLI entry points are importable and callable
# ===========================================================================
#
# If pyproject.toml lists `delfin-foo = "module:func"` but `module` can't
# be imported or `func` isn't callable, `pip install -e .` still appears
# to work — the break only surfaces when the user runs the command.

_ENTRY_POINTS = {
    "delfin": ("delfin.main", "main"),
    "delfin_ESD": ("delfin.cli_esd_report", "main"),
    "delfin_IR": ("delfin.cli_ir_report", "main"),
    "delfin_NMR": ("delfin.cli_nmr_report", "main"),
    "delfin-json": ("delfin.cli_delfin_collect", "main"),
    "delfin-build": ("delfin.build_up_complex", "main"),
    "delfin-guppy": ("delfin.guppy_sampling", "main"),
    "delfin-voila": ("delfin.cli_voila", "main"),
    "delfin-step": ("delfin.cli_step", "main"),
    "delfin-pipeline": ("delfin.cli_pipeline", "main"),
    "delfin-docs-index": ("delfin.doc_server.indexer", "main"),
}


@pytest.mark.parametrize(
    "cli_name,target",
    sorted(_ENTRY_POINTS.items()),
    ids=[name for name in sorted(_ENTRY_POINTS.keys())],
)
def test_cli_entry_point_is_importable_and_callable(cli_name, target):
    """Every CLI entry point declared in pyproject.toml must resolve.

    This catches:
    - `pyproject.toml` was updated to point at a renamed module
    - `module:main` was renamed to `module:run` without updating pyproject
    - A refactor deleted the module entirely
    """
    module_name, func_name = target
    try:
        module = importlib.import_module(module_name)
    except Exception as exc:
        pytest.fail(
            f"CLI entry point '{cli_name}' targets {module_name}:{func_name}, "
            f"but module failed to import: {type(exc).__name__}: {exc}"
        )
    func = getattr(module, func_name, None)
    assert func is not None, (
        f"CLI entry point '{cli_name}' targets {module_name}:{func_name}, "
        f"but {func_name} is missing from {module_name}."
    )
    assert callable(func), (
        f"{module_name}.{func_name} exists but is not callable. "
        f"CLI entry point '{cli_name}' is broken."
    )


# ===========================================================================
# 2. Critical function signatures
# ===========================================================================
#
# A signature contract is stronger than "the name still exists": it pins
# the parameter names and whether they are keyword-only / positional. If
# a refactor renames a parameter from `config=` to `settings=`, every
# kwargs-call at the boundary breaks. This section locks those in.

_SIGNATURE_CONTRACTS = [
    # (module, function, required_params_in_order, optional_params)
    (
        "delfin.orca",
        "run_orca",
        ["input_file_path", "output_log"],
        ["timeout", "scratch_subdir", "working_dir", "isolate", "copy_files"],
    ),
    (
        "delfin.orca",
        "run_orca_with_intelligent_recovery",
        ["input_file_path", "output_log"],
        ["timeout", "scratch_subdir", "working_dir", "isolate", "copy_files", "config"],
    ),
    (
        "delfin.user_settings",
        "load_settings",
        [],
        ["settings_path"],
    ),
    (
        "delfin.smart_recalc",
        "should_skip",
        ["inp_path", "out_path"],
        ["extra_deps", "required_outputs"],
    ),
    (
        "delfin.common.logging",
        "get_logger",
        ["name"],
        [],
    ),
    (
        "delfin.workflows.scheduling.manager",
        "get_global_manager",
        [],
        [],
    ),
]


@pytest.mark.parametrize(
    "module_name,func_name,required,optional",
    _SIGNATURE_CONTRACTS,
    ids=[f"{m}.{f}" for m, f, _, _ in _SIGNATURE_CONTRACTS],
)
def test_function_signature_contract(module_name, func_name, required, optional):
    """Pinned signatures for functions other code / users call.

    If a refactor renames a parameter, this test fails BEFORE a user
    gets a TypeError. Example: run_orca's first arg must stay called
    `input_file_path`, not `input_path` or `inp_file`, because callers
    and documentation use that name.
    """
    module = importlib.import_module(module_name)
    func = getattr(module, func_name)
    sig = inspect.signature(func)

    params = list(sig.parameters.values())
    # Drop `self` / `cls` if present (shouldn't be — these are module-level)
    params = [p for p in params if p.name not in ("self", "cls")]

    # Required params: the first N named parameters must match names & order
    for i, want in enumerate(required):
        assert i < len(params), (
            f"{module_name}.{func_name}: expected {len(required)} required "
            f"params starting with {required}, but got only {len(params)}: "
            f"{[p.name for p in params]}"
        )
        assert params[i].name == want, (
            f"{module_name}.{func_name}: required param #{i} is "
            f"'{params[i].name}', contract says '{want}'. "
            f"Full signature: {sig}"
        )

    # Optional params: just require they exist somewhere in the signature
    param_names = {p.name for p in params}
    missing = [o for o in optional if o not in param_names]
    assert not missing, (
        f"{module_name}.{func_name}: optional params {missing} are gone. "
        f"Either add back or update the contract in "
        f"tests/test_functional_contracts.py. Full signature: {sig}"
    )


# ===========================================================================
# 3. WorkflowJob has the fields callers rely on
# ===========================================================================
#
# Many call sites construct WorkflowJob(...) with specific kwargs. If a
# field is renamed, only the unit tests that construct a WorkflowJob
# catch it — but most construction sites are inside pipeline code that
# only runs under a real CONTROL.txt. Pin the dataclass fields here.

_WORKFLOWJOB_FIELDS = {
    "job_id",
    "work",
    "description",
    "dependencies",
    "cores_min",
    "cores_optimal",
    "cores_max",
    "priority",
    "memory_mb",
    "estimated_duration",
    "inline",
    "preserve_cores_optimal",
    "working_dir",
    "precomplete_check",
}


def test_workflowjob_fields_are_stable():
    """WorkflowJob's fields must not silently disappear.

    Pipelines all over the code base construct WorkflowJob(
    job_id=..., work=..., description=..., cores_min=..., ...
    ). Losing a field means every caller fails at runtime when that
    kwarg is supplied.
    """
    from delfin.workflows.engine.classic import WorkflowJob
    import dataclasses

    actual_fields = {f.name for f in dataclasses.fields(WorkflowJob)}
    missing = _WORKFLOWJOB_FIELDS - actual_fields
    assert not missing, (
        f"WorkflowJob is missing contract fields {sorted(missing)}. "
        f"Actual: {sorted(actual_fields)}. Either restore the fields "
        f"or update _WORKFLOWJOB_FIELDS after a deliberate API change."
    )


def test_workflowjob_accepts_all_contract_kwargs():
    """Construct a WorkflowJob with every contract kwarg — catches kwarg
    renames that dataclass default-handling hides."""
    from delfin.workflows.engine.classic import WorkflowJob
    from delfin.workflows.scheduling.pool import JobPriority

    # Minimal viable constructor, exercising every public kwarg
    job = WorkflowJob(
        job_id="contract_test",
        work=lambda cores: None,
        description="contract test",
        dependencies=set(),
        cores_min=1,
        cores_optimal=1,
        cores_max=1,
        priority=JobPriority.NORMAL,
        memory_mb=1000,
        estimated_duration=10.0,
        inline=False,
        preserve_cores_optimal=False,
        working_dir=None,
        precomplete_check=None,
    )
    assert job.job_id == "contract_test"


# ===========================================================================
# 4. Settings loader tolerates legacy shapes
# ===========================================================================
#
# user_settings.py:382 was the KeyError: 'docs' bug — load_settings
# crashed on a settings file that didn't have the (newer) docs section.
# This contract guards the general case: legacy / minimal settings
# files must still load cleanly.

def test_load_settings_accepts_minimal_settings_file(tmp_path):
    """A settings file with only one top-level key must load cleanly.

    Reproduces the 39b1919 regression where load_settings crashed on
    KeyError: 'docs' when a legacy settings file lacked the docs section.
    """
    import json
    from delfin.user_settings import load_settings

    path = tmp_path / "settings.json"
    path.write_text(json.dumps({"transfer": {"host": "example.org"}}))
    loaded = load_settings(path)
    assert loaded["transfer"]["host"] == "example.org"
    # Defaults were filled in for missing sections
    assert "runtime" in loaded
    assert "docs" in loaded


def test_load_settings_accepts_empty_settings_file(tmp_path):
    """An empty {} settings file must load cleanly (all defaults)."""
    import json
    from delfin.user_settings import load_settings

    path = tmp_path / "settings.json"
    path.write_text(json.dumps({}))
    loaded = load_settings(path)
    assert isinstance(loaded, dict)
    # Core sections always present
    for section in ("runtime", "features", "docs"):
        assert section in loaded, (
            f"load_settings dropped default section '{section}'. "
            f"Legacy settings files relying on this will crash."
        )


def test_load_settings_accepts_nonexistent_path(tmp_path):
    """Settings loader must tolerate missing file (fresh install)."""
    from delfin.user_settings import load_settings

    path = tmp_path / "does_not_exist.json"
    loaded = load_settings(path)
    assert isinstance(loaded, dict)


# ===========================================================================
# 5. pyproject.toml entry points stay in sync with reality
# ===========================================================================

_PYPROJECT = Path(__file__).resolve().parents[1] / "pyproject.toml"


def test_pyproject_entry_points_all_declared_in_test():
    """Any new CLI entry point added to pyproject.toml must also be
    added to _ENTRY_POINTS above, so it gets covered by the per-CLI
    contract test.
    """
    if not _PYPROJECT.exists():
        pytest.skip("pyproject.toml not found")
    try:
        import tomllib  # Python 3.11+
    except ImportError:
        try:
            import tomli as tomllib  # type: ignore
        except ImportError:
            pytest.skip("tomllib/tomli not available")

    data = tomllib.loads(_PYPROJECT.read_text())
    declared = set((data.get("project") or {}).get("scripts", {}).keys())
    tested = set(_ENTRY_POINTS.keys())

    new_untested = declared - tested
    assert not new_untested, (
        f"pyproject.toml declares CLI entry point(s) {sorted(new_untested)} "
        f"that are NOT covered by test_cli_entry_point_is_importable_and_callable. "
        f"Add them to _ENTRY_POINTS in tests/test_functional_contracts.py."
    )
    stale_in_test = tested - declared
    assert not stale_in_test, (
        f"_ENTRY_POINTS references {sorted(stale_in_test)} but pyproject.toml "
        f"no longer declares them. Remove from the test."
    )


# ===========================================================================
# 6. Smart recalc: skip behavior contracts
# ===========================================================================
#
# Smart recalc is the feature that lets users re-submit a DELFIN job
# and skip already-finished ORCA steps. If should_skip() semantics
# change silently, users end up either re-running the whole pipeline
# (wasting hours) or skipping real work (corrupting outputs).

def test_should_skip_returns_false_when_output_missing(tmp_path):
    """If the output file doesn't exist, we must NOT skip."""
    from delfin.smart_recalc import should_skip

    inp = tmp_path / "job.inp"
    inp.write_text("! B3LYP def2-SVP\n* xyz 0 1\nH 0 0 0\n*\n")
    out = tmp_path / "job.out"  # does NOT exist
    assert should_skip(inp, out) is False, (
        "should_skip returned True but the output file doesn't exist — "
        "users would lose their calculation output."
    )


def test_should_skip_returns_false_when_input_missing(tmp_path):
    """Missing input shouldn't claim success."""
    from delfin.smart_recalc import should_skip

    inp = tmp_path / "never.inp"  # does NOT exist
    out = tmp_path / "never.out"
    assert should_skip(inp, out) is False


# ===========================================================================
# 7. SLURM wrapper script still has the post-fix shape
# ===========================================================================
#
# These are belt-and-suspenders for the set -e / _postprocess_nmr fixes.
# test_submit_wrapper.py already tests the behaviour — here we also
# check that the script on disk still has the fix lines, so a bad
# merge can't quietly revert it.

_SUBMIT_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "submit_templates"
    / "submit_delfin.sh"
)


def test_submit_script_still_has_set_plus_e_guard():
    """The `set +e` before `wait $_DELFIN_PID` protects the post-ORCA
    sync from being killed by errexit when ORCA exits non-zero. If this
    line vanishes, every failed-ORCA job loses its output sync again.
    """
    text = _SUBMIT_SCRIPT.read_text()
    assert "set +e" in text and 'wait "$_DELFIN_PID"' in text, (
        "set +e / wait combination missing from submit_delfin.sh — "
        "the post-ORCA cleanup can die silently again."
    )
    # Also check order: set +e must appear before the wait
    assert text.index("set +e") < text.index('wait "$_DELFIN_PID"'), (
        "set +e appears AFTER wait — the guard is in the wrong place."
    )


def test_submit_script_has_completion_marker():
    """The wrapper must print a completion marker at the very end.
    Without it, an aborted wrapper looks identical to a normal exit
    and users can't tell why their files are missing."""
    text = _SUBMIT_SCRIPT.read_text()
    assert "DELFIN WRAPPER COMPLETED" in text, (
        "Completion marker line is gone — silent wrapper deaths will be "
        "indistinguishable from successful runs in delfin_*.out."
    )


def test_submit_script_sets_pythonunbuffered():
    """Python stdout must stay unbuffered on SLURM so errors appear
    immediately in the job log rather than being lost when the process
    dies before flushing."""
    text = _SUBMIT_SCRIPT.read_text()
    assert "PYTHONUNBUFFERED" in text, (
        "PYTHONUNBUFFERED export is missing — Python crashes in "
        "local_runner will be invisible in delfin_*.out (Tilman's "
        "original symptom)."
    )


def test_submit_script_rsync_failures_are_logged():
    """rsync failures must print a WARNING, not be swallowed silently."""
    text = _SUBMIT_SCRIPT.read_text()
    # The fix uses a captured rsync_rc and prints a warning
    assert "rsync_rc" in text, (
        "rsync exit code is no longer captured — a full home partition "
        "or NFS hiccup will cause silent partial syncs."
    )
    assert "WARNING: rsync exited" in text, (
        "rsync failure warning is gone — silent data loss possible."
    )
