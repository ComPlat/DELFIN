"""The session budget rides on a real submit_calculation through the gate.

Package BG1 wave-28, wiring deliverable. The pre-submit hook lives in
``calc_policy.judge_delfin_submission`` (budget/job/client/run_fn, built by
the operator from the approved .gate/px pair); this file is the RED CONTROL
for the call site: ``api_client._mcp_calculation_verdict`` called it WITHOUT
budget/job, so a real ops ``submit_calculation`` (mode=delfin) never reached
the over-budget verdict (operator, real-run fail PR #622). The fix passes
the session Budget (from the user's settings, opt-in) and the job dict
(from the folder's CONTROL.txt + geometry) at that one call site.

The path is the real one: ``_doc_executor._gate_mcp_tool`` ->
``_gate_side_effecting_mcp_tool`` -> ``_mcp_calculation_verdict`` ->
``calc_policy.judge_delfin_submission``. Over budget, the reason reaches
the model through the dialog path: the confirm preview carries
"Not started unasked: <calc_reason>" (api_client.py _gate_side_effecting_
mcp_tool), so these tests run with a confirm_callback that REFUSES and
assert the budget reason in the refusal — the head-less generic refusal
would mask it (operator defect report, 2026-11-11).
"""
from __future__ import annotations

import re
import json

import pytest

from delfin.agent import api_client as A
from delfin.agent import calc_policy as C

WATER = "O 0.0 0.0 0.0\nH 0.0 0.76 0.59\nH 0.0 -0.76 0.59\n"


def _control(pal=2):
    from delfin.define import TEMPLATE
    t = TEMPLATE.replace("charge=[CHARGE]", "charge=0")
    t = t.replace("method=[classic|manually|OCCUPIER]", "method=classic")
    t = t.replace("[SOLVENT]", "water")
    t = re.sub(r"^smiles_converter=.*$", "smiles_converter=QUICK", t, flags=re.M)
    t = re.sub(r"^PAL=.*$", f"PAL={pal}", t, flags=re.M)
    return t


@pytest.fixture
def cluster(monkeypatch):
    monkeypatch.setattr(C, "host_kind", lambda env=None: "login_node")
    monkeypatch.setattr(C, "partition_limits", lambda name: {
        "PartitionName": "cpu", "MaxTime": "3-00:00:00",
        "MaxCPUsPerNode": "64"})
    monkeypatch.setattr(C, "queued_jobs", lambda: 0)


@pytest.fixture
def no_classifier(monkeypatch):
    executor = A._doc_executor
    monkeypatch.setattr(executor, "_classifier_second_opinion",
                        lambda perms, name, args, mode: None)


@pytest.fixture
def ws(tmp_path, monkeypatch):
    """The session's workspace: a project dir under the redirected home.

    ``KitToolPermissions`` refuses $HOME and system roots as a workspace,
    so the workspace is a project sub-directory of the tmp home, and the
    calculation folder a sub-directory of that (matching the real shape:
    submit_calculation gets a folder under the workspace).
    """
    home = tmp_path / "home"
    (home / ".delfin").mkdir(parents=True)
    (home / ".delfin" / "settings.json").write_text("{}")
    monkeypatch.setenv("HOME", str(home))
    ws = home / "agent_workspace" / "proj"
    ws.mkdir(parents=True)
    return ws


def _folder(ws):
    folder = ws / "water_job"
    folder.mkdir(exist_ok=True)
    (folder / "CONTROL.txt").write_text(_control())
    (folder / "input.txt").write_text(WATER)
    return folder


def _set_budget(monkeypatch, tmp_path, budget=None, payload=None):
    """Write the budget key into the redirected home's settings file.

    DELFIN's user settings live in ~/.delfin_settings.json
    (user_settings.py:16 SETTINGS_FILE_NAME, :518 get_settings_path =
    Path.home()/SETTINGS_FILE_NAME) — NOT ~/.delfin/settings.json, a
    path v3 invented and load_settings() never read (operator v3 run:
    budget never seen, verdict silently (None, None)).
    """
    home = tmp_path / "home"
    home.mkdir(parents=True, exist_ok=True)
    if payload is None:
        payload = {}
        if budget is not None:
            payload = {"agent": {"calculation_budget_core_hours": budget}}
    from delfin import user_settings
    monkeypatch.setattr(
        user_settings, "get_settings_path", lambda base=None: home /
        user_settings.SETTINGS_FILE_NAME)
    user_settings.save_settings(payload, home /
                                user_settings.SETTINGS_FILE_NAME)


def _gate(ws, args):
    """The dialog path: a confirm_callback that refuses, as a user who
    reads the preview and says no. The budget reason must be IN the
    refusal the callback path produces; the head-less generic refusal
    (confirm_callback=None) hides it."""
    perms = A.KitToolPermissions(workspace=ws, mode="default",
                                 agent_role="solo_agent")
    seen = {}
    perms.confirm_callback = lambda name, args, preview="": (
        seen.update(preview=preview) or False)
    out = A._doc_executor._gate_mcp_tool(
        "mcp__delfin-ops__submit_calculation", args, perms)
    return out, seen


def test_over_budget_submission_is_refused_with_the_estimate(
        ws, tmp_path, monkeypatch, cluster, no_classifier):
    """RED until the wiring: the call site passes the session Budget."""
    folder = _folder(ws)
    _set_budget(monkeypatch, tmp_path, budget=0.001)
    out, seen = _gate(ws, {
        "folder": str(folder), "mode": "delfin", "pal": 2,
        "maxcore": 2000, "time_limit": "12:00:00"})
    # The budget reason lives in the PREVIEW the callback is shown
    # (api_client.py:24278-24280); the refusal text itself is the
    # callback's own "user denied" answer (24295-24298). Assert on the
    # preview (operator defect 2). Red (unpatched): the verdict admits
    # before the dialog, seen stays {} and the failure says so.
    reasons = f"{seen.get('preview', '')} {out}"
    assert "over budget -> ask the person" in reasons, (
        f"over-budget submission was admitted unasked (no dialog with the "
        f"budget estimate was shown: preview={seen.get('preview')!r}, "
        f"out={out!r})")
    assert "core-hours" in reasons, reasons   # the estimate is shown


def test_without_a_budget_setting_the_behaviour_is_unchanged(
        ws, tmp_path, monkeypatch, cluster, no_classifier):
    """Opt-in: no budget key -> today's admission path, no over-budget text."""
    folder = _folder(ws)
    _set_budget(monkeypatch, tmp_path, budget=None)
    out, seen = _gate(ws, {
        "folder": str(folder), "mode": "delfin", "pal": 2,
        "maxcore": 2000, "time_limit": "12:00:00"})
    assert out is None, out
    assert "over budget" not in str(seen.get("preview", ""))


def test_within_a_generous_budget_the_job_is_admitted(
        ws, tmp_path, monkeypatch, cluster, no_classifier):
    folder = _folder(ws)
    _set_budget(monkeypatch, tmp_path, budget=1.0e6)
    out, seen = _gate(ws, {
        "folder": str(folder), "mode": "delfin", "pal": 2,
        "maxcore": 2000, "time_limit": "12:00:00"})
    assert out is None, out
    assert "over budget" not in str(seen.get("preview", ""))


def test_the_budget_never_silently_shrinks_the_request(
        ws, tmp_path, monkeypatch, cluster, no_classifier):
    """Over budget refuses; it does not rewrite PAL/maxcore/time."""
    folder = _folder(ws)
    _set_budget(monkeypatch, tmp_path, budget=0.001)
    out, seen = _gate(ws, {
        "folder": str(folder), "mode": "delfin", "pal": 2,
        "maxcore": 2000, "time_limit": "12:00:00"})
    reasons = f"{seen.get('preview', '')} {out}"
    assert "over budget -> ask the person" in reasons, (
        f"over-budget submission was admitted unasked (no dialog with the "
        f"budget estimate was shown: preview={seen.get('preview')!r}, "
        f"out={out!r})")
    assert "core-hours" in reasons, reasons
    assert "is not on the auto-allow list" not in reasons, reasons


# -- operator review 2026-11-10: three defects, red controls first ----------


def test_basis_functions_count_from_the_basis_definition(
        ws, tmp_path, monkeypatch, cluster, no_classifier):
    """H2O def2-SVP must count 24 basis functions (ORCA-printed).

    RED: the pre-review helper estimated electrons/2 = 5. A budget fed a
    5x-too-small size admits jobs it must refuse.
    """
    from delfin.agent.api_client import _session_budget_verdict_args
    _set_budget(monkeypatch, tmp_path, budget=0.001)
    folder = _folder(ws)
    _budget_args, job = _session_budget_verdict_args(folder)
    assert job is not None, "water def2-SVP must be sizeable"
    assert job["basis_functions"] == 24, job


def test_basis_functions_count_triple_zeta(
        ws, tmp_path, monkeypatch, cluster, no_classifier):
    """H2O def2-TZVP counts 43 (ORCA-printed)."""
    from delfin.agent.api_client import _session_budget_verdict_args
    _set_budget(monkeypatch, tmp_path, budget=0.001)
    folder = ws / "water_job_tz"
    folder.mkdir(exist_ok=True)
    (folder / "CONTROL.txt").write_text(
        _control().replace("main_basisset=def2-SVP",
                           "main_basisset=def2-TZVP"))
    (folder / "input.txt").write_text(WATER)
    _budget_args, job = _session_budget_verdict_args(folder)
    assert job is not None, "water def2-TZVP must be sizeable"
    assert job["basis_functions"] == 43, job


def test_basis_functions_count_methane_svp(
        ws, tmp_path, monkeypatch, cluster, no_classifier):
    """CH4 def2-SVP counts 34 (ORCA-printed)."""
    from delfin.agent.api_client import _session_budget_verdict_args
    _set_budget(monkeypatch, tmp_path, budget=0.001)
    folder = ws / "ch4_job"
    folder.mkdir(exist_ok=True)
    (folder / "CONTROL.txt").write_text(_control())
    (folder / "input.txt").write_text(
        "C 0.0 0.0 0.0\nH 0.63 0.63 0.63\nH -0.63 -0.63 0.63\n"
        "H -0.63 0.63 -0.63\nH 0.63 -0.63 -0.63\n")
    _budget_args, job = _session_budget_verdict_args(folder)
    assert job is not None, "methane def2-SVP must be sizeable"
    assert job["basis_functions"] == 34, job


def test_an_unknown_element_or_basis_gives_no_budget_not_a_guess(
        ws, tmp_path, monkeypatch, cluster, no_classifier):
    """An element/basis-function table miss -> no budget."""
    from delfin.agent.api_client import _session_budget_verdict_args
    _set_budget(monkeypatch, tmp_path, budget=0.001)
    folder = ws / "ru_job"
    folder.mkdir(exist_ok=True)
    (folder / "CONTROL.txt").write_text(_control())
    (folder / "input.txt").write_text(
        "Ru 0.0 0.0 0.0\nH 1.7 0.0 0.0\n")
    _budget_args, job = _session_budget_verdict_args(folder)
    assert job is None, "unknown element/basis must not be guessed"


def test_partition_comes_from_the_submit_path_s_source(
        ws, tmp_path, monkeypatch, cluster, no_classifier):
    """The job's partition is what DELFIN's submit path would use
    (runtime.slurm.partitions in the settings), never a hard-coded 'cpu'."""
    from delfin.agent.api_client import _session_budget_verdict_args
    _set_budget(monkeypatch, tmp_path, payload={
        "agent": {"calculation_budget_core_hours": 0.001},
        "runtime": {"slurm": {"partitions": "cpu_short, cpu"}}})
    folder = _folder(ws)
    _budget_args, job = _session_budget_verdict_args(folder)
    assert job is not None
    assert job["partition"] == "cpu_short", job


def test_a_skipped_budget_is_logged_not_silent(
        ws, tmp_path, monkeypatch, cluster, no_classifier, caplog):
    """A user who set a budget sees one line when it was not applied."""
    import logging
    from delfin.agent.api_client import _session_budget_verdict_args
    _set_budget(monkeypatch, tmp_path, budget=0.001)
    folder = ws / "broken_job"
    folder.mkdir(exist_ok=True)
    (folder / "CONTROL.txt").write_text(_control())
    (folder / "input.txt").write_text("Ru 0.0 0.0 0.0\n")  # unsizable
    with caplog.at_level(logging.WARNING, logger="delfin.agent.api_client"):
        _budget_args, _job = _session_budget_verdict_args(folder)
    assert _budget_args is None
    assert any("calculation_budget" in r.message for r in caplog.records), \
        caplog.text


def test_def2_variant_prefixes_are_not_billed_at_the_family_counts(
        ws, tmp_path, monkeypatch, cluster, no_classifier):
    """F-A: def2-SVPD/-TZVPD/-TZVPP/-QZVP carry MORE functions than the
    family their name starts with; billing them at the SVP/TZVP counts
    UNDER-estimates a runaway calculation. Exact family match: they get
    NO budget, never a guess. RED pre-patch: ImportError.
    """
    from delfin.agent.api_client import _session_budget_verdict_args
    _set_budget(monkeypatch, tmp_path, budget=0.001)
    for basis in ("def2-SVPD", "def2-TZVPD", "def2-TZVPP", "def2-QZVP"):
        folder = ws / f"var_job_{basis}"
        folder.mkdir(exist_ok=True)
        (folder / "CONTROL.txt").write_text(
            _control().replace("main_basisset=def2-SVP",
                               f"main_basisset={basis}"))
        (folder / "input.txt").write_text(WATER)
        _budget_args, job = _session_budget_verdict_args(folder)
        assert job is None, (
            f"{basis} starts with a smaller family's name; it must get "
            f"no budget, not the family's counts (got {job})")


def test_smiles_control_gets_no_budget_and_the_skip_is_logged(
        ws, tmp_path, monkeypatch, cluster, no_classifier, caplog):
    """F-B pinned: a CONTROL with SMILES= and no geometry file feeds the
    SMILES string to geometry_text_symbols -> no symbols -> no budget.
    Fail-safe direction (no budget rather than a wrong one), but it must
    be VISIBLE: logged with the reason, never silent. RED pre-patch:
    ImportError.
    """
    import logging
    from delfin.agent.api_client import _session_budget_verdict_args
    _set_budget(monkeypatch, tmp_path, budget=0.001)
    folder = ws / "smiles_job"
    folder.mkdir(exist_ok=True)
    (folder / "CONTROL.txt").write_text(
        _control().replace("input.txt", "input.txt")
        + "SMILES=O\n")
    with caplog.at_level(logging.WARNING, logger="delfin.agent.api_client"):
        _budget_args, job = _session_budget_verdict_args(folder)
    assert job is None, (
        "a SMILES string is not xyz text; the budget must be skipped "
        f"(got {job})")
    assert any("calculation_budget" in r.message for r in caplog.records), \
        caplog.text
