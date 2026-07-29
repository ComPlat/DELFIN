"""Tests for the aggregate environment health report (agent doctor)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli
from delfin.agent import doctor


def _ctx(**over) -> dict:
    ctx = {"workspace": "", "settings": {}, "fast": True}
    ctx.update(over)
    return ctx


# ---------------------------------------------------------------------------
# run_doctor aggregation — each check monkeypatched to a chosen status
# ---------------------------------------------------------------------------


def test_run_doctor_aggregates_stubbed_checks(monkeypatch):
    cycle = ["PASS", "WARN", "FAIL"]
    statuses = [cycle[i % 3] for i in range(len(doctor._CHECK_ATTRS))]
    for (group, attr), status in zip(doctor._CHECK_ATTRS, statuses):
        def make(g, s):
            return lambda ctx: [doctor._row(g, s, f"stub {s}", "do X")]
        monkeypatch.setattr(doctor, attr, make(group, status))
    results = doctor.run_doctor()
    assert [r["status"] for r in results] == statuses
    assert [r["check"] for r in results] == [
        g for g, _ in doctor._CHECK_ATTRS
    ]
    for r in results:
        assert set(r) == {"check", "status", "detail", "fix"}


def test_run_doctor_single_check_flips_status(monkeypatch):
    """Flipping ONE probe changes only that row (checks are isolated)."""
    for _, attr in doctor._CHECK_ATTRS:
        monkeypatch.setattr(
            doctor, attr,
            lambda ctx, a=attr: [doctor._row(a, "PASS", "ok")],
        )
    monkeypatch.setattr(
        doctor, "_check_credentials",
        lambda ctx: [doctor._row("credentials", "FAIL", "no keys", "set")],
    )
    results = doctor.run_doctor()
    fails = [r for r in results if r["status"] == "FAIL"]
    assert len(fails) == 1
    assert fails[0]["check"] == "credentials"


def test_run_doctor_contains_crash_to_one_row(monkeypatch):
    """A probe that raises becomes one FAIL row; the rest still run."""
    def boom(ctx):
        raise RuntimeError("kaput")
    monkeypatch.setattr(doctor, "_check_doc_index", boom)
    for _, attr in doctor._CHECK_ATTRS[1:]:
        monkeypatch.setattr(
            doctor, attr, lambda ctx: [doctor._row("x", "PASS", "ok")],
        )
    results = doctor.run_doctor()
    assert results[0]["status"] == "FAIL"
    assert "kaput" in results[0]["detail"]
    assert all(r["status"] == "PASS" for r in results[1:])


def test_run_doctor_never_raises_when_every_probe_explodes(monkeypatch):
    def boom(ctx):
        raise OSError("disk on fire")
    for _, attr in doctor._CHECK_ATTRS:
        monkeypatch.setattr(doctor, attr, boom)
    results = doctor.run_doctor()          # must not raise
    assert len(results) == len(doctor._CHECK_ATTRS)
    assert all(r["status"] == "FAIL" for r in results)
    assert all("disk on fire" in r["detail"] for r in results)
    text = doctor.format_doctor(results)   # must not raise either
    assert f"0 pass, 0 warn, {len(results)} fail" in text


def test_run_doctor_normalises_bad_check_returns(monkeypatch):
    """Weird return shapes become FAIL rows, never exceptions."""
    monkeypatch.setattr(doctor, "_check_doc_index", lambda ctx: "nonsense")
    monkeypatch.setattr(doctor, "_check_credentials", lambda ctx: [])
    for _, attr in doctor._CHECK_ATTRS[2:]:
        monkeypatch.setattr(
            doctor, attr, lambda ctx: [doctor._row("x", "PASS", "ok")],
        )
    results = doctor.run_doctor()
    assert results[0]["status"] == "FAIL"
    assert results[1]["status"] == "FAIL"


# ---------------------------------------------------------------------------
# Individual probes
# ---------------------------------------------------------------------------


def test_missing_binary_warns_with_fix(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: None)
    rows = doctor._check_binaries(_ctx())
    assert [r["check"] for r in rows] == ["binary: xtb", "binary: orca"]
    for r in rows:
        assert r["status"] == "WARN"
        assert r["fix"]                       # install hint present
    orca = rows[1]
    assert "ORCA" in orca["fix"] or "PATH" in orca["fix"]


def test_present_binary_passes(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: f"/usr/bin/{name}")
    rows = doctor._check_binaries(_ctx())
    assert all(r["status"] == "PASS" for r in rows)
    assert "/usr/bin/xtb" in rows[0]["detail"]


def test_credentials_missing_key_never_leaks_value(monkeypatch):
    secret = "sk-super-secret-value-12345"
    from delfin.agent import credentials as cred_mod

    def fake_load(name, **kw):
        return secret if name == "ANTHROPIC_API_KEY" else ""
    monkeypatch.setattr(cred_mod, "load_credential", fake_load)
    (row,) = doctor._check_credentials(_ctx())
    assert row["status"] == "PASS"
    assert "Anthropic" in row["detail"]
    assert "KIT" in row["detail"] and "OpenAI" in row["detail"]
    assert secret not in row["detail"]
    assert secret not in row["fix"]
    assert secret[:6] not in row["detail"]    # not even a prefix


def test_credentials_none_configured_fails_with_fix(monkeypatch):
    from delfin.agent import credentials as cred_mod
    monkeypatch.setattr(cred_mod, "load_credential", lambda n, **kw: "")
    (row,) = doctor._check_credentials(_ctx())
    assert row["status"] == "FAIL"
    assert "credentials set" in row["fix"]


def test_doc_index_missing_warns(monkeypatch, tmp_path):
    from delfin.doc_server import indexer
    monkeypatch.setattr(
        indexer, "get_default_index_path",
        lambda: tmp_path / "doc_index.json",
    )
    (row,) = doctor._check_doc_index(_ctx())
    assert row["status"] == "WARN"
    assert "delfin-docs-index" in row["fix"]


def test_doc_index_present_reports_count(monkeypatch, tmp_path):
    from delfin.doc_server import indexer
    idx = tmp_path / "doc_index.json"
    idx.write_text(
        '{"documents": {"a": {}, "b": {}}}', encoding="utf-8",
    )
    monkeypatch.setattr(indexer, "get_default_index_path", lambda: idx)
    (row,) = doctor._check_doc_index(_ctx())
    assert row["status"] == "PASS"
    assert "2 document(s)" in row["detail"]


def test_doc_index_corrupt_fails(monkeypatch, tmp_path):
    from delfin.doc_server import indexer
    idx = tmp_path / "doc_index.json"
    idx.write_text("{not json", encoding="utf-8")
    monkeypatch.setattr(indexer, "get_default_index_path", lambda: idx)
    (row,) = doctor._check_doc_index(_ctx())
    assert row["status"] == "FAIL"


def test_scheduler_not_running_with_entries_warns(monkeypatch):
    from delfin.agent import scheduler_daemon as sd
    monkeypatch.setattr(
        sd, "daemon_status",
        lambda *a, **k: {"running": False, "pid": 0,
                         "entries": 3, "disabled": 0},
    )
    (row,) = doctor._check_scheduler(_ctx())
    assert row["status"] == "WARN"
    assert "scheduler start" in row["fix"]


def test_scheduler_idle_without_entries_passes(monkeypatch):
    from delfin.agent import scheduler_daemon as sd
    monkeypatch.setattr(
        sd, "daemon_status",
        lambda *a, **k: {"running": False, "pid": 0,
                         "entries": 0, "disabled": 0},
    )
    (row,) = doctor._check_scheduler(_ctx())
    assert row["status"] == "PASS"


def test_attention_pending_warns(monkeypatch):
    from delfin.agent import attention
    monkeypatch.setattr(
        attention, "list_pending", lambda *a, **k: [{"id": "1"}, {"id": "2"}],
    )
    (row,) = doctor._check_attention(_ctx())
    assert row["status"] == "WARN"
    assert "2 pending" in row["detail"]


def test_benchmark_summary_reuses_optimize_check(monkeypatch):
    from delfin.agent import optimize_check as oc
    monkeypatch.setattr(
        oc, "run_checks",
        lambda **kw: [oc.Issue("warn", "tasks", "todo placeholder")],
    )
    (row,) = doctor._check_benchmark(_ctx())
    assert row["status"] == "WARN"
    assert "1 warning(s)" in row["detail"]
    monkeypatch.setattr(
        oc, "run_checks",
        lambda **kw: [oc.Issue("error", "ground-truth", "broken")],
    )
    (row,) = doctor._check_benchmark(_ctx())
    assert row["status"] == "FAIL"


def test_mcp_fast_lists_configured_without_network(monkeypatch):
    from delfin.agent import mcp_client
    calls = {"registry": 0}

    class NoRegistry:  # instantiating it would mean a network/process probe
        def __init__(self):
            calls["registry"] += 1
    monkeypatch.setattr(mcp_client, "MCPRegistry", NoRegistry)
    monkeypatch.setattr(
        mcp_client, "_load_configs",
        lambda ws: {"delfin-tools": {}, "extra": {}},
    )
    (row,) = doctor._check_mcp(_ctx(fast=True))
    assert row["status"] == "PASS"
    assert "2 configured" in row["detail"]
    assert calls["registry"] == 0


def test_memory_store_unwritable_fails(monkeypatch, tmp_path):
    home = tmp_path / "home"
    delfin_dir = home / ".delfin"
    delfin_dir.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", lambda: home)
    delfin_dir.chmod(0o500)
    try:
        (row,) = doctor._check_memory_store(_ctx())
    finally:
        delfin_dir.chmod(0o700)
    assert row["status"] == "FAIL"
    assert "chmod" in row["fix"]


def test_memory_store_writable_passes(monkeypatch, tmp_path):
    home = tmp_path / "home"
    (home / ".delfin").mkdir(parents=True)
    monkeypatch.setattr(Path, "home", lambda: home)
    (row,) = doctor._check_memory_store(_ctx())
    assert row["status"] == "PASS"
    assert not (home / ".delfin" / ".doctor_probe").exists()


def test_disk_space_below_threshold_warns(monkeypatch, tmp_path):
    import collections
    Usage = collections.namedtuple("usage", "total used free")
    monkeypatch.setattr(
        shutil, "disk_usage",
        lambda p: Usage(10 * 1024**3, 10 * 1024**3, 512 * 1024**2),
    )
    (row,) = doctor._check_disk(_ctx(workspace=str(tmp_path)))
    assert row["status"] == "WARN"
    assert "free up disk space" in row["fix"]


# ---------------------------------------------------------------------------
# format_doctor
# ---------------------------------------------------------------------------


def test_format_doctor_structure():
    results = [
        doctor._row("alpha", "PASS", "all good"),
        doctor._row("beta-longer-name", "WARN", "meh", "try this"),
        doctor._row("gamma", "FAIL", "broken", "fix me"),
    ]
    text = doctor.format_doctor(results)
    lines = text.splitlines()
    assert lines[0].startswith("✅ alpha")
    assert any(line.startswith("⚠️ beta-longer-name") for line in lines)
    assert any(line.startswith("❌ gamma") for line in lines)
    assert sum("fix:" in line for line in lines) == 2   # not for PASS
    assert "fix: try this" in text and "fix: fix me" in text
    assert lines[-1] == "1 pass, 1 warn, 1 fail"


def test_format_doctor_empty():
    text = doctor.format_doctor([])
    assert "0 pass, 0 warn, 0 fail" in text


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_cli_parser_accepts_doctor_and_workspace():
    p = agent_cli.build_parser()
    args = p.parse_args(["doctor"])
    assert args.cmd == "doctor"
    assert args.workspace == ""
    args = p.parse_args(["doctor", "--workspace", "/tmp/ws"])
    assert args.workspace == "/tmp/ws"


def test_cli_doctor_exit_0_without_fail(monkeypatch, capsys):
    monkeypatch.setattr(
        doctor, "run_doctor",
        lambda ws=None, **kw: [
            doctor._row("a", "PASS", "ok"),
            doctor._row("b", "WARN", "meh", "hint"),
        ],
    )
    rc = agent_cli.main(["doctor"])
    assert rc == 0
    out, _ = capsys.readouterr()
    assert "1 pass, 1 warn, 0 fail" in out


def test_cli_doctor_exit_1_on_fail(monkeypatch, capsys):
    monkeypatch.setattr(
        doctor, "run_doctor",
        lambda ws=None, **kw: [doctor._row("a", "FAIL", "broken", "fix")],
    )
    rc = agent_cli.main(["doctor"])
    assert rc == 1
    out, _ = capsys.readouterr()
    assert "0 pass, 0 warn, 1 fail" in out


def test_cli_doctor_passes_workspace_through(monkeypatch):
    seen = {}

    def fake(ws=None, **kw):
        seen["ws"] = ws
        return [doctor._row("a", "PASS", "ok")]
    monkeypatch.setattr(doctor, "run_doctor", fake)
    rc = agent_cli.main(["doctor", "--workspace", "/tmp/ws"])
    assert rc == 0
    assert seen["ws"] == "/tmp/ws"


def test_cli_doctor_survives_exploding_probes(monkeypatch, capsys):
    def boom(ctx):
        raise RuntimeError("nope")
    for _, attr in doctor._CHECK_ATTRS:
        monkeypatch.setattr(doctor, attr, boom)
    rc = agent_cli.main(["doctor"])       # must not raise
    assert rc == 1
    out, _ = capsys.readouterr()
    assert f"{len(doctor._CHECK_ATTRS)} fail" in out


# ---------------------------------------------------------------------------
# Filesystem isolation for shell commands must be visible, not assumed
# ---------------------------------------------------------------------------


def _isolation_row(monkeypatch, mode, usable):
    import delfin.agent.doctor as doc
    import delfin.user_settings as us
    monkeypatch.setattr(
        us, "load_settings", lambda: {"agent": {"bash_isolation": mode}})
    import delfin.agent.api_client as ac
    monkeypatch.setattr(ac, "_bwrap_functional", lambda: usable)
    return doc._check_bash_isolation({})[0]


def test_isolation_active_is_a_pass(monkeypatch):
    row = _isolation_row(monkeypatch, "bwrap", True)
    assert row["status"] == "PASS" and "active" in row["detail"]


def test_isolation_configured_but_unusable_is_a_failure(monkeypatch):
    row = _isolation_row(monkeypatch, "bwrap", False)
    assert row["status"] == "FAIL" and "does not work here" in row["detail"]


def test_auto_mode_says_when_it_does_not_isolate(monkeypatch):
    row = _isolation_row(monkeypatch, "auto", True)
    assert row["status"] == "WARN"
    assert "bypassPermissions only" in row["detail"]
    assert "bash_isolation" in row["fix"]


def test_auto_without_working_bwrap_says_never_isolated(monkeypatch):
    row = _isolation_row(monkeypatch, "auto", False)
    assert "never isolated" in row["detail"]


def test_isolation_off_is_reported_with_the_remaining_protection(monkeypatch):
    row = _isolation_row(monkeypatch, "off", True)
    assert row["status"] == "WARN" and "write-target gate" in row["detail"]


def test_isolation_check_is_registered_and_never_raises(monkeypatch):
    import delfin.agent.doctor as doc
    assert ("bash isolation", "_check_bash_isolation") in doc._CHECK_ATTRS
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings",
                        lambda: (_ for _ in ()).throw(RuntimeError("boom")))
    assert doc._check_bash_isolation({})[0]["check"] == "bash isolation"
