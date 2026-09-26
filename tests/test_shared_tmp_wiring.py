"""Phase 4 (LH): the WIRING, as xfail(strict) over the real handlers.

The cage and api_client are NOT my files in this assignment -- these
tests describe the wiring the OPERATOR will do, driven through the
real public call paths:

1. ``_bash_isolation_argv`` (api_client) soll statt ``--tmpfs /tmp``
   den Sitzungs-Ordner als /tmp binden (``shared_tmp.bwrap_substitution``).
2. ``_execute_read_file`` (api_client) soll ``/tmp/...`` über
   ``shared_tmp.map_to_host`` auf den Sitzungs-Ordner abbilden, BEVOR
   der Read-Gate läuft.
3. Ein Symlink IM Sitzungs-Ordner auf eine geheime Datei (~/.ssh/...)
   muss trotzdem abgelehnt werden: die Abbildung ersetzt die
   Pfad-Prüfungen nicht, sie läuft davor.

xfail(strict=True): while the wiring is missing these tests MUST fail
(XPASS would be a bug). When the operator wires it up they turn green,
and strict guarantees no wiring slips past silently.
"""

import json
from pathlib import Path

import pytest

from delfin.agent import shared_tmp


@pytest.fixture
def env(tmp_path, monkeypatch):
    ws = tmp_path / "ws"
    ws.mkdir()
    monkeypatch.setenv("DELFIN_TMP_ROOT", str(tmp_path / "rt"))
    monkeypatch.setattr("delfin.agent.api_client._process_cage_enabled",
                        lambda: False, raising=False)
    monkeypatch.setattr("delfin.agent.api_client._process_cage_options",
                        lambda: [], raising=False)
    from delfin.agent.api_client import KitToolPermissions
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "acceptEdits"
    perms.task_session_id = "wire-test"
    sd = shared_tmp.session_dir("wire-test")
    shared_tmp.ensure_session_dir(sd)
    return perms, sd


@pytest.mark.xfail(strict=True, reason=(
    "operator wiring pending: _bash_isolation_argv still uses "
    "--tmpfs /tmp; the swap is bwrap_substitution(session_dir)"))
def test_isolation_argv_binds_the_session_dir_as_tmp(env, monkeypatch):
    from delfin.agent import api_client
    perms, sd = env
    monkeypatch.setattr("delfin.agent.api_client._git_metadata_roots",
                        lambda ws: [], raising=False)
    argv = api_client._bash_isolation_argv("true", perms.workspace, perms,
                                           mode="bwrap")
    old, new = shared_tmp.bwrap_substitution(sd)
    assert new == ["--bind", str(sd), "/tmp"]
    # the bind must be CONSECUTIVE in argv (list `<=` is an element
    # subset test and passed vacuously against --tmpfs /tmp)
    triples = [tuple(argv[i:i + 3]) for i in range(len(argv) - 2)]
    pairs = [tuple(argv[i:i + 2]) for i in range(len(argv) - 1)]
    assert tuple(new) in triples
    assert tuple(old) not in pairs  # --tmpfs /tmp must be GONE


@pytest.mark.xfail(strict=True, reason=(
    "operator wiring pending: read_file does not map /tmp/... through "
    "map_to_host before the read gate"))
def test_read_file_maps_cage_tmp_through_the_real_handler(env):
    from delfin.agent.api_client import _DocToolExecutor
    perms, sd = env
    (sd / "probe.txt").write_text("from the cage\n", encoding="utf-8")
    ex = _DocToolExecutor()
    out = ex._execute_read_file({"path": "/tmp/probe.txt"}, perms)
    assert "from the cage" in out
    assert not out.startswith('{"error"')


@pytest.mark.xfail(strict=True, reason=(
    "operator wiring pending: the mapping must run BEFORE the gate, so "
    "a symlink inside the session dir to ~/.ssh/... stays denied"))
def test_symlink_in_session_dir_to_a_secret_stays_denied(env, tmp_path,
                                                         monkeypatch):
    from delfin.agent.api_client import _DocToolExecutor
    perms, sd = env
    # a LEGITIMATE file first: the mapping must make it readable, so
    # the denial below cannot pass vacuously via "file not found"
    (sd / "plain.txt").write_text("ordinary scratch", encoding="utf-8")
    fake_home = tmp_path / "home"
    (fake_home / ".ssh").mkdir(parents=True)
    secret = fake_home / ".ssh" / "id_rsa"
    secret.write_text("SECRET KEY MATERIAL", encoding="utf-8")
    monkeypatch.setenv("HOME", str(fake_home))
    (sd / "leak").symlink_to(secret)
    ex = _DocToolExecutor()
    ok = ex._execute_read_file({"path": "/tmp/plain.txt"}, perms)
    assert "ordinary scratch" in ok  # mapping works for real files
    out = ex._execute_read_file({"path": "/tmp/leak"}, perms)
    parsed = json.loads(out)
    assert "error" in parsed
    assert "SECRET" not in out
