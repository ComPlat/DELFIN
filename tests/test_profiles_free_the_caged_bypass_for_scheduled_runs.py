"""The caged bypass rule reaches scheduled runs and the per-workspace
profiles that grant it.

Swarm side (merged #493): in Bypass, workspace work that runs inside the
bash cage runs free; push, network, installers and anything uncaged still
asks. Here: a scheduled run may run in bypass only when the host can cage
AND the schedule's workspace has the profile "free"; per-workspace
profiles are the user's, masked and confirmed like schedules; the doctor
says whether this host can cage.
"""

import json

import pytest


def _settings(bypass=True, unattended=True):
    if not bypass:
        return {}
    return {"agent": {"scheduler": {"permission_mode": "bypassPermissions",
                                    "unattended": unattended}}}


@pytest.fixture()
def cage(monkeypatch):
    """Patch the host-cage decision the runner reads: 'bwrap' by default."""
    from delfin.agent import workspace_profiles
    monkeypatch.setattr(workspace_profiles, "host_cage_decision",
                        lambda: "bwrap")
    return workspace_profiles


@pytest.fixture()
def profiles(tmp_path, monkeypatch, cage):
    """A redirected, masked profile store."""
    from delfin.agent import trust_store, workspace_profiles
    path = tmp_path / "user_schedules"
    monkeypatch.setattr(workspace_profiles, "_DIR", path)
    monkeypatch.setattr(trust_store, "store_dir",
                        lambda: tmp_path / ".delfin" / "trust_store")
    trust_store.forget_cached_keys()
    return path


def test_a_scheduled_run_with_the_free_profile_runs_in_bypass_on_a_caged_host(
        tmp_path, profiles, cage):
    from delfin.agent.service import schedule_runner
    ws = tmp_path / "proj"
    ws.mkdir()
    cage.set_profile(str(ws), "free")
    entry = {"id": 4, "task": "x", "when": "alle 1h", "workspace": str(ws)}
    params = schedule_runner.session_params(entry, _settings())
    assert params["permissionMode"] == "bypassPermissions"
    assert params["unattended"] is True and params["attach"] is False


def test_a_scheduled_run_stays_capped_without_the_free_profile(
        tmp_path, profiles, cage):
    from delfin.agent.service import schedule_runner
    ws = tmp_path / "proj"
    ws.mkdir()
    entry = {"id": 4, "task": "x", "when": "alle 1h", "workspace": str(ws)}
    assert schedule_runner.session_params(entry, _settings())[
        "permissionMode"] == "acceptEdits"


def test_a_scheduled_run_stays_capped_on_an_uncaged_host(
        tmp_path, profiles, monkeypatch):
    from delfin.agent import workspace_profiles
    from delfin.agent.service import schedule_runner
    monkeypatch.setattr(workspace_profiles, "host_cage_decision",
                        lambda: "none")
    ws = tmp_path / "proj"
    ws.mkdir()
    with pytest.raises(Exception):
        # No cage: the profile "free" is refused with the reason.
        workspace_profiles.set_profile(str(ws), "free")
    assert workspace_profiles.get_profile(str(ws)) == "default"
    entry = {"id": 4, "task": "x", "when": "alle 1h", "workspace": str(ws)}
    assert schedule_runner.session_params(entry, _settings())[
        "permissionMode"] == "acceptEdits"


def test_a_scheduled_run_stays_capped_without_the_bypass_opt_in(
        tmp_path, profiles, cage):
    from delfin.agent.service import schedule_runner
    ws = tmp_path / "proj"
    ws.mkdir()
    cage.set_profile(str(ws), "free")
    entry = {"id": 4, "task": "x", "when": "alle 1h", "workspace": str(ws)}
    assert schedule_runner.session_params(entry, _settings(unattended=False))[
        "permissionMode"] == "acceptEdits"


def test_the_profile_store_masks_and_confirms(tmp_path, profiles, cage):
    from delfin.agent import workspace_profiles
    ws = tmp_path / "proj"
    ws.mkdir()
    result = cage.set_profile(str(ws), "free")
    assert result["profile"] == "free"
    assert result["confirmed"] is True
    # Masked: 0700 directory, 0600 file -- the agent's command cage and the
    # deny globs cover the folder name (shared with the schedules store).
    mode = profiles.stat().st_mode & 0o777
    file_mode = (profiles / "workspace_profiles.json").stat().st_mode & 0o777
    assert mode == 0o700 and file_mode == 0o600
    # Confirmed exactly: the HMAC under the account key.
    data = json.loads((profiles / "workspace_profiles.json").read_text())
    assert data["profiles"][str(ws)]["confirmed_sha"].startswith("sha")
    assert workspace_profiles.confirmed(str(ws)) is True
    # A later content change without a fresh confirmation is unconfirmed.
    data["profiles"][str(ws)]["profile"] = "free"
    data["profiles"][str(ws)]["confirmed_sha"] = "sha:0"
    (profiles / "workspace_profiles.json").write_text(json.dumps(data))
    assert workspace_profiles.confirmed(str(ws)) is False
    assert workspace_profiles.get_profile(str(ws)) == "default"


def test_the_profile_store_refuses_an_unknown_profile(tmp_path, profiles, cage):
    ws = tmp_path / "proj"
    ws.mkdir()
    with pytest.raises(ValueError):
        cage.set_profile(str(ws), "bypass-everything")
    assert cage.get_profile(str(ws)) == "default"


def test_set_from_the_messenger_words(tmp_path, profiles, cage):
    from delfin.agent import workspace_profiles
    ws = tmp_path / "proj"
    ws.mkdir()
    assert workspace_profiles.set_from_word(str(ws), "frei")["profile"] == "free"
    assert workspace_profiles.get_profile(str(ws)) == "free"
    assert workspace_profiles.set_from_word(str(ws), "standard")[
        "profile"] == "default"
    assert workspace_profiles.get_profile(str(ws)) == "default"


def test_the_bridge_sets_the_profile_for_a_session_workspace(
        tmp_path, profiles, cage, monkeypatch):
    from delfin.agent.service import bridge
    ws = tmp_path / "proj"
    ws.mkdir()
    b = bridge.Bridge.__new__(bridge.Bridge)   # no sockets
    async def fake_request(method, params=None):
        assert method == "sessions.list"
        return {"sessions": [{"sessionId": "s1", "workspace": str(ws)}]}
    monkeypatch.setattr(b, "_request", fake_request)
    import asyncio
    reply = asyncio.run(b._cmd_profile({"session": "s1", "text": "frei"}))
    assert reply.get("ok") is True
    assert cage.get_profile(str(ws)) == "free"
    reply = asyncio.run(b._cmd_profile({"session": "s1", "text": "blubb"}))
    assert reply.get("ok") is not True


def test_doctor_says_whether_this_host_can_cage(monkeypatch):
    from delfin.agent import doctor, workspace_profiles
    monkeypatch.setattr(workspace_profiles, "host_cage_decision",
                        lambda: "bwrap")
    rows = doctor._check_bash_cage({})
    assert rows and rows[0]["check"] == "bash cage"
    assert rows[0]["status"] == "PASS"
    monkeypatch.setattr(workspace_profiles, "host_cage_decision",
                        lambda: "none")
    rows = doctor._check_bash_cage({})
    assert rows[0]["status"] == "WARN"
    assert "free" in rows[0]["detail"]


def test_caged_workspace_work_runs_free_but_push_and_network_ask(monkeypatch):
    from delfin.agent import api_client

    class Perms:
        scope_locked = False

    monkeypatch.setattr(api_client, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    assert api_client._bash_runs_caged(Perms()) is True

    class Action:
        tool = "bash"

        def __init__(self, command):
            self.action = {"tool": "bash", "command": command}

    write = Action("python -c \"open('out.txt','w').write('x')\"")
    assert api_client._leaves_workspace(Perms(), write.action, True) is False
    push = Action("git push origin main")
    assert api_client._leaves_workspace(Perms(), push.action, True) is True
    net = Action("curl https://example.org/data.tar.gz -o data.tar.gz")
    assert api_client._leaves_workspace(Perms(), net.action, True) is True
    install = Action("pip install numpy")
    assert api_client._leaves_workspace(Perms(), install.action, True) is True
