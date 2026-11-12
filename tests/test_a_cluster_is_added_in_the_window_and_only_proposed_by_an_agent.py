"""A cluster is added and removed in the window; an agent can only propose one.

Input: ``delfin.remote.setup`` with DELFIN's settings and ctun's config in
the test's directory (ctun's real ``add_cluster``/``remove_cluster`` write
it; its host-key calls, which would reach a server, are replaced), the
attention inbox (the suite redirects it), the agent's ``cluster_propose``,
and ``scripts/orca.sh`` run with bash in a fake home. Output: the two
config files, the results the window gets, the inbox records.

Semantics:

* an add writes ctun's entry (a preset's login prompts and budget guard,
  ``guard_scope: every-command``; none for a cluster given by host) and
  DELFIN's (partition, directories); a preset's published key is pinned
  when the server presents exactly it, any other key is only shown, and
  pinned by ``pin`` with the fingerprint the person confirmed;
* a bad field, a name either config has, or settings that cannot be
  written leave neither entry behind;
* a remove takes both out, and waits while DELFIN still follows a job there;
* only a window may add, pin or remove (the service's ACTING gate);
  ``cluster_propose`` records a NEEDS YOU item with the form's fields and
  writes nothing else; a text answer to it is refused; the add resolves it;
* ORCA is looked for where the installer looks (an ORCA directory, else a
  tarball), never taken for /usr/bin/orca, and its absence says what to do.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from cluster_tunnel import api as ctun_api
from delfin.agent import attention
from delfin.agent.service import attention_feed, protocol, supervisor
from delfin.remote import setup

PUBLISHED = "SHA256:5mZYEpKigwK5ibBMHRrh3WIkOtCqomJW6H7OMbPk3ec"
OTHER = "SHA256:not-the-published-key-0000000000000000000000"


@pytest.fixture
def where(tmp_path, monkeypatch):
    """DELFIN's settings and ctun's config in tmp_path; the host-key calls
    answered as a server presenting ``where.offered`` would."""
    w = SimpleNamespace(settings=tmp_path / "delfin_settings.json",
                        ctun=tmp_path / "ctun" / "config.yaml", offered=PUBLISHED, pinned=[])
    w.ctun.parent.mkdir()

    def pin_host_key(name, fingerprint, *, replace=False, config_path=None):
        if fingerprint != w.offered:
            raise ctun_api.HostKeyMismatch(name, w.offered, fingerprint)
        w.pinned.append((name, fingerprint))

    def host_key(name, *, config=None, config_path=None):
        return SimpleNamespace(offered=SimpleNamespace(fingerprint=w.offered, kind="ssh-ed25519"),
                               strict=False, matches=False)

    monkeypatch.setattr(ctun_api, "pin_host_key", pin_host_key)
    monkeypatch.setattr(ctun_api, "host_key", host_key)
    monkeypatch.setattr(ctun_api, "logout",
                        lambda *a, **k: [SimpleNamespace(was_live=False)])
    real_add, real_remove = setup.add, setup.remove
    w.add = lambda **p: real_add(p, settings_path=str(w.settings), ctun_config=str(w.ctun))
    w.remove = lambda name: real_remove(name, settings_path=str(w.settings),
                                        ctun_config=str(w.ctun))
    w.delfin = lambda: (json.loads(w.settings.read_text()).get("remote") or {}).get("clusters") or {} \
        if w.settings.exists() else {}
    w.entries = lambda: ((yaml.safe_load(w.ctun.read_text()) or {}).get("clusters") or {}) \
        if w.ctun.exists() else {}
    return w


def test_a_preset_whose_server_presents_the_published_key_is_pinned_at_once(where):
    got = where.add(name="bwuni3", preset="bwuni3", user="ka_ab1234", partition="cpu")
    assert got["ok"] and got["hostKey"]["pinned"] is True
    assert got["hostKey"]["fingerprint"] == PUBLISHED
    assert "bwhpc" in (got["hostKey"]["source"] or "")
    assert where.pinned == [("bwuni3", PUBLISHED)]
    entry = where.entries()["bwuni3"]
    assert (entry["host"], entry["user"], entry["requires_otp"]) == ("uc3.scc.kit.edu", "ka_ab1234", True)
    assert entry["budget"]["guard_scope"] == "every-command"
    assert where.delfin() == {"bwuni3": {"default_partition": "cpu"}}


def test_another_key_is_only_shown_and_pinned_once_the_person_confirms_it(where):
    where.offered = OTHER
    got = where.add(name="bwuni3", preset="bwuni3", user="ka_ab1234")
    key = got["hostKey"]
    assert got["ok"] and key["pinned"] is False
    assert (key["offered"], key["expected"]) == (OTHER, PUBLISHED)
    assert where.pinned == []
    assert setup.pin("bwuni3", PUBLISHED, ctun_config=str(where.ctun))["code"] == "host_key_mismatch"
    assert setup.pin("bwuni3", OTHER, ctun_config=str(where.ctun))["ok"] is True
    assert where.pinned == [("bwuni3", OTHER)]


def test_a_cluster_given_by_its_host_has_no_budget_guard_and_shows_its_key(where):
    where.offered = OTHER
    got = where.add(name="lab", host="login.lab.example", user="me",
                    workDir="~/runs", installRoot="/scratch/me/delfin")
    assert got["ok"] and got["hostKey"]["offered"] == OTHER and not got["hostKey"]["pinned"]
    assert "budget" not in where.entries()["lab"]
    assert where.delfin() == {"lab": {"work_dir": "~/runs", "install_root": "/scratch/me/delfin"}}


@pytest.mark.parametrize("fields, word", [
    ({"partition": "bad name"}, "partition"),
    ({"workDir": "~"}, "work_dir"),
    ({"installRoot": "~/a/../b"}, "install_root"),
    ({"user": ""}, "user"),
    ({"preset": None, "host": None}, "host"),
    ({"preset": "nowhere"}, "preset"),
])
def test_a_bad_field_writes_neither_entry(where, fields, word):
    params = {"name": "bwuni3", "preset": "bwuni3", "user": "ka_ab1234", **fields}
    got = where.add(**params)
    assert got["ok"] is False and got["code"] == "invalid" and word in got["error"]
    assert where.entries() == {} and where.delfin() == {}


def test_a_name_either_config_has_is_refused(where):
    assert where.add(name="kcist", preset="kcist", user="me")["ok"]
    again = where.add(name="kcist", preset="kcist", user="me")
    assert (again["ok"], again["code"]) == (False, "exists")
    ctun_api.add_cluster("horeka", preset="horeka", user="me", config_path=str(where.ctun))
    got = where.add(name="horeka", preset="horeka", user="me")
    assert (got["ok"], got["code"]) == (False, "exists") and "ctun remove" in got["error"]
    assert "horeka" not in where.delfin()


def test_settings_that_cannot_be_written_take_ctuns_entry_out_again(where, monkeypatch):
    import delfin.user_settings as user_settings

    def broken(*_a, **_k):
        raise OSError("read-only file system")
    monkeypatch.setattr(user_settings, "save_settings", broken)
    got = where.add(name="kcist", preset="kcist", user="me")
    assert (got["ok"], got["code"]) == (False, "settings_unwritable")
    assert "kcist" not in where.entries()


def test_a_remove_takes_both_entries_out(where):
    assert where.add(name="kcist", preset="kcist", user="me", partition="small")["ok"]
    assert where.add(name="lab", host="login.lab.example", user="me")["ok"]
    got = where.remove("kcist")
    assert got["ok"] is True
    assert set(where.entries()) == {"lab"} and set(where.delfin()) == {"lab"}
    assert where.remove("kcist")["code"] == "unknown"


def test_a_remove_waits_while_delfin_follows_a_job_there(where, monkeypatch):
    assert where.add(name="kcist", preset="kcist", user="me")["ok"]
    monkeypatch.setattr(setup, "_followed", lambda name: 2 if name == "kcist" else 0)
    got = where.remove("kcist")
    assert (got["ok"], got["code"]) == (False, "jobs_followed") and "2 jobs" in got["error"]
    assert "kcist" in where.entries() and "kcist" in where.delfin()


def test_only_a_window_adds_pins_or_removes_a_cluster():
    writes = {"cluster.add", "cluster.hostkey", "cluster.pin", "cluster.remove"}
    assert writes <= supervisor.ACTING
    assert writes | {"clusters.presets"} <= set(protocol.SERVICE_REQUESTS)
    assert supervisor.SERVICE_VERSION >= 7


def test_the_presets_say_their_published_keys_and_where_from():
    got = {p["name"]: p for p in setup.presets()["presets"]}
    assert got["bwuni3"]["fingerprint"] == PUBLISHED and got["bwuni3"]["partition"] == "cpu"
    assert got["kcist"]["fingerprint"] is None
    assert got["bwuni3"]["userHint"] and got["bwuni3"]["budget"]["unit"]


def _proposals(name: str) -> list[dict]:
    return [e for e in attention.list_pending("cluster_add") if e.get("cluster") == name]


def test_an_agent_proposes_a_cluster_as_an_item_and_writes_nothing_else(where):
    from delfin.agent import cluster_tools
    perms = SimpleNamespace(task_session_id="s-1", workspace=None, mode="default")
    got = json.loads(cluster_tools.propose(None, perms, name="bwuni3", preset="bwuni3",
                                           user="ka_ab1234", partition="cpu"))
    assert got["ok"] and "nothing is added" in got["next"]
    (event,) = _proposals("bwuni3")
    assert event["proposal"] == {"name": "bwuni3", "preset": "bwuni3", "user": "ka_ab1234",
                                 "host": "", "partition": "cpu"}
    assert where.entries() == {} and where.delfin() == {}
    item = attention_feed.item_of(event, findings={}, presence={}, titles=attention_feed.Titles())
    assert (item["kind"], item["group"], item["cluster"]) == ("addCluster", "needs", "bwuni3")
    assert item["proposal"]["user"] == "ka_ab1234"
    refused = attention.answer_item(event["id"], "yes add it")
    assert refused["ok"] is False and "form" in refused["error"]
    # Proposed again: one item, with the newer values.
    cluster_tools.propose(None, perms, name="bwuni3", preset="bwuni3", user="ka_cd5678")
    (event,) = _proposals("bwuni3")
    assert event["proposal"]["user"] == "ka_cd5678"
    # The person's add closes it, with nothing for a model.
    assert where.add(name="bwuni3", preset="bwuni3", user="ka_cd5678")["ok"]
    assert _proposals("bwuni3") == []


def test_a_proposal_names_a_known_preset_or_a_host():
    from delfin.agent import cluster_tools
    perms = SimpleNamespace(task_session_id="", workspace=None, mode="default")
    assert '"not_proposed"' in cluster_tools.propose(None, perms, name="x", preset="nowhere")
    assert '"not_proposed"' in cluster_tools.propose(None, perms, name="x")
    assert '"not_proposed"' in cluster_tools.propose(None, perms, preset="kcist")
    assert _proposals("x") == []


def _orca(home: Path, *extra: str) -> str:
    script = Path(setup.ORCA_SCRIPT)
    env = {"HOME": str(home), "PATH": f"{home / 'bin'}:/usr/bin:/bin"}
    return subprocess.run(["bash", "-s", "--", *extra], input=script.read_text(), text=True,
                          capture_output=True, env=env, check=True).stdout.strip()


def _orca_dir(path: Path) -> Path:
    path.mkdir(parents=True)
    for exe in ("orca", "orca_plot"):
        (path / exe).write_text("#!/bin/sh\n")
        (path / exe).chmod(0o755)
    return path


def test_orca_is_looked_for_where_the_installer_looks(tmp_path):
    home = tmp_path / "home"
    (home / "bin").mkdir(parents=True)
    assert _orca(home) == "none"
    screen_reader = home / "bin" / "orca"            # an `orca` that is no ORCA directory
    screen_reader.write_text("#!/bin/sh\n")
    screen_reader.chmod(0o755)
    assert _orca(home) == "none"
    (home / "orca_6_1_1_linux_x86-64_shared_openmpi418.tar.xz").write_text("x")
    assert _orca(home) == f"tarball {home}/orca_6_1_1_linux_x86-64_shared_openmpi418.tar.xz"
    shared = _orca_dir(home / "delfin" / "software" / "orca_6_0_1")
    assert _orca(home, "~/delfin/software") == f"dir {shared}"
    _orca_dir(home / "orca_6_0_1")
    newer = _orca_dir(home / "orca_6_1_1")                 # the newest, as the installer sorts
    assert _orca(home, "~/delfin/software") == f"dir {newer}"


def test_a_missing_orca_says_what_to_do_and_never_copies():
    class Conn:
        cluster = SimpleNamespace(software_dir="~/delfin/software")

        def __init__(self, out):
            self.out = out

        def feed_script(self, script, args, timeout=None):
            assert Path(script).name == "orca.sh" and args == ["~/delfin/software"]
            return SimpleNamespace(ok=True, stdout=self.out)

    missing = setup.orca_status(Conn("none\n"))
    assert missing["found"] is False and "home" in missing["hint"] and "never copies" in missing["hint"]
    found = setup.orca_status(Conn("dir /home/me/orca_6_1_1\n"))
    assert found == {"found": True, "path": "/home/me/orca_6_1_1", "how": "dir", "hint": ""}
    packed = setup.orca_status(Conn("tarball /home/me/orca.tar.xz\n"))
    assert packed["found"] and "unpacks" in packed["hint"]
