"""Settings in two layers: shared defaults under a personal file.

DELFIN_SETTINGS_DEFAULTS names a defaults file (the launcher's people mode:
the account's ~/.delfin_settings.json). Effective settings are
DEFAULT_SETTINGS < defaults file < personal file; save_settings writes to the
personal file only what differs from DEFAULT_SETTINGS < defaults file.
Without the variable the single-file behaviour is unchanged, byte for byte.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from delfin import user_settings as us

ROOT = Path(__file__).resolve().parents[1]
ENV = us.SETTINGS_DEFAULTS_ENV


def _write(path, payload):
    path.write_text(json.dumps(payload), encoding="utf-8")


@pytest.fixture
def layered(tmp_path, monkeypatch):
    """An account file (the defaults) and two persons' files under it."""
    monkeypatch.setenv("HOME", str(tmp_path / "account"))
    account = tmp_path / "account" / us.SETTINGS_FILE_NAME
    account.parent.mkdir()
    anna = tmp_path / "anna" / us.SETTINGS_FILE_NAME
    bert = tmp_path / "bert" / us.SETTINGS_FILE_NAME
    anna.parent.mkdir()
    bert.parent.mkdir()
    monkeypatch.setenv(ENV, str(account))
    return account, anna, bert


def test_the_module_under_test_is_this_checkout():
    assert Path(us.__file__).resolve().parents[1] == ROOT


# ---------------------------------------------------------------------------
# Without the variable nothing changes
# ---------------------------------------------------------------------------

def _single_file_bytes(payload):
    """What save_settings wrote before the layer existed: the normalised
    settings with every default filled in, indent 2, sorted keys."""
    full = us._merge_missing_defaults(us._normalized_settings_dict(payload),
                                      us.DEFAULT_SETTINGS)
    return json.dumps(full, indent=2, sort_keys=True).encode("utf-8")


_CUSTOM = {
    "agent": {"model": "kit.example", "effort": "high",
              "permissions": {"deny": ["Bash(rm *)"]}},
    "ui": {"tabs": {"order": ["agent", "home"], "hidden": ["submit"]}},
    "runtime": {"backend": "slurm"},
    "transfer": {"host": "example.org", "user": "u", "remote_path": "/r",
                 "port": 22},
}


@pytest.mark.parametrize("variable", ["unset", "empty", "self", "link"])
def test_without_a_layer_the_file_round_trips_byte_for_byte(tmp_path, monkeypatch, variable):
    path = tmp_path / us.SETTINGS_FILE_NAME
    us.save_settings(_CUSTOM, path)
    written = path.read_bytes()
    assert written == _single_file_bytes(_CUSTOM)
    if variable == "empty":
        monkeypatch.setenv(ENV, "  ")
    elif variable == "self":
        monkeypatch.setenv(ENV, str(path))
    elif variable == "link":
        # A person file that still links to the account's is one file: no
        # layer, read through the link. (A write replaces the link by a
        # file, as _write_json_atomic always has; the bootstrap removes such
        # links before DELFIN starts.)
        us.load_settings(path)
        target = tmp_path / "account.json"
        target.write_bytes(path.read_bytes())
        path.unlink()
        path.symlink_to(target)
        monkeypatch.setenv(ENV, str(target))
        assert us.get_settings_defaults_path(path) is None
        loaded = us.load_settings(path)
        assert path.is_symlink()
        assert target.read_bytes() == _single_file_bytes(loaded)
        return
    assert us.get_settings_defaults_path(path) is None
    loaded = us.load_settings(path)
    # load completes the file with the one-time update records, as before:
    # the whole effective dict, in the same format.
    after_load = path.read_bytes()
    assert after_load == json.dumps(loaded, indent=2, sort_keys=True).encode()
    assert after_load == _single_file_bytes(loaded)
    us.save_settings(loaded, path)
    assert path.read_bytes() == after_load
    us.load_settings(path)
    assert path.read_bytes() == after_load
    assert us.load_declared_settings(path) == json.loads(after_load)


def test_without_a_layer_a_partial_file_is_completed_as_before(tmp_path):
    path = tmp_path / us.SETTINGS_FILE_NAME
    _write(path, {"agent": {"model": "kit.example"}})
    loaded = us.load_settings(path)
    on_disk = json.loads(path.read_text())
    assert on_disk == loaded
    assert on_disk["runtime"]["backend"] == us.DEFAULT_SETTINGS["runtime"]["backend"]
    assert us.settings_layers(path) == [path]


# ---------------------------------------------------------------------------
# The layer
# ---------------------------------------------------------------------------

def test_layering_order_is_shipped_then_defaults_then_personal(layered):
    account, anna, _bert = layered
    _write(account, {"agent": {"model": "kit.default", "effort": "low"},
                     "runtime": {"backend": "slurm"}})
    _write(anna, {"agent": {"effort": "high"}})
    eff = us.load_settings(anna)
    assert eff["agent"]["effort"] == "high"                  # personal
    assert eff["agent"]["model"] == "kit.default"            # defaults file
    assert eff["runtime"]["backend"] == "slurm"              # defaults file
    shipped = us.DEFAULT_SETTINGS["runtime"]["local"]["max_cores"]
    assert eff["runtime"]["local"]["max_cores"] == shipped   # DEFAULT_SETTINGS
    assert us.settings_layers(anna) == [anna, account]
    declared = us.load_declared_settings(anna)
    assert declared["agent"] == {"model": "kit.default", "effort": "high"}


def test_save_writes_only_the_difference(layered):
    account, anna, _bert = layered
    _write(account, {"agent": {"model": "kit.default"}})
    _write(anna, {})
    eff = us.load_settings(anna)
    assert json.loads(anna.read_text()) == {}                # load wrote nothing new
    eff["agent"]["effort"] = "high"
    eff["ui"]["tabs"] = {"order": ["agent"], "hidden": ["home"]}
    returned = us.save_settings(eff, anna)
    assert json.loads(anna.read_text()) == {
        "agent": {"effort": "high"},
        "ui": {"tabs": {"order": ["agent"], "hidden": ["home"]}},
    }
    assert returned["agent"]["model"] == "kit.default"
    assert oct(anna.stat().st_mode & 0o777) == "0o600"
    assert json.loads(account.read_text()) == {"agent": {"model": "kit.default"}}


def test_a_partial_save_keeps_the_default_not_the_shipped_value(layered):
    account, anna, _bert = layered
    _write(account, {"runtime": {"backend": "slurm"}})
    us.save_settings({"agent": {"effort": "high"}}, anna)
    assert json.loads(anna.read_text()) == {"agent": {"effort": "high"}}
    assert us.load_settings(anna)["runtime"]["backend"] == "slurm"


def test_a_changed_default_reaches_who_did_not_override_it(layered):
    account, anna, bert = layered
    _write(account, {"agent": {"model": "kit.old"}})
    for person in (anna, bert):
        _write(person, {})
    eff = us.load_settings(anna)
    eff["agent"]["model"] = "kit.mine"
    us.save_settings(eff, anna)
    us.save_settings(us.load_settings(bert), bert)   # bert saves something else
    _write(account, {"agent": {"model": "kit.new"}})
    assert us.load_settings(anna)["agent"]["model"] == "kit.mine"
    assert us.load_settings(bert)["agent"]["model"] == "kit.new"
    assert json.loads(bert.read_text()) == {}


def test_an_old_full_copy_shrinks_to_its_own_choices(layered):
    # A person file written whole (before the layer) holds every default as
    # it was then; on load it keeps only what differs, and follows the
    # defaults afterwards.
    account, anna, _bert = layered
    _write(account, {"agent": {"model": "kit.default"}})
    full = us._merge_missing_defaults(
        us._normalized_settings_dict({"agent": {"model": "kit.default",
                                                "effort": "high"}}),
        us.DEFAULT_SETTINGS)
    anna.write_text(json.dumps(full, indent=2, sort_keys=True))
    eff = us.load_settings(anna)
    assert eff["agent"]["effort"] == "high"
    own = json.loads(anna.read_text())
    assert own["agent"] == {"effort": "high"}
    _write(account, {"agent": {"model": "kit.new"}})
    assert us.load_settings(anna)["agent"]["model"] == "kit.new"


def test_a_missing_personal_file_is_not_created_by_loading(layered):
    account, anna, _bert = layered
    _write(account, {"agent": {"model": "kit.default"}})
    assert us.load_settings(anna)["agent"]["model"] == "kit.default"
    assert not anna.exists()


def test_a_broken_defaults_file_is_an_empty_layer(layered):
    account, anna, _bert = layered
    account.write_text("{not json")
    _write(anna, {"agent": {"effort": "high"}})
    eff = us.load_settings(anna)
    assert eff["agent"]["effort"] == "high"
    assert eff["runtime"]["backend"] == us.DEFAULT_SETTINGS["runtime"]["backend"]
    assert account.read_text() == "{not json"


def test_security_notices_of_the_defaults_are_cleared_in_the_personal_file(layered):
    account, anna, _bert = layered
    _write(account, {"security_notices": ["shown once"]})
    _write(anna, {})
    assert us.take_security_notices(anna) == ["shown once"]
    assert us.take_security_notices(anna) == []
    assert json.loads(account.read_text()) == {"security_notices": ["shown once"]}


def test_the_security_update_of_the_defaults_does_not_override_a_person(layered):
    # The defaults file says bash_isolation off and was never migrated: the
    # one-time update applies to everyone's view, and a person who then sets
    # "off" on purpose is respected.
    account, anna, _bert = layered
    _write(account, {"agent": {"bash_isolation": "off"}})
    _write(anna, {})
    assert us.load_settings(anna)["agent"]["bash_isolation"] == "auto"
    eff = us.load_settings(anna)
    eff["agent"]["bash_isolation"] = "off"
    us.save_settings(eff, anna)
    assert us.load_settings(anna)["agent"]["bash_isolation"] == "off"


def test_readers_of_raw_keys_see_the_defaults(layered):
    account, anna, _bert = layered
    _write(account, {"agent": {
        "model_overrides": {"kit.example": {"max_tool_rounds": 7}},
        "permissions": {"deny": ["Bash(shutdown *)"]}}})
    _write(anna, {"agent": {"permissions": {"deny": ["Bash(reboot *)"]}}})
    from delfin.agent import model_profiles as mp
    mp._OVERRIDES_CACHE = None
    try:
        assert mp._load_overrides() == {"kit.example": {"max_tool_rounds": 7}}
    finally:
        mp._OVERRIDES_CACHE = None


def test_the_defaults_deny_rules_reach_a_person_with_rules_of_their_own(layered, monkeypatch):
    account, anna, _bert = layered
    _write(account, {"permissions": {"deny": ["Bash(shutdown *)"],
                                     "allow": ["Read"]}})
    _write(anna, {"permissions": {"deny": ["Bash(reboot *)"]}})
    from delfin.agent import permission_rules as pr
    monkeypatch.setattr(us, "get_settings_path", lambda *a, **k: anna)
    pr._CACHE.clear()
    try:
        rs = pr._load_local(None)
        values = [str(r.value) for r in rs.deny]
        assert any("shutdown" in v for v in values), values
        assert any("reboot" in v for v in values), values
        assert any(str(r.value) == "Read" for r in rs.allow)
        assert not rs.errors
        assert rs.user_path == str(anna)
    finally:
        pr._CACHE.clear()


# ---------------------------------------------------------------------------
# pi_install.sh writes the defaults, never a person's file
# ---------------------------------------------------------------------------

_PI_INSTALL = ROOT / "scripts" / "pi" / "pi_install.sh"


def _heredoc(marker):
    """The python program pi_install.sh feeds to "$VENV_PY" after ``marker``,
    and the command line that runs it."""
    text = _PI_INSTALL.read_text(encoding="utf-8")
    m = re.search(r"^([^\n]*" + re.escape(marker) + r"[^\n]*<<'EOF'[^\n]*)\n(.*?)\nEOF\n",
                  text, re.M | re.S)
    assert m, marker
    return m.group(1), m.group(2)


def _run_apply(program, args, home):
    env = {k: v for k, v in os.environ.items()
           if not k.startswith("DELFIN_") and k not in ("PYTHONPATH",)}
    env.update(HOME=str(home), PYTHONPATH=str(ROOT))
    return subprocess.run([sys.executable, "-", *args], input=program,
                          env=env, cwd=str(home), text=True,
                          capture_output=True, timeout=120)


@pytest.mark.parametrize("marker", ['"$agent_conf"', '"$tabs_conf"'])
def test_the_apply_blocks_clear_the_layer_variable(marker):
    line, _program = _heredoc(marker)
    assert "env -u DELFIN_SETTINGS_DEFAULTS" in line, line


def test_the_apply_blocks_leave_a_personal_override_alone(tmp_path, monkeypatch):
    account_home = tmp_path / "account"
    account_home.mkdir()
    account = account_home / us.SETTINGS_FILE_NAME
    anna = tmp_path / "anna" / us.SETTINGS_FILE_NAME
    bert = tmp_path / "bert" / us.SETTINGS_FILE_NAME
    anna.parent.mkdir()
    bert.parent.mkdir()
    _write(anna, {"agent": {"model": "kit.mine"},
                  "ui": {"tabs": {"order": ["home"], "hidden": []}}})
    _write(bert, {})
    before = {p: p.read_bytes() for p in (anna, bert)}

    agent_conf = tmp_path / "pi-agent.conf"
    agent_conf.write_text("mode=solo\nprovider=kit\nmodel=kit.glm-5.3-flash\n"
                          "permission_profile=all_free\neffort=\n")
    tabs_conf = tmp_path / "pi-tabs.conf"
    tabs_conf.write_text("agent home\n")
    _line, agent_program = _heredoc('"$agent_conf"')
    _line, tabs_program = _heredoc('"$tabs_conf"')
    for program, args in (
            (agent_program, [str(agent_conf)]),
            (tabs_program, [str(tabs_conf), str(ROOT / "delfin" / "dashboard" / "__init__.py")])):
        result = _run_apply(program, args, account_home)
        assert result.returncode == 0, result.stderr
    assert {p: p.read_bytes() for p in (anna, bert)} == before

    written = json.loads(account.read_text())
    assert written["agent"]["model"] == "kit.glm-5.3-flash"
    assert written["agent"]["permission_profile"] == "all_free"
    assert written["ui"]["tabs"]["order"][:2] == ["agent", "home"]

    monkeypatch.setenv("HOME", str(account_home))
    monkeypatch.setenv(ENV, str(account))
    a, b = us.load_settings(anna), us.load_settings(bert)
    assert a["agent"]["model"] == "kit.mine"
    assert a["ui"]["tabs"] == {"order": ["home"], "hidden": []}
    assert a["agent"]["permission_profile"] == "all_free"
    assert b["agent"]["model"] == "kit.glm-5.3-flash"
    assert b["ui"]["tabs"]["order"][:2] == ["agent", "home"]

    # The next deploy changes a default: bert follows, anna keeps hers.
    agent_conf.write_text("model=kit.next\n")
    assert _run_apply(agent_program, [str(agent_conf)], account_home).returncode == 0
    assert us.load_settings(bert)["agent"]["model"] == "kit.next"
    assert us.load_settings(anna)["agent"]["model"] == "kit.mine"
