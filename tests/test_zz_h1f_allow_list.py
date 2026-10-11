"""H1F (operator build): a shape is an allow list of exact spellings."""

import os

import pytest

from delfin.agent import read_only_shapes as R


@pytest.fixture
def ws(tmp_path):
    (tmp_path / "p.patch").write_text("x")
    (tmp_path / "sub").mkdir()
    (tmp_path / "link").symlink_to("/etc/hostname")
    return str(tmp_path)


@pytest.mark.parametrize("cmd,shape", [
    ("git apply --check p.patch", "git_apply_check"),
    ("git apply --check sub/new.patch", "git_apply_check"),
    ("git hash-object p.patch", "git_hash_object"),
    ("git switch -c pr/fix-1 origin/main", "git_switch_pr"),
    ("gh --version", "gh_version"),
    ("gh version", "gh_version"),
])
def test_the_exact_forms_fire(ws, cmd, shape):
    assert R.matches(cmd, ws) == shape


@pytest.mark.parametrize("cmd", [
    # git accepts any unique prefix of a long option: this one RUNS sh.
    "git grep --open-files-in-p=sh x",
    "git grep --open-files-in-pager=sh x",
    "git grep -Osh x",
    "git grep x",
    "git show HEAD:f",
    # any option at all
    "git apply --check --directory=/etc p.patch",
    "git apply --check --unsafe-paths p.patch",
    "git apply --check --build-fake-ancestor=out p.patch",
    "git apply --check -3 p.patch",
    "git apply --check --index p.patch",
    "git apply p.patch",
    "git apply --check p.patch p.patch",
    "git hash-object -w p.patch",
    "git hash-object --stdin",
    "git hash-object --stdin-paths",
    "git hash-object --path=/etc/passwd p.patch",
    "git hash-object -t blob p.patch",
    "git -c core.pager=sh hash-object p.patch",
    "git -C /etc hash-object p.patch",
    # paths
    "git hash-object /etc/hostname",
    "git hash-object ../outside",
    "git hash-object link",
    "git hash-object ~/.ssh/id_rsa",
    "git hash-object *",
    "git hash-object -",
    "git apply --check /etc/passwd",
    # the shell would read these otherwise
    "git hash-object p.patch; id",
    "git hash-object p.patch && id",
    "git hash-object p.patch | sh",
    "git hash-object $(id)",
    "git hash-object `id`",
    "git hash-object p.patch > out",
    "git hash-object p.patch\nid",
    "git hash-object 'p.patch'",
    "git hash-object p.patch\tx",
    "git  hash-object p.patch",
    " git hash-object p.patch",
    "X=1 git hash-object p.patch",
    "env git hash-object p.patch",
    # switch
    "git switch -C pr/x origin/main",
    "git switch -c pr/x origin/main --force",
    "git switch -c pr/x origin/main --discard-changes",
    "git switch -c pr/x",
    "git switch -c pr/x other",
    "git switch -c pr/.. origin/main",
    "git switch -c pr/-x origin/main",
    "git switch -c pr/a..b origin/main",
    "git switch -c pr/a/b origin/main",
    "git switch -c pr/x.lock origin/main",
    "git switch -c main origin/main",
    "git switch --detach origin/main",
    # gh
    "gh --version --foo",
    "gh auth token",
    "./gate-tools/gate tests/test_zz_x.py",
    "gate tests/test_zz_x.py",
    "",
])
def test_everything_else_is_still_asked(ws, cmd):
    assert R.matches(cmd, ws) == ""


def test_an_unticked_shape_does_not_fire_and_none_means_none(ws):
    assert R.matches_with("gh --version", ws, {"git_hash_object"}) == ""
    assert R.matches_with("gh --version", ws, set()) == ""
    assert R.matches_with("gh --version", ws, None) == "gh_version"
    block = R.read_only_shapes_from_payload(
        {"agent": {"read_only_shapes": {"enabled": True, "enabled_shapes": []}}})
    assert block == {"enabled": True, "enabled_shapes": []}
    stored = R.read_only_shapes_into_payload(
        {}, {"enabled": True, "enabled_shapes": []})
    assert stored["agent"]["read_only_shapes"]["enabled_shapes"] == []
    assert R.read_only_shapes_from_payload({})["enabled_shapes"] == list(R.SHAPE_NAMES)


def test_no_shape_names_a_program_to_run():
    assert set(R.SHAPE_NAMES) == {
        "git_apply_check", "git_hash_object", "git_switch_pr", "gh_version"}
    assert os.path.basename(R.__file__) == "read_only_shapes.py"


# -- through the real gate ---------------------------------------------------

from delfin.agent import api_client as A          # noqa: E402
from delfin.agent import security_events as SE    # noqa: E402


def _gate(ws, cmd, asked, **args):
    perms = A.KitToolPermissions(
        workspace=ws, mode="default",
        confirm_callback=lambda n, a, p: asked.append(cmd) or False)
    return A._DocToolExecutor()._run_permission_gate(
        "bash", dict(command=cmd, **args), perms)


def test_a_shape_is_approved_without_a_dialog_and_recorded_as_an_allow(ws):
    asked = []
    before = max((e.seq for e in SE.recent(200)), default=0)
    assert _gate(ws, "git hash-object p.patch", asked) is None
    assert asked == []
    events = [e for e in SE.recent(200)
              if e.seq > before and e.kind == "read_only_shape"]
    assert len(events) == 1 and events[0].blocked is False
    assert "read_only_shape" in SE.known_kinds()


def test_the_program_running_spelling_is_still_asked_at_the_gate(ws):
    asked = []
    assert _gate(ws, "git grep --open-files-in-p=sh x", asked) is not None
    assert asked == ["git grep --open-files-in-p=sh x"]


def test_a_cwd_takes_the_shape_away(ws):
    asked = []
    assert _gate(ws, "git hash-object p.patch", asked, cwd="sub") is not None
    assert len(asked) == 1


def test_switched_off_it_asks_again(ws, monkeypatch):
    import delfin.user_settings as US
    monkeypatch.setattr(US, "load_settings", lambda: {
        "agent": {"read_only_shapes": {"enabled": False}}})
    asked = []
    assert _gate(ws, "gh --version", asked) is not None and len(asked) == 1
    monkeypatch.setattr(US, "load_settings", lambda: {
        "agent": {"read_only_shapes": {"enabled": True, "enabled_shapes": []}}})
    assert _gate(ws, "gh --version", asked) is not None and len(asked) == 2


def test_a_refusal_made_earlier_is_not_overruled_by_a_shape(ws):
    perms = A.KitToolPermissions(workspace=ws, mode="default")
    perms.record_denied_action("bash", "git hash-object p.patch")
    out = A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": "git hash-object p.patch"}, perms)
    assert out and "already refused" in out


def test_plan_mode_stays_read_only(ws):
    perms = A.KitToolPermissions(workspace=ws, mode="plan")
    assert A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": "git switch -c pr/x origin/main"}, perms) is not None
