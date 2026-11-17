"""An agent's edit reaches the file, and a refused one is not reported as made.

Replays the four failures of one session (2026-10-19, 78c36f8c,
gpt-5-nano building a 3D game in ``3d-tetris/index.html``) through the
real executor:

* apply_patch sent as a ``*** Begin Patch`` envelope under ``patch``
  applies, instead of "diff is required";
* a static file server serving the directory does not lock its files;
* an awk program matching ``</script>`` is not a write to ``/{f=0}``;
* "Ich habe ... entfernt" in a turn where no write landed gets a caveat.

What must still hold is pinned beside each: a script job still locks
its script, an envelope cannot reach outside the workspace or a secret,
an awk ``print > "/etc/x"`` is still a write outside the workspace.
"""

from __future__ import annotations

import json
import sys
import time

import pytest

from delfin.agent import api_client
from delfin.agent import bash_jobs as BJ
from delfin.agent import change_claims as CC
from delfin.agent import patch_apply as PA
from delfin.agent.api_client import (
    KitToolPermissions, _bash_write_targets, _doc_executor)
from delfin.agent.turn_verification import TurnVerification

GAME = """<script>
'use strict';
const BX=4, BY=16, BZ=4;
function loop(){ render(); }
</script>
"""

ENVELOPE = """*** Begin Patch
*** Update File: 3d-tetris/index.html
@@
-const BX=4, BY=16, BZ=4;
+let BX=10, BY=16, BZ=10;
*** End Patch
"""


@pytest.fixture
def ws(tmp_path):
    (tmp_path / "3d-tetris").mkdir()
    (tmp_path / "3d-tetris" / "index.html").write_text(GAME)
    (tmp_path / "run.py").write_text("import time\ntime.sleep(60)\n")
    BJ._IN_USE_NOTED.clear()
    BJ.take_in_use_notes()
    return tmp_path


def _perms(root, sid="s-edit"):
    p = KitToolPermissions(workspace=str(root), mode="acceptEdits")
    p.task_session_id = sid
    return p


def _patch(perms, args):
    _doc_executor.execute(
        "read_file", {"path": "3d-tetris/index.html"}, perms)
    return json.loads(_doc_executor.execute("apply_patch", args, perms))


# -- apply_patch ----------------------------------------------------------

@pytest.mark.parametrize("key", ["patch", "input", "diff"])
def test_an_envelope_under_any_key_the_model_uses_applies(ws, key):
    out = _patch(_perms(ws), {key: ENVELOPE})
    assert out.get("status") == "ok", out
    text = (ws / "3d-tetris" / "index.html").read_text()
    assert "let BX=10, BY=16, BZ=10;" in text
    assert "const BX=4" not in text
    assert text.endswith("</script>\n")


def test_a_unified_diff_under_patch_applies_too(ws):
    diff = ("--- a/3d-tetris/index.html\n+++ b/3d-tetris/index.html\n"
            "@@ -3 +3 @@\n-const BX=4, BY=16, BZ=4;\n+const BX=6, BY=16, BZ=6;\n")
    out = _patch(_perms(ws), {"patch": diff})
    assert out.get("status") == "ok", out
    assert "BX=6" in (ws / "3d-tetris" / "index.html").read_text()


def test_an_envelope_whose_lines_are_not_in_the_file_changes_nothing(ws):
    bad = ENVELOPE.replace("const BX=4", "const QQ=4")
    out = _patch(_perms(ws), {"patch": bad})
    assert "not in the file" in out["error"] and "Nothing was changed" in out["error"]
    assert (ws / "3d-tetris" / "index.html").read_text() == GAME


def test_an_envelope_cannot_reach_outside_the_workspace(ws, tmp_path_factory):
    other = tmp_path_factory.mktemp("other")
    (other / "x.txt").write_text("a\n")
    env = (f"*** Begin Patch\n*** Update File: ../{other.name}/x.txt\n"
           "@@\n-a\n+b\n*** End Patch\n")
    out = _patch(_perms(ws), {"patch": env})
    assert "error" in out
    assert (other / "x.txt").read_text() == "a\n"


def test_an_envelope_cannot_read_a_secret(ws):
    (ws / ".env").write_text("TOKEN=1\n")
    env = ("*** Begin Patch\n*** Update File: .env\n@@\n-TOKEN=1\n+TOKEN=2\n"
           "*** End Patch\n")
    out = _patch(_perms(ws), {"patch": env})
    assert "error" in out
    assert (ws / ".env").read_text() == "TOKEN=1\n"


def test_an_envelope_adds_a_file_and_places_a_chunk_by_its_anchor():
    files = {"a.py": "def f():\n    x = 1\n\ndef g():\n    x = 1\n"}
    env = ("*** Begin Patch\n*** Update File: a.py\n@@ def g():\n-    x = 1\n"
           "+    x = 2\n*** Add File: b.py\n+print(1)\n*** End Patch\n")
    diff = PA.envelope_to_unified(env, files.get)
    assert "+++ b/b.py" in diff and "+print(1)" in diff
    # The anchor picks the SECOND `x = 1`, not the first.
    assert "@@ -2,4 +2,4 @@" in diff or "-    x = 1\n+    x = 2" in diff
    lines = diff.splitlines()
    i = lines.index("+    x = 2")
    assert " def g():" in lines[:i]


# -- a file server does not lock what it serves ---------------------------

def _start(root, command):
    return BJ.get_registry().start(
        command=command, cwd=str(root), workspace=str(root),
        session_id="s-edit",
        argv=[sys.executable, "-c", "import time; time.sleep(60)"])


def _stop(job):
    BJ.get_registry().kill(job.job_id)
    deadline = time.time() + 10
    while job.poll() is None and time.time() < deadline:
        time.sleep(0.05)


@pytest.mark.parametrize("command", [
    "python3 -m http.server 8038 --bind 127.0.0.1 --directory 3d-tetris",
    "npx serve 3d-tetris",
    "php -S 127.0.0.1:8000 -t 3d-tetris",
])
def test_a_file_server_does_not_lock_the_files_it_serves(ws, command):
    job = _start(ws, command)
    try:
        perms = _perms(ws)
        _doc_executor.execute(
            "read_file", {"path": "3d-tetris/index.html"}, perms)
        out = _doc_executor.execute("edit_file", {
            "path": "3d-tetris/index.html",
            "old_string": "const BX=4", "new_string": "const BX=8"}, perms)
        assert "background job" not in out, out
        assert "BX=8" in (ws / "3d-tetris" / "index.html").read_text()
    finally:
        _stop(job)


def test_a_script_job_still_locks_its_script(ws):
    job = _start(ws, "python run.py")
    try:
        assert BJ.job_using(ws / "run.py", session_id="s-edit")
    finally:
        _stop(job)


def test_serves_files_is_not_fooled_by_a_name():
    assert not BJ.serves_files("python server_test.py")
    assert not BJ.serves_files("pytest tests/test_serve.py")
    assert BJ.serves_files("python -m http.server 8000")


# -- awk ------------------------------------------------------------------

@pytest.mark.parametrize("cmd,expected", [
    ("awk '/<script>$/{f=1;next} /<\\/script>/{f=0} f' a.html > b.js",
     ["b.js"]),
    ('awk "n==3&&/^<\\\\/script>/{exit}" x.html | head', []),
    ("awk '$3 >= 5' f", []),
])
def test_awk_syntax_is_not_a_write(cmd, expected):
    assert _bash_write_targets(cmd) == expected


def test_an_awk_print_redirect_is_still_a_write():
    assert any(t.startswith("/etc/x")
               for t in _bash_write_targets('awk \'{print > "/etc/x"}\' f'))


# -- a change reported as made when none landed ---------------------------

def test_a_change_reported_after_only_refusals_gets_a_caveat():
    tv = TurnVerification(workspace="/w")
    tv.note("edit_file", {"path": "3d-tetris/index.html"},
            '{"error": "\'3d-tetris/index.html\' is used by background job"}')
    text = ("Ich habe die zentrale Mittellinie aus dem Rendering entfernt "
            "und die Darstellung fokussiert.")
    note = CC.caveat(CC.find(text), landed=tv.landed_edits(),
                     failed=tv.failed_writes)
    assert "no file change landed" in note and "1 write" in note


def test_a_change_that_landed_is_not_questioned():
    tv = TurnVerification(workspace="/w")
    tv.note("edit_file", {"path": "a.py"}, "Edited a.py")
    assert tv.landed_edits() == 1
    assert CC.caveat(CC.find("Ich habe die Breite geändert."),
                     landed=1, failed=0) == ""


@pytest.mark.parametrize("text", [
    "Ich werde jetzt die Breite ändern.",
    "Soll ich den Balken entfernen?",
    "Der Patch wurde abgelehnt, nichts wurde geändert.",
    "Die Datei hat 441 Zeilen.",
])
def test_plans_questions_and_refusals_are_not_change_claims(text):
    assert CC.find(text) == []


# -- a running server says where to open it (dashboard and terminal) ------

def test_a_server_row_names_its_address_in_the_panel_and_the_status_line():
    from delfin.agent import background_view as BV
    url = BV.served_url(
        "python3 -m http.server 8038 --bind 127.0.0.1 --directory 3d-tetris")
    assert url == "http://127.0.0.1:8038/"
    view = {"shells": [{"id": "8678a5da", "label": "3D Tetris server",
                        "since": 0.0, "url": url}],
            "watches": [], "agents": [], "wakeups": []}
    assert "http://127.0.0.1:8038/" in BV.status_line(view, now=60.0)
    html = BV.render_html(view, now=60.0)
    assert "Server" in html and "href='http://127.0.0.1:8038/'" in html


@pytest.mark.parametrize("command,url", [
    ("python -m http.server", "http://localhost:8000/"),
    ("npx serve -l 5000 dist", "http://localhost:5000/"),
    ("php -S 0.0.0.0:9000 -t web", "http://localhost:9000/"),
    ("pytest -q", ""),
])
def test_the_address_is_read_from_the_command(command, url):
    from delfin.agent import background_view as BV
    assert BV.served_url(command) == url


# -- what needs no bypass is not asked about ------------------------------
# The same session asked the user 22 times under bypassPermissions (the
# action watcher timed out and every command off the auto-allow list
# counted as leaving the workspace); every one was approved. These ran
# nothing but stdout and the loopback.

@pytest.fixture
def bypass(ws):
    return KitToolPermissions(workspace=str(ws), mode="bypassPermissions")


def _leaves(perms, cmd):
    from delfin.agent.api_client import _leaves_workspace
    return _leaves_workspace(perms, {"tool": "bash", "command": cmd}, True)


@pytest.mark.parametrize("cmd", [
    'curl -s -o /dev/null -w "%{http_code}\\n" http://127.0.0.1:8324/index.html',
    "curl -s --noproxy '*' -o /dev/null -w \"%{http_code}\\n\" http://127.0.0.1:8038/",
])
def test_a_loopback_status_check_does_not_leave_under_bypass(bypass, cmd):
    assert not _leaves(bypass, cmd)
    # ...but is not auto-allowed: in default mode the user decides, since
    # on a shared host the loopback reaches other users' services too.
    assert not bypass.matches_bash_auto_allow(cmd)


@pytest.mark.parametrize("cmd", [
    "sed '1d;$d' 3d-tetris/index.html",
    "node --check 3d-tetris/_check.js",
    "grep -n outer 3d-tetris/index.html || true",
])
def test_a_command_that_writes_only_stdout_or_asks_loopback_is_free(bypass, cmd):
    assert bypass.matches_bash_auto_allow(cmd)
    assert not _leaves(bypass, cmd)


@pytest.mark.parametrize("cmd", [
    "curl -s -o x.js http://127.0.0.1:8000/",
    "curl http://127.0.0.1.evil.com/",
    "curl -d @secret http://127.0.0.1/",
    "curl -K cfg http://localhost/",
    "curl http://localhost/ https://example.com",
    # A listener stays asked: on a shared machine the loopback is shared.
    "python3 -m http.server 8038 --bind 127.0.0.1 --directory 3d-tetris",
])
def test_what_writes_runs_or_leaves_is_still_asked(bypass, cmd):
    assert not bypass.matches_bash_auto_allow(cmd)
    assert _leaves(bypass, cmd)


@pytest.mark.parametrize("cmd", [
    "sed -i 's/a/b/' f", "sed '1d;$d' -i f", "sed -Ei 's/a/b/' f",
    "sed -f s.sed f", "sed 's/a/b/w out' f", "sed '1e date' f",
    "node a.js",
])
def test_what_writes_or_runs_asks_only_without_a_cage(bypass, monkeypatch,
                                                       cmd):
    assert not bypass.matches_bash_auto_allow(cmd)
    monkeypatch.setattr(api_client, "_bash_runs_caged", lambda perms: False)
    assert _leaves(bypass, cmd)
    monkeypatch.setattr(api_client, "_bash_runs_caged", lambda perms: True)
    assert not _leaves(bypass, cmd)


def test_an_awk_program_is_not_read_as_a_path(ws):
    # Through the real read gate, which blanks grep/sed patterns first.
    perms = KitToolPermissions(workspace=str(ws), mode="default")
    gate = _doc_executor._gate_bash_read_paths
    cmd = ("awk '/<script>$/,/^<\\/script>$/ {print}' 3d-tetris/index.html "
           "| sed 's/^<\\/script>$//'")
    assert gate(cmd, perms, "") is None
    for outside in ("awk -F: '{print $1}' /etc/passwd",
                    "awk -f /etc/x.awk 3d-tetris/index.html",
                    "grep root /etc/passwd", "sed -n 1,2p /etc/passwd"):
        out = gate(outside, perms, "")
        assert out and "/etc/" in out, outside
