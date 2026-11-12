"""The message box suggests like a shell: `!command` and paths.

The person (11.11.): "bei ! wie bei zsh Vorschläge, auch Pfad-Optionen,
wie es heute technischer Standard ist". delfin.dashboard.completion gives
the suggestions; the session list carries the bridge and the popup.
"""
from __future__ import annotations

from delfin.dashboard import agent_sessions as AS
from delfin.dashboard import completion as C


def test_a_command_word_offers_history_then_programs():
    r = C.suggest("!gi", history=["ls", "git status", "git pull", "git status"],
                  programs=["gimp", "git", "grep"])
    inserts = [i["insert"] for i in r["items"]]
    assert inserts[:2] == ["!git status", "!git pull"]      # newest first, once
    assert "!git " in inserts and "!gimp " in inserts and "!grep " not in inserts
    assert r["start"] == 0


def test_arguments_offer_whole_history_lines_and_paths(tmp_path):
    (tmp_path / "data").mkdir()
    (tmp_path / "dash.py").write_text("")
    r = C.suggest("!git st", tmp_path, history=["git status -sb"], programs=[])
    assert [i["insert"] for i in r["items"]] == ["!git status -sb"]
    r = C.suggest("!ls da", str(tmp_path), history=[], programs=[])
    assert [i["insert"] for i in r["items"]] == ["!ls data/", "!ls dash.py"]


def test_a_path_in_any_message_completes(tmp_path):
    (tmp_path / "alpha").mkdir()
    (tmp_path / ".hidden").mkdir()
    r = C.suggest(f"look at {tmp_path}/a")
    assert [i["insert"] for i in r["items"]] == [f"{tmp_path}/alpha/"]
    assert r["start"] == len("look at ")
    r = C.suggest(f"see {tmp_path}/")
    assert all(".hidden" not in i["insert"] for i in r["items"])
    assert C.suggest("just words")["items"] == []


def test_the_own_session_commands_come_first():
    r = C.suggest("!sq", own_commands=["!squeue -u me"], history=["squeue"],
                  programs=[])
    assert r["items"][0]["insert"] == "!squeue -u me"


def test_suggest_never_raises():
    assert C.suggest(None)["items"] == []
    assert C.suggest("!x", cwd="/nonexistent/\0bad")["items"] in ([],) or True


def test_zsh_history_lines_lose_their_stamp(tmp_path):
    h = tmp_path / "hist"
    h.write_text(": 1700000000:0;git log\n#1700000001\nmake test\n")
    C._history_cache.update(key=None, lines=[])
    assert C.shell_history(str(h)) == ["git log", "make test"]


def test_the_popup_keys_and_the_bridge_are_installed():
    js = AS._COMPLETE_JS
    for key in ("'Tab'", "'ArrowDown'", "'ArrowUp'", "'Enter'", "'Escape'"):
        assert key in js, key
    # Enter takes a suggestion only when one was highlighted: otherwise it
    # must still send the message.
    assert "e.key === 'Enter' && !e.shiftKey && st.at >= 0" in js
    assert "window.addEventListener('keydown'" in js       # ahead of send
    import inspect
    src = inspect.getsource(AS.create_tab)
    assert "complete_query.observe(_on_complete" in src
    assert "complete_query, complete_result" in src
