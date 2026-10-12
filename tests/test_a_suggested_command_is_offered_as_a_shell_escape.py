"""A shell command the answer gives the user is offered as `! <command>`.

Input: the final text of a turn. A fenced ``bash``/``sh``/``shell``/
``zsh`` block, or a ``console`` block's ``$ `` lines, holding ONE logical
command. Output: ``! <command>`` as the message box's offer -- the
placeholder Tab takes and the button beside the box -- in the style every
agent offer has. Nothing runs until the user sends it, and sending it is
the same `!` path, through the same gate, as typing it.

Field case: the answer said "Du kannst dir die Partitionen selbst so
anzeigen lassen" over a ``bash`` block holding ``sinfo -o "%P" | sort -u``,
and the user had to retype it with a ``!`` in front.
"""

from __future__ import annotations

import json
import threading

import pytest

pytest.importorskip("ipywidgets")

import ipywidgets as widgets                                 # noqa: E402

from delfin.agent import scheduler as S                      # noqa: E402
from delfin.dashboard import tab_agent                       # noqa: E402
from delfin.dashboard.context import DashboardContext        # noqa: E402
from delfin.dashboard.tab_agent import proposed_shell_command  # noqa: E402

FIELD_ANSWER = (
    "Die Partitionen konnte ich nicht abfragen. Du kannst dir die "
    "Partitionen selbst so anzeigen lassen (im Terminal auf dem Cluster):\n"
    "```bash\nsinfo -o \"%P\" | sort -u\n```\n")
FIELD_OFFER = '! sinfo -o "%P" | sort -u'


# ---------------------------------------------------------------------------
# What is offered
# ---------------------------------------------------------------------------

def test_the_field_case_is_offered_as_a_shell_escape():
    assert proposed_shell_command(FIELD_ANSWER) == FIELD_OFFER


@pytest.mark.parametrize("lang", ["bash", "sh", "shell", "zsh", "BASH"])
def test_every_shell_tag_is_read(lang):
    assert proposed_shell_command(f"x\n```{lang}\nsqueue -u $USER\n```") == (
        "! squeue -u $USER")


def test_comments_blank_lines_and_a_prompt_sign_are_not_the_command():
    text = "```bash\n# the partitions\n\n$ sinfo -s\n```"
    assert proposed_shell_command(text) == "! sinfo -s"


def test_a_console_block_offers_its_prompt_line_not_its_output():
    text = "```console\n$ sinfo -s\nPARTITION AVAIL\ncpu up\n```"
    assert proposed_shell_command(text) == "! sinfo -s"


def test_a_continued_line_is_one_command():
    text = "```bash\nsqueue -u $USER \\\n  -o '%i %T'\n```"
    assert proposed_shell_command(text) == "! squeue -u $USER -o '%i %T'"


@pytest.mark.parametrize("body", [
    "cd /tmp\nls",                                   # two commands
    "cat > f <<'EOF'\nx\nEOF",                        # a heredoc
    "for f in *; do\n  echo $f\ndone",                # a loop over lines
    "$ ls\n$ pwd",                                    # two prompt lines
])
def test_a_block_of_several_commands_is_not_offered(body):
    assert proposed_shell_command(f"```bash\n{body}\n```") is None


def test_the_first_single_command_block_is_offered():
    text = ("```bash\ncd /tmp\nls\n```\nund dann\n"
            "```bash\nsqueue -u $USER\n```\n```bash\nsinfo\n```")
    assert proposed_shell_command(text) == "! squeue -u $USER"


@pytest.mark.parametrize("text", [
    "```python\nprint('x')\n```",
    "```\nsinfo\n```",                                # untagged: not a claim
    "Run `sinfo -s` yourself.",                       # inline code
    "", None,
])
def test_what_is_not_a_shell_block_offers_nothing(text):
    assert proposed_shell_command(text) is None


def test_a_command_longer_than_the_box_shows_is_not_offered():
    """Tab copies the placeholder verbatim, and a long placeholder is cut
    with an ellipsis: a cut command must never be the thing Tab sends."""
    long_cmd = "echo " + "a" * 200
    assert proposed_shell_command(f"```bash\n{long_cmd}\n```") is None


def test_a_command_already_written_with_the_bang_keeps_one():
    assert proposed_shell_command("```bash\n! sinfo -s\n```") == "! sinfo -s"


def test_control_characters_are_never_offered():
    assert proposed_shell_command("```bash\nls\x1b[2J\n```") is None


# ---------------------------------------------------------------------------
# The tab: offered, not run; sent through the same `!` path
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _scheduler(tmp_path, monkeypatch):
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))


class _Engine:
    """Records every command the gate is asked to run."""

    kit_permissions = None

    def __init__(self, result=None):
        self.calls: list = []
        self.typed: list = []
        self._result = result or json.dumps(
            {"exit_code": 0, "stdout": "cpu\ngpu\n", "stderr": ""})

    def run_gated_bash(self, command, **kw):
        self.calls.append(command)
        return self._result

    def note_typed(self, text):
        self.typed.append(text)


def _tab(tmp_path):
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    ctx.agent_workspace = str(tmp_path / "ws")
    tab, refs = tab_agent.create_tab(ctx)
    return tab, refs


def _walk(root):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(getattr(node, "children", ()) or ())


def _input(tab):
    return next(w for w in _walk(tab) if isinstance(w, widgets.Textarea)
                and "delfin-agent-input" in w._dom_classes)


def _send(tab):
    return next(w for w in _walk(tab) if isinstance(w, widgets.Button)
                and w.description == "Send")


def _offer_buttons(tab):
    return [w for w in _walk(tab) if isinstance(w, widgets.Button)
            and "delfin-next-step" in w._dom_classes]


def _wait_escape():
    for t in list(threading.enumerate()):
        if t.name == "dashboard-shell-escape":
            t.join(10)


def _chat(tab):
    return "\n".join(w.value for w in _walk(tab) if isinstance(w, widgets.HTML)
                     and "delfin-chat" in str(w.value))


def test_the_answer_fills_the_offer_and_nothing_runs(tmp_path):
    tab, refs = _tab(tmp_path)
    try:
        eng = _Engine()
        refs["state"]["engine"] = eng
        refs["state"]["_last_answer_text"] = FIELD_ANSWER
        refs["refresh_next_steps"]()
        box = _input(tab)
        assert box.placeholder == FIELD_OFFER
        assert "delfin-agent-input-proposed" in box._dom_classes
        assert box.value == "", "an offer is never typed for the user"
        buttons = [b for b in _offer_buttons(tab) if b.tooltip == FIELD_OFFER]
        assert buttons, "the offer has its button beside the box"
        buttons[0].click()
        assert box.value == FIELD_OFFER
        assert eng.calls == [], "filling the box runs nothing"
    finally:
        refs["shutdown"]()


def test_sending_the_offer_is_the_shell_escape_and_feeds_no_user_model(
        tmp_path, monkeypatch):
    noted: list = []
    monkeypatch.setattr("delfin.agent.user_model.note_typed_message",
                        lambda text: noted.append(text))
    tab, refs = _tab(tmp_path)
    try:
        eng = _Engine()
        refs["state"]["engine"] = eng
        refs["state"]["_last_answer_text"] = FIELD_ANSWER
        refs["refresh_next_steps"]()
        next(b for b in _offer_buttons(tab) if b.tooltip == FIELD_OFFER).click()
        _send(tab).click()
        _wait_escape()
        assert eng.calls == ['sinfo -o "%P" | sort -u']
        assert noted == [] and eng.typed == [], (
            "a shell escape is not the user describing themselves")
        assert _input(tab).value == ""
    finally:
        refs["shutdown"]()


# ---------------------------------------------------------------------------
# The output: all of it, as a terminal shows it
# ---------------------------------------------------------------------------

def _special_result(n_lines=250):
    out = "".join(f"  {i:04d} %P *a* _b_ <b>x</b> & |c|d|\n"
                  for i in range(n_lines))
    return json.dumps({"exit_code": 3, "stdout": out,
                       "stderr": "warn: <tag> 50% *x*\n"})


def test_the_whole_output_reaches_the_chat_escaped_with_its_exit_code(tmp_path):
    tab, refs = _tab(tmp_path)
    try:
        eng = _Engine(_special_result())
        refs["state"]["engine"] = eng
        _input(tab).value = "! printf stuff"
        _send(tab).click()
        _wait_escape()
        html = _chat(tab)
        for i in range(250):
            assert (f"  {i:04d} %P *a* _b_ &lt;b&gt;x&lt;/b&gt; &amp; |c|d|"
                    in html), i
        assert "<b>x</b>" not in html and "<i>a</i>" not in html
        assert "warn: &lt;tag&gt; 50% *x*" in html
        assert "exit 3" in html
        assert "delfin-shell-out" in html and "<pre" in html
        assert "more lines" not in html, "nothing is cut below the chat cap"
    finally:
        refs["shutdown"]()


def test_output_beyond_the_chat_cap_is_saved_and_its_path_shown(tmp_path):
    tab, refs = _tab(tmp_path)
    try:
        big = "".join(f"line {i}\n" for i in range(60_000))
        eng = _Engine(json.dumps({"exit_code": 0, "stdout": big, "stderr": ""}))
        refs["state"]["engine"] = eng
        _input(tab).value = "! seq 60000"
        _send(tab).click()
        _wait_escape()
        html = _chat(tab)
        saved = list((tmp_path / "ws").rglob("*.txt"))
        assert len(saved) == 1, saved
        assert saved[0].read_text() == big
        assert str(saved[0]) in html
        assert f"{len(big):,}" in html or str(len(big)) in html
    finally:
        refs["shutdown"]()


def test_the_html_renderer_escapes_the_command_too():
    html = tab_agent.shell_output_html({
        "command": "echo '<script>x</script>'", "stdout": "", "stderr": "",
        "exit_code": 0, "refused": "", "error": ""})
    assert "<script>" not in html and "&lt;script&gt;" in html


def test_a_refusal_is_shown_as_one():
    html = tab_agent.shell_output_html({
        "command": "rm -rf /", "stdout": "", "stderr": "", "exit_code": None,
        "refused": "blocked: <deny-list>", "error": ""})
    assert "refused" in html and "&lt;deny-list&gt;" in html


# ---------------------------------------------------------------------------
# The gate's own output cap: lifted for the person, kept for the model
# ---------------------------------------------------------------------------

def _engine_with_real_gate(tmp_path, monkeypatch):
    from unittest.mock import MagicMock, patch

    from delfin.agent import engine as E
    from delfin.agent.api_client import KitToolPermissions

    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                            model="m", mode="solo")
    perms = KitToolPermissions(mode="default", workspace=str(tmp_path))
    monkeypatch.setattr(type(eng), "kit_permissions",
                        property(lambda self: perms))
    return eng, perms


def test_a_person_gets_the_whole_output_and_the_model_keeps_its_cap(
        tmp_path, monkeypatch):
    from delfin.agent import repl as R
    from delfin.agent.api_client import _doc_executor

    eng, perms = _engine_with_real_gate(tmp_path, monkeypatch)
    run = R.shell_escape_result(eng, "seq 1 20000")
    assert run["refused"] == "" and run["error"] == "", run
    assert run["exit_code"] == 0
    assert run["stdout"].splitlines() == [str(i) for i in range(1, 20001)]

    model_view = json.loads(_doc_executor.execute(
        "bash", {"command": "seq 1 20000"}, perms))
    assert "truncated" in model_view["stdout"], (
        "the model's cap is unchanged by a person's escape")
    assert len(model_view["stdout"]) < len(run["stdout"])
