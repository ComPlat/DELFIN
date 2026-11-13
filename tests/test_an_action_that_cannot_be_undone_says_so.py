"""An approval for what cannot be undone says so, before the answer keys.

Input: a call the gate puts to a person -- a shell command (a delete of a
file nothing keeps, a forced or plain push, a rebase, a job cancel, an
upload), a cluster cancel -- and the surfaces it reaches: the terminal
dialog (terminal_confirm + repl), the dashboard row (kit_confirm), the
service's description and the phone card (messenger.render, German texts in
messenger.lang.DE).

Output and semantics:

- ``irreversible.classify`` names the level (reversible, recoverable with
  effort, irreversible) and, for the last two, the concrete facts: the file
  and its size, the branch, remote and number of commits, the host, the job;
  and the reversible way where there is one (the trash, a new branch);
- ``_ask_confirm`` puts that statement at the head of every preview it hands
  a broker, so each surface reads the same sentences;
- an irreversible request offers no grant that outlives it (no ``always``,
  no ``commands like this``, no ``stop asking about writes``, no ``Allow +
  Permanent``), Enter in the terminal answers no, and an allow pattern, an
  allow rule or a grant does not answer it in the person's place;
- the phone says it in the person's language: every sentence has a German
  text with the same placeholders;
- the rule that the agent itself names the consequence before it asks is in
  the built prompt of every role that can ask.
"""

from __future__ import annotations

import io
import json
import os
import subprocess
import threading
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import irreversible as IR
from delfin.agent import terminal_confirm as tc

BIG = 2 * 1024 * 1024          # over the undo journal's per-file pre-image


def _git(cwd, *args):
    env = {**os.environ, "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t",
           "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t",
           "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}
    return subprocess.run(["git", *args], cwd=cwd, env=env, check=True,
                          capture_output=True, text=True).stdout


@pytest.fixture
def ws(tmp_path):
    w = tmp_path / "ws"
    w.mkdir()
    (w / "big.dat").write_text("E -1234.5678\n" * (BIG // 13 + 1))
    (w / "small.txt").write_text("small\n")
    return w


@pytest.fixture
def repo(tmp_path):
    """A clone whose main is one commit ahead of, and one commit behind,
    origin/main (fetched): a forced push drops exactly one commit."""
    bare = tmp_path / "remote.git"
    _git(tmp_path, "init", "-q", "--bare", "-b", "main", str(bare))
    a = tmp_path / "a"
    _git(tmp_path, "clone", "-q", str(bare), str(a))
    (a / "f").write_text("1\n")
    _git(a, "add", "f")
    _git(a, "commit", "-q", "-m", "one")
    _git(a, "push", "-q", "origin", "HEAD:main")
    b = tmp_path / "b"
    _git(tmp_path, "clone", "-q", str(bare), str(b))
    (b / "g").write_text("2\n")
    _git(b, "add", "g")
    _git(b, "commit", "-q", "-m", "theirs")
    _git(b, "push", "-q", "origin", "main")
    (a / "h").write_text("3\n")
    _git(a, "add", "h")
    _git(a, "commit", "-q", "-m", "mine")
    _git(a, "fetch", "-q", "origin")
    return a


# ---------------------------------------------------------------------------
# One classification
# ---------------------------------------------------------------------------

def test_a_delete_nothing_keeps_names_the_file_its_size_and_the_trash(ws):
    cons = IR.classify("bash", {"command": "rm big.dat"}, cwd=ws)
    assert cons.level == IR.IRREVERSIBLE
    said = IR.statement(cons)
    assert said.startswith(IR.HEAD_IRREVERSIBLE)
    assert "Deletes big.dat (2.0 MiB)" in said
    assert "mv big.dat" in said and ".delfin-trash/" in said


def test_what_the_undo_journal_holds_is_reversible(ws):
    assert IR.classify("bash", {"command": "rm small.txt"}, cwd=ws).level \
        == IR.REVERSIBLE
    assert IR.classify("bash", {"command": "ls -la"}, cwd=ws).level \
        == IR.REVERSIBLE
    assert IR.statement(IR.classify("bash", {"command": "cat big.dat"},
                                    cwd=ws)) == ""


def test_a_forced_push_names_branch_remote_lost_commits_and_a_new_branch(repo):
    cons = IR.classify(
        "bash", {"command": "git push --force-with-lease origin main"},
        cwd=repo)
    assert cons.level == IR.IRREVERSIBLE
    said = IR.statement(cons)
    assert "Replaces the history of branch main on origin" in said
    assert "remote.git" in said
    assert "1 commit(s) there that your branch does not contain" in said
    assert "git push origin HEAD:main-2" in said


def test_a_plain_push_publishes_and_counts_the_commits(repo):
    cons = IR.classify("bash", {"command": "git push origin HEAD:topic"},
                       cwd=repo)
    assert cons.level == IR.IRREVERSIBLE
    assert "Publishes the new branch topic on origin" in IR.statement(cons)


def test_a_rebase_is_recoverable_with_effort(repo):
    cons = IR.classify("bash", {"command": "git rebase origin/main"}, cwd=repo)
    assert cons.level == IR.RECOVERABLE
    assert IR.statement(cons).startswith(IR.HEAD_RECOVERABLE)
    assert "reflog" in IR.statement(cons)


def test_a_job_cancel_and_an_upload_are_irreversible(ws):
    assert "Cancels job 4711" in IR.statement(
        IR.classify("bash", {"command": "scancel 4711"}, cwd=ws))
    said = IR.statement(IR.classify(
        "bash", {"command": "scp big.dat user@hpc.example.org:run/"}, cwd=ws))
    assert "Sends 1 file(s), 2.0 MiB (big.dat), to hpc.example.org" in said
    said = IR.statement(IR.classify(
        "cluster_cancel", {"cluster": "horeka", "job_ids": ["12", "13"]}))
    assert "Cancels job 12, 13 on horeka" in said


def test_the_statement_reads_back_off_a_preview(ws):
    cons = IR.classify("bash", {"command": "rm big.dat"}, cwd=ws)
    preview = IR.with_statement("$ rm big.dat", cons)
    assert IR.with_statement(preview, cons) == preview      # once
    found, rest = IR.split(preview)
    assert found is not None and found.irreversible
    assert rest == "$ rm big.dat"
    assert found.render() == IR.statement(cons).splitlines()
    assert IR.split("$ ls")[0] is None


# ---------------------------------------------------------------------------
# The gate hands every broker the statement
# ---------------------------------------------------------------------------

class Asker:
    def __init__(self, answer=False):
        self.answer = answer
        self.asked = []

    def __call__(self, name, args, preview):
        self.asked.append((name, dict(args or {}), preview))
        return self.answer


def _perms(ws, answer=False, **over):
    perms = A.KitToolPermissions(workspace=ws, mode="default", **over)
    perms.confirm_callback = Asker(answer)
    return perms


def test_the_approval_states_the_consequence_before_the_action(ws):
    perms = _perms(ws)
    A._ask_confirm(perms, "bash", {"command": "rm big.dat"}, "$ rm big.dat")
    (_, _, preview), = perms.confirm_callback.asked
    assert preview.startswith(IR.HEAD_IRREVERSIBLE)
    assert preview.index("Deletes big.dat") < preview.index("$ rm big.dat")


def test_an_allow_pattern_does_not_answer_an_irreversible_command(ws, monkeypatch):
    monkeypatch.setattr(A, "_caged_command_may_run_unasked",
                        lambda *a, **k: False)
    perms = _perms(ws, answer=False)
    perms.bash_auto_allow_patterns = (*perms.bash_auto_allow_patterns,
                                      r"^rm\b")
    # What the journal holds stays answered by the pattern.
    assert A._doc_executor._run_permission_gate(
        "bash", {"command": "rm small.txt"}, perms) is None
    assert perms.confirm_callback.asked == []
    # What nothing keeps is put to the person, with the statement.
    refused = A._doc_executor._run_permission_gate(
        "bash", {"command": "rm big.dat"}, perms)
    assert refused
    (_, _, preview), = perms.confirm_callback.asked
    assert preview.startswith(IR.HEAD_IRREVERSIBLE)


# ---------------------------------------------------------------------------
# The terminal
# ---------------------------------------------------------------------------

def _irreversible_request(ws, command="rm big.dat"):
    cons = IR.classify("bash", {"command": command}, cwd=ws)
    return tc.ConfirmRequest(kind=tc.CONFIRM, tool="bash",
                             args={"command": command},
                             preview=IR.with_statement(f"$ {command}", cons))


def test_the_terminal_says_it_first_and_offers_no_lasting_grant(ws):
    from delfin.agent import repl_render as rr
    req = _irreversible_request(ws)
    keys = {o.key for o in tc.options_for(req, suggestion=r"^rm\b")}
    assert not keys & {"A", "k", "e", "R"}
    assert {"y", "n"} <= keys
    screen = tc.render_request(req, theme=rr.Theme(enabled=False), width=100)
    assert screen.index("This cannot be undone") < screen.index("$ rm big.dat")
    assert "Deletes big.dat (2.0 MiB)" in screen
    assert "Enter" in screen


class _Keys:
    active = True

    def __init__(self, chunks):
        self._chunks = list(chunks)

    def read_ready(self, timeout):
        return self._chunks.pop(0) if self._chunks else ""


class _Engine:
    token_usage = {"input": 0, "output": 0}

    def get_status(self):
        return {}

    def request_stop(self):
        pass


def test_enter_answers_no_to_what_cannot_be_undone(ws, tmp_path, monkeypatch):
    from delfin.agent import repl
    from delfin.agent import repl_render as rr
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
    broker = tc.TerminalConfirmBroker(timeout_s=5.0)
    agent = repl.TerminalAgent(_Engine(), out=io.StringIO(),
                               err=io.StringIO(), broker=broker)
    agent.transcript.theme = rr.Theme(enabled=False)
    req = _irreversible_request(ws)
    t = threading.Thread(target=agent._answer_request,
                         args=(req, _Keys(["\r"])), daemon=True)
    t.start()
    t.join(5)
    if t.is_alive():
        broker.resolve(req, True)
        t.join(5)
    assert req.resolved and req.decision is False


# ---------------------------------------------------------------------------
# The dashboard
# ---------------------------------------------------------------------------

def test_the_dashboard_row_says_it_first_and_has_no_permanent_button(ws, tmp_path, monkeypatch):
    import ipywidgets as widgets
    from delfin.agent.kit_confirm import KitConfirmBroker, _ConfirmRequest
    broker = KitConfirmBroker(default_timeout_s=1.0,
                              persist_callback=lambda kind, pat: (True, pat))
    req = _irreversible_request(ws)
    row = broker._build_request_row(
        _ConfirmRequest(seq=1, tool_name="bash", args=dict(req.args),
                        preview=req.preview), widgets)

    def walk(box):
        for child in getattr(box, "children", ()):
            yield child
            yield from walk(child)

    html = "".join(c.value for c in walk(row) if isinstance(c, widgets.HTML))
    assert html.index("This cannot be undone") < html.index("$ rm big.dat")
    labels = [c.description for c in walk(row) if isinstance(c, widgets.Button)]
    assert not any("Permanent" in label for label in labels)


# ---------------------------------------------------------------------------
# The service and the phone
# ---------------------------------------------------------------------------

def test_the_service_marks_it_and_offers_no_lasting_grant(ws):
    from delfin.agent.service.approvals import ProtocolBroker
    broker = ProtocolBroker.__new__(ProtocolBroker)
    broker.suggest_pattern = lambda command: r"^rm\b"
    broker.request_id = lambda req: "r1"
    fields = broker.describe(_irreversible_request(ws))
    assert fields["flags"]["irreversible"] is True
    assert not {o["key"] for o in fields["options"]} & {"A", "k", "e"}


def _card(ws, monkeypatch, lang):
    from delfin.messenger import render
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", lang)
    req = _irreversible_request(ws)
    return render.approval_card({
        "kind": "confirm", "tool": "bash", "action": "shell",
        "text": req.preview, "phoneApproval": True,
        "options": [{"key": "y", "label": "yes"}, {"key": "n", "label": "no"}]})


def test_the_phone_card_says_it_first_in_english(ws, monkeypatch):
    card = _card(ws, monkeypatch, "en")
    assert card.index("This cannot be undone") < card.index("$ rm big.dat")
    assert card.index("Deletes big.dat (2.0 MiB)") < card.index("Reply")


def test_the_phone_card_says_it_first_in_german(ws, monkeypatch):
    from delfin.messenger import lang
    card = _card(ws, monkeypatch, "de")
    head = lang.DE[IR.HEAD_IRREVERSIBLE]
    deletes = lang.DE[IR.DELETE].split("{path}")[0]
    assert card.index(head) < card.index("$ rm big.dat")
    assert deletes and card.index(deletes) < card.index("$ rm big.dat")


def test_every_sentence_has_a_german_text_with_the_same_placeholders():
    import string
    from delfin.messenger import lang

    def fields(s):
        return sorted(f for _, f, _, _ in string.Formatter().parse(s) if f)

    missing = [t[:50] for t in IR.TEMPLATES if t not in lang.DE]
    assert not missing
    for t in IR.TEMPLATES:
        assert fields(lang.DE[t]) == fields(t), t[:50]


# ---------------------------------------------------------------------------
# The agent says it too
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("role_id,mode_id", [
    ("solo_agent", "solo"), ("dashboard_agent", "dashboard"),
    ("office_agent", "office"), ("builder_agent", "quick"),
    ("critic_agent", "quick"), ("session_manager", "quick"),
])
def test_the_rule_is_in_the_prompt_of_every_role_that_can_ask(role_id, mode_id):
    from delfin.agent.prompt_loader import PromptLoader
    prompt = " ".join(PromptLoader().build_system_prompt(
        role_id=role_id, mode_id=mode_id, task_text="tidy the workspace").split())
    assert "cannot be undone" in prompt
    assert "say what will be lost or sent" in prompt
