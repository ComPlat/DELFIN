"""The land gate refuses a branch that carries no proof, and fixes the signers.

Three things a landed branch was taken on trust for:

1. PROOF. ``land.check`` read ``scratch/PR.md`` for a home path and for the
   word WIP, and for nothing else. A package whose test had never been shown
   red on the base, and a package with no real run at all, were both landed on
   the author's word. The gate now RUNS the test the package's own decision
   table names -- on the base (it must fail there) and on the branch (it must
   pass) -- and reads the real-run row the stage recorded outside the member's
   cage.

2. SIGNERS. Only the AUTHOR was compared with the configured identity. A
   committer is a second signature on the same commit and was unchecked, and
   ``Co-Authored-By:`` lines -- which name an account the forge then credits as
   a contributor -- were read by nobody at all.

3. THE SQUASH MESSAGE. The message a squash merge publishes on the base now
   carries exactly one ``Co-Authored-By: DELFIN-Agent <agent email>`` and no
   model, provider or site name.

Every test runs the SAME entry point the product runs -- ``land.main``, which
``delfin-agent swarm wave land`` calls (cli_swarm.py:227-233) -- so none of
them measures a helper that only exists here. Each states the input that would
make it fail. Two of them are CONTROLS that pass on the base and on the
branch: they are what shows the refusals are not a blanket.
"""

from __future__ import annotations

import json
import os
import socket
import subprocess
import sys

import pytest

from delfin.agent.swarm import land as pr_bot
from delfin.agent.swarm import wave

#: The two identities of the fixture: the person the commits are signed as,
#: and the agent account the trailer must name.
PERSON = "person@example.org"
AGENT = "424242+delfin-agent@users.noreply.github.com"
#: The same account NAME with a different NUMBER: the address a forge would
#: attach to whoever owns that number, which is not the agent.
AGENT_OTHER_NUMBER = "999999+delfin-agent@users.noreply.github.com"

RUN, TID = "w54t", "g1a1"

#: A test that passes wherever it exists: the base's own, so "green on the
#: base" is a fact and not an accident of a missing file.
_OLD_TEST = "def test_already_passing():\n    assert True\n"
_NEW_TEST = "def test_the_change():\n    assert 1 + 1 == 2\n"


def _git(*args, cwd):
    r = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


@pytest.fixture
def swarm(tmp_path, monkeypatch):
    """A real git worktree with a base and a branch commit, a run directory,
    and every path redirected.

    The two trees are real because the gate now runs a test on each of them.
    ``review_gate`` is stubbed to a pass for EVERY test: without it the gate's
    own review refuses the branch on the base as well, so a refusal test would
    be green on the base for the review's reason and would have measured
    nothing about proof.
    """
    wt = tmp_path / "wt"
    (wt / "scratch").mkdir(parents=True)
    (wt / "tests").mkdir()
    env = {"GIT_AUTHOR_NAME": "person", "GIT_AUTHOR_EMAIL": PERSON,
           "GIT_COMMITTER_NAME": "person", "GIT_COMMITTER_EMAIL": PERSON,
           "GIT_CONFIG_GLOBAL": str(tmp_path / "gitconfig"), "HOME": str(tmp_path),
           "PATH": os.environ["PATH"]}
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    _git("init", "-q", "-b", "main", cwd=wt)
    (wt / "scratch" / ".gitignore").write_text("*\n")
    (wt / "tests" / "test_old.py").write_text(_OLD_TEST)
    _git("add", "-A", cwd=wt)
    _git("commit", "-qm", "base", cwd=wt)
    _git("update-ref", "refs/remotes/origin/main", "HEAD", cwd=wt)
    (wt / "tests" / "test_new.py").write_text(_NEW_TEST)
    _git("add", "-A", cwd=wt)
    _git("commit", "-qm", "add the change", cwd=wt)

    run = tmp_path / RUN
    run.mkdir()
    (run / "state.json").write_text(json.dumps({"agents": {TID: {
        "worktree": str(wt), "branch": f"wave/{RUN}/{TID}",
        "task": {"id": TID, "base": "origin/main", "model": "alias-code"}}}}))

    monkeypatch.setattr(pr_bot, "ROOT", tmp_path)
    monkeypatch.setattr(pr_bot, "identity_email", lambda: PERSON)
    monkeypatch.setattr(pr_bot, "PYCI", sys.executable)
    if hasattr(pr_bot, "coauthor_email"):       # absent on the base
        monkeypatch.setattr(pr_bot, "coauthor_email", lambda: AGENT)
    monkeypatch.setattr(pr_bot, "review_gate",
                        lambda wt, run, tid, agent: (True, "stubbed pass", {"stub": True}))
    monkeypatch.setattr(wave, "root", lambda: tmp_path)
    return {"root": tmp_path, "wt": wt, "run": RUN, "tid": TID}


def _record_real_run(swarm, probe="inline-xyz", *, side="branch", error="", tools=3,
                     touched=True):
    """The rows the stage records outside the member's cage, through the
    product's own writer, plus the ``<tid>.realrun`` marker.

    ``realrun._one`` writes the probe's log next to those rows, so a faithful
    fixture writes one too: ``touched`` says whether the log shows a call on a
    file the branch changed, which is the difference between a run that
    reached the change and one that answered in chat.
    """
    rows = [{"probe": probe, "side": "base", "stop": "end_turn", "tools": tools, "error": ""},
            {"probe": probe, "side": side, "stop": "end_turn", "tools": tools, "error": error}]
    wave.record(swarm["run"], "realrun", swarm["tid"], rows=rows)
    (swarm["root"] / swarm["run"] / f"{swarm['tid']}.realrun").write_text("0")
    if error:
        return                              # a run that failed wrote no answer
    d = swarm["wt"] / "scratch" / "realrun"
    d.mkdir(parents=True, exist_ok=True)
    ev = ({"type": "tool_use", "name": "edit_file",
           "input": {"path": "tests/test_new.py", "content": "x"}} if touched
          else {"type": "text", "text": "I read the code and answered the question."})
    (d / f"{probe}-{side}.jsonl").write_text(json.dumps(ev) + "\n")


def _proof_table(*rows) -> str:
    head = "| package | test | realrun | decision |\n|---|---|---|---|\n"
    return head + "".join(f"| {p} | `{t}` | {r} | {d} |\n" for p, t, r, d in rows)


def _write_pr(swarm, body: str) -> None:
    (swarm["wt"] / "scratch" / "PR.md").write_text("# a title\n\n" + body + "\n")


def _gate(swarm):
    """The product entry point: what ``delfin-agent swarm wave land`` calls."""
    return pr_bot.main(["--dry-run", swarm["run"]])


def _landing(swarm) -> str:
    return (swarm["wt"] / "scratch" / "LANDING.txt").read_text()


def _commit_as(swarm, env_extra: dict, message: str):
    """One more commit on the branch, signed as ``env_extra`` says."""
    (swarm["wt"] / "extra.txt").write_text("x\n")
    _git("add", "-A", cwd=swarm["wt"])
    r = subprocess.run(["git", "commit", "-qm", message], cwd=swarm["wt"],
                       capture_output=True, text=True, env={**os.environ, **env_extra})
    assert r.returncode == 0, r.stderr


def _ready_proof(swarm):
    """Everything the gate needs from the branch itself: a table naming the
    test that is red on the base and green here, and a real run that reached it."""
    _record_real_run(swarm)
    _write_pr(swarm, _proof_table(("g1", "tests/test_new.py", "inline-xyz", "done")))


# --- 0. The control: a branch that meets every rule passes on BOTH trees ---


def test_a_branch_with_red_on_the_base_a_green_test_and_a_real_run_passes(swarm):
    """CONTROL. Passes on the base (which checks none of it) and on the branch
    (which checks all of it). If this ever goes red on the branch, the new
    refusals are a blanket rather than a gate."""
    _ready_proof(swarm)
    _gate(swarm)
    assert _landing(swarm).startswith("ready:"), _landing(swarm)


# --- 1. PROOF: the gate runs the test, it does not trust the text ----------


def test_a_branch_with_no_proof_table_is_refused_and_says_so(swarm):
    """The base accepts a body with no proof at all and prints it as ready, so
    this is red on the base for the reason it names."""
    _write_pr(swarm, "Just a body, no table. Nothing here was measured.")
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "no proof table" in _landing(swarm), _landing(swarm)


def test_a_package_row_without_a_test_is_refused_and_named(swarm):
    """A row with an empty test cell names no evidence. The refusal names the
    package, so the member knows which row to fix."""
    _record_real_run(swarm)
    _write_pr(swarm, _proof_table(("g1", "", "inline-xyz", "done")))
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "g1" in _landing(swarm) and "no test" in _landing(swarm), _landing(swarm)


def test_a_test_that_is_green_on_the_base_proves_nothing_and_is_refused(swarm):
    """THE falsifying case for this check. ``tests/test_old.py`` exists on the
    base and passes there, so naming it measures nothing -- exactly the arm
    that reports "green" for a test that never distinguished the two trees.
    A text matcher would accept it; the gate has to RUN it to refuse."""
    _record_real_run(swarm)
    _write_pr(swarm, _proof_table(("g1", "tests/test_old.py", "inline-xyz", "done")))
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "g1" in _landing(swarm) and "origin/main" in _landing(swarm), _landing(swarm)


def test_proof_without_a_real_run_is_refused_and_named(swarm):
    """A green test is not a run: the change was never used the way a person
    would use it. Red on the base, which never reads the real-run record."""
    _write_pr(swarm, _proof_table(("g1", "tests/test_new.py", "inline-xyz", "done")))
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "g1" in _landing(swarm) and "real run" in _landing(swarm), _landing(swarm)


def test_a_real_run_whose_probe_never_reached_the_change_is_refused(swarm):
    """A probe that started and errored recorded a row and reached nothing
    (the import check refused it -- six real runs were refused that way on
    2026-11-24). A row is not evidence; the row has to be clean."""
    _record_real_run(swarm, error="the run imported DELFIN from elsewhere")
    _write_pr(swarm, _proof_table(("g1", "tests/test_new.py", "inline-xyz", "done")))
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "g1" in _landing(swarm) and "reached nothing" in _landing(swarm), _landing(swarm)


# --- 2. SIGNERS: both signatures, and the account the trailer names --------


def test_a_commit_whose_committer_is_someone_else_is_refused(swarm):
    """Red on the base, which compares the AUTHOR only: a commit authored as
    the configured identity and committed as somebody else passed."""
    _ready_proof(swarm)
    _commit_as(swarm, {"GIT_COMMITTER_NAME": "other",
                       "GIT_COMMITTER_EMAIL": "other@example.org"},
               "a commit by another committer")
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "committer" in _landing(swarm), _landing(swarm)


def test_a_coauthor_line_naming_a_stranger_is_refused(swarm):
    """Red on the base, which reads no Co-Authored-By line at all. That line
    attaches an account to the commit in the forge's history."""
    _ready_proof(swarm)
    _commit_as(swarm, {}, "a commit\n\nCo-Authored-By: Someone <someone@example.org>\n")
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "someone@example.org" in _landing(swarm), _landing(swarm)


def test_the_agent_account_with_a_different_number_is_refused(swarm):
    """The NUMBER is what names the account, not the spelling. A line written
    for the agent's name with another account's number renders in the forge as
    that other account's contribution, so it is refused and the number named."""
    _ready_proof(swarm)
    _commit_as(swarm, {}, f"a commit\n\nCo-Authored-By: DELFIN-Agent <{AGENT_OTHER_NUMBER}>\n")
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "account number 999999" in _landing(swarm), _landing(swarm)


def test_a_model_suffixed_trailer_is_the_same_account_and_its_model_stays_out(swarm, capsys):
    """``DELFIN-Agent (model)`` is what the executor writes (prompt_loader.py:970),
    so it must be ACCEPTED as the agent account -- the account is its address,
    and a name suffix does not change which account that is. Refusing it would
    refuse every branch in a deployment that puts the model in the trailer.

    What the package forbids is the model NAME reaching the published message,
    and that is the second half: the squash message carries no model word.
    """
    _ready_proof(swarm)
    _commit_as(swarm, {}, f"a commit\n\nCo-Authored-By: DELFIN-Agent (a-model) <{AGENT}>\n")
    _gate(swarm)
    assert _landing(swarm).startswith("ready:"), _landing(swarm)
    out = capsys.readouterr().out
    assert "a-model" not in out.split("squash message")[-1], out
    assert out.count("Co-Authored-By:") == 1, out


def test_the_identity_with_one_agent_trailer_passes(swarm):
    """CONTROL for the signers half: author and committer are the identity and
    the trailer is the agent account. Green on the base and on the branch."""
    _ready_proof(swarm)
    _commit_as(swarm, {}, f"a commit\n\nCo-Authored-By: DELFIN-Agent <{AGENT}>\n")
    _gate(swarm)
    assert _landing(swarm).startswith("ready:"), _landing(swarm)


# --- 3. THE SQUASH MESSAGE the lander proposes -----------------------------


def test_the_proposed_squash_message_carries_one_agent_line(swarm, capsys):
    """The message published on the base names the agent once, at the
    configured address. Red on the base, which composes no message at all."""
    _ready_proof(swarm)
    _gate(swarm)
    out = capsys.readouterr().out
    assert f"Co-Authored-By: DELFIN-Agent <{AGENT}>" in out, out
    assert out.count("Co-Authored-By:") == 1, out
    assert "squash message" in out, out


def test_the_squash_message_drops_the_trailers_the_commits_carry(swarm, capsys):
    """A branch whose commits carry their own trailer publishes ONE line, the
    one the gate writes -- not two contributors of one branch. Red on the base,
    which composes no message."""
    _ready_proof(swarm)
    _commit_as(swarm, {}, f"a commit\n\nCo-Authored-By: DELFIN-Agent <{AGENT}>\n")
    _gate(swarm)
    out = capsys.readouterr().out
    assert out.count("Co-Authored-By:") == 1, out


def test_a_squash_message_naming_the_model_is_refused(swarm):
    """A commit subject naming the model reaches the published message. Red on
    the base, which never looks at the text of a message."""
    _ready_proof(swarm)
    _commit_as(swarm, {}, "kit.glm-5.3-flash fixed the gate")
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "model, provider or site" in _landing(swarm), _landing(swarm)


def test_a_squash_message_naming_the_site_is_refused(swarm, monkeypatch):
    """A site name in the published message names an installation, which rule
    16 forbids. The host is pinned here, so this measures the gate and not the
    machine it happens to run on (the base composes no message at all)."""
    monkeypatch.setattr(socket, "gethostname", lambda: "compile-node-7.example")
    _ready_proof(swarm)
    _commit_as(swarm, {}, "runs on compile-node-7 only")
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "names this site" in _landing(swarm), _landing(swarm)


# --- 4. The reach test the first version only named --------------------------
# Measured on this branch with scratch/measure_arms.py: a stage real-run row
# whose probe log never touched a file the branch changed was ACCEPTED, and a
# branch whose decision table named a row but whose probe answered in chat
# alone passed. "The run happened" is not "the run reached the change", and the
# package's first version measured only the first.


def _probe_log(swarm, name="inline-xyz", *, touched=True, side="branch"):
    """The probe's own log, in the shape ``--output-format stream-json`` writes
    (delfin/agent/cli_stream.py:20): one JSON object per line, a tool call
    carrying ``name`` and ``input``."""
    d = swarm["wt"] / "scratch" / "realrun"
    d.mkdir(parents=True, exist_ok=True)
    if touched:
        ev = {"type": "tool_use", "name": "edit_file",
              "input": {"path": "tests/test_new.py", "content": "x"}}
    else:
        ev = {"type": "text", "text": "I read the code and answered the question."}
    (d / f"{name}-{side}.jsonl").write_text(json.dumps(ev) + "\n")


def test_a_real_run_that_never_touched_the_change_is_refused(swarm):
    """The stage recorded a real run on the branch, so "the run happened" --
    but its log shows a chat and no file the branch changed. Red before this
    test: the gate accepted it (measured, six-case probe, case B)."""
    _ready_proof(swarm)
    _probe_log(swarm, touched=False)
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "reached" in _landing(swarm), _landing(swarm)


def test_a_real_run_that_touched_a_changed_file_is_accepted(swarm):
    """CONTROL for the reach test, on the same commit: the same stage row with
    a log that names a file the branch changed passes. Without this the
    refusal above could be a blanket."""
    _ready_proof(swarm)
    _probe_log(swarm, touched=True)
    _gate(swarm)
    assert _landing(swarm).startswith("ready:"), _landing(swarm)


# --- 5. The squash message is the PR BODY, which is what is published --------
# Measured on this branch: a model name in the PR body was ACCEPTED, because
# the message was composed from the title and the commit subjects only. The
# body is most of what a squash merge publishes, so a model, provider or site
# name in it reaches the base exactly like one in the title.


def test_a_model_name_in_the_pr_body_is_refused(swarm):
    """Red before this test: the body was never scanned (six-case probe, D)."""
    _ready_proof(swarm)
    _write_pr(swarm, _proof_table(("g1", "tests/test_new.py", "inline-xyz", "done"))
              + "\nthis bundle was built with kit.glm-5.3-flash\n")
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "model, provider or site" in _landing(swarm), _landing(swarm)


def test_a_provider_name_in_the_pr_body_is_refused(swarm):
    """Same rule, the provider word: a name no model family maps to."""
    _ready_proof(swarm)
    _write_pr(swarm, _proof_table(("g1", "tests/test_new.py", "inline-xyz", "done"))
              + "\nserved through blablador\n")
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)


def test_a_body_without_such_a_name_is_accepted(swarm):
    """CONTROL: the same body shape with ordinary prose passes, so the two
    refusals above measure the name and not the body's presence."""
    _ready_proof(swarm)
    _write_pr(swarm, _proof_table(("g1", "tests/test_new.py", "inline-xyz", "done"))
              + "\nthe gate runs the test on both trees\n")
    _gate(swarm)
    assert _landing(swarm).startswith("ready:"), _landing(swarm)


# --- 6. The table the wave itself tells members to write ---------------------
# DELIVERY_RULES 17 says "one table row per item (item, decision, commit,
# reason)". The proof check looked for a column called "package" and refused
# anything else, so a branch written to the wave's own rule was refused for
# carrying no proof table (measured: scratch/probe_house_format.py).


_HOUSE_TABLE = ("| item | decision | test | real run | commit | reason |\n"
                "|---|---|---|---|---|---|\n"
                "| g1 | landed | tests/test_new.py | inline-xyz | abc1234 | measured |\n")


def test_the_table_rule_17_describes_is_accepted(swarm):
    """The wave's own four columns plus the two the package asks for, with the
    row word "item". Red before this fix: refused as "no proof table"."""
    _ready_proof(swarm)
    _write_pr(swarm, _HOUSE_TABLE)
    _gate(swarm)
    assert _landing(swarm).startswith("ready:"), _landing(swarm)


def test_a_table_with_no_test_column_is_still_refused(swarm):
    """CONTROL: naming the row "item" does not make proof optional. The same
    table without a test column is refused, and it names the row."""
    _ready_proof(swarm)
    _write_pr(swarm, "| item | decision | reason |\n|---|---|---|\n"
                     "| g1 | landed | it looks fine |\n")
    _gate(swarm)
    assert _landing(swarm).startswith("refused:"), _landing(swarm)
    assert "g1" in _landing(swarm), _landing(swarm)
