"""The agent improvised around a missing prerequisite. Now it proposes one.

A field report: with no pytest in the interpreter, the agent built a
runner of its own -- a venv, or a wrapper script in the home directory.
It could already DETECT the gap; what it did with the detection was work
around it, because "no report file produced" reads like an obstacle
rather than a fact about the environment.

The loop is now: the check detects, the check DECLARES its remedy, the
user approves that remedy verbatim, and only then does anything run.

Most of what is asserted here is what the mechanism REFUSES, because it
is a mechanism that executes:

* the action comes from the repository, never from model text -- a check
  states its own command or setting, and nothing is parsed out of prose;
* approval is verbatim: ``apply_proposal`` compares what it is given
  against the action it would perform, so approving one thing cannot run
  another;
* advice is never "applied". Most prerequisites have no runnable
  remedy -- installing a system package, logging a credential helper in,
  reordering a library path -- and saying so is the answer rather than a
  gap;
* a settings change is offered only where it can be honoured: proposing
  bwrap isolation on a host with no bwrap would turn a warning into a
  refusal of every shell command.

Universal: every case is driven on rows built in the test, so no host
decides the outcome.
"""

from __future__ import annotations

import types

import pytest

from delfin.agent import doctor as D
from delfin.agent import prerequisites as P


def _row(**kw):
    base = {"check": "a check", "status": D.WARN, "detail": "d", "fix": "f"}
    base.update(kw)
    return base


# ---------------------------------------------------------------------------
# What becomes a proposal
# ---------------------------------------------------------------------------

def test_a_passing_check_is_not_a_proposal():
    rows = [_row(status=D.PASS), _row(check="bad", status=D.WARN)]
    assert [p.check for p in P.proposals(rows=rows)] == ["bad"]


def test_the_id_is_derived_from_the_name_so_it_can_be_typed_back():
    p = P.proposals(rows=[_row(check="documents: OCR")])[0]
    assert p.pid == "documents-ocr"


def test_two_checks_that_slugify_alike_stay_reachable():
    rows = [_row(check="a b", detail="one"), _row(check="a-b", detail="two")]
    pids = [p.pid for p in P.proposals(rows=rows)]
    assert len(set(pids)) == 2, pids


def test_the_doctor_order_is_kept():
    rows = [_row(check="z"), _row(check="a"), _row(check="m")]
    assert [p.check for p in P.proposals(rows=rows)] == ["z", "a", "m"]


def test_a_unique_prefix_finds_one():
    rows = [_row(check="bash isolation"), _row(check="documents: OCR")]
    assert P.find("bash", rows=rows).check == "bash isolation"
    assert P.find("nope", rows=rows) is None


def test_an_ambiguous_prefix_finds_none():
    rows = [_row(check="bash one"), _row(check="bash two")]
    assert P.find("bash", rows=rows) is None


def test_a_malformed_report_does_not_raise():
    assert P.proposals(rows=[None, 5, {}, _row(check="ok")]) [0].check == "ok"


# ---------------------------------------------------------------------------
# Applicable versus advice
# ---------------------------------------------------------------------------

def test_a_row_with_no_remedy_is_advice():
    p = P.proposals(rows=[_row()])[0]
    assert p.kind == P.ADVICE
    assert p.action == ""
    assert "yours" in P.render(p)


def test_a_declared_command_makes_it_applicable():
    p = P.proposals(rows=[_row(command="pip install x")])[0]
    assert p.kind == P.APPLICABLE
    assert p.action == "run: pip install x"
    assert f"/fix {p.pid} run" in P.render(p)


def test_a_declared_setting_makes_it_applicable():
    p = P.proposals(rows=[_row(setting=("agent.x", "y"))])[0]
    assert p.kind == P.APPLICABLE
    assert p.action == "set agent.x = 'y'"


def test_advice_is_never_applied():
    p = P.proposals(rows=[_row()])[0]
    out = P.apply_proposal(p, p.action, run_command=_forbidden)
    assert out["applied"] is False
    assert "advice" in out["refused"]


# ---------------------------------------------------------------------------
# The approval is verbatim
# ---------------------------------------------------------------------------

def _forbidden(*a, **k):
    raise AssertionError("nothing may run here")


def test_a_mismatched_approval_runs_nothing():
    p = P.proposals(rows=[_row(command="pip install x")])[0]
    for wrong in ("", "yes", "run: pip install y", "pip install x"):
        out = P.apply_proposal(p, wrong, run_command=_forbidden)
        assert out["applied"] is False, wrong
        assert "does not repeat the action" in out["refused"], wrong


def test_the_action_shown_is_the_action_performed():
    ran: list[str] = []
    p = P.proposals(rows=[_row(command="pip install x")])[0]
    out = P.apply_proposal(
        p, p.action,
        run_command=lambda c, w: (ran.append(c), {"ok": True, "output": ""})[1])
    assert out["applied"] is True
    assert ran == ["pip install x"], (
        "the command run must be the one the action named")


def test_a_failing_command_is_reported_not_swallowed():
    p = P.proposals(rows=[_row(command="false")])[0]
    out = P.apply_proposal(
        p, p.action, run_command=lambda c, w: {"ok": False, "output": "nope"})
    assert out["applied"] is False
    assert out["refused"]
    assert "nope" in out["output"]


def test_a_command_that_cannot_start_is_a_sentence_not_a_traceback():
    p = P.proposals(rows=[_row(command="x")])[0]
    out = P.apply_proposal(p, p.action, run_command=_forbidden)
    assert out["applied"] is False
    assert "could not be started" in out["refused"]


def test_a_setting_is_written_through_one_dotted_key():
    seen: list[tuple] = []
    p = P.proposals(rows=[_row(setting=("agent.deep.key", 3))])[0]
    out = P.apply_proposal(
        p, p.action, save_setting=lambda k, v: seen.append((k, v)))
    assert out["applied"] is True
    assert seen == [("agent.deep.key", 3)]


def test_a_setting_that_cannot_be_written_is_reported():
    p = P.proposals(rows=[_row(setting=("agent.x", 1))])[0]
    out = P.apply_proposal(
        p, p.action, save_setting=_forbidden)
    assert out["applied"] is False
    assert "could not be written" in out["refused"]


# ---------------------------------------------------------------------------
# The action comes from the repository, not from text
# ---------------------------------------------------------------------------

def test_a_remedy_is_not_parsed_out_of_the_prose():
    """The prose says a command; the row declares none. Nothing runnable
    may be inferred from a sentence -- a wording change would move it."""
    p = P.proposals(rows=[_row(fix="run: pip install something")])[0]
    assert p.kind == P.ADVICE
    assert p.command == ""


def test_normalisation_refuses_a_command_that_is_not_one_line():
    """_normalise coerces whatever a check returned, including a
    monkeypatched one, and what it returns is offered for approval."""
    row = D._normalise({"check": "c", "status": D.WARN,
                        "command": "a\nrm -rf /"}, "g")
    assert row.get("command", "") == ""


def test_normalisation_refuses_a_malformed_setting():
    for bad in ("x", ["only-one"], ["", 1], [1, 2], {"a": 1}):
        row = D._normalise({"check": "c", "status": D.WARN, "setting": bad},
                           "g")
        assert "setting" not in row, bad


def test_normalisation_keeps_a_well_formed_pair():
    row = D._normalise({"check": "c", "status": D.WARN,
                        "setting": ("agent.a", False)}, "g")
    assert row["setting"] == ["agent.a", False]


# ---------------------------------------------------------------------------
# The two checks that declare a remedy
# ---------------------------------------------------------------------------

def test_the_test_runner_check_declares_the_install(monkeypatch):
    import importlib.util as iu

    monkeypatch.setattr(iu, "find_spec",
                        lambda name, *a, **k: None if name == "pytest" else 1)
    row = D._check_test_runner({})[0]
    assert row["status"] != D.PASS
    assert "pip install" in row.get("command", "")
    assert "delfin-complat[test]" in row["command"]
    assert "must say so rather than build a runner of its own" in row["fix"]


def test_the_test_runner_check_passes_when_pytest_is_there():
    row = D._check_test_runner({})[0]
    assert row["status"] == D.PASS
    assert "command" not in row


def test_isolation_offers_bwrap_only_where_something_can_hold_it():
    """Proposing bwrap with no bwrap turns a warning into a refusal of
    every shell command -- a worse state than the one being fixed."""
    import inspect

    body = inspect.getsource(D._check_bash_isolation)
    i = body.index('setting=("agent.bash_isolation", "bwrap")')
    assert "if held_by else None" in body[i:i + 120]


def test_the_registry_includes_the_new_check():
    assert any(a == "_check_test_runner" for _n, a in D._CHECK_ATTRS)


# ---------------------------------------------------------------------------
# The command surface
# ---------------------------------------------------------------------------

def _ctx():
    return types.SimpleNamespace(workspace=".")


def test_fix_lists_without_acting(monkeypatch):
    from delfin.agent import repl_commands as RC

    monkeypatch.setattr(P, "apply_proposal", _forbidden)
    out = RC._fix(_ctx(), "").output
    assert out


def test_fix_with_an_id_shows_but_does_not_act(monkeypatch):
    from delfin.agent import repl_commands as RC

    real = P.proposals
    monkeypatch.setattr(
        P, "proposals",
        lambda *a, **k: real(rows=[_row(check="c", command="pip install x")]))
    monkeypatch.setattr(P, "apply_proposal", _forbidden)
    out = RC._fix(_ctx(), "c").output
    assert "run: pip install x" in out
    assert "/fix c run" in out


def test_only_the_word_run_acts(monkeypatch):
    from delfin.agent import repl_commands as RC

    calls: list[str] = []
    real = P.proposals
    monkeypatch.setattr(
        P, "proposals",
        lambda *a, **k: real(rows=[_row(check="c", command="pip install x")]))
    monkeypatch.setattr(
        P, "apply_proposal",
        lambda prop, approved, **k: calls.append(approved) or
        {"applied": True, "action": approved, "refused": ""})
    for word in ("", "yes", "go", "--force"):
        RC._fix(_ctx(), f"c {word}".strip())
    assert calls == [], f"acted on {word!r}"
    RC._fix(_ctx(), "c run")
    assert calls == ["run: pip install x"]


def test_an_unknown_id_says_where_the_list_is(monkeypatch):
    from delfin.agent import repl_commands as RC

    monkeypatch.setattr(P, "proposals", lambda *a, **k: [])
    out = RC._fix(_ctx(), "nope run").output
    assert "/fix" in out


def test_fix_never_raises(monkeypatch):
    from delfin.agent import repl_commands as RC

    monkeypatch.setattr(P, "proposals", _forbidden)
    assert "fix failed" in RC._fix(_ctx(), "").output


# ---------------------------------------------------------------------------
# And the agent is told not to improvise
# ---------------------------------------------------------------------------

def test_the_shared_addendum_forbids_routing_around_a_prerequisite():
    from pathlib import Path

    text = (Path(D.__file__).resolve().parent / "pack" / "shared"
            / "refusal_addendum.md").read_text(encoding="utf-8")
    assert "missing" in text and "prerequisite" in text
    # The behavioural contract, not the command names: those live in
    # /doctor's own output and the command help, which is where a reader
    # looks for them, and every token here is paid on every turn.
    for forbidden in ("second interpreter", "wrapper",
                      "install into a home"):
        assert forbidden in text, forbidden
    assert "/doctor" in text
    assert "tell the user what" in text
