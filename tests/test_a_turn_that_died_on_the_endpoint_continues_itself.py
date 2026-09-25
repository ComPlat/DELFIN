"""A turn that died on a transient endpoint failure continues itself.

On 2026-09-25 the kit.glm-5.3 endpoint was down from 09:45 to 10:50
("Model 'glm-5.3' not found", "No deployments available", "Response
payload is not completed").  Every running terminal session ended its
turn with ``[error] ...`` and then stood at the prompt until the
operator nudged it with a message -- 31 times in one day.  The
in-turn retry (``api_client._is_transient_api_error``) covers only the
first ~3 attempts inside one turn; an outage longer than that still
ends the turn, and the standstill after it is what this file judges.

What is judged here, and the instrument:

  retry            a scripted turn whose engine raises a transient
                   endpoint error: after a visible pause line and a
                   growing wait the terminal sends the continuation
                   prompt ("The endpoint failed mid-turn; continue
                   exactly where you were.") as a NEW turn

  give-up          the pause budget (60 min) is not exceeded one
                   retry at a time: after the budget the session
                   stands at the prompt like today and the operator
                   hears it over session_message

  any key          a key during the wait cancels it: no further
                   turn is started

  classification   a 400 context-length error, an auth error and a
                   "model not found" from a model that never answered
                   in this session are NOT retried; a "model not
                   found" after the same model already answered IS
                   (that is the outage shape, not a typo)

  timing           the wait delays are read from an injectable
                   clock/sleep so no test sleeps real seconds
"""

import io

import pytest


# -----------------------------------------------------------------------
# The scripted engine
# -----------------------------------------------------------------------

class _EndpointEngine:
    """Turns scripted as ("answer"|"raise", exception_or_None).

    ``stream_response`` either returns an answer (a normal end) or
    raises, which is the path the real engine takes out of a turn
    when the endpoint fails mid-turn.  ``first_model_answers``
    decides whether the model answered at least once BEFORE the
    error -- the fact a "model not found" retry is judged on.
    """

    session_id = "endpoint-check-0001"
    token_usage = {"input": 0, "output": 0}
    last_turn_stop_reason = "end_turn"

    def __init__(self, script, first_model_answers=True):
        self.script = list(script)
        self.messages = []
        self.turns = 0
        self.prompts = []
        self.first_model_answers = first_model_answers
        self.notified = []
        self.model = "glm-5.3"

    def stream_response(self, user_message="", max_tokens=0, **kw):
        self.turns += 1
        self.prompts.append(user_message)
        kind, payload = self.script.pop(0) if self.script \
            else ("answer", None)
        if self.turns > 1 or self.first_model_answers:
            # An answer from THIS model happened at least once.
            pass
        if kind == "answer":
            self.messages.append({"role": "user", "content": user_message})
            self.messages.append({"role": "assistant", "content": "done"})
            return "done"
        # The failure path: the engine classifies and re-raises; the
        # classification itself is the api_client function under judge,
        # imported by the code under test, so replay its behaviour.
        raise payload

    def _notify_operator_of_stall(self, what):
        self.notified.append(what)

    def get_status(self):
        return {}

    def export_state(self):
        return {"engine_messages": list(self.messages),
                "token_usage": dict(self.token_usage)}


class _Clock:
    """Injectable clock + sleep; records the waits it was asked for."""

    def __init__(self):
        self.now = 1000.0
        self.waits = []

    def time(self):
        return self.now

    def sleep(self, seconds):
        self.waits.append(seconds)
        self.now += seconds


class _Err(Exception):
    pass


def _transient(msg="Response payload is not completed"):
    e = _Err(msg)
    e.status_code = 503
    return e


def _context_error():
    e = _Err("This model supports a maximum context length of 32768 "
             "tokens. However, your messages resulted in 50000 tokens")
    e.status_code = 400
    return e


def _auth_error():
    e = _Err("Incorrect API key provided")
    e.status_code = 401
    return e


def _model_not_found():
    e = _Err("Model 'glm-5.3' not found")
    e.status_code = 404
    return e


def _run(script, tmp_path, first_model_answers=True, read_key=None,
         user_lines=("do the work", "go on", "and again")):
    """Drive TerminalAgent.run() over the scripted engine.

    ``read_key``: called with the pause length during each wait;
    return the byte chunk a keypress would deliver ("" = no key).
    Tests use it to press a key mid-wait.
    """
    from delfin.agent import repl as R
    engine = _EndpointEngine(script,
                             first_model_answers=first_model_answers)
    out, err = io.StringIO(), io.StringIO()
    out.isatty = lambda: False
    err.isatty = lambda: False
    lines = list(user_lines)

    def _read_line(_prompt=""):
        # The script IS the session: when it is spent, the user stops
        # typing (the loop leaves on the next interrupt).
        if lines and engine.script:
            return lines.pop(0)
        raise KeyboardInterrupt

    agent = R.TerminalAgent(
        engine, opts=R.ReplOptions(cwd=tmp_path, max_tokens=0),
        out=out, err=err, read_line=_read_line)
    agent._stdin = type("S", (), {"isatty": lambda self: False})()
    agent._idle_interrupts = 1           # one interrupt leaves (130)
    clock = _Clock()
    agent._endpoint_clock = clock
    agent._endpoint_sleep = clock.sleep
    engine.clock = clock                 # the schedule test reads the waits
    if read_key is not None:
        agent._endpoint_read_key = read_key
    try:
        agent.run()
    except KeyboardInterrupt:
        pass
    return engine, err


# -----------------------------------------------------------------------
# retry: error -> wait -> continuation
# -----------------------------------------------------------------------

def test_a_transient_error_waits_and_continues(tmp_path):
    engine, err = _run([("raise", _transient()), ("answer", None)],
                       tmp_path)
    assert engine.turns == 2, (
        "a turn that died on a transient endpoint error must continue "
        "itself as a new turn; the session ran "
        f"{engine.turns} turns instead of 2")
    assert "endpoint failed mid-turn; continue exactly where you were" \
        in engine.prompts[1], (
            "the continuation must tell the model the endpoint failed "
            "mid-turn and to continue exactly where it was")


def test_the_pause_is_visible(tmp_path):
    engine, err = _run([("raise", _transient()), ("answer", None)],
                       tmp_path)
    assert "endpoint unavailable" in err.getvalue(), (
        "the user must see a line saying the endpoint is unavailable "
        "and the session is retrying")


def test_the_wait_grows_60_120_then_300(tmp_path):
    """The waits grow 60, 120, 300 and then stay at 300 (see the
    schedule test below, which drives the same rule directly)."""
    script = [("raise", _transient()) for _ in range(4)] \
        + [("answer", None)]
    engine, err = _run(script, tmp_path)
    assert engine.turns == 5


def test_wait_schedule_60_120_300_300(tmp_path):
    """The waits grow 60, 120, 300 and then stay at 300."""
    script = [("raise", _transient()) for _ in range(4)] \
        + [("answer", None)]
    engine, err = _run(script, tmp_path)
    assert engine.clock.waits == [60.0, 120.0, 300.0, 300.0], (
        f"the wait must grow 60, 120, 300 and then stay at 300; it was "
        f"{engine.clock.waits}")


# -----------------------------------------------------------------------
# give-up: budget and operator note
# -----------------------------------------------------------------------

def test_the_budget_ends_the_retrying(tmp_path):
    """60 min of waiting is the ceiling: after it, no further turn."""
    script = [("raise", _transient()) for _ in range(20)]
    engine, err = _run(script, tmp_path,
                       user_lines=("do the work",))
    assert engine.turns <= 14, (
        f"with a 60-minute budget the session must stop retrying; it "
        f"ran {engine.turns} turns")
    assert engine.notified, (
        "when the budget is spent the operator must hear it over "
        "session_message -- the channel that reaches a person not "
        "watching six panes")


# -----------------------------------------------------------------------
# any key cancels
# -----------------------------------------------------------------------

def test_a_key_cancels_the_wait(tmp_path):
    engine, err = _run(
        [("raise", _transient())],
        tmp_path, read_key=lambda timeout: "x",
        user_lines=("do the work",))
    assert engine.turns == 1, (
        "a key during the wait must cancel it: no further turn, the "
        "session stands at the prompt like today")


# -----------------------------------------------------------------------
# classification
# -----------------------------------------------------------------------

def test_a_context_error_is_never_retried(tmp_path):
    engine, err = _run([("raise", _context_error())], tmp_path)
    assert engine.turns == 1, (
        "a 400 context-length error is deterministic: retrying it just "
        "fails again and burns the budget")


def test_an_auth_error_is_never_retried(tmp_path):
    engine, err = _run([("raise", _auth_error())], tmp_path)
    assert engine.turns == 1, (
        "a 401 is a credential problem, not an outage: no retry")


def test_model_not_found_without_prior_answer_is_not_retried(tmp_path):
    engine, err = _run([("raise", _model_not_found())], tmp_path,
                       first_model_answers=False)
    assert engine.turns == 1, (
        "a 'model not found' from a model that never answered in this "
        "session is a typo, not an outage: no retry")


def test_model_not_found_after_a_prior_answer_is_retried(tmp_path):
    engine, err = _run(
        [("answer", None), ("raise", _model_not_found()),
         ("answer", None)],
        tmp_path)
    assert engine.turns == 3, (
        "a 'model not found' after the same model already answered is "
        "the outage shape (2026-09-25): it must be retried")


def test_a_normal_answer_ends_the_retry_state(tmp_path):
    """An answered turn clears the pause budget so the next outage
    starts from 60 s again."""
    engine, err = _run(
        [("answer", None), ("raise", _transient()), ("answer", None),
         ("raise", _transient()), ("answer", None)],
        tmp_path)
    assert engine.turns == 5


# -----------------------------------------------------------------------
# review fixes (operator, 2026-09-25)
# -----------------------------------------------------------------------

def test_a_terminal_wait_is_not_followed_by_a_second_sleep(tmp_path):
    # On a terminal the key wait itself takes the whole pause; sleeping
    # again afterwards doubled every pause.
    engine, _ = _run([("raise", _transient()), ("answer", None)],
                     tmp_path, read_key=lambda timeout: "")
    assert engine.turns == 2
    assert engine.clock.waits == []


def test_ctrl_c_during_the_pause_cancels_it(tmp_path):
    def interrupted(timeout):
        raise KeyboardInterrupt

    engine, _ = _run([("raise", _transient()), ("answer", None)],
                     tmp_path, read_key=interrupted)
    # No automatic continuation: whatever runs next is the person's own
    # line, and run() did not leave through the interrupt.
    assert not any("endpoint failed mid-turn" in p
                   for p in engine.prompts[1:])
    assert engine.turns == 2      # the session went on to the next line
