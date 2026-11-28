"""Measure the background-job wake package on DELFIN's own entry points.

Package p2. The claim under test is a behaviour a user has: a command started in
the background reports back when it ends, so nobody has to poll it or chain
sleeps, and waiting for a condition is one bounded call.

HOW A VARIANT IS DRIVEN. Through the product's front door, ``cli.main([...])``:
``delfin-agent chat -p ... --output-format stream-json`` (the door a script or
supervisor drives, and the one the package names), ``delfin-agent run``, or
``delfin-agent chat`` -- with the real engine, the real OpenAI-compatible tool
loop and the real permission gate. Only the model's HTTP replies are scripted,
in the stream-chunk shape the provider sends (the same scripted behaviour
``tests/agent_terminal_probes/fake_model.py`` answers with, restated here
because these variants need a per-round script rather than a marker in the
typed prompt). Where a variant's whole answer is argparse's, the engine builder
is stopped at construction with a sentinel so nothing reaches the network; that
is still the front door, because the refusal happens before it.

THE ONE SUBSTITUTION, and why it is not a stand-in. A real ``bash_background``
job IS started through the real dispatch and the cage refuses it on this host
(exit 126, "sandbox could not be applied ([Errno 16] seccomp filter: Device or
resource busy)"), so no variant can wait on a job that takes measurable time.
Loosening the cage would be loosening a security check, so the scenario does
not touch it -- and ``notice_sandbox_refused`` measures that refusal instead.
For the timed variants the job enters the registry through the product's OWN
demote path, ``bash_jobs.get_registry().adopt(...)``: what a foreground command
that hit its timeout goes through. It shares ``_register`` with
``bash_background`` -- the same record fields, the same persistent registry
file, the same watchdog thread, the same session scoping, the same completion
event. What is substituted is how a live OS process joins the registry, never
what the registry or the wake path does with it. The control variant
``notice_sandbox_refused`` proves the real dispatch route is exercised too.

WHAT EACH VARIANT RECORDS.

  REACHED  the capability's own code ran, proven by an observable it produces
           (the wait loop names the bound it was given on stderr; the notice
           renderer emitted a line), never by a patched flag;
  USED     the wake happened where it should -- a second turn ran whose prompt
           carried the job id, or the wait loop said it was waiting;
  CORRECT  the outcome is the one the contract promises -- the right number of
           extra model calls, the bound honoured, the right line, the right
           refusal. ``expected()`` states that promise independently of what
           this tree does, so no variant can be made green by moving the goal.

Plus the cost the package asks to be measured: ``model_calls``, ``wake_turns``,
``wall_s``.

Each variant runs in its OWN process: ``bash_jobs`` fixes its index path from
``Path.home()`` at import, so one shared process would carry a variant's state
into the next.

Usage (from the repository root):
    PYTHONPATH=. python experiments/parity/p2-bg-job-wake/scenario.py \
        --all --label before \
        --out experiments/parity/p2-bg-job-wake/before.json
"""
from __future__ import annotations

import argparse
import io
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
#: Where the command was launched, snapshotted before any variant chdirs into
#: its scratch workspace. `Path.absolute()` does not normalise `..`, so a
#: relative --out is resolved against this instead.
LAUNCH_CWD = Path.cwd()

#: How long a "still running" job stays alive: well past one 30 s look, so a
#: variant that must catch it running really does.
LONG_LIVE_S = 60


class EngineRefused(RuntimeError):
    """Stops an engine from being built once argparse has answered."""


# ---------------------------------------------------------------------------
# the scripted model: stream chunks in the shape the provider sends
# ---------------------------------------------------------------------------

class _Fn:
    def __init__(self, name=None, arguments=None):
        self.name, self.arguments = name, arguments


class _TCDelta:
    def __init__(self, index, call_id=None, name=None, arguments=None):
        self.index, self.id = index, call_id
        self.function = _Fn(name, arguments)


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content, self.tool_calls = content, tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta, self.finish_reason = delta, finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices, self.usage = choices, usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


def call_stream(name: str, args: dict, call_id: str = "c1") -> _Stream:
    return _Stream([_Chunk([_Choice(_Delta(tool_calls=[_TCDelta(
        0, call_id=call_id, name=name, arguments=json.dumps(args))]),
        finish="tool_calls")], usage=_Usage())])


def answer_stream(text: str = "Started it; I will be told when it is done."
                  ) -> _Stream:
    return _Stream([_Chunk([_Choice(_Delta(content=text), finish="stop")],
                           usage=_Usage())])


# ---------------------------------------------------------------------------
# the scratch world one variant runs in
# ---------------------------------------------------------------------------

def _world(variant_id: str) -> dict:
    """A scratch HOME / workspace / state root for ONE variant."""
    root = HERE / ".runs" / variant_id
    if root.exists():
        shutil.rmtree(root)
    home, ws, scratch = root / "home", root / "ws", root / "scratch"
    (home / ".delfin").mkdir(parents=True)
    ws.mkdir()
    (scratch / ".delfin").mkdir(parents=True)
    os.environ["HOME"] = str(home)
    os.environ["DELFIN_SCRATCH_STATE"] = str(scratch)
    os.environ["NO_COLOR"] = "1"
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    from delfin.agent import state_paths
    state_paths.scratch_state_from_environment()
    return {"root": root, "home": home, "ws": ws, "scratch": scratch}


def _adopt(ws: Path, session_id: str, live_s: float, code: str = "",
           description: str = "the suite") -> str:
    """Register a real live process as a background job (see the header).

    ``start_new_session=True`` is what the product itself does
    (``bash_jobs.py:1432``, ``preexec_fn=os.setsid``). It matters here: the
    registry ends a job with ``os.killpg``, so a child left in the scenario's
    own process group would take the scenario down with it at exit.
    """
    import delfin.agent.bash_jobs as bj
    body = code or f"import time; time.sleep({live_s})"
    proc = subprocess.Popen([sys.executable, "-c", body],
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL,
                            start_new_session=True)
    job = bj.get_registry().adopt(f"python -c {body!r}", str(ws), proc,
                                 description=description, workspace=str(ws),
                                 session_id=session_id)
    return job.job_id


def _adopt_pid(ws: Path, session_id: str, live_s: float) -> tuple[str, int]:
    """As :func:`_adopt`, and the process id the variant must report.

    Whether ``--job-wait 0`` really stops the run's shells is decided by an
    atexit handler, which runs when the variant process ends -- so the
    variant cannot see it. It reports the pid; the parent looks for it after
    the variant exited.
    """
    import delfin.agent.bash_jobs as bj
    body = f"import time; time.sleep({live_s})"
    proc = subprocess.Popen([sys.executable, "-c", body],
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL,
                            start_new_session=True)
    job = bj.get_registry().adopt(f"python -c {body!r}", str(ws), proc,
                                  description="the suite",
                                  workspace=str(ws), session_id=session_id)
    return job.job_id, proc.pid


def _pid_alive(pid: int) -> bool:
    """Whether a process is still there -- the product's own judgement.

    Deliberately not re-implemented here: `bash_jobs._pid_alive` is what the
    registry itself uses to decide a job is alive (and it additionally rules
    out a recycled pid). Measuring "did the run stop its shell" with a
    different predicate than the one that stops it would ask a different
    question than the product answers.
    """
    from delfin.agent.bash_jobs import _pid_alive as _product_alive
    try:
        return bool(_product_alive(int(pid)))
    except Exception:
        return False


def _start_with_dispatch(ws: Path, session_id: str, command: str) -> str:
    """Start a job through the model's own tool dispatch (the real route)."""
    from delfin.agent.api_client import KitToolPermissions, _doc_executor
    perms = KitToolPermissions(workspace=ws, mode="acceptEdits",
                               task_session_id=session_id)
    out = json.loads(_doc_executor.execute(
        "bash_background", {"command": command, "description": "the suite"},
        perms))
    if "job_id" not in out:
        raise AssertionError(f"dispatch refused the job: {out}")
    return out["job_id"]


def _quiet_stdin(tty: bool = False) -> None:
    sys.stdin = io.StringIO("")
    sys.stdin.isatty = lambda: tty      # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# the front door, with the answer captured
# ---------------------------------------------------------------------------

def _drive_cli(argv: list) -> dict:
    """``cli.main(argv)`` with stderr captured and SystemExit recorded."""
    from delfin.agent import cli

    err = io.StringIO()
    real_err, sys.stderr = sys.stderr, err
    code, stopped = None, False
    try:
        code = cli.main(argv)
    except SystemExit as exc:
        code = exc.code
    except EngineRefused:
        # The door accepted everything and went on to start a session: the
        # sentinel is the evidence that nothing was refused before it.
        code, stopped = "engine-stopped", True
    finally:
        sys.stderr = real_err
    text = err.getvalue()
    return {"rc": code, "stderr": text, "engine_stopped": stopped,
            "argparse_rejected": "unrecognized arguments" in text,
            "argparse_message": _line_holding(text, "unrecognized arguments")}


#: `-p` as its own token: the -p inside --permission-mode must not count.
_DASH_P_TOKEN = re.compile(r"(?<![\w-])-p(?![\w-])")


def _names_dash_p(text: str) -> bool:
    """Whether one line names the --job-wait flag AND `-p` as a word.

    Both on the same line, because that is the sentence a caller needs: what
    was refused, and where it does work. Matching `-p` loosely can never
    fail -- `--permission-mode` contains it -- so every run would answer
    "yes, it named -p".
    """
    for line in (text or "").splitlines():
        if "--job-wait" in line and _DASH_P_TOKEN.search(line):
            return True
    return False


def _line_holding(text: str, needle: str) -> str:
    for line in text.splitlines():
        if needle in line:
            return line.strip()[:200]
    return ""


# ---------------------------------------------------------------------------
# variant groups
# ---------------------------------------------------------------------------

def run_variant(variant: str) -> dict:
    """One variant -> its record. Runs in a process of its own."""
    if variant.startswith("parser_"):
        return _variant_parser(variant)
    if variant.startswith("notice_"):
        return _variant_notice(variant)
    if variant in _HEADLESS:
        return _variant_headless(variant)
    raise AssertionError(f"unknown variant {variant}")


def _variant_parser(variant: str) -> dict:
    """The door's answer to ``--job-wait``, driven through ``cli.main``.

    The engine builder raises a sentinel once argparse has answered, so the
    run cannot reach the network: what is under test here IS argparse's answer
    plus the value that arrives on the namespace the door hands on.
    """
    world = _world(variant)
    ws = Path(world["ws"])
    os.chdir(ws)
    from delfin.agent import cli

    cfg = _PARSER[variant]
    _quiet_stdin(tty=cfg["tty"])
    seen: dict = {}

    def _build(args):
        # The sentinel proves the door accepted the flag and went on to start
        # a session. cmd_run / _start_interactive_session catch it as an
        # engine-init failure, so what a variant asserts is never this
        # exception's path but the observable the door owes the caller:
        # `job_wait_on_namespace` (the value that reached cmd_run) or the
        # refusal note on stderr.
        seen["job_wait"] = getattr(args, "job_wait", "__absent__")
        raise EngineRefused("stop-here")

    cli._build_engine = _build
    argv = [cfg["door"]] + (["-p"] if cfg["dash_p"] else [])
    argv += [cfg["prompt"]] if cfg["prompt"] else []
    argv += ["--cwd", str(ws), "--provider", "ollama", "--model", "m",
             "--permission-mode", "acceptEdits", "--job-wait", "0"]
    got = _drive_cli(argv)
    got["job_wait_on_namespace"] = seen.get("job_wait", "__never_built__")
    # A line that names BOTH the flag and `-p` as a standalone word. Matching
    # `-p` loosely can never fail: `--permission-mode` contains it, so every
    # run would answer "yes, it named -p".
    got["names_dash_p"] = _names_dash_p(got["stderr"])
    got["stderr"] = got["stderr"].strip().splitlines()[-1][:180] \
        if got["stderr"].strip() else ""
    return got


def _variant_notice(variant: str) -> dict:
    """Real registry rows through the real renderer, as the terminal calls it."""
    from delfin.agent import job_wake

    world = _world(variant)
    ws = Path(world["ws"])
    events = _notice_events(variant, world)
    text = job_wake.wake_prompt(events)
    return {"text": text,
            "emitted": bool(text),
            "names_job": bool(events) and all(
                str(e.get("job_id")) in text for e in events),
            "state_words": [e.get("state") for e in events],
            "names_output": "bash_output" in text,
            "n_rows": len(events),
            "_ws": str(ws)}


def _watch_with_tool(ws: Path, job_id: str, description: str,
                     session_id: str) -> dict:
    """Register a watch through the model's own ``watch_job`` call."""
    from delfin.agent.api_client import KitToolPermissions, _doc_executor
    perms = KitToolPermissions(workspace=ws, mode="acceptEdits",
                               task_session_id=session_id)
    out = json.loads(_doc_executor.execute(
        "watch_job", {"job_id": job_id, "description": description}, perms))
    if out.get("status") != "watching":
        raise AssertionError(f"watch_job did not watch {job_id}: {out}")
    return out


def _notice_events(variant: str, world: dict) -> list:
    """The rows the real producers built for a real job."""
    import delfin.agent.bash_jobs as bj
    from delfin.agent import job_wake

    ws = Path(world["ws"])
    if variant == "notice_shell_ok":
        jid = _adopt(ws, "s-1", 0.1)
        return _shells(bj, jid)
    if variant == "notice_shell_failed":
        # Sleeps first: adopt() refuses a process that already exited.
        jid = _adopt(ws, "s-1", 0.0,
                     code="import time, sys; time.sleep(0.3); sys.exit(3)")
        return _shells(bj, jid)
    if variant == "notice_slurm":
        # Registered and reported the way the model does it: the watch_job
        # tool answers "end your turn instead of polling", and the idle look
        # is what brings the job back. The scheduler itself is shadowed by
        # the run_fn the product already takes -- squeue is not on this host.
        _watch_with_tool(ws, "4711", "the chain", "s-1")
        seen: set = set()
        return job_wake.finished_watched_jobs(
            ws, seen, session_id="s-1",
            run_fn=lambda argv: "4711 COMPLETED")
    if variant == "notice_sandbox_refused":
        # The real dispatch on THIS host: the cage refuses (exit 126) and the
        # notice must say the command never ran. The control: true on base.
        jid = _start_with_dispatch(ws, "s-1", "echo nothing-is-run-here")
        return _shells(bj, jid)
    raise AssertionError(variant)


def _shells(bj, job_id: str, limit_s: float = 25.0) -> list:
    from delfin.agent import job_wake
    seen: set = set()
    deadline = time.monotonic() + limit_s
    while time.monotonic() < deadline:
        got = job_wake.finished_shells(seen, session_id="s-1")
        if any(str(e.get("job_id")) == job_id for e in got):
            return got
        time.sleep(0.05)
    raise AssertionError(f"job {job_id} never reported")


def _variant_headless(variant: str) -> dict:
    """The one-shot run, end to end, with a job its session owns."""
    world = _world(variant)
    ws = Path(world["ws"])
    os.chdir(ws)
    from delfin.agent import cli

    cfg = _HEADLESS[variant]
    _quiet_stdin()
    job: dict = {}
    built: dict = {}
    model_calls = {"n": 0}
    turn_prompts: list = []
    # The clock for what this package is about: how long after the job ENDED
    # the session noticed. The run's wall clock is not that number -- it also
    # carries the engine's own start-up -- so the answer's turn and the wake's
    # turn are stamped separately.
    clock = {"answer": None, "wake": None, "job_end": None}

    streams = [call_stream("bash_background",
                           {"command": "python -m pytest -q",
                            "description": "the suite"}, "b1"),
               answer_stream(cfg["answer"])]

    def _create(**kwargs):
        # The job joins the registry at the moment the model asks for it --
        # the scripted first reply IS the bash_background call. That timing
        # matters: _resume_or_create has already run by then, so the job's
        # session id is the one the engine will keep, exactly as a real tool
        # call's is. Adopting at engine-build time instead would file the job
        # under the pre-resume id and the wait loop, which scopes by session,
        # would never see it.
        if model_calls["n"] == 0 and cfg["live_s"] is not None:
            sid = str(getattr(built["engine"], "session_id", "") or "")
            job["job_id"], job["pid"] = _adopt_pid(ws, sid, cfg["live_s"])
            job["expected_end"] = time.monotonic() + cfg["live_s"]
        msgs = kwargs.get("messages") or []
        if msgs:
            prompt = str(msgs[-1].get("content") or "")
            turn_prompts.append(prompt)
            if "not the user" in prompt:
                clock["wake"] = time.monotonic()
        model_calls["n"] += 1
        return streams[min(model_calls["n"] - 1, len(streams) - 1)]

    def _build(args):
        engine = _real_build(args)
        engine.client.client.chat.completions.create = _create
        built["engine"] = engine
        return engine

    _real_build = cli._build_engine
    cli._build_engine = _build

    door = "run" if variant.startswith("run_") else "chat"
    argv = [door] + (["-p"] if door == "chat" else [])
    argv += ["do the thing in the background", "--cwd", str(ws),
             "--provider", "ollama", "--model", "m",
             "--permission-mode", "acceptEdits"]
    # Each door with the machine-readable flag IT declares: chat has
    # --output-format (cli.py:6033), run has --json (cli.py:6187). Asking the
    # run door for --output-format is refused by argparse, which is a real
    # asymmetry between the doors but not what this variant measures.
    argv += (["--output-format", "stream-json"] if door == "chat"
             else ["--json"])
    if cfg["flag"] is not None:
        argv += ["--job-wait", str(cfg["flag"])]

    t0 = time.monotonic()
    got = _drive_cli(argv)
    wall = time.monotonic() - t0

    out = {"rc": got["rc"], "wall_s": round(wall, 2),
           "model_calls": model_calls["n"],
           "wake_turns": sum(1 for p in turn_prompts
                             if "not the user" in p),
           "waited": "waiting for it" in got["stderr"],
           "gave_up": "the run ends now" in got["stderr"],
           "bound_named": _bound_named(got["stderr"]),
           "argparse_rejected": got["argparse_rejected"],
           "stderr_tail": got["stderr"].strip().splitlines()[-1][:180]
           if got["stderr"].strip() else ""}
    if "job_id" in job:
        out["job_id"] = job["job_id"]
        if "pid" in job:
            out["job_pid"] = job["pid"]
        # Only a WAKE prompt counts. Every later request also carries the job
        # id -- inside the tool result of the call that started it -- so
        # searching all prompts answers "did the run ever mention the id",
        # not "did the wake report the job".
        wake_prompts = [p for p in turn_prompts if "not the user" in p]
        out["wake_names_job"] = any(job["job_id"] in p for p in wake_prompts)
        # Latency from the job's own end to the turn that reported it. None
        # when nothing ever reported it -- a job nobody noticed is the
        # failure the package names, and it must not read as a fast one.
        out["notice_latency_s"] = (
            round(clock["wake"] - job["expected_end"], 2)
            if clock["wake"] is not None else None)
    return out


def _bound_named(stderr_text: str):
    """The bound the wait loop itself printed, or None."""
    import re
    match = re.search(r"up to (\d+) s", stderr_text or "")
    return int(match.group(1)) if match else None


#: door, -p, prompt, a tty stdin. The engine stops at construction.
_PARSER = {
    "parser_run_takes_it": {"door": "run", "dash_p": False, "tty": False,
                            "prompt": "hi"},
    "parser_chat_rejects": {"door": "chat", "dash_p": True, "tty": False,
                            "prompt": "hi"},
    "parser_chat_interactive": {"door": "chat", "dash_p": False, "tty": True,
                                "prompt": ""},
}

#: door comes from the variant name; flag = what the caller passed.
_HEADLESS = {
    # Controls: the `run` door already has the flag on the base tree.
    "run_job_wait_zero": {"flag": 0, "live_s": LONG_LIVE_S,
                          "answer": "Started it; I will be told."},
    "run_job_wait_bounded": {"flag": 45, "live_s": 1.0,
                             "answer": "Started it; I will be told."},
    # The documented headless door: chat -p, which IS run (cli.py:2581).
    "headless_job_wait_zero": {"flag": 0, "live_s": LONG_LIVE_S,
                               "answer": "Started it; I will be told."},
    "headless_job_wait_bounded": {"flag": 45, "live_s": 1.0,
                                  "answer": "Started it; I will be told."},
    # The default on that door: a job that ends just after the answer.
    "headless_default_wakes": {"flag": None, "live_s": 1.0,
                               "answer": "Started it; I will be told."},
    # A bound the caller asks for must be the bound it gets.
    "headless_bound_is_honoured": {"flag": 5, "live_s": LONG_LIVE_S,
                                   "answer": "Started it; I will be told."},
}

#: Every variant, in report order.
VARIANTS = ["headless_job_wait_zero", "headless_job_wait_bounded",
            "headless_default_wakes", "headless_bound_is_honoured",
            "run_job_wait_zero", "run_job_wait_bounded",
            "parser_chat_rejects", "parser_chat_interactive",
            "parser_run_takes_it",
            "notice_shell_ok", "notice_shell_failed", "notice_slurm",
            "notice_sandbox_refused"]


# ---------------------------------------------------------------------------
# the promise each variant is held to, stated independently of this tree
# ---------------------------------------------------------------------------

def expected(variant: str) -> dict:
    """What a CORRECT product answers for this variant.

    Written from the flag's own help text (``cli.py:6190``: "0: end with the
    answer and stop them"), from ``_cmd_chat_one_shot`` handing the prompt to
    ``cmd_run`` (so the two doors answer one question), and from the package's
    wording ("reports when it ends (exit code, where its output is)"). Nothing
    here is read off what this tree happens to do.
    """
    if variant == "headless_job_wait_zero":
        return {"rc": 0, "wake_turns": 0, "waited": False, "job_stopped": True,
                "reason": "--job-wait 0 ends with the answer and stops the "
                          "shells"}
    if variant == "headless_job_wait_bounded":
        return {"rc": 0, "wake_turns": 1, "waited": True, "bound_named": 45,
                "wake_names_job": True, "latency_max_s": 8.0,
                "reason": "a job that ends 1 s after the answer wakes exactly "
                          "one turn that names it, noticed within seconds of "
                          "the job ending"}
    if variant == "headless_default_wakes":
        return {"rc": 0, "wake_turns": 1, "waited": True,
                "wake_names_job": True,
                "reason": "the default keeps the run up and wakes on the job"}
    if variant == "headless_bound_is_honoured":
        return {"rc": 0, "waited": True, "bound_named": 5, "gave_up": True,
                "wall_max_s": 20.0,
                "reason": "--job-wait 5 holds the run about 5 s and says it "
                          "ended the wait"}
    if variant == "run_job_wait_zero":
        return {"rc": 0, "wake_turns": 0, "waited": False, "job_stopped": True,
                "reason": "control: the run door already answers this"}
    if variant == "run_job_wait_bounded":
        return {"rc": 0, "wake_turns": 1, "waited": True, "bound_named": 45,
                "wake_names_job": True, "latency_max_s": 8.0,
                "reason": "the same promise on the run door: one wake turn "
                          "that names the job, noticed within seconds"}
    if variant == "parser_chat_rejects":
        return {"argparse_rejected": False, "job_wait_on_namespace": 0,
                "reason": "chat -p IS run, so the bound it forwards must "
                          "parse on this door too"}
    if variant == "parser_chat_interactive":
        return {"argparse_rejected": False, "names_dash_p": True,
                "reason": "an interactive session has no one answer to bound "
                          "(its jobs wake at the idle prompt), so the flag "
                          "must parse and SAY it cannot deliver here, naming "
                          "-p -- the way --output-format does at cli.py:2807"}
    if variant == "parser_run_takes_it":
        return {"argparse_rejected": False, "job_wait_on_namespace": 0,
                "reason": "control: already true on the base tree"}
    if variant == "notice_shell_ok":
        return {"emitted": True, "n_rows": 1, "names_job": True,
                "state_has": "ok", "names_output": True,
                "reason": "the notice says how it ended AND where its output "
                          "is, as the in-turn notice already does"}
    if variant == "notice_shell_failed":
        return {"emitted": True, "n_rows": 1, "names_job": True,
                "state_has": "exit 3", "names_output": True,
                "reason": "a failed job says which code it ended with and "
                          "where to read it"}
    if variant == "notice_slurm":
        return {"emitted": True, "n_rows": 1, "names_job": True,
                "state_has": "COMPLETED",
                "reason": "a SLURM row says what the scheduler answered"}
    if variant == "notice_sandbox_refused":
        return {"emitted": True, "n_rows": 1, "names_job": True,
                "state_has": "not run",
                "reason": "control: a cage refusal must not read as a failing "
                          "command (job_wake.py:89)"}
    raise AssertionError(variant)


def grade(variant: str, got: dict) -> dict:
    """REACHED / USED / CORRECT for one variant, from its measurement."""
    want = expected(variant)
    problems = []
    for key, value in want.items():
        if key == "reason":
            continue
        if key == "wall_max_s":
            if got.get("wall_s", 1e9) > value:
                problems.append(f"held the run {got['wall_s']} s, promised "
                                f"<= {value} s")
            continue
        if key == "latency_max_s":
            late = got.get("notice_latency_s")
            if late is None:
                problems.append("nothing ever reported the finished job")
            elif late > value:
                problems.append(f"noticed the job {late} s after it ended, "
                                f"promised <= {value} s")
            continue
        if key == "state_has":
            joined = " ".join(str(s) for s in got.get("state_words", []))
            if str(value).lower() not in joined.lower():
                problems.append(f"state {joined!r} does not say {value!r}")
            continue
        have = got.get(key, "__missing__")
        if have != value:
            problems.append(f"{key}={have!r}, promised {value!r}")
    reached, used = _reached(variant, got), _used(variant, got)
    return {"REACHED": reached, "USED": used,
            "CORRECT": bool(reached and used and not problems),
            "problems": problems, "promised": want.get("reason", "")}


def _reached(variant: str, got: dict) -> bool:
    """The capability's own code ran -- by an observable, not a patched flag."""
    if variant.startswith("parser_"):
        # The step that decides everything after it: did the door take the
        # caller's bound at all. The engine behind it is stopped by a
        # sentinel here, so its exit status is not the capability's answer.
        return not got.get("argparse_rejected")
    if got.get("argparse_rejected"):
        # argparse refused before the caller's value reached anything: the
        # waiting code never saw the request.
        return False
    if variant.startswith("notice_"):
        return bool(got.get("emitted"))
    return got.get("rc") == 0


def _used(variant: str, got: dict) -> bool:
    """The capability was exercised, not merely reachable."""
    if variant.startswith("parser_"):
        # USED here means the caller's value reached the object that acts on
        # it: the namespace cmd_run reads, or -- for an interactive session,
        # where there is nothing to bound -- the note the person is shown.
        # A parser that accepts a flag into nothing is the exact defect
        # test_the_flags_a_run_can_be_bounded_by.py was written for.
        if variant == "parser_chat_interactive":
            return bool(got.get("names_dash_p"))
        return got.get("job_wait_on_namespace") == 0
    if variant.startswith("notice_"):
        return bool(got.get("names_job"))
    return bool(got.get("waited")) or got.get("wake_turns", 0) > 0 \
        or bool(got.get("gave_up"))


# ---------------------------------------------------------------------------
# parent: spawn every variant, aggregate, write the file
# ---------------------------------------------------------------------------

def run_all(out_path: Path, label: str) -> int:
    rows = []
    for vid in VARIANTS:
        dest = HERE / ".runs" / f"{vid}.json"
        cmd = [sys.executable, str(HERE / "scenario.py"), "--variant", vid,
               "--out", str(dest)]
        t0 = time.monotonic()
        try:
            proc = subprocess.run(cmd, cwd=str(REPO), capture_output=True,
                                  text=True, timeout=900,
                                  env={**os.environ, "PYTHONPATH": str(REPO)})
        except subprocess.TimeoutExpired as exc:
            rows.append({"variant": vid, "REACHED": False, "USED": False,
                         "CORRECT": False,
                         "problems": [f"the variant timed out: "
                                      f"{str(exc.stderr)[-200:]}"]})
            continue
        if not dest.exists():
            rows.append({"variant": vid, "REACHED": False, "USED": False,
                         "CORRECT": False,
                         "problems": [f"the variant process died (rc "
                                      f"{proc.returncode}): "
                                      f"{proc.stderr.strip()[-300:]}"]})
            continue
        got = json.loads(dest.read_text())
        got.pop("_ws", None)
        # A job the run was supposed to end with: was it really gone once the
        # run's own process had exited? (The stopping is an atexit handler.)
        pid = got.get("job_pid")
        if pid:
            time.sleep(0.5)
            got["job_alive_after_run"] = _pid_alive(int(pid))
            got["job_stopped"] = not got["job_alive_after_run"]
        row = {"variant": vid, **grade(vid, got), "measured": got}
        rows.append(row)
        print(f"{vid:28} "
              f"R={'Y' if row['REACHED'] else 'N'} "
              f"U={'Y' if row['USED'] else 'N'} "
              f"C={'Y' if row['CORRECT'] else 'N'}  "
              f"{'; '.join(row['problems'])[:100]:100} "
              f"[{time.monotonic() - t0:.0f}s]", flush=True)

    tally = {k: sum(1 for r in rows if r[k])
             for k in ("REACHED", "USED", "CORRECT")}
    timed = [r for r in rows if "wall_s" in r.get("measured", {})]
    doc = {"package": "p2-bg-job-wake", "label": label,
           "command": ("PYTHONPATH=. python "
                       "experiments/parity/p2-bg-job-wake/scenario.py "
                       "--all --label <label> --out <this file>"),
           "entry_point": ("cli.main([...]): delfin-agent chat -p ... "
                           "--output-format stream-json, delfin-agent run, "
                           "delfin-agent chat"),
           "n_variants": len(rows), "tally": tally, "rows": rows,
           "waiting_cost": {
               "population": "the 6 one-shot-run variants only; the parser "
                             "and notice variants make no model call",
               "model_calls": sum(r["measured"]["model_calls"]
                                  for r in timed),
               "wake_turns": sum(r["measured"]["wake_turns"] for r in timed),
               "wall_s": round(sum(r["measured"]["wall_s"] for r in timed),
                               2)},
           "git": _git()}
    out_path.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\n{label}: REACHED {tally['REACHED']}/{len(rows)}  "
          f"USED {tally['USED']}/{len(rows)}  "
          f"CORRECT {tally['CORRECT']}/{len(rows)}  -> {out_path}")
    return 0


def _git() -> dict:
    def _run(*a):
        try:
            return subprocess.run(["git", *a], cwd=str(REPO),
                                  capture_output=True, text=True,
                                  timeout=20).stdout.strip()
        except Exception:
            return ""
    return {"rev": _run("rev-parse", "HEAD"),
            "branch": _run("rev-parse", "--abbrev-ref", "HEAD"),
            "dirty": bool(_run("status", "--porcelain"))}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--label", default="measured")
    args = ap.parse_args()
    # Resolved before anything else: a variant chdirs into its scratch
    # workspace, and a relative --out would then be written in there.
    raw = Path(args.out) if args.out else None
    out = None if raw is None else (
        raw if raw.is_absolute() else LAUNCH_CWD / raw)
    if args.all:
        return run_all(out, args.label)
    got = run_variant(args.variant)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(got, indent=2, default=str) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
