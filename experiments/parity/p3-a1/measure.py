"""Parity measurement for PACKAGE p3 — the review as a flow.

Drives DELFIN's REAL entry point (``delfin/agent/cli.py:main``, the same
front door ``delfin-agent chat -p ... --output-format stream-json`` uses)
with the model's HTTP stream scripted. Everything between the entry point
and the finding is real: the engine's turn, the tool dispatch, the
``subagent`` executor, ``run_subagent``, the preset registry, the schema
validation. Only the model's words are scripted.

Per variant three endpoints are recorded:

    REACHED  the findings-contract code executed — proven by a ``sys.settrace``
             call-hook that fires only for frames whose code object lives in
             the contract module. On a tree without that module the hook
             cannot fire, which is a fact about the tree, not a guess.
    USED     the product path carried the contract's answer to the parent:
             the ``subagent`` tool's payload holds the normalised review.
    CORRECT  the outcome is right: every verifiable finding survives with its
             ``file:line`` and its failure scenario, every unverifiable one
             (missing path, line past EOF, no scenario, verified=false) is
             dropped AND counted, and the survivors are severity-ranked.

What this CANNOT measure, and does not report: recall and precision of a
review against real defects. Those are properties of the model reading the
code; scripting the model would make the number a measurement of this
script. They are named as probes in scratch/realrun.json instead.

Usage:
    python experiments/parity/p3-a1/measure.py --out experiments/parity/p3-a1/before.json
"""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
_TESTS = REPO / "tests"
if str(_TESTS) not in sys.path:
    sys.path.insert(0, str(_TESTS))
# The script's OWN directory is sys.path[0] when it is run as a file, and the
# worktree root is not on the path at all -- so `import delfin` would resolve to
# whatever copy is installed. A measurement of another tree is a measurement of
# nothing, so the worktree goes first and the summary records what was loaded.
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

CONTRACT_MODULE = "review_findings"


# ---------------------------------------------------------------------------
# The scripted model stream (parent turn + the delegate's turn)
# ---------------------------------------------------------------------------

class _FnDelta:
    def __init__(self, name=None, arguments=None):
        self.name = name
        self.arguments = arguments


class _TCDelta:
    def __init__(self, index, id=None, name=None, arguments=None):
        self.index = index
        self.id = id
        self.function = _FnDelta(name, arguments)


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


def _tool_stream(name: str, args: dict, call_id: str = "c1") -> _Stream:
    return _Stream([_Chunk([_Choice(_Delta(tool_calls=[_TCDelta(
        0, id=call_id, name=name, arguments=json.dumps(args))]),
        finish="tool_calls")], usage=_Usage())])


def _text_stream(text: str) -> _Stream:
    return _Stream([_Chunk([_Choice(_Delta(content=text), finish="stop")],
                           usage=_Usage())])


# ---------------------------------------------------------------------------
# The corpus: real DELFIN code with a defect planted in it, per variant
# ---------------------------------------------------------------------------

#: A copy of a real, small DELFIN module. Each variant plants ONE defect
#: into it at the line the scripted reviewer will cite, so the location the
#: finding names is a location that really exists (or, in the adversarial
#: variants, one that really does not).
_CORPUS_SOURCE = '''\
"""A small module, copied into the workspace so a finding has a real file."""

from __future__ import annotations

import json


DEFAULT_LIMIT = 10


def read_rows(path, limit=DEFAULT_LIMIT):
    """Read at most ``limit`` rows from a JSON-lines file."""
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))
            if len(rows) >= limit:
                break
    return rows


def total(rows, field):
    return sum(row[field] for row in rows)


def summarise(path):
    rows = read_rows(path)
    return {"count": len(rows), "total": total(rows, "value")}
'''

#: (variant id, defect description, planted line, extra findings the scripted
#: reviewer reports). ``planted line`` is the 1-based line of ``DEFAULT_LIMIT``
#: in the corpus above; the adversarial values deliberately point at things
#: that do not exist.
_CORPUS_LINES = {
    "limit": _CORPUS_SOURCE.splitlines().index("DEFAULT_LIMIT = 10") + 1,
}
_PLANTED = _CORPUS_LINES["limit"]
_FILE = "delfin/agent/inventory.py"


def _finding(severity: str, line: int, scenario: str, *,
             verified: bool = True, path: str = _FILE) -> dict:
    return {
        "severity": severity,
        "file": path,
        "line": line,
        "what": "the default limit is mutable state the caller can mutate",
        "failure_scenario": scenario,
        "verified": verified,
    }


#: 12 variants. Each carries a MIX: findings that can be checked against the
#: workspace and findings that cannot. The product must keep the first kind
#: and drop-and-count the second.
_VARIANTS: tuple[dict, ...] = (
    {"id": "v01_verified_only",
     "findings": [_finding("major", _PLANTED,
                           "read_rows() called twice after a caller sets "
                           "DEFAULT_LIMIT to 0 returns [] both times")],
     "expect_kept": 1, "expect_dropped": 0},
    {"id": "v02_missing_file",
     "findings": [
         _finding("major", _PLANTED,
                  "read_rows() called twice after a caller sets "
                  "DEFAULT_LIMIT to 0 returns [] both times"),
         _finding("critical", 12, "the parser drops every row after a blank "
                  "line", path="delfin/agent/does_not_exist.py")],
     "expect_kept": 1, "expect_dropped": 1},
    {"id": "v03_line_past_eof",
     "findings": [
         _finding("major", _PLANTED, "a caller that sets DEFAULT_LIMIT to 0 "
                  "silently gets an empty list"),
         _finding("critical", 9000, "the loop never terminates")],
     "expect_kept": 1, "expect_dropped": 1},
    {"id": "v04_no_scenario",
     "findings": [
         _finding("major", _PLANTED, "a caller that sets DEFAULT_LIMIT to 0 "
                  "silently gets an empty list"),
         _finding("minor", _PLANTED, "")],
     "expect_kept": 1, "expect_dropped": 1},
    {"id": "v05_unverified",
     "findings": [
         _finding("major", _PLANTED, "a caller that sets DEFAULT_LIMIT to 0 "
                  "silently gets an empty list"),
         _finding("critical", _PLANTED,
                  "I think this is probably wrong",
                  verified=False)],
     "expect_kept": 1, "expect_dropped": 1},
    {"id": "v06_ranking",
     "findings": [
         _finding("minor", _PLANTED, "a caller that sets DEFAULT_LIMIT to 0 "
                  "silently gets an empty list"),
         _finding("critical", _PLANTED, "two concurrent callers mutate "
                  "DEFAULT_LIMIT and both lose rows")],
     "expect_kept": 2, "expect_dropped": 0},
    {"id": "v07_all_unverifiable",
     "findings": [
         _finding("major", 9000, "the loop never terminates"),
         _finding("critical", 12, "the parser drops rows",
                  path="delfin/agent/not_here.py")],
     "expect_kept": 0, "expect_dropped": 2},
    {"id": "v08_three_good_one_bad",
     "findings": [
         _finding("major", _PLANTED, "a caller that sets DEFAULT_LIMIT to 0 "
                  "silently gets an empty list"),
         _finding("minor", _PLANTED, "total() raises KeyError on a row "
                  "without 'value'"),
         _finding("critical", _PLANTED, "summarise() on an empty file returns "
                  "count 0 and crashes on total"),
         _finding("major", _PLANTED, "guess", verified=False)],
     "expect_kept": 3, "expect_dropped": 1},
    {"id": "v09_path_outside_workspace",
     "findings": [
         _finding("major", _PLANTED, "a caller that sets DEFAULT_LIMIT to 0 "
                  "silently gets an empty list"),
         _finding("critical", 1, "the whole file is wrong",
                  path="../../etc/passwd")],
     "expect_kept": 1, "expect_dropped": 1},
    {"id": "v10_empty_findings",
     "findings": [],
     "expect_kept": 0, "expect_dropped": 0},
    {"id": "v11_line_zero",
     "findings": [
         _finding("major", 0, "the file is unreadable"),
         _finding("major", _PLANTED, "a caller that sets DEFAULT_LIMIT to 0 "
                  "silently gets an empty list")],
     "expect_kept": 1, "expect_dropped": 1},
    {"id": "v12_severity_unknown",
     "findings": [
         _finding("blocker", _PLANTED, "a caller that sets DEFAULT_LIMIT to 0 "
                  "silently gets an empty list"),
         _finding("major", _PLANTED, "total() raises KeyError on a row "
                  "without 'value'")],
     "expect_kept": 2, "expect_dropped": 0},
)


# ---------------------------------------------------------------------------
# The trace hook: did the contract module run at all?
# ---------------------------------------------------------------------------

class _TraceHit:
    """Records a call into the findings contract, however deeply reached.

    Matched on the FUNCTION, not on the file. ``subagents.py`` is entered on
    any delegation at all and exists on the untouched tree, so a file match
    would report the capability as reached on a run that has no findings
    contract whatsoever -- a set membership standing in for the call.
    """

    NAMES = frozenset({"verify_review_findings", "attach_review_findings"})

    def __init__(self) -> None:
        self.hit = False
        self.seen: set[str] = set()

    def __call__(self, frame, event, arg):
        if event == "call":
            code = frame.f_code
            if code.co_name in self.NAMES and \
                    code.co_filename.endswith("subagents.py"):
                self.hit = True
                self.seen.add(code.co_name)
        # The tracer must return ITSELF, not None. A local trace function
        # returning None switches tracing off for that frame, so the call is
        # recorded and every function that frame then calls is never seen --
        # which is how the first version of this hook reported REACHED=False
        # on a run whose contract code demonstrably ran (USED was 12/12).
        return self


# ---------------------------------------------------------------------------
# One variant, end to end through cli.main
# ---------------------------------------------------------------------------

def _run_variant(variant: dict, root: Path) -> dict:
    """Drive ``cli.main`` for one variant. Returns the recorded endpoints."""
    import io
    from unittest import mock

    from delfin.agent import api_client as _A
    from delfin.agent import cli as agent_cli
    from delfin.agent import mcp_client as M
    from delfin.agent import model_capabilities as mc
    from delfin.agent import provider_profile as PP

    home = root / "home"
    ws = home / "ws"
    (ws / "delfin" / "agent").mkdir(parents=True, exist_ok=True)
    (home / ".delfin").mkdir(parents=True, exist_ok=True)
    (ws / _FILE).write_text(_CORPUS_SOURCE, encoding="utf-8")

    review_payload = {"findings": variant["findings"]}
    # The marker the DELEGATE's prompt carries and the parent's never does. The
    # parent's prompt names the variant too (it is what the run is about), so a
    # variant id used as the marker made the parent's first round answer with
    # the review JSON and no tool call was ever made -- the run "passed" while
    # measuring nothing.
    delegate_marker = f"DELEGATE-{variant['id']}"

    class _EmptyReg:
        def discover_all(self): return []
        def discover_resources(self): return []
        def discover_prompts(self): return []

    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)

    seen: dict = {"delegate_payload": None, "dispatch": 0}

    def _create(**kwargs):
        """Answer by WHO is asking: the delegate's prompt carries the marker,
        the parent's does not."""
        blob = json.dumps(kwargs.get("messages") or [])
        if delegate_marker in blob:
            # The delegate's own turn: its final message is the review JSON.
            return _text_stream(json.dumps(review_payload))
        if seen["dispatch"] == 0:
            seen["dispatch"] += 1
            return _tool_stream("subagent", {
                "subagent_type": "code-reviewer",
                "description": "review the diff",
                "prompt": (f"{delegate_marker}: review {_FILE} for defects. "
                           "Report only what you verified."),
            })
        return _text_stream("review complete.")

    def _prepare(engine):
        engine.client.client.chat.completions.create = _create
        real = engine.kit_permissions.subagent_runner

        def _wrapped(**kw):
            payload = real(**kw)
            seen["delegate_payload"] = payload
            return payload

        engine.kit_permissions.subagent_runner = _wrapped
        engine.client.set_permissions(engine.kit_permissions)
        # set_permissions re-binds the runner from the client; put ours back.
        engine.kit_permissions.subagent_runner = _wrapped

    out = io.StringIO()
    trace = _TraceHit()
    result: dict = {
        "variant": variant["id"],
        "reachable": False, "used": False, "correct": False,
        "detail": "", "kept": None, "dropped": None,
    }
    cwd = Path.cwd()
    try:
        with mock.patch.object(Path, "home", classmethod(lambda cls: home)), \
             mock.patch.object(M, "get_registry", lambda *a, **k: _EmptyReg()), \
             mock.patch.object(mc, "resolve", lambda *a, **k: caps), \
             mock.patch.object(agent_cli, "_save_session",
                               lambda e, r, **kw: "sid-p3"), \
             mock.patch.object(agent_cli, "_resume_or_create",
                               lambda e, a: ""), \
             mock.patch.object(_A, "_caged_command_may_run_unasked",
                               lambda *a, **k: False), \
             mock.patch.object(PP, "_DEFAULT_PATH",
                               root / "learned_profiles.json"), \
             mock.patch("sys.stdout", out), \
             mock.patch("sys.stdin", _fake_stdin()):
            (root / "learned_profiles.json").write_bytes(
                (REPO / "delfin" / "agent" / "learned_profiles.json").read_bytes())
            import os
            os.environ["HOME"] = str(home)
            os.chdir(ws)
            real_build = agent_cli._build_engine

            def _build(args):
                engine = real_build(args)
                _prepare(engine)
                return engine

            with mock.patch.object(agent_cli, "_build_engine", _build):
                sys.settrace(trace)
                try:
                    agent_cli.main(["-p", f"review the diff ({variant['id']})",
                                    "--cwd", str(ws), "--provider", "ollama",
                                    "--model", "m",
                                    "--permission-mode", "default",
                                    "--output-format", "json"])
                finally:
                    sys.settrace(None)
    except Exception:
        result["detail"] = "run raised: " + traceback.format_exc(limit=3)
        return result
    finally:
        try:
            import os
            os.chdir(cwd)
        except Exception:
            pass

    payload = seen["delegate_payload"] or {}
    result["reachable"] = bool(trace.hit)
    review = payload.get("review") if isinstance(payload, dict) else None
    result["used"] = isinstance(review, dict)
    if not isinstance(review, dict):
        result["detail"] = (
            "the parent's subagent payload carries no normalised review "
            f"(keys: {sorted(payload) if isinstance(payload, dict) else type(payload).__name__})")
        return result

    findings = review.get("findings") or []
    dropped = review.get("dropped") or []
    result["kept"] = len(findings)
    result["dropped"] = len(dropped)

    problems = []
    if len(findings) != variant["expect_kept"]:
        problems.append(f"kept {len(findings)} of {variant['expect_kept']}")
    if len(dropped) != variant["expect_dropped"]:
        problems.append(f"dropped {len(dropped)} of {variant['expect_dropped']}")
    for f in findings:
        if not isinstance(f.get("line"), int) or f["line"] < 1:
            problems.append(f"kept a finding with no usable line: {f}")
        if not str(f.get("failure_scenario") or "").strip():
            problems.append(f"kept a finding with no scenario: {f}")
    # Ranking: severity descending by the contract's own order.
    order = {"critical": 0, "major": 1, "minor": 2}
    ranks = [order.get(str(f.get("severity") or "").lower(), 9)
             for f in findings]
    if ranks != sorted(ranks):
        problems.append(f"not severity-ranked: {ranks}")
    # Nothing that could not be verified may appear in the report.
    for f in findings:
        if "verification" in f or "verified" in f:
            if not f.get("verification") and f.get("verified") is not True:
                problems.append(f"kept an unverified finding: {f}")
    result["correct"] = not problems
    result["detail"] = "; ".join(problems) if problems else "ok"
    return result


def _fake_stdin():
    import io
    s = io.StringIO("")
    s.isatty = lambda: False
    return s


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", default="")
    args = ap.parse_args(argv)

    try:
        import delfin.agent.review_findings as _rf
        contract = Path(_rf.__file__).name
    except Exception:
        contract = ""

    import delfin as _delfin
    loaded_from = str(Path(_delfin.__file__).resolve().parent.parent)
    if loaded_from != str(REPO):
        print(f"REFUSING TO MEASURE: `delfin` loaded from {loaded_from}, "
              f"not the tree under test {REPO}. Re-run with "
              f"PYTHONPATH={REPO}.", file=sys.stderr)
        return 2

    results = []
    for variant in _VARIANTS:
        if args.only and args.only not in variant["id"]:
            continue
        import tempfile
        with tempfile.TemporaryDirectory(prefix="p3-") as td:
            results.append(_run_variant(variant, Path(td)))

    summary = {
        "package": "p3-a1",
        "contract_module_present": bool(contract),
        "commit": _git_head(),
        "measured_tree": loaded_from,
        "command": ("python experiments/parity/p3-a1/measure.py --out "
                    f"{args.out}"),
        "variants": len(results),
        "reached": sum(1 for r in results if r["reachable"]),
        "used": sum(1 for r in results if r["used"]),
        "correct": sum(1 for r in results if r["correct"]),
        "rows": results,
    }
    Path(args.out).write_text(json.dumps(summary, indent=2) + "\n",
                              encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"},
                     indent=2))
    for r in results:
        print(f"  {r['variant']:<28} reached={r['reachable']!s:<5} "
              f"used={r['used']!s:<5} correct={r['correct']!s:<5} {r['detail']}")
    return 0


def _git_head() -> str:
    """The measured tree, through the product's own reader.

    Not a second implementation: ``commit_metrics.git_commit`` already
    answers this and appends ``-dirty``, which matters here — a run that
    describes a commit while measuring an edited tree is a row about
    neither.
    """
    from delfin.agent import commit_metrics
    return commit_metrics.git_commit(REPO)


if __name__ == "__main__":
    raise SystemExit(main())
