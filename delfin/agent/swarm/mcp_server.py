"""MCP server for running and watching waves of DELFIN agents.

An orchestrator -- a person's own coding agent, another DELFIN session, an
IDE -- starts a wave, waits for what happens in it, reads what the members
produced, steers them and stops them, all through these tools instead of
scripts of its own around DELFIN's internals:

  wave_compose                                        the wave's task list from packages
  wave_runs, wave_status, wave_review, wave_log      read
  wave_watch                                          wait for the next events
  wave_start, wave_add_tasks, wave_resume, wave_stop  act
  wave_tick                                           one pass now
  wave_emergency_stop, wave_clear_stop                the brake
  wave_friction                                       what the members ran into
  wave_land                                           checked, reviewed, opened as a PR
  wave_guard_install, wave_guard_remove               the tick as a service or schedule

``delfin-agent swarm wave mcp`` serves them on stdio. The handlers are
plain functions returning JSON text, testable without the ``mcp``
package; ``run_server`` only wires them up, as ``delfin/tools/mcp_server.py``
does.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any

from . import friction, guard, wave


def _dumps(obj: Any) -> str:
    return json.dumps(obj, indent=1, default=str)


def _answer(fn, *a, **kw) -> str:
    try:
        return _dumps(fn(*a, **kw))
    except (wave.WaveError, ValueError) as exc:
        return _dumps({"error": str(exc)})
    except subprocess.CalledProcessError as exc:
        return _dumps({"error": f"git: {(exc.stderr or exc.stdout or '').strip()[:400]}"})


# -- handlers -------------------------------------------------------------------------

def h_runs() -> str:
    def runs():
        out = []
        for run in wave.runs():
            members = wave.status(run)["members"]
            states = {}
            for m in members:
                states[m["state"]] = states.get(m["state"], 0) + 1
            out.append({"run": run, "members": len(members), "states": states,
                        "final": sum("final" in m["marks"] for m in members)})
        return {"root": str(wave.root()), "runs": out,
                "stopped": guard._flag().read_text().strip() if guard._flag().exists() else ""}
    return _answer(runs)


def h_status(run: str) -> str:
    return _answer(wave.status, run)


def h_start(run: str, tasks: list, repo: str, stages_on: bool = True,
            parallel: int = 0) -> str:
    def start():
        if guard._flag().exists():
            raise wave.WaveError("the swarm is stopped (EMERGENCY_STOPPED): clear the stop first")
        res = wave.start(run, tasks, repo, parallel=parallel or None, wait=False)
        with wave._state_lock(run):
            st = wave.load_state(run)
            st["stages"] = bool(stages_on)
            st["driven"] = True
            wave._save_state(run, st)
        return res
    return _answer(start)


def h_compose(run: str, packages: list, arms: list, integrator: dict, lander: dict,
              reviewer: dict, fixer: dict, base: str = "origin/main",
              stages: list | None = None, arm_refs: dict | None = None) -> str:
    """The whole task list of a wave from its packages; starts nothing."""
    def composed():
        from . import compose
        kw = {"base": base, "arm_refs": arm_refs}
        if stages:
            kw["stages"] = tuple(stages)
        tasks = compose.compose(run, packages, arms=arms, integrator=integrator,
                                lander=lander, reviewer=reviewer, fixer=fixer, **kw)
        return {"run": run, "tasks": tasks,
                "summary": [{"id": t["id"], "model": t.get("model"),
                             "after": t.get("after", [])} for t in tasks]}
    return _answer(composed)


def h_add_tasks(run: str, tasks: list) -> str:
    return _answer(wave.add_tasks, run, tasks)


def h_watch(runs: list, timeout_s: float = 600, cursors: dict | None = None) -> str:
    return _answer(wave.watch, list(runs or wave.runs()), timeout_s=min(float(timeout_s), 3600),
                   cursors=cursors or {})


def h_resume(run: str, member: str, note: str = "", fresh: bool = False,
             after_stop: bool = False) -> str:
    """Continue an ended member. ``after_stop`` is the person lifting their
    own stop: without it a stopped member is not restarted."""
    return _answer(wave.resume, run, member, note, fresh=fresh, wait_for_room=False,
                   after_stop=after_stop)


def h_stop(run: str, member: str = "") -> str:
    return _answer(wave.stop, run, member)


def h_tick() -> str:
    return _answer(guard.tick)


def h_emergency_stop(why: str) -> str:
    return _answer(guard.emergency_stop, why or "stopped by the orchestrator")


def h_clear_stop() -> str:
    return _answer(lambda: {"cleared": guard.clear_stop()})


def h_review(run: str, member: str) -> str:
    """What a member produced: its write-up, its commits and diff size
    against its base, how its last session ended, its real runs."""
    def review():
        st = wave.load_state(run)
        a = st["agents"].get(member)
        if a is None:
            raise wave.WaveError(f"{run}/{member} was never started")
        wt = Path(a["worktree"])
        base = str((a.get("task") or {}).get("base") or "origin/main")

        def git(*args):
            from delfin.agent import safe_git
            return safe_git.out(wt, *args, timeout=60).strip()

        def text(p: Path, limit: int = 20000) -> str:
            try:
                return p.read_text(errors="replace")[:limit]
            except OSError:
                return ""

        res = wave._last_result(Path(a["log"]))
        return {"run": run, "member": member, "branch": a["branch"], "base": base,
                "commits": git("log", "--oneline", f"{base}..HEAD").splitlines()[:50],
                "diffstat": git("diff", "--stat", f"{base}...HEAD").splitlines()[-40:],
                "uncommitted": git("status", "--porcelain").splitlines()[:40],
                "pr_md": text(wt / "scratch" / "PR.md"),
                "realrun_summary": text(wt / "scratch" / "realrun" / "summary.md", 8000),
                "last_session": {"stop_reason": res.get("stop_reason"),
                                 "error": (res.get("error") or "")[:400],
                                 "answer_tail": (res.get("text") or "")[-1500:]}}
    return _answer(review)


def h_log(run: str, member: str, tail: int = 80) -> str:
    """The member's recent steps: tool calls, notices, refusals and its
    text, newest last (redacted)."""
    def log():
        a = wave.load_state(run)["agents"].get(member)
        if a is None:
            raise wave.WaveError(f"{run}/{member} was never started")
        out, text = [], []
        for line in Path(a["log"]).read_text(errors="replace").splitlines():
            if not line.startswith("{"):
                continue
            try:
                e = json.loads(line)
            except ValueError:
                continue
            t = e.get("type")
            if t == "text":
                text.append(e.get("text") or "")
                continue
            if text:
                out.append({"type": "text", "text": "".join(text)[-600:]})
                text = []
            if t == "tool_use":
                out.append({"type": "tool", "name": e.get("name"),
                            "input": json.dumps(e.get("input"))[:300]})
            elif t in ("notice", "permission_denied", "stall", "result"):
                out.append({"type": t, **{k: (str(v)[:300] if isinstance(v, str) else v)
                                          for k, v in e.items()
                                          if k in ("text", "tool", "reason_code", "class",
                                                   "stop_reason", "error")}})
        if text:
            out.append({"type": "text", "text": "".join(text)[-600:]})
        rows = out[-max(1, min(int(tail), 400)):]
        for r in rows:
            for k in ("text", "input", "error"):
                if isinstance(r.get(k), str):
                    r[k] = friction._redact(r[k])
        return {"run": run, "member": member, "steps": rows}
    return _answer(log)


def h_friction(runs: list | None = None, top: int = 30) -> str:
    return _answer(friction.report, list(runs) if runs else None, top=int(top))


def h_land(runs: list | None = None, member: str = "", push: bool = False) -> str:
    """Check, independently review and (with ``push``) open the PRs of the
    members that wrote scratch/PR.md. Without ``push`` nothing leaves this
    machine: the checks and the review run, and each record says what would
    happen."""
    def land_them():
        from . import land
        lines: list[str] = []
        results = land.land(list(runs) if runs else None, dry_run=not push,
                            only={member} if member else None, log=lines.append)
        return {"results": results, "log": lines}
    return _answer(land_them)


def h_guard_install(every_minutes: int = 1, systemd: bool = False) -> str:
    if systemd:
        return _answer(guard.install_unit, int(every_minutes) * 60)
    return _answer(guard.install_entry, int(every_minutes))


def h_guard_remove() -> str:
    return _answer(lambda: {"removed": guard.remove_entry(), "unit_removed": guard.remove_unit()})


# -- the stdio server ----------------------------------------------------------------

def run_server(argv: list[str] | None = None) -> None:
    from delfin.mcp_compat import load_server_class

    mcp = load_server_class()(
        "delfin-swarm",
        instructions=(
            "Run and watch waves of DELFIN agents. A wave is a task list; each task is one "
            "caged agent session in its own git worktree and branch, started when what it "
            "waits for has ended. Start with wave_start, then call wave_watch in a loop: it "
            "returns when something happened (a member ended, went silent, crashed, became "
            "final). Read a member's work with wave_review and wave_log, steer it with "
            "wave_resume and a note, stop it with wave_stop. wave_emergency_stop stops "
            "everything; only a person should clear it. Members commit but cannot push: "
            "review and land their branches yourself."))

    @mcp.tool()
    def wave_runs() -> str:
        """Every run under the waves root: members by state, how many are final, and
        whether the swarm is stopped."""
        return h_runs()

    @mcp.tool()
    def wave_status(run: str) -> str:
        """One row per task of a run: state, commits, uncommitted files, PR.md written,
        how its last session ended, stage marks, minutes since its log grew."""
        return h_status(run)

    @mcp.tool()
    def wave_start(run: str, tasks: list, repo: str, stages: bool = True,
                   parallel: int = 0) -> str:
        """Start a wave: tasks = [{"id", "brief", "base"?, "provider"?, "model"?,
        "seconds"?, "turns"?, "after"?, "after_min"?, "after_wait_s"?, "collect"?,
        "attach"?, "hide"?, "no_realrun"?}], repo = the git repository. Starts what is
        ready now; the tick starts the rest as their dependencies end. With stages,
        every member that ends gets its real runs and a resume to read them."""
        return h_start(run, tasks, repo, stages, parallel)

    @mcp.tool()
    def wave_compose(run: str, packages: list, arms: list, integrator: dict, lander: dict,
                     reviewer: dict, fixer: dict, base: str = "origin/main",
                     stages: list | None = None, arm_refs: dict | None = None) -> str:
        """The wave as it has been run, as a task list: packages [{id, title, body}];
        arms [{provider, model}, ...] of different families; integrator, lander,
        reviewer (another family than the lander), fixer {provider, model}. Returns the
        tasks for wave_start; starts nothing. arm_refs {pkg: ["run/tid", ...]} reuses arms
        that already ran. docs/ORCHESTRATING_WAVES.md says how."""
        return h_compose(run, packages, arms, integrator, lander, reviewer, fixer, base, stages,
                         arm_refs)

    @mcp.tool()
    def wave_add_tasks(run: str, tasks: list) -> str:
        """Append tasks to a running wave (a later stage: reviewers after a lander, a
        fixer after them); they start when what they wait for has ended."""
        return h_add_tasks(run, tasks)

    @mcp.tool()
    def wave_watch(runs: list | None = None, timeout_s: float = 600,
                   cursors: dict | None = None) -> str:
        """Wait (at most timeout_s, up to 3600) until something happens in the runs
        (default: all); returns the events and the cursors to pass to the next call."""
        return h_watch(runs or [], timeout_s, cursors)

    @mcp.tool()
    def wave_review(run: str, member: str) -> str:
        """A member's write-up (scratch/PR.md), commits and diffstat against its base,
        uncommitted files, how its last session ended, and its real-run summary."""
        return h_review(run, member)

    @mcp.tool()
    def wave_log(run: str, member: str, tail: int = 80) -> str:
        """A member's recent steps (tool calls, notices, refusals, its text), redacted."""
        return h_log(run, member, tail)

    @mcp.tool()
    def wave_resume(run: str, member: str, note: str = "", fresh: bool = False) -> str:
        """Continue an ended member's session with a note (your review) ahead of its
        task; fresh starts a new session told to read the branch first.

        A member a person stopped is not resumed: the answer is
        ``why: "stopped by a person"``. Lifting a stop is that person's next
        decision, so it is not a parameter here -- they run
        ``delfin-agent swarm wave resume --after-stop``."""
        return h_resume(run, member, note, fresh)

    @mcp.tool()
    def wave_stop(run: str, member: str = "") -> str:
        """End the members of a run, or one member. Worktrees and commits stay."""
        return h_stop(run, member)

    @mcp.tool()
    def wave_tick() -> str:
        """Run one pass now: start what became ready, advance the stages, record events,
        check the host and the members' listeners."""
        return h_tick()

    @mcp.tool()
    def wave_emergency_stop(why: str) -> str:
        """Stop every member of every run and hold the swarm until a person clears it."""
        return h_emergency_stop(why)

    @mcp.tool()
    def wave_clear_stop() -> str:
        """Lift an emergency stop. A person's decision: do not call it on your own."""
        return h_clear_stop()

    @mcp.tool()
    def wave_friction(runs: list | None = None, top: int = 30) -> str:
        """What the members ran into (refusals, retries, stalls, hints, crashes),
        grouped by reason and per model as rates per 100 tool calls."""
        return h_friction(runs, top)

    @mcp.tool()
    def wave_land(runs: list | None = None, member: str = "", push: bool = False) -> str:
        """Land members that wrote scratch/PR.md: branch checks, the tests they add,
        an independent review by a calibrated pair of another model family; with
        push=true the branch is pushed and its PR opened. Default: a dry run."""
        return h_land(runs, member, push)

    @mcp.tool()
    def wave_guard_install(every_minutes: int = 1, systemd: bool = False) -> str:
        """Run the tick every N minutes (idempotent): with systemd=true as a user service
        that survives a reboot, otherwise as an entry of the scheduler daemon."""
        return h_guard_install(every_minutes, systemd)

    @mcp.tool()
    def wave_guard_remove() -> str:
        """Remove the tick: the scheduler entry and the user service."""
        return h_guard_remove()

    mcp.run()


def main(argv: list[str] | None = None) -> int:
    run_server(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
