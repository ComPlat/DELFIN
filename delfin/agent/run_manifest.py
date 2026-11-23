"""One artifact per run that says what produced a result.

A session report answers "what did this session do". It cannot answer
"could someone else arrive at this result again", because the facts that
decide that are not activity: which commit ran, which model answered
under which sampling settings, which prompt the model was given, which
engine binaries were on the machine, and whether the work the result
claims is present in the workspace at all.

Those facts exist in this tree already, in six places. This module
joins them into one file beside the session report and adds nothing it
can measure elsewhere:

    delfin.utils.get_git_commit_info      -> the commit, and whether it was dirty
    session_report.collect_session_report -> tools, files, denials, tokens, cost
    run_evidence.collect / .contract      -> the execution contract over the artifacts
    change_journal.list_changes           -> the event ledger of this session
    change_journal.sha256_file            -> the digest of each evidence-bearing file
    prompt_loader                         -> the pack the model actually read

What is NOT known is named. Every field this module could not resolve is
``None`` AND its name is listed in ``unresolved``, because a manifest
whose gaps look like values is worse than no manifest: it invites a
reader to believe a run was reproducible when the record cannot say so.
``complete()`` is therefore a question about the manifest, never about
the run.

Never raises. A manifest is written at the end of a run, after the
answer; a failure here must not change the run's outcome, so every
collector is guarded on its own and a collector that fails leaves its
field unresolved rather than taking the file down with it.
"""
from __future__ import annotations

import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Optional

#: A version in an engine's own words: "Program Version 6.1.1",
#: "xtb version 6.7.1 (edcfbbe)", "Version 3.0.2, Sun, 25 August".
#: Anchored on the word, so a banner rule cannot match.
_VERSION_RE = re.compile(r"[Vv]ersion\s+v?([0-9][0-9A-Za-z._+-]{0,30})")

#: The manifest's own version. A reader that does not know this number
#: must refuse to interpret the file rather than guess at a field.
SCHEMA = 1

#: At most this many evidence-bearing files are digested. A run that
#: writes thousands of files is not made more reproducible by listing
#: them all, and an unbounded list would make the manifest the largest
#: artifact of the run. When the cap bites it is reported, not hidden.
_MAX_FILES = 200

#: Files worth digesting: the inputs a calculation was run from and the
#: outputs a claim is read out of. Deliberately not "every file" --
#: scratch, integrals and temporary directories are not evidence and
#: they differ between two correct runs of the same work.
_EVIDENCE_SUFFIXES = (
    ".inp", ".out", ".xyz", ".gbw", ".hess", ".log", ".json", ".md",
    ".csv", ".property.txt",
)

#: Engines whose version belongs in the record: a result is only
#: comparable against another run of the same program.
_ENGINES = ("orca", "xtb", "crest")


def _unresolved(bag: list[str], name: str) -> None:
    if name not in bag:
        bag.append(name)


def _code(unresolved: list[str]) -> dict[str, Any]:
    """The commit that ran, and the version that names it."""
    out: dict[str, Any] = {"delfin_version": None, "commit": None,
                           "dirty": None}
    try:
        from delfin import __version__
        out["delfin_version"] = str(__version__)
    except Exception:
        _unresolved(unresolved, "code.delfin_version")
    try:
        from delfin.utils import get_git_commit_info
        raw = get_git_commit_info()
        if raw:
            out["commit"] = str(raw).replace("-dirty", "")
            out["dirty"] = str(raw).endswith("-dirty")
        else:
            _unresolved(unresolved, "code.commit")
    except Exception:
        _unresolved(unresolved, "code.commit")
    return out


def _environment(unresolved: list[str]) -> dict[str, Any]:
    """The machine, and the engine versions a result depends on."""
    out: dict[str, Any] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "engines": {},
    }
    for name in _ENGINES:
        found = shutil.which(name)
        if not found:
            out["engines"][name] = None
            continue
        out["engines"][name] = {"path": found, "version": None}
        # The version is read from the program, not guessed from the
        # path: two installations under one name are the ordinary case
        # on a cluster.
        #
        # Matched out of the whole output, not taken as its first line.
        # Measured 2026-11-23 on a real run: all three of these print a
        # banner first, so the first line gave "*****************" for
        # ORCA and a box-drawing rule for CREST. A version field holding
        # decoration is worse than an empty one -- it reads as a version
        # to anyone comparing two manifests.
        for args in ([found, "--version"], [found, "-v"]):
            try:
                res = subprocess.run(args, capture_output=True, text=True,
                                     timeout=10, check=False)
            except Exception:
                continue
            found_version = _VERSION_RE.search(
                (res.stdout or "") + (res.stderr or ""))
            if found_version:
                out["engines"][name]["version"] = found_version.group(1)
                break
        if out["engines"][name]["version"] is None:
            _unresolved(unresolved, f"environment.engines.{name}.version")
    return out


def _prompt(unresolved: list[str], role: str) -> dict[str, Any]:
    """A digest of the pack the model read.

    The digest covers the pack's markdown as it sits on disk, in sorted
    path order, so a changed instruction changes the manifest. It does
    not cover the BUILT prompt: that carries the session's own memory
    and workspace, which differ between two runs of the same configured
    agent and would make every manifest unique for no gain.
    """
    out: dict[str, Any] = {"role": role or None, "pack_digest": None,
                           "pack_files": None}
    if not str(role or "").strip():
        # Which role prompt the model read decides what it was told it
        # may do, so an unknown role is a gap in the record and says so
        # rather than sitting there as a null.
        _unresolved(unresolved, "prompt.role")
    try:
        from .change_journal import sha256_text
        from .prompt_loader import PromptLoader
        pack = PromptLoader().agent_dir
        files = sorted(p for p in pack.rglob("*.md") if p.is_file())
        if not files:
            _unresolved(unresolved, "prompt.pack_digest")
            return out
        joined = "\n".join(
            f"{p.relative_to(pack).as_posix()}:"
            f"{sha256_text(p.read_text(encoding='utf-8', errors='replace'))}"
            for p in files)
        out["pack_digest"] = sha256_text(joined)
        out["pack_files"] = len(files)
    except Exception:
        _unresolved(unresolved, "prompt.pack_digest")
    return out


def _model(unresolved: list[str], report: Any,
           settings: Optional[dict]) -> dict[str, Any]:
    """Which model answered, and under which sampling settings.

    The id is taken from the session's own per-turn metrics first, not
    from configuration: a run that switched model mid-session did not
    run under the configured one, and the record should say which model
    actually answered.

    A short run may have no per-turn metrics at all, and a manifest that
    then names no model is useless where it is needed most. So the
    caller's configured id is the fallback, and ``id_source`` says which
    of the two this is -- the distinction matters to a reader comparing
    two runs, and collapsing it would make a configured id look like a
    measured one.

    The sampling settings are the caller's to pass: nothing in the
    session record keeps them, so an absent temperature is reported
    absent rather than assumed.
    """
    out: dict[str, Any] = {"id": None, "id_source": None, "provider": None,
                           "temperature": None, "seed": None,
                           "switches": None}
    measured = str(getattr(report, "model", "") or "")
    configured = str((settings or {}).get("model") or "")
    if measured:
        out["id"], out["id_source"] = measured, "session_metrics"
    elif configured:
        out["id"], out["id_source"] = configured, "configured"
    else:
        _unresolved(unresolved, "model.id")
    for key in ("provider", "temperature", "seed"):
        value = (settings or {}).get(key)
        if value is None:
            _unresolved(unresolved, f"model.{key}")
        else:
            out[key] = value
    switches = (settings or {}).get("model_switches")
    out["switches"] = list(switches) if switches else []
    return out


def _activity(unresolved: list[str], report: Any) -> dict[str, Any]:
    """What the session did, as the session report already measured it."""
    if report is None:
        _unresolved(unresolved, "activity")
        return {}
    return {
        "started_at": getattr(report, "started_at", 0.0) or None,
        "ended_at": getattr(report, "ended_at", 0.0) or None,
        "tool_calls": len(getattr(report, "tool_calls", []) or []),
        "files_changed": len(getattr(report, "files_changed", []) or []),
        "commands_run": len(getattr(report, "commands_run", []) or []),
        "tests_run": len(getattr(report, "tests_run", []) or []),
        "denials": len(getattr(report, "denials", []) or []),
        "input_tokens": int(getattr(report, "input_tokens", 0) or 0),
        "output_tokens": int(getattr(report, "output_tokens", 0) or 0),
        "cost_usd": float(getattr(report, "cost_usd", 0.0) or 0.0),
    }


def _ledger(unresolved: list[str], session_id: str) -> dict[str, Any]:
    """The event ledger's own index range for this session.

    A reader checking a replay compares event identity and order, so the
    range is what it needs; the records themselves stay in the journal,
    which is the one place that can restore a file.
    """
    out: dict[str, Any] = {"changes": None, "first_seq": None,
                           "last_seq": None}
    try:
        from .change_journal import list_changes
        rows = list_changes(session_id) or []
        seqs = [int(r.get("seq")) for r in rows
                if str(r.get("seq", "")).isdigit()]
        out["changes"] = len(rows)
        out["first_seq"] = min(seqs) if seqs else None
        out["last_seq"] = max(seqs) if seqs else None
    except Exception:
        _unresolved(unresolved, "ledger")
    return out


def _contract(unresolved: list[str], workspace: Optional[Path],
              expected_steps: Optional[int],
              window: Optional[tuple[float, float]] = None,
              ) -> dict[str, Any]:
    """The execution contract over the files the run left behind.

    This is the manifest's only judgement, and it is the one the field
    asks for: an oracle that reads the artifacts rather than the
    model's account of them. An unjudged flag is not a pass, so the
    flags are carried as they are -- True, False or None.

    ``window`` is the run's own ``(started_at, ended_at)``, already
    established above from the session's activity. The egress flag is
    judged from the audit log over exactly that window: without it the
    flag was None for EVERY run DELFIN has ever made, because nothing
    in the product filled ``egress_events`` -- the parameter existed and
    had no caller. With no window the flag stays None, which is the
    honest answer and not a pass.
    """
    out: dict[str, Any] = {"flags": None, "certified": None, "failed": None,
                           "unjudged": None, "problems": None,
                           "steps_found": None,
                           "expected_steps": expected_steps}
    if workspace is None:
        _unresolved(unresolved, "contract")
        return out
    try:
        from . import run_evidence
        ev = run_evidence.collect(workspace)
        egress = None
        if window and window[0] and window[1] and window[1] >= window[0]:
            # Returns None wherever absence cannot be demonstrated (no
            # log, or a log that begins after the run did), so a rotated
            # log leaves the flag unjudged instead of quiet.
            egress = run_evidence.egress_in_window(window[0], window[1])
        con = run_evidence.contract(ev, expected_steps=expected_steps,
                                    egress_events=egress)
        out.update({
            "flags": dict(con.flags),
            "certified": con.certified(),
            "failed": con.failed(),
            "unjudged": con.unjudged(),
            "problems": list(con.problems),
            "steps_found": len(ev.steps),
        })
    except Exception:
        _unresolved(unresolved, "contract")
    return out


def _workspace(unresolved: list[str],
               workspace: Optional[Path]) -> dict[str, Any]:
    """A digest per evidence-bearing file, so a claim can be traced.

    Bounded and sorted. When the cap bites, the count that was left out
    is reported: a truncated list that does not say it was truncated
    reads as a complete one.
    """
    out: dict[str, Any] = {"root": None, "files": None, "omitted": None}
    if workspace is None:
        _unresolved(unresolved, "workspace")
        return out
    try:
        from .change_journal import sha256_file
        root = Path(workspace)
        out["root"] = str(root)
        found = sorted(
            p for p in root.rglob("*")
            if p.is_file() and not p.is_symlink()
            and p.name.endswith(_EVIDENCE_SUFFIXES)
            and ".delfin-trash" not in p.parts)
        kept = found[:_MAX_FILES]
        out["files"] = [
            {"path": p.relative_to(root).as_posix(),
             "sha256": sha256_file(p),
             "bytes": p.stat().st_size}
            for p in kept]
        out["omitted"] = len(found) - len(kept)
    except Exception:
        _unresolved(unresolved, "workspace")
    return out


def build(*, session_id: str = "", workspace: Any = None,
          expected_steps: Optional[int] = None,
          termination: str = "", settings: Optional[dict] = None,
          role: str = "") -> dict[str, Any]:
    """Join what this tree already measures into one record.

    Input: the session id, the workspace a result was produced in, and
    the settings only the caller knows (provider, temperature, seed,
    model switches). Output: the manifest, with every field it could not
    resolve named in ``unresolved``.

    ``termination`` is the run's own stop reason, passed in rather than
    inferred: the record must say whether a run ended because it was
    finished, because it hit a ceiling, or because something stopped it,
    and only the caller knows which.
    """
    unresolved: list[str] = []
    sid = str(session_id or "").strip()
    if not sid:
        _unresolved(unresolved, "run.session_id")
    ws: Optional[Path] = None
    if workspace is not None:
        try:
            candidate = Path(workspace)
            ws = candidate if candidate.is_dir() else None
            if ws is None:
                _unresolved(unresolved, "workspace.root")
        except Exception:
            _unresolved(unresolved, "workspace.root")

    report = None
    try:
        from .session_report import collect_session_report
        report = collect_session_report(sid) if sid else None
    except Exception:
        report = None

    activity = _activity(unresolved, report)
    started = (activity or {}).get("started_at")
    ended = (activity or {}).get("ended_at")
    manifest: dict[str, Any] = {
        "schema": SCHEMA,
        "written_at": time.time(),
        "run": {
            "session_id": sid or None,
            "termination": str(termination or "") or None,
            "wall_clock_s": (round(ended - started, 1)
                             if started and ended and ended >= started
                             else None),
            "label": os.environ.get("DELFIN_AGENT_LABEL") or None,
        },
        "code": _code(unresolved),
        "model": _model(unresolved, report, settings),
        "prompt": _prompt(unresolved, role),
        "environment": _environment(unresolved),
        "activity": activity,
        "ledger": _ledger(unresolved, sid) if sid else {},
        "contract": _contract(unresolved, ws, expected_steps,
                              window=((started, ended)
                                      if started and ended else None)),
        "workspace": _workspace(unresolved, ws),
    }
    if not str(termination or "").strip():
        _unresolved(unresolved, "run.termination")
    if manifest["run"]["wall_clock_s"] is None:
        _unresolved(unresolved, "run.wall_clock_s")
    manifest["unresolved"] = sorted(unresolved)
    return manifest


#: What a reader needs before it may call a run reconstructable. Named
#: here rather than left to each reader's judgement, so "complete" means
#: one thing across the tree.
REQUIRED = (
    "code.commit",
    "model.id",
    "prompt.pack_digest",
    "run.termination",
    "contract",
)


def missing(manifest: dict) -> list[str]:
    """The REQUIRED fields this manifest does not carry.

    Reads the manifest's own ``unresolved`` list and the values, so a
    field that was written as None by a caller that never knew it is
    reported as missing even if its collector did not complain.
    """
    gaps: list[str] = []
    unresolved = set(manifest.get("unresolved") or ())
    for name in REQUIRED:
        if name in unresolved:
            gaps.append(name)
            continue
        node: Any = manifest
        for part in name.split("."):
            node = (node or {}).get(part) if isinstance(node, dict) else None
        if node is None or node == {}:
            gaps.append(name)
    return gaps


def complete(manifest: dict) -> bool:
    """Whether the MANIFEST is complete -- never whether the run passed.

    A complete manifest may record a failed contract; an incomplete one
    says nothing either way, which is the distinction the name protects.
    """
    return not missing(manifest)


def write(manifest: dict, path: Any) -> Optional[Path]:
    """Write the manifest where a reader will look, atomically.

    Returns the path, or None when it could not be written: a manifest
    is a record of a run that already finished, and failing to store it
    must not change that run's outcome.
    """
    try:
        target = Path(path)
        target.parent.mkdir(parents=True, mode=0o700, exist_ok=True)
        tmp = target.with_suffix(target.suffix + ".tmp")
        tmp.write_text(json.dumps(manifest, indent=1, sort_keys=True,
                                  ensure_ascii=False) + "\n",
                       encoding="utf-8")
        # Narrowed before it takes the final name: a chmod afterwards
        # leaves a window in which the record is readable by others, and
        # this file names every evidence file of the run.
        try:
            os.chmod(tmp, 0o600)
        except OSError:
            pass
        tmp.replace(target)
        return target
    except Exception:
        return None


def read(path: Any) -> Optional[dict]:
    """Read a manifest, or None. A file of a schema this code does not
    know is refused rather than interpreted: an unknown schema may have
    moved a field this reader would then read as absent."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return None
    if not isinstance(data, dict) or data.get("schema") != SCHEMA:
        return None
    return data


def path_for(session_id: str) -> Optional[Path]:
    """Where this session's manifest belongs: beside its report."""
    try:
        from .session_report import _report_dir, _safe_session_id
        sid = str(session_id or "").strip()
        if not sid:
            return None
        return _report_dir() / (_safe_session_id(sid) + ".manifest.json")
    except Exception:
        return None


def write_run_manifest(session_id: str, *, workspace: Any = None,
                       expected_steps: Optional[int] = None,
                       termination: str = "",
                       settings: Optional[dict] = None,
                       role: str = "") -> Optional[Path]:
    """Build and store one run's manifest; never raises.

    The session-end call. Returns the written path, or None when there
    was nothing to write or it could not be written.
    """
    try:
        sid = str(session_id or "").strip()
        if not sid:
            return None
        manifest = build(session_id=sid, workspace=workspace,
                         expected_steps=expected_steps,
                         termination=termination, settings=settings,
                         role=role)
        target = path_for(sid)
        return write(manifest, target) if target else None
    except Exception:
        return None
