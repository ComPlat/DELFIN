"""Team skill registry — share chemistry skills over the remote archive.

ClawHub-style CONCEPT (own code, no external code): the team's shared
skills live in ``<transfer.remote_path>/AGENT_SKILLS/`` on the same SSH
remote the dashboard already uses for transfers and bug reports.

- ``publish_skill(name, evidence=[...])`` pushes a local skill file there
  (rsync/SSH via the proven ``ssh_transfer_jobs`` builders) — but only
  with an evidence record that re-verifies (learning wave): no valid
  evidence, no publication; domain "chemistry" additionally needs a
  calc-or-test evidence entry.
- ``pull_skills()`` downloads shared skills as PROPOSALS (source
  "team", status "pending") via the proposal pipeline (learning wave:
  a pulled skill goes through the same safety check and human
  acceptance as a self-written one; it never lands directly in the
  active ``~/.delfin/skills/`` directory).

Safety: pulls NEVER overwrite an existing local skill — on conflict the
incoming file is written as ``<name>-shared.md`` (user no-delete rule).
Skills are plain markdown (no code execution risk). No paths hardcoded —
host/user/remote_path come from the transfer settings.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Callable, Optional


_SKILLS_SUBDIR = "AGENT_SKILLS"
_LOCAL_SKILLS_DIR = Path.home() / ".delfin" / "skills"
# Built-in pack skills are publishable too (they are the curated set).
_PACK_SKILLS_DIR = Path(__file__).resolve().parent / "pack" / "skills"

RunFn = Callable[[list[str]], tuple[int, str, str]]


def _default_run(cmd: list[str]) -> tuple[int, str, str]:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        return out.returncode, out.stdout, out.stderr
    except Exception as exc:
        return 1, "", str(exc)


def find_skill_file(name: str, workspace: str | Path | None = None) -> Optional[Path]:
    """Resolve a skill name to its markdown file (user/project/pack)."""
    name = (name or "").strip().lstrip("/")
    if not name:
        return None
    try:
        from delfin.agent.skills import discover_skills
        for sk in discover_skills(workspace):
            if sk.name == name and getattr(sk, "source", None):
                return Path(sk.source)
    except Exception:
        pass
    pack = _PACK_SKILLS_DIR / f"{name}.md"
    if pack.is_file():
        return pack
    return None


def _remote_dir(remote_path: str) -> str:
    return f"{(remote_path or '').rstrip('/')}/{_SKILLS_SUBDIR}"


# ---------------------------------------------------------------------------
# Evidence gate (learning wave, package 8): the team archive only ever
# receives a skill whose evidence record re-verifies NOW. Reuses
# delfin.agent.evidence (which itself reuses evidence_freshness,
# result_critic, verify_recipe, the throttled job listing) -- nothing
# is re-implemented here.
# ---------------------------------------------------------------------------

def _gate_evidence(evidence, *, domain=None, workspace=None, runs=None,
                   list_jobs=None) -> tuple[bool, str]:
    """Check *evidence* (list of entries) for publication.

    Rules: at least one entry must re-verify as valid; for domain
    "chemistry" at least one VALID entry must be of kind "calc" or
    "test". Refusal messages are English and name the reason(s).
    """
    from delfin.agent.evidence import verify_evidence
    entries = list(evidence or [])
    if not entries:
        return False, ("no evidence record given: a skill may only be "
                       "published with evidence that re-verifies")
    valid: list[str] = []
    invalid: list[str] = []
    valid_kinds: set[str] = set()
    for ev in entries:
        ok, detail = verify_evidence(ev, workspace=workspace, runs=runs,
                                     list_jobs=list_jobs)
        ref = str((ev.get("ref") if isinstance(ev, dict)
                   else getattr(ev, "ref", "")) or "")
        if ok:
            valid.append(ref)
            kind = str((ev.get("kind") if isinstance(ev, dict)
                        else getattr(ev, "kind", "")) or "")
            valid_kinds.add(kind.lower())
        else:
            invalid.append(f"{ref}: {detail}")
    if not valid:
        return False, ("evidence does not verify: "
                       + "; ".join(invalid))
    if str(domain or "").strip().lower() == "chemistry" \
            and not (valid_kinds & {"calc", "test"}):
        return False, ("chemistry skills need at least one valid 'calc' or "
                       "'test' evidence entry")
    return True, f"evidence verified: {', '.join(valid)}"


def _propose_pulled_skill(name: str, text: str) -> tuple[str, str]:
    """Route one downloaded team skill into the PROPOSAL pipeline.

    Pulled skills are proposals (status "pending", source "team"),
    never directly active -- they still go through the safety check and
    human acceptance. Delegates to skill_proposals (package 1) when that
    module is available; without it the skill is NOT installed and the
    status says so. The monkeypatch point for tests.
    """
    try:
        from delfin.agent import skill_proposals as sp
        proposal = sp.propose(name, text, evidence=[], source="team")
        return "proposed", str(getattr(proposal, "name", name))
    except ImportError:
        return "proposed", (f"{name} (proposal pipeline not yet available; "
                            f"kept out of the active skill dir)")
    except Exception as exc:
        return "proposed-failed", f"{name}: {exc}"


def publish_skill(
    name: str,
    *,
    host: str,
    user: str,
    remote_path: str,
    port: int = 22,
    workspace: str | Path | None = None,
    evidence: Optional[list] = None,
    domain: str = "",
    runs: Optional[list] = None,
    list_jobs: Optional[Callable] = None,
    run_fn: RunFn = _default_run,
) -> tuple[bool, str]:
    """Push one local skill to the shared registry. Returns (ok, message).

    Evidence gate (learning wave): *evidence* is a list of evidence
    entries (dicts or objects with kind/ref); each is re-verified via
    delfin.agent.evidence. Refused -- with the reason -- when no entry
    verifies, or when *domain* is "chemistry" and no valid entry is of
    kind "calc" or "test". Refusal happens BEFORE any command is built.
    """
    host, user, remote_path = (host or "").strip(), (user or "").strip(), (remote_path or "").strip()
    if not (host and user and remote_path):
        return False, "no remote configured (host/user/remote path missing)"
    ok, msg = _gate_evidence(evidence, domain=domain, workspace=workspace,
                             runs=runs, list_jobs=list_jobs)
    if not ok:
        return False, msg
    src = find_skill_file(name, workspace)
    if src is None:
        return False, f"skill '{name}' not found (neither local nor in the pack)"

    from delfin.ssh_transfer_jobs import (
        build_ssh_mkdir_command,
        build_rsync_transfer_command,
    )
    target = _remote_dir(remote_path)
    rc, _, err = run_fn(build_ssh_mkdir_command(host, user, target, port))
    if rc != 0:
        return False, f"remote mkdir failed: {err.strip()[:200]}"
    rc, _, err = run_fn(
        build_rsync_transfer_command([str(src)], host, user, target, port))
    if rc != 0:
        return False, f"rsync failed: {err.strip()[:200]}"
    return True, f"{target}/{src.name}"


def list_shared(
    *,
    host: str,
    user: str,
    remote_path: str,
    port: int = 22,
    run_fn: RunFn = _default_run,
) -> tuple[bool, list[str]]:
    """List skill files in the shared registry."""
    from delfin.ssh_transfer_jobs import build_ssh_command_base
    cmd = build_ssh_command_base(host, user, port)
    cmd.append(f"ls -1 -- {_remote_dir(remote_path)} 2>/dev/null")
    rc, out, _ = run_fn(cmd)
    if rc != 0:
        return False, []
    return True, [l.strip() for l in out.splitlines()
                  if l.strip().endswith(".md")]


def install_skill_text(name: str, text: str,
                       dest_dir: Path | None = None) -> tuple[str, Path]:
    """Install skill content locally WITHOUT overwriting.

    Returns (status, path): status in {installed, identical, conflict}.
    On conflict the file is written as ``<name>-shared.md`` instead.
    """
    d = dest_dir or _LOCAL_SKILLS_DIR
    d.mkdir(parents=True, exist_ok=True)
    target = d / f"{name}.md"
    if not target.exists():
        target.write_text(text, encoding="utf-8")
        return "installed", target
    if target.read_text(encoding="utf-8") == text:
        return "identical", target
    alt = d / f"{name}-shared.md"
    alt.write_text(text, encoding="utf-8")
    return "conflict", alt


def pull_skills(
    name: str = "",
    *,
    host: str,
    user: str,
    remote_path: str,
    port: int = 22,
    dest_dir: Path | None = None,
    run_fn: RunFn = _default_run,
) -> tuple[bool, list[tuple[str, str]]]:
    """Download shared skills into the local user skill dir.

    ``name`` empty = all shared skills. Returns (ok, [(skill, status)]),
    status from :func:`install_skill_text`. Never overwrites local files.
    """
    host, user, remote_path = (host or "").strip(), (user or "").strip(), (remote_path or "").strip()
    if not (host and user and remote_path):
        return False, [("", "no remote configured")]

    from delfin.ssh_transfer_jobs import build_rsync_download_command
    with tempfile.TemporaryDirectory(prefix="delfin_skills_") as tmp:
        sources = [f"{_SKILLS_SUBDIR}/{name}.md" if name
                   else f"{_SKILLS_SUBDIR}/"]
        rc, _, err = run_fn(build_rsync_download_command(
            sources, host, user, remote_path, port, tmp))
        if rc != 0:
            return False, [("", f"rsync failed: {err.strip()[:200]}")]
        results: list[tuple[str, str]] = []
        for f in sorted(Path(tmp).rglob("*.md")):
            status, info = _propose_pulled_skill(
                f.stem, f.read_text(encoding="utf-8"))
            results.append((f.stem, status))
        if not results:
            return True, [("", "no skills found in the remote registry")]
        return True, results
