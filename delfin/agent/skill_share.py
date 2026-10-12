"""Share an approved skill as a local JSON bundle; import it as a PROPOSAL.

Package L6. The problem: what one user's DELFIN learned stays with that user.
This module lets an *accepted* skill (one that lives in ``~/.delfin/skills/
<name>/``) be exported to a shared folder on disk and imported by another
DELFIN as a PROPOSAL (source "shared") -- never active until a human accepts
it. Part (a) of the sharing plan only: no PR path, no network.

The export is opt-in per skill (every item is named explicitly), and the
bundle is one JSON file carrying ``{name, text, version, author, sha256}``
where ``sha256`` covers the canonical payload so any tamper is detected on
import. Nothing is silently rewritten:

* ``export_bundle`` runs the scrub BEFORE writing and BLOCKS -- naming the
  line -- when the text carries something machine-private or unsafe. It
  REUSES DELFIN's own ``skill_safety.check`` (the command + secret scan, the
  same one ``skill_proposals.propose`` runs) on top of a small targeted scan
  for an absolute path, a user@host, or a job id -- those shapes are not
  covered by the generic scan. A finding on either side blocks the export.

* ``import_bundle`` does a bounded read (regular file only -- no symlink, no
  FIFO -- under a size cap) and verifies the sha256 of the bundle payload,
  then runs the same scrub. A text that trips git-push / curl / an absolute
  path / a secret is RECORDED as status "blocked", naming the line -- that
  blocked outcome is the package's "stored blocked" result, never a live
  skill. A clean, hash-verified bundle returns status "verified" but does NOT
  create a proposal: it becomes one only when a human attaches REAL evidence
  via ``skill_proposals.propose(source="shared")``. No test/calc/job ref is
  ever invented, accept's evidence gate (skill_proposals.py:315) is
  untouched, and nothing is active without a human ``accept``.

Paths are resolved per call from ``Path.home()`` (never at import time), so a
redirected HOME redirects the live skills tree and the proposal store and no
test ever touches the real ``~/.delfin``.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import stat
from pathlib import Path

__all__ = ["export_bundle", "import_bundle"]

#: An import reads at most this many bytes; larger bundles are refused.
DEFAULT_MAX_BUNDLE_BYTES = 1_048_576  # 1 MiB
_MAX_BUNDLE_BYTES = DEFAULT_MAX_BUNDLE_BYTES


def _skills_root() -> Path:
    """The live skills tree, resolved per call (HOME may be redirected)."""
    return Path.home() / ".delfin" / "skills"


# --------------------------------------------------------------------------
# bundle building / verification
# --------------------------------------------------------------------------

def _payload(name: str, text: str, version: str, author: str) -> dict:
    """The canonical payload the sha256 covers. Export and import MUST use
    identical serialisation so a round-trip verifies byte-for-byte."""
    return {"name": name, "text": text, "version": version, "author": author}


def _sha_of(payload: dict) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, ensure_ascii=False).encode()
    ).hexdigest()


def _bundle_blob(payload: dict, sha: str) -> str:
    return json.dumps({"sha256": sha, **payload}, ensure_ascii=False)


# --------------------------------------------------------------------------
# the export scrub
# --------------------------------------------------------------------------

#: The machine-private shapes a shared skill text must not carry. These are
#: NOT covered by skill_safety.check (which guards commands + secrets while
#: a bare path, a user@host or a bare long-number job id are not commands).
#: A path is a token that BEGINS with "/" or "~/" (start of line or after
#: whitespace, a quote, a bracket, "=" or ":"), never a "/" inside a word.
#: That keeps "B3LYP/def2-TZVP", "and/or" and "kcal/mol" (method strings,
#: chemistry) exportable while "/home/alice/x" and "~/calc" (machine paths)
#: still block (operator F4).
_ABS_PATH_RE = re.compile(
    r"(?<![A-Za-z0-9_.~-])(?:~/|/)[A-Za-z0-9_@.~/-]*")
_USERHOST_RE = re.compile(r"(?<![A-Za-z0-9_.-])[A-Za-z0-9_.-]+@[A-Za-z0-9.-]+")
#: A job id is only that when it sits next to a scheduler word; a bare
#: 6+ digit number in a chemistry skill (e.g. a mass or a year) must not
#: block the export (operator F4).
_JOB_ID_RE = re.compile(
    r"(?i)(?:job|slurm|sbatch|squeue|jobid)\W{0,3}(\d{3,})")


def _machine_ref_snippets(text: str) -> list[tuple[int, str, str]]:
    """One ``(line_no, snippet, kind)`` per machine-private shape in text."""
    hits: list[tuple[int, str, str]] = []
    for i, line in enumerate(text.splitlines(), start=1):
        for pattern, kind in ((_ABS_PATH_RE, "absolute path"),
                              (_USERHOST_RE, "user@host"),
                              (_JOB_ID_RE, "job id")):
            m = pattern.search(line)
            if m:
                hits.append((i, m.group(0).strip()[:80], kind))
    return hits


def _export_scrub(text: str, field: str = "text") -> str:
    """Return a blocking reason naming the offending line, or "" when clean.

    Reuses ``skill_safety.check`` (commands + secrets DELFIN already scans)
    and adds the targeted machine-reference scan (absolute path / user@host
    / job id). ``field`` is the bundle field being scanned, so a hit names
    where the machine identity sits. Nothing is rewritten.
    """
    reasons: list[str] = []
    try:
        from . import skill_safety
        for finding in skill_safety.check(text):
            reasons.append(f"{field}: skill_safety: {finding}")
    except Exception as exc:
        # the reused scan must never pass what it could not read
        reasons.append(f"{field}: skill_safety unavailable: {exc!r}")
    for line_no, snippet, kind in _machine_ref_snippets(text):
        reasons.append(
            f"{kind} on line {line_no} of field {field!r}: "
            f"{snippet!r} must not leave this machine in a shared skill")
    return "; ".join(reasons)


def _scrub_leaving(payload: dict) -> str:
    """Block reason for every bundle field that leaves: text, version, author.

    The export's promise is that no machine identity -- an absolute path, a
    user/host name, a job id or a secret -- leaves in the shared bundle. The
    bundle carries name, text, version and author; ``name`` is validated at
    the door (no slash, no leading dot) and the rest are each scrubbed here,
    field-labelled, so a hit names the field and the line.
    """
    reasons: list[str] = []
    for field in ("text", "version", "author"):
        value = str(payload.get(field) or "")
        if not value:
            continue
        r = _export_scrub(value, field=field)
        if r:
            reasons.append(r)
    return "; ".join(reasons)


# --------------------------------------------------------------------------
# public API
# --------------------------------------------------------------------------

def export_bundle(name: str, target_dir: str | os.PathLike, *,
                  author: str = "") -> str:
    """Export one ACCEPTED skill into ``target_dir`` as a JSON bundle.

    Only an accepted skill (a live ``~/.delfin/skills/<name>/SKILL.md``)
    may be exported, and only by an explicit call -- there is no bulk path.
    The scrub runs BEFORE anything is written; a machine-private shape or a
    skill_safety finding BLOCKS the export and names the line.

    Returns exactly the skill text that leaves this machine.
    """
    name = (name or "").strip()
    if not name or "/" in name or name.startswith("."):
        raise ValueError(f"invalid skill name for export: {name!r}")
    skill_dir = _skills_root() / name
    text_path = skill_dir / "SKILL.md"
    if not text_path.is_file():
        raise ValueError(
            f"skill {name!r} is not an accepted skill: no live "
            f"~/.delfin/skills/{name}/SKILL.md")
    text = text_path.read_text(encoding="utf-8")

    version = ""
    meta_path = skill_dir / "proposal.json"
    if meta_path.is_file():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            version = str(meta.get("base_version") or "")
        except (OSError, ValueError, TypeError):
            version = ""

    payload = _payload(name, text, version, (author or "").strip())
    # scrub EVERY field that leaves (text, version, author), not just text:
    # an author label "alice@srv" or a version "/home/x/base" is also machine
    # identity and must not reach the shared bundle (s17/s29 finding).
    reason = _scrub_leaving(payload)
    if reason:
        raise ValueError(f"export of {name!r} blocked: {reason}")

    sha = _sha_of(payload)
    out = Path(target_dir)
    out.mkdir(parents=True, exist_ok=True)
    json_path = out / f"{name}.skill.json"
    # F3: never follow a planted symlink. The final path must not be a
    # symlink, and the temp is created with O_EXCL|O_NOFOLLOW (0600) so an
    # attacker cannot point it somewhere else and os.replace land outside.
    try:
        if os.path.islink(json_path):
            raise ValueError(
                f"refusing to overwrite a symlink at the export target: "
                f"{json_path}")
        # append .tmp to the full filename (with_suffix on a name with
        # embedded dots, e.g. <name>.skill.json, would shift the real dot
        # beyond the planted temp name -- F3 must hit exactly the temp).
        tmp = json_path.with_name(json_path.name + ".tmp")
        flags = (os.O_WRONLY | os.O_CREAT | os.O_EXCL
                 | getattr(os, "O_NOFOLLOW", 0))
        fd = os.open(str(tmp), flags, 0o600)
        try:
            st = os.fstat(fd)
            if not stat.S_ISREG(st.st_mode):
                raise ValueError(
                    f"refusing to write a non-regular export temp: {tmp}")
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                fh.write(_bundle_blob(payload, sha))
            fd = -1
        finally:
            if fd >= 0:
                os.close(fd)
        os.replace(tmp, json_path)
    except (FileExistsError, IsADirectoryError):
        raise ValueError(
            f"refusing to export: a file/dir named "
            f"{json_path.name}.tmp or the final "
            f"bundle already exists")
    return text


def import_bundle(path: str | os.PathLike, *, source: str = "shared") -> dict:
    """Read a shared skill bundle, verify it, and scrub it.

    import_bundle does NOT create a proposal (QS ruling, corrected). It is a
    bound + integrity + safety gate:

    * bounded read: regular file only (no symlink, no FIFO) under a strict
      size cap; the bundle sha256 must match the payload -- a tampered or
      truncated bundle raises and proposes nothing;
    * scrub: the same scan as export (skill_safety.check reused for
      commands/secrets + the targeted machine-private net). A text that
      trips git-push / curl / an absolute path / a secret is RECORDED here as
      status "blocked", naming the line -- that blocked outcome is the
      package's "stored blocked" result, never a live skill;
    * a clean, hash-verified bundle returns status "verified" (with the
      bundle fields) and becomes a proposal ONLY when a human attaches real
      evidence via ``skill_proposals.propose(source=...)`` -- no fabricated
      ref, accept's evidence gate untouched, nothing active without accept.

    Returns ``{"status", ...}`` where status is "verified" or "blocked".
    """
    # F2: decide on the OPEN descriptor, not on the path -- what was a file a
    # moment ago can be swapped for a FIFO or a symlink by the time it is
    # opened. Open once O_RDONLY|O_NOFOLLOW|O_NONBLOCK, fstat the descriptor,
    # check a bounded regular file, read at most cap+1 bytes from that fd.
    flags = (os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
             | getattr(os, "O_NONBLOCK", 0) | getattr(os, "O_CLOEXEC", 0))
    try:
        fd = os.open(path, flags)
    except OSError as exc:
        raise ValueError(f"cannot open bundle {path}: {exc}") from exc
    try:
        st = os.fstat(fd)
        if not stat.S_ISREG(st.st_mode):
            raise ValueError(f"not a regular file: {path}")
        if st.st_size > _MAX_BUNDLE_BYTES:
            raise ValueError(
                f"bundle too large: {st.st_size} bytes "
                f"(cap {_MAX_BUNDLE_BYTES})")
        blob = os.read(fd, _MAX_BUNDLE_BYTES + 1)
    finally:
        os.close(fd)
    if len(blob) > _MAX_BUNDLE_BYTES:
        raise ValueError(
            f"bundle too large (read cap {_MAX_BUNDLE_BYTES} bytes): {path}")
    try:
        text_data = blob.decode("utf-8")
        data = json.loads(text_data)
    except (OSError, ValueError) as exc:
        raise ValueError(f"invalid bundle JSON: {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError("invalid bundle: not a JSON object")

    name = str(data.get("name") or "")
    text = str(data.get("text") or "")
    version = str(data.get("version") or "")
    author = str(data.get("author") or "")
    stored_sha = str(data.get("sha256") or "")
    payload = _payload(name, text, version, author)
    if not stored_sha or _sha_of(payload) != stored_sha:
        raise ValueError(
            f"bundle hash mismatch: {path!r} was tampered with or truncated; "
            "nothing is proposed")

    # F4 (operator): the SAME machine-reference scan also runs on import.
    # A bundle hand-made to carry /home/<user>/ (or a user@host / job id) is
    # "blocked" here, naming the line, and never reaches a proposal. This is
    # separate from propose's own command+secret scan below, which stores
    # git-push / curl text as a "blocked" proposal. Both are never-active.
    net_reason = _scrub_leaving(payload)
    if net_reason:
        return {"status": "blocked", "source": source, "name": name,
                "version": version, "author": author,
                "reason": net_reason}

    # F1 (operator overrules the earlier QS ruling): an import ARRIVES AS A
    # PROPOSAL. propose() with ONE honest provenance evidence entry stores it
    # (never a claim of a test): kind "shared", ref the bundle sha256, detail
    # saying it is not verified on this machine. The kind matters: "recipe"
    # verifies in any workspace that has a verify recipe, whatever its ref,
    # so a bundle would have been acceptable with no check at all. Unsafe text (git-push /
    # curl / a secret) is STORED "blocked" by propose's own safety scan; clean
    # text is STORED "pending"; accept() still refuses at the evidence gate
    # until local evidence verifies (skill_proposals docstring :236 names this
    # pulled-proposal case). Nothing is ever active without accept.
    from . import skill_proposals
    from .skill_proposals import Evidence
    provenance = Evidence(
        kind="shared",
        ref=stored_sha,
        detail=f"imported bundle, author label {author!r}, "
               "not verified on this machine",
    )
    proposal = skill_proposals.propose(
        name,
        text,
        evidence=[provenance],
        source=source,
        # No base_version: with one, propose() treats the proposal as a PATCH
        # of the live skill of the same name, and accepting it would replace
        # this user's own skill. An import is always a new proposal under a
        # free name; the exporter's version stays in the returned record.
    )
    return {"status": proposal.status,
            "source": proposal.source,
            "name": proposal.name,
            "text": text,
            "version": version,
            "author": author,
            "reason": "; ".join(proposal.findings) if proposal.findings
            else None}
