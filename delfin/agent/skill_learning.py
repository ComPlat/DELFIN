"""Learn a reusable skill from one finished DELFIN session.

Hermes-style trigger: a session that earned a skill proposal produces at
most ONE proposal — and only when it also carries evidence (a green test
run, a verified calculation/job, or a verify_recipe). A session without
evidence never becomes a rule; an unsafe success stays a suggestion at
most (skill_proposals.propose runs the safety check).

Message forms handled here are the ones DELFIN transcripts actually
contain: ``chat_messages`` dicts with roles user/system/thinking/assistant
and ``tool`` messages whose content embeds an HTML tool chip
(``<span class="tool-name">NAME</span>``) and a result block
(``<details><summary> -> {json}</details>``). Parsing is tolerant: any
message dict with a ``name``/``tool``/``tool_name`` key or an embedded
tool chip counts as a tool call.
"""

from __future__ import annotations

import re

_TOOL_CHIP = re.compile(r'tool-name">([^<]+)<')
_RESULT_BLOCK = re.compile(r"<summary>\s*(?:&rarr;|->|→)\s*(.*?)</details>",
                           re.S)

# Tool calls inside one working stretch (between user turns) that mark
# deep, non-obvious work.
_DEEP_STRETCH_CALLS = 5


def _content(msg: dict) -> str:
    """Best-effort text of one message; never raises."""
    try:
        c = msg.get("content")
        return c if isinstance(c, str) else str(c or "")
    except Exception:
        return ""


def _tool_name(msg: dict) -> str | None:
    """Name of the tool this message represents, or None."""
    for key in ("name", "tool", "tool_name"):
        v = msg.get(key)
        if isinstance(v, str) and v:
            return v
    if msg.get("role") in ("tool", "function"):
        m = _TOOL_CHIP.search(_content(msg))
        if m:
            return m.group(1)
        return "tool"
    m = _TOOL_CHIP.search(_content(msg))
    if m and msg.get("role") in (None, "tool"):
        return m.group(1)
    return None


def _result_text(msg: dict) -> str:
    m = _RESULT_BLOCK.search(_content(msg))
    return m.group(1) if m else ""


def _is_tool_failed(msg: dict, name: str | None) -> bool:
    """Did this tool result report failure? Tolerant, best-effort."""
    text = (_content(msg) + " " + _result_text(msg)).lower()
    m = re.search(r'"exit_code":\s*(-?\d+)', text)
    if m and m.group(1) not in ("0", ""):
        return True
    if re.search(r'"(status|state)":\s*"(error|failed|denied)"', text):
        return True
    if re.search(r'\b(failed|error)\b', text) and "passed" not in text:
        # "1 failed" style summaries; a summary that also says "passed"
        # (e.g. "12 passed, 1 failed") is handled above only via counts —
        # keep it conservative: count-based check first.
        if not re.search(r'"failed":\s*0\b', text):
            return True
    if "traceback" in text or "modulenotfounderror" in text:
        return True
    if name and name.lower() in ("read", "edit", "write") and "refused" in text:
        return True
    return False


def _is_tool_green(msg: dict, name: str | None) -> bool:
    """Did this tool result report clear success?"""
    text = (_content(msg) + " " + _result_text(msg)).lower()
    m = re.search(r'"exit_code":\s*(\d+)', text)
    if m and m.group(1) == "0":
        return True
    if re.search(r'"failed":\s*0\b.*"passed":\s*[1-9]', text):
        return True
    if re.search(r'"(?:passed|updated|created|written)"\s*:\s*[1-9]', text):
        return True
    if name and name.lower() == "toolrunner" and '"failed": 0' in text:
        return True
    return False


def _looks_like_correction(text: str) -> bool:
    """Does a user message read as a course correction?"""
    t = text.lower()
    markers = (
        "no, ", "nein", "falsch", "wrong", "not what i meant",
        "i meant", "stattdessen", "instead", "revert", "undo",
        "please revert", "undo that", "rollback", "that's not",
        "das ist falsch", "das war falsch", "ich meinte",
        "nicht so", "anders", "stop, ", "halt,",
    )
    return any(m in t for m in markers)


def qualifies(messages) -> list[str]:
    """Named reasons why this session earned a skill proposal.

    Empty list means: no proposal. One entry per trigger, each a short
    named English reason. Never raises on odd input.
    """
    try:
        msgs = [m for m in (messages or []) if isinstance(m, dict)]
    except Exception:
        return []
    if len(msgs) < 3:
        return []

    reasons: list[str] = []
    user_msgs = [m for m in msgs if m.get("role") == "user"]

    # 1) Deep work: >= 5 tool calls in one working stretch (between two
    #    user turns, or from the first user turn to the end).
    stretches: list[list[dict]] = [[]]
    for m in msgs:
        if m.get("role") == "user" and user_msgs and m is not user_msgs[0]:
            stretches.append([])
        stretches[-1].append(m)
    for stretch in stretches:
        n = sum(1 for m in stretch if _tool_name(m))
        if n >= _DEEP_STRETCH_CALLS:
            reasons.append(
                f"deep-work: {n} tool calls in one working stretch")
            break

    # 2) Error followed by recovery: a tool result reporting failure,
    #    then a later tool result of the same kind going green.
    seen_failed: set[str] = set()
    for m in msgs:
        name = _tool_name(m)
        if not name:
            continue
        key = name.lower()
        if _is_tool_failed(m, name):
            seen_failed.add(key)
        elif key in seen_failed and _is_tool_green(m, name):
            reasons.append(
                f"error-recovery: a {name} call failed, then a later "
                f"{name} call succeeded")
            seen_failed.discard(key)

    # 3) A course correction by the user.
    for m in user_msgs:
        if _looks_like_correction(_content(m)):
            reasons.append("user-correction: the user corrected the "
                           "agent's course mid-session")
            break

    # 4) A non-obvious route: several different tools used (>= 4 distinct)
    #    hints the solution needed exploration beyond the obvious path.
    distinct = {(_tool_name(m) or "").lower() for m in msgs} - {""}
    if len(distinct) >= 4:
        reasons.append(
            f"non-obvious-route: {len(distinct)} distinct tools used "
            f"to reach the result")

    return reasons


# ---------------------------------------------------------------------------
# Evidence extraction (Phase 2)
# ---------------------------------------------------------------------------

# Ranking: a green test run is the strongest proof, then a verified
# calculation/job, then a verify_recipe rendered in text.
_KIND_RANK = {"test": 0, "calc": 1, "job": 2, "recipe": 3}

_TEST_REF = re.compile(r"(tests/[A-Za-z0-9_./]+\.py(?:::[A-Za-z0-9_]+)?)")
_JOB_ID = re.compile(r'"job_id":\s*"?(\d+)"?')


def _evidence_class():
    """The contract's Evidence type when package 1 is merged in, else a
    local duck-typed stand-in with the same fields."""
    try:
        from delfin.agent.skill_proposals import Evidence  # type: ignore
        return Evidence
    except Exception:
        from dataclasses import dataclass

        @dataclass
        class Evidence:  # type: ignore[no-redef]
            kind: str
            ref: str
            detail: str = ""
            verified_at: str = ""
            runs: list = None  # type: ignore[assignment]

        return Evidence


def _all_text(msg: dict) -> str:
    return _content(msg) + "\n" + _result_text(msg)


def _green_test_result(text: str) -> bool:
    """A pytest-style summary that is fully green."""
    if re.search(r'"failed":\s*[1-9]', text) or re.search(
            r'"errors":\s*[1-9]', text):
        return False
    if re.search(r'"passed":\s*[1-9]', text):
        return True
    # "7 passed in 2.4s", "12 passed" — but not "2 failed, 5 passed".
    m = re.search(r"(\d+)\s+passed", text)
    if m and int(m.group(1)) >= 1:
        return not re.search(r"\d+\s+failed", text)
    return False


def _test_refs(text: str) -> list[str]:
    return _TEST_REF.findall(text)


def extract_evidence(messages):
    """The strongest verifiable proof in this session, or None.

    Kinds (contract): "test" — a green pytest/Gate run with test-IDs;
    "calc"/"job" — a submitted/verified calculation or job; "recipe" —
    a verify_recipe rendering in assistant text. No evidence means the
    session NEVER becomes a skill proposal ("no proof, no proposal").
    """
    try:
        msgs = [m for m in (messages or []) if isinstance(m, dict)]
    except Exception:
        return None
    Evidence = _evidence_class()
    best = None
    best_rank = 99

    def _offer(kind, ref, detail, at):
        nonlocal best, best_rank
        rank = _KIND_RANK.get(kind, 9)
        if rank < best_rank:
            best = Evidence(kind=kind, ref=ref, detail=detail,
                            verified_at=at)
            best_rank = rank

    import datetime as _dt
    now = _dt.datetime.now(_dt.timezone.utc).isoformat(
        timespec="seconds")

    for msg in msgs:
        text = _all_text(msg)
        name = (_tool_name(msg) or "").lower()

        is_tool_msg = msg.get("role") == "tool" or bool(name)

        if is_tool_msg and _green_test_result(text):
            refs = _test_refs(text) or _test_refs(_content(msg))
            if refs:
                _offer("test", refs[0],
                       f"green run covering {len(refs)} test id(s)", now)
                continue

        if is_tool_msg and (name in ("submit_calculation", "run_calculation",
                                    "submit_application")
                            or '"job_id"' in text):
            m = _JOB_ID.search(text)
            if m:
                _offer("job", f"job:{m.group(1)}",
                       "calculation/job submitted through the sanctioned "
                       "tools", now)
                continue

        if msg.get("role") == "assistant" and "gate" in text \
                and _green_test_result(text):
            refs = _test_refs(text)
            if refs:
                _offer("recipe", refs[0],
                       "verify_recipe rendering in assistant text", now)

    return best


# ---------------------------------------------------------------------------
# learn_from_session (Phase 3): draft one skill, propose it once
# ---------------------------------------------------------------------------

_DRAFT_SYSTEM = """You draft DELFIN skills. A skill is a short, reusable
playbook in Markdown with a YAML front-matter (name, description)
followed by numbered steps. Rules:
- English only; every string a DELFIN user could see is English.
- Steps must be actionable and specific to the session's lesson.
- Never include: approvals, sandbox or deny exceptions, secrets,
  network bypasses, writes outside the workspace, or actions that were
  previously denied.
- Name: lowercase-hyphenated, at most 40 chars.
Output ONLY the SKILL.md text, starting with the front-matter '---'."""

_FM_NAME = re.compile(r"^name:\s*(\S[^\n]*)$", re.M)
_MAX_DRAFT_CHARS = 4000


def _is_skill_draft(text: str) -> bool:
    """Does the draft look like a SKILL.md (frontmatter + name)?"""
    t = (text or "").strip()
    if not t.startswith("---") or len(t) > _MAX_DRAFT_CHARS:
        return False
    return bool(_FM_NAME.search(t))


def _skill_name(text: str) -> str:
    m = _FM_NAME.search(text or "")
    name = (m.group(1).strip() if m else "")
    name = name.strip("`\"' ")
    return name or "learned-skill"


def _default_propose(name, text, *, evidence, source, base_version=""):
    """The real propose() once package 1 is merged; absent until then."""
    from delfin.agent.skill_proposals import propose  # type: ignore
    return propose(name, text, evidence=evidence, source=source,
                   base_version=base_version)


def _attach_runs(evidence, runs) -> None:
    """Give test evidence the session's green runs that cover it.

    ``runs`` is the session's own test ledger (``evidence["tests"]`` of
    the saved session: command, exit code, tree fingerprint). Only the
    entries that ran the cited file green are kept -- selected by the
    evidence module's own rule, so propose and accept agree on what
    covers what. Never raises.
    """
    try:
        if getattr(evidence, "kind", "") != "test" or not runs:
            return
        from delfin.agent.evidence import _green_runs_for
        node = str(evidence.ref).split("::")[0].replace("\\", "/").lstrip("./")
        evidence.runs = [dict(r) for r in _green_runs_for(runs, node)]
    except Exception:
        pass


def learn_from_session(messages, *, settings=None, llm=None,
                       _propose=None, runs=None) -> object | None:
    """At most ONE skill proposal from a finished session, or None.

    Order of the gates (each one may end the learning silently):
    qualification (qualifies) -> evidence (extract_evidence) -> ONE
    cheap LLM call that drafts a SKILL.md from the transcript excerpt,
    the reasons and the evidence -> propose() through skill_proposals.
    A draft without frontmatter is rejected, not "repaired". Never
    raises; any failure returns None. ``llm`` mirrors the
    ``llm_fn(prompt, system, settings) -> str`` injection used by
    memory_distill; without it one real cheap-tier client call runs.
    """
    try:
        reasons = qualifies(messages)
        if not reasons:
            return None
        evidence = extract_evidence(messages)
        if evidence is None:
            return None  # no proof, no proposal

        from delfin.agent.memory_distill import _transcript_excerpt
        excerpt = _transcript_excerpt(
            [m for m in (messages or []) if isinstance(m, dict)])
        ev_block = (f"kind: {evidence.kind}\nref: {evidence.ref}\n"
                    f"detail: {getattr(evidence, 'detail', '')}")
        prompt = (
            "Session excerpt:\n" + (excerpt or "(empty)") +
            "\n\nWhy this session earned a skill:\n- " +
            "\n- ".join(reasons) +
            "\n\nEvidence:\n" + ev_block +
            "\n\nDraft ONE reusable skill from this session.")

        if llm is not None:
            draft = llm(prompt, _DRAFT_SYSTEM, settings)
        else:
            from delfin.agent.memory_distill import _default_llm
            draft = _default_llm(prompt, _DRAFT_SYSTEM, settings)
        if not _is_skill_draft(draft):
            return None

        _attach_runs(evidence, runs)
        propose_fn = _propose or _default_propose
        return propose_fn(_skill_name(draft), draft.strip() + "\n",
                          evidence=[evidence],
                          source="skill_learning")
    except Exception:
        return None
