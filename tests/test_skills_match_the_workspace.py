"""A skill is only offered where its tools exist.

The skill catalogue is pasted into the `skill` tool's description for
every role alike, so both modes were told about playbooks they cannot
follow:

  * Office lost the chemistry TOOLS when its surface was cut to 36, but
    kept being offered the chemistry PLAYBOOKS -- nine of the twelve.
    Invoking one returns a procedure whose every step calls a tool the
    role may not reach, which costs a round and reads like the agent's
    own mistake.
  * The reverse holds too: the coding role has no `read_document`,
    `compare_tables` or `fill_pdf_form`, and was offered all three
    office skills.

The fix is the vocabulary the memory store already uses -- office vs
code, derived from the workspace folder -- applied to skills. A skill
that names no domain stays universal, so nothing can silently lose a
skill it has today.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import skills as S


PACK = Path(S.__file__).resolve().parent / "pack" / "skills"

# Tools that only exist on one side of the split. A skill that calls one
# of these belongs to that domain -- and this is checked against the file
# body, so a new skill cannot be added into the wrong half quietly.
_OFFICE_TOOLS = ("read_document", "compare_tables", "edit_sheet",
                 "fill_pdf_form", "fill_docx_template", "create_docx")
_CODE_TOOLS = ("search_docs", "search_calcs")


# ---------------------------------------------------------------------------
# The packaged skills declare where they belong
# ---------------------------------------------------------------------------

def test_every_packaged_skill_declares_a_domain():
    """An undeclared skill is offered everywhere -- fine for a user's own
    playbook, wrong for the curated set, which is split down the middle."""
    undeclared = [s.name for s in S.discover_skills() if not s.domains]
    assert not undeclared, (
        f"packaged skills without a domain: {undeclared}")


def test_a_skill_is_declared_where_its_tools_are():
    for sk in S.discover_skills():
        body = sk.source.read_text(encoding="utf-8")
        if any(t in body for t in _OFFICE_TOOLS):
            assert "office" in sk.domains, (
                f"{sk.name} calls office tools but is not an office skill")
        if any(t in body for t in _CODE_TOOLS):
            assert "code" in sk.domains, (
                f"{sk.name} calls technical tools but is not a code skill")


def test_both_halves_are_non_empty():
    """A split that empties one side would be a rename, not a split."""
    office = S.discover_skills(domain="office")
    code = S.discover_skills(domain="code")
    assert office and code
    assert {s.name for s in office}.isdisjoint({s.name for s in code})


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def test_office_is_not_offered_chemistry(tmp_path):
    names = {s.name for s in S.discover_skills(domain="office")}
    assert "casscf-setup" not in names
    assert "tddft-excited-states" not in names
    assert "reconcile-tables" in names


def test_the_coding_role_is_not_offered_office_playbooks():
    names = {s.name for s in S.discover_skills(domain="code")}
    assert "fill-form" not in names
    assert "casscf-setup" in names


def test_no_domain_asked_for_means_everything(tmp_path):
    """Callers that never learn about domains keep what they had."""
    assert len(S.discover_skills()) >= len(S.discover_skills(domain="office"))
    assert {s.name for s in S.discover_skills()} >= {
        s.name for s in S.discover_skills(domain="code")}


def test_an_undeclared_skill_is_offered_everywhere(tmp_path):
    ws = tmp_path / "ws"
    (ws / ".delfin" / "skills").mkdir(parents=True)
    (ws / ".delfin" / "skills" / "my-own.md").write_text(
        "# My own\n> Something I wrote myself.\n", encoding="utf-8")
    for domain in ("office", "code", ""):
        assert "my-own" in {s.name for s in S.discover_skills(ws, domain=domain)}


def test_a_declared_domain_survives_the_frontmatter(tmp_path):
    ws = tmp_path / "ws"
    (ws / ".delfin" / "skills").mkdir(parents=True)
    (ws / ".delfin" / "skills" / "mine.md").write_text(
        "---\nname: mine\ndomains: office\n---\n# Mine\n> Body.\n",
        encoding="utf-8")
    got = S.get_skill("mine", ws)
    assert got is not None and got.domains == ("office",)
    assert "mine" not in {s.name for s in S.discover_skills(ws, domain="code")}


def test_a_list_of_domains_is_accepted(tmp_path):
    ws = tmp_path / "ws"
    (ws / ".delfin" / "skills").mkdir(parents=True)
    (ws / ".delfin" / "skills" / "both.md").write_text(
        "---\nname: both\ndomains: [office, code]\n---\n# Both\n> Body.\n",
        encoding="utf-8")
    for domain in ("office", "code"):
        assert "both" in {s.name for s in S.discover_skills(ws, domain=domain)}


# ---------------------------------------------------------------------------
# The refusal has to say why
# ---------------------------------------------------------------------------

def test_calling_an_out_of_domain_skill_is_refused_by_name(tmp_path):
    """Not "not found" -- the skill exists, it just does not apply here."""
    from delfin.agent import api_client as A

    perms = A.KitToolPermissions(workspace=tmp_path, agent_role="office_agent")
    client = A._DocToolExecutor.__new__(A._DocToolExecutor)
    out = json.loads(client._execute_skill({"name": "casscf-setup"}, perms))
    assert "error" in out
    assert "casscf-setup" in out["error"]
    assert "office" in out["error"].lower()
    # It must not read as a typo: the available list is the office half.
    assert "reconcile-tables" in out.get("available", [])
    assert "casscf-setup" not in out.get("available", [])


def test_an_in_domain_skill_still_runs(tmp_path):
    from delfin.agent import api_client as A

    perms = A.KitToolPermissions(workspace=tmp_path, agent_role="office_agent")
    client = A._DocToolExecutor.__new__(A._DocToolExecutor)
    out = json.loads(client._execute_skill({"name": "reconcile-tables"}, perms))
    assert out.get("status") == "ok"
    assert out.get("skill") == "reconcile-tables"


def test_a_missing_skill_still_reads_as_missing(tmp_path):
    from delfin.agent import api_client as A

    perms = A.KitToolPermissions(workspace=tmp_path, agent_role="office_agent")
    client = A._DocToolExecutor.__new__(A._DocToolExecutor)
    out = json.loads(client._execute_skill({"name": "no-such-skill"}, perms))
    assert "not found" in out.get("error", "")
