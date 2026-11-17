"""LS1 phase 3 red test: the cite-sources skill ships and stays in step.

The skill `cite-sources` is the human/agent instruction half of LS1: it
tells the answering agent to cite every method choice and reference value
to the local index or mark it "(recalled, no source)". Two contracts are
pinned so the skill cannot silently drift away from the code it documents:

* The skill must exist as a *packaged* (built-in) skill and declare the
  `code` domain — the domain whose tools (`search_docs` / cite) it uses.
* Its body must actually name the LS1 functions (`cite`,
  `check_answer_citations`) — otherwise the playbook teaches an interface
  that does not exist.

RED before the skill ships: discover_skills() knows no built-in
`cite-sources`, so both assertions fail.
"""

from delfin.agent import skills as S


def _find(name):
    for sk in S.discover_skills():
        if sk.name == name:
            return sk
    return None


def test_cite_sources_skill_is_packaged_and_declares_a_domain():
    sk = _find("cite-sources")
    assert sk is not None, "the cite-sources skill is not packaged"
    assert "code" in sk.domains, (
        f"cite-sources uses the code-side tools and must declare the "
        f"code domain, got {sk.domains}")


def test_skill_body_names_the_ls1_functions_it_documents():
    sk = _find("cite-sources")
    assert sk is not None
    body = sk.source.read_text(encoding="utf-8")
    # cite() is how an answer obtains a citable passage with a locator;
    # check_answer_citations is the audit it tells the agent to run.
    assert "cite(" in body, "skill must name cite(query)"
    assert "check_answer_citations" in body, (
        "skill must name the citation audit it asks the agent to run")
