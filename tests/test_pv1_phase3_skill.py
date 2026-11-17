"""PV1 phase 3 (red): the packaged provenance skill exists and teaches citation.

Target 5: skill text that tells the agent how to cite provenance in its
answer.  The packaged skill ``delfin/agent/pack/skills/provenance.md`` must
declare a domain (the CI gate ``test_every_packaged_skill_declares_a_domain``
fails on an undeclared skill), teach provenance citation against the published
API (``PROVENANCE.json``, ``provenance_record``), and carry no absolute home
path, user name, host name or ``~`` (PR-public rule, same class the C1
skill-text lint guards).

These cases are red on the current tree: the skill file does not exist.
"""

import re
from pathlib import Path

ACCENT = Path("delfin/agent/pack/skills/provenance.md")


def _body() -> str:
    return ACCENT.read_text(encoding="utf-8")


class TestSkillExistsAndDeclaresDomain:
    def test_skill_file_exists(self):
        assert ACCENT.is_file(), "packaged skill delfin/.../provenance.md is missing"

    def test_skill_declares_a_domain(self):
        text = _body()
        assert text.startswith("---"), f"{ACCENT.name} must open with YAML front matter"
        front, _, _ = text[3:].partition("---")
        m = re.search(r"(?m)^\s*domains\s*:\s*(\S.*)$", front)
        assert m and m.group(1).strip(), (
            f"{ACCENT.name} must declare a non-empty 'domains:' domain"
        )


class TestSkillTeachesCitation:
    def test_mentions_the_record_file(self):
        assert "PROVENANCE.json" in _body()

    def test_mentions_the_module(self):
        assert "provenance_record" in _body()

    def test_instructs_citing_program_and_version(self):
        text = _body().lower()
        assert "program" in text and "version" in text

    def test_instructs_citing_delfin_commit(self):
        text = _body().lower()
        assert "commit" in text

    def test_instructs_naming_the_chain(self):
        text = _body().lower()
        assert "chain" in text or "step" in text


class TestSkillCarriesNoPrivateContent:
    def test_no_home_path_or_tilde(self):
        text = _body()
        assert "~" not in text, "skill must not reference a home path via ~"

    def test_no_absolute_home_path(self):
        text = _body()
        assert "/home/" not in text, "skill must not name an absolute home path"

    def test_no_user_or_host_name(self):
        text = _body()
        assert "@" not in text, "skill must not carry an email/account-like token"
