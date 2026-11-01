"""Adversarial tests for delfin.agent.git_role (reviewer nacht-s2, GH1).

Attacks the URL normalisation and the role-resolution table: URL shapes a
user may write, fail-closed behaviour, entry-wins-in-both-directions, and
a non-mapping / invalid git_roles table.
"""
from __future__ import annotations

import delfin.agent.git_role as gr


class TestNormaliseRemote:
    def test_scp_style_case_insensitive_host(self):
        assert gr.normalise_remote("git@GitHub.com:Me/Repo.git") == \
            "github.com/Me/Repo"

    def test_ssh_url_with_port(self):
        assert gr.normalise_remote("ssh://git@github.com:22/Me/Repo") == \
            "github.com/Me/Repo"

    def test_https_with_user_and_trailing_slash(self):
        assert gr.normalise_remote("https://user@github.com/Me/Repo/") == \
            "github.com/Me/Repo"

    def test_git_protocol(self):
        assert gr.normalise_remote("git://github.com/Me/Repo.git") == \
            "github.com/Me/Repo"

    def test_ipv6_bracketed_host(self):
        # MEASURED (f61acf19b): the brackets are KEPT in the key
        # ('[2001:db8::1]/Me/Repo'). Deterministic — the same URL always
        # normalises to the same key, so settings entries still match —
        # but a settings entry written without brackets will not match.
        # Cosmetic inconsistency, sent to the builder as a question, not
        # a fail-open defect.
        assert gr.normalise_remote(
            "ssh://git@[2001:db8::1]/Me/Repo.git") == "[2001:db8::1]/Me/Repo"

    def test_file_url_is_path_key(self, tmp_path):
        repo = tmp_path / "repo"
        repo.mkdir()
        assert gr.normalise_remote(f"file://{tmp_path}/repo") == str(
            (tmp_path / "repo").resolve())

    def test_unknown_scheme_rejected(self):
        assert gr.normalise_remote("ftp://github.com/Me/Repo") is None

    def test_whitespace_rejected(self):
        assert gr.normalise_remote("git@github.com:Me/Repo .git") is None

    def test_traversal_in_path_rejected(self):
        assert gr.normalise_remote("git@github.com:Me/../Other") is None

    def test_empty_and_none(self):
        assert gr.normalise_remote("") is None
        assert gr.normalise_remote(None) is None

    def test_unicode_host_rejected(self):
        # str.lower can map look-alike hosts onto ASCII; the module
        # rejects non-ASCII hosts outright. Pin that decision.
        assert gr.normalise_remote("git@gïthub.com:Me/Repo") is None

    def test_relative_scp_host_with_slash_rejected_as_path(self):
        # 'host:with/slash' with a slash in the host part is not a
        # scp-style remote; it falls through to a path key.
        assert gr.normalise_remote("bad:re/po", base="/tmp") is not None


class TestEntryKey:
    def test_bare_host_path_kept_case_of_path(self):
        assert gr.entry_key("github.com/Me/Repo") == "github.com/Me/Repo"
        assert gr.entry_key("GitHub.com/Me/Repo") == "github.com/Me/Repo"

    def test_url_through_normalise(self):
        assert gr.entry_key("git@GitHub.com:Me/Repo.git") == \
            "github.com/Me/Repo"

    def test_absolute_path_rstrips_slash(self):
        assert gr.entry_key("/some/repo/") == "/some/repo"

    def test_invalid_forms(self):
        assert gr.entry_key("") is None
        assert gr.entry_key(None) is None
        assert gr.entry_key("host/") is None
        assert gr.entry_key("/Me/Repo with space") == "/Me/Repo with space"


class TestRoleForKeys:
    def _agent(self, table, role="contributor"):
        return {"git_role": role, "git_roles": table}

    def _fallback(self):
        return "contributor"

    def test_entry_wins_up(self):
        role, src = gr.role_for_keys(
            ["github.com/Me/Repo"], fallback=self._fallback,
            agent=self._agent({"github.com/Me/Repo": "maintainer"}))
        assert (role, src) == ("maintainer", "repository")

    def test_entry_wins_down_contributor_under_global_maintainer(self):
        role, src = gr.role_for_keys(
            ["github.com/Me/Repo"], fallback=lambda: "maintainer",
            agent=self._agent({"github.com/Me/Repo": "contributor"}))
        assert (role, src) == ("contributor", "repository")

    def test_no_entry_uses_fallback(self):
        role, src = gr.role_for_keys(
            ["github.com/Other/Repo"], fallback=lambda: "maintainer",
            agent=self._agent({}))
        assert (role, src) == ("maintainer", "installation")

    def test_multiple_keys_all_maintainer_required(self):
        role, _ = gr.role_for_keys(
            ["github.com/A/B", "github.com/C/D"],
            fallback=lambda: "maintainer",
            agent=self._agent({"github.com/C/D": "contributor"}))
        assert role == "contributor"

    def test_empty_keys_unresolved_contributor(self):
        assert gr.role_for_keys([], fallback=self._fallback) == \
            ("contributor", "unresolved")
        assert gr.role_for_keys(None, fallback=self._fallback) == \
            ("contributor", "unresolved")

    def test_non_mapping_table_fails_closed(self):
        role, src = gr.role_for_keys(
            ["github.com/A/B"], fallback=lambda: "maintainer",
            agent={"git_roles": "maintainer-everything"})
        assert (role, src) == ("contributor", "unresolved")

    def test_invalid_role_value_fails_closed(self):
        role, src = gr.role_for_keys(
            ["github.com/A/B"], fallback=self._fallback,
            agent=self._agent({"github.com/A/B": "root"}))
        assert (role, src) == ("contributor", "unresolved")

    def test_unnormalisable_entry_key_ignored(self):
        # A key that cannot be normalised never matches anything, so the
        # fallback decides - never an invalid entry.
        role, _ = gr.role_for_keys(
            ["github.com/A/B"], fallback=self._fallback,
            agent=self._agent({"": "maintainer", "bad host": "maintainer"}))
        assert role == "contributor"
