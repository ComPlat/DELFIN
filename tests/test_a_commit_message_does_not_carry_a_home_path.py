"""A commit message may not carry a home path, account or host name.

Every brief this project hands an agent repeats the rule -- no home
paths, no account names, no machine names in code, comments or commit
texts -- and nothing enforced it. The repository is public. A rule that
lives only in a prompt is a hope; the model that forgets it is the same
model the rule is addressed to.

The gate already knows what a commit message is: ``_prose_blanked``
blanks it so the PATH scan does not read it, because a message naming
files it skips was refused as if it were reaching for them. This is the
mirror of that -- the same prose, read for what it must not contain.

  a home path in the message      refused, and the message says what to
                                  write instead
  the account or host name        refused for the same reason: together
                                  with a public repository they name a
                                  person and a machine
  a relative path                 allowed; that is the fix
  the same text NOT in a message  not this check's business -- a command
                                  that touches a path is a path question
"""

from __future__ import annotations

import pytest

WHO = {"home": "/pfs/data6/home/xx/yy_gr0042/zz_ab1234",
       "account": "zz_ab1234", "host": "node0991"}


def home_paths_in(*a, **kw):
    # Imported per call, so a missing name fails each case on its own
    # instead of taking the whole file out at collection.
    from delfin.agent.output_guard import home_paths_in as _f
    return _f(*a, **kw)


def _find(text):
    return home_paths_in(text, **WHO)


class TestTheDetector:
    def test_the_exact_home_path_is_found(self):
        assert _find("see /pfs/data6/home/xx/yy_gr0042/zz_ab1234/src/a.py")

    def test_a_foreign_home_path_is_found_too(self):
        assert _find("copied from /home/alice/project/run.sh")

    def test_a_macos_home_path_is_found(self):
        assert _find("/Users/bob/Code/thing.py was the model")

    def test_the_account_name_alone_is_found(self):
        assert _find("reported by zz_ab1234 on Tuesday")

    def test_the_host_name_alone_is_found(self):
        assert _find("only reproduces on node0991")

    def test_a_relative_path_is_clean(self):
        assert _find("see delfin/agent/api_client.py:2199") == []

    def test_a_tilde_path_is_clean(self):
        assert _find("the store lives in ~/.delfin/credentials.json") == []

    def test_ordinary_prose_is_clean(self):
        assert _find("Read why it is there before removing it") == []

    def test_a_short_token_does_not_count_as_an_account(self):
        assert home_paths_in("run it", home="/home/ab", account="ab",
                             host="n1") == []

    def test_it_returns_what_it_found(self):
        hits = _find("from /home/alice/x.py on node0991")
        assert any("alice" in h for h in hits)
        assert "node0991" in hits


class TestTheGate:
    @pytest.fixture
    def gate(self, monkeypatch):
        from delfin.agent import api_client as ac
        monkeypatch.setattr(ac, "_identity_for_commit_texts",
                            lambda: dict(WHO))
        return ac._commit_text_identity_hit

    def test_a_message_with_a_home_path_is_refused(self, gate):
        assert gate('git commit -m "fix /home/alice/x.py"')

    def test_a_clean_message_passes(self, gate):
        assert gate('git commit -m "Read why it is there"') is None

    def test_git_tag_is_covered(self, gate):
        assert gate('git tag -a v1 -m "built on node0991"')

    def test_a_path_outside_a_message_is_not_this_check(self, gate):
        assert gate("cat /home/alice/x.py") is None

    def test_a_non_git_command_is_not_this_check(self, gate):
        assert gate('echo "/home/alice/x.py"') is None

    def test_the_refusal_says_what_to_write_instead(self, gate):
        err = gate('git commit -m "see /home/alice/x.py"')
        assert "relative" in err.lower()

    def test_the_refusal_does_not_repeat_the_account_name(self, gate):
        err = gate('git commit -m "by zz_ab1234"')
        assert "zz_ab1234" not in err, "the refusal must not restate it"
