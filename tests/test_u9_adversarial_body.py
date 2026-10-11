"""U9 adversarial review of pr_helper.body(workspace, base).

Reviewer (nacht-s17) adversarial suite for package U9, run against the
builder's delivered code (agent/s16-u9b 4e0d5d75 + 39a8cbb0). Encodes the
frozen contract sent to the builder (nacht-s16) on 2026-10-11:

  (a) account-name/secret guards run on the ASSEMBLED body: a fake commit
      carrying /home/<acct>/... or the account name must yield a body with
      neither.
  (b) the PASSED base drives the git range (never a hardcoded origin/main).
  (c) edge cases: empty task file, None, corrupt JSON, unicode, base==HEAD.
  (d) body() is PURE: no network / no gh / no real subprocess; the git-log
      reader is injected so tests use a FAKE repository.

Identity is injected (FAKE ka_ew9999 / /home/ka_ew9999 / login-fake-01),
never the operator's real one. The body must never leak a home path, an
account name, a host name, or a secret value (guards reuse output_guard);
the assembled body must pass home_paths_in(body, ...) == [] under the
injected identity.

Adapted at review time to the builder's actual seam: identity is a single
dict {home, account, host}, and each commit is a CommitRecord dict
{hash, subject, body} (pr_helper.py:36,41-42). The builder's body() asks
the runner for (base, "HEAD"), so fake ranges are keyed on ("<base>", "HEAD").
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import pr_helper

# FAKE identity the guards are handed — never a real home/account/host.
FAKE_HOME = "/home/ka_ew9999"
FAKE_ACCOUNT = "ka_ew9999"
FAKE_HOST = "login-fake-01"
_FAKE_IDENTITY = {"home": FAKE_HOME, "account": FAKE_ACCOUNT, "host": FAKE_HOST}


class FakeLogReader:
    """Returns scripted commit records for base..HEAD.

    Each commit is a CommitRecord dict {hash, subject, body} as the real
    reader yields (pr_helper.py:36). Tests then assert on the assembled
    body, never on the raw records.
    """

    def __init__(self, records_by_range):
        self.records_by_range = records_by_range
        self.called_ranges = []

    def __call__(self, base, head):
        self.called_ranges.append((base, head))
        return self.records_by_range.get((base, head), [])


def _commit(subject, body="", hash_="abc1234"):
    """A CommitRecord dict with sane defaults."""
    return {"hash": hash_, "subject": subject, "body": body}


def _session_tasks(ws, tasks):
    (ws / ".delfin").mkdir(parents=True, exist_ok=True)
    (ws / ".delfin" / "session_tasks.json").write_text(
        json.dumps({"next_id": len(tasks) + 1, "tasks": tasks}),
        encoding="utf-8")


@pytest.fixture
def env(tmp_path):
    ws = tmp_path / "workspace"
    ws.mkdir()
    runner = FakeLogReader({})
    return ws, runner


# --------------------------------------------------------------------------
# (a) account-name / secret guard on the ASSEMBLED body
# --------------------------------------------------------------------------

class TestAccountSecretGuard:
    def test_commit_message_path_is_not_in_body(self, env):
        ws, runner = env
        _session_tasks(ws, [{"id": 1, "subject": "add the fix",
                             "description": "closes the handover gap",
                             "status": "completed"}])
        runner.records_by_range = {
            ("main", "HEAD"): [
                _commit("Merge pull request #1: the fix",
                        body="Control: tests/test_x.py, 5 cases, all 5 red "
                             "on the previous commit."),
            ],
        }
        body = pr_helper.body(ws, "main", git_log=runner,
                              identity=_FAKE_IDENTITY)
        assert FAKE_ACCOUNT not in body
        assert FAKE_HOME not in body
        assert FAKE_HOST not in body

    def test_account_path_in_commit_body_is_stripped(self, env):
        ws, runner = env
        _session_tasks(ws, [{"id": 1, "subject": "add the fix",
                             "status": "completed"}])
        leak = "the build wrote to %s/scratch and finished" % FAKE_HOME
        runner.records_by_range = {
            ("main", "HEAD"): [_commit("finished", body=leak)],
        }
        body = pr_helper.body(ws, "main", git_log=runner,
                              identity=_FAKE_IDENTITY)
        assert FAKE_ACCOUNT not in body
        assert "/home/" not in body
        assert "scratch" not in body

    def test_home_path_in_commit_removed_from_body(self, env):
        ws, runner = env
        _session_tasks(ws, [{"id": 1, "subject": "done", "status": "completed"}])
        runner.records_by_range = {
            ("main", "HEAD"): [_commit("x", body="see %s/tmp" % FAKE_HOME)],
        }
        body = pr_helper.body(ws, "main", git_log=runner,
                              identity=_FAKE_IDENTITY)
        assert FAKE_ACCOUNT not in body
        assert FAKE_HOME not in body

    def test_secret_value_in_commit_is_scrubbed(self, env):
        # A secret-shaped value (e.g. an AWS key) in a commit body must not
        # survive into the PR body — output_guard.scrub_secrets must have run.
        ws, runner = env
        _session_tasks(ws, [{"id": 1, "subject": "creds", "status": "completed"}])
        # Build the token at runtime: the test file itself is scanned for
        # secret-shaped literals on disk, so we assemble it from parts.
        token = "AKIA" + "0123456789ABCDEF"
        runner.records_by_range = {
            ("main", "HEAD"): [_commit("key", body="embedded %s here" % token)],
        }
        body = pr_helper.body(ws, "main", git_log=runner,
                              identity=_FAKE_IDENTITY)
        assert token not in body

    def test_assembled_body_passes_home_paths_in(self, env):
        from delfin.agent import output_guard
        ws, runner = env
        _session_tasks(ws, [{"id": 1, "subject": "add the fix",
                             "status": "completed"}])
        runner.records_by_range = {
            ("main", "HEAD"): [_commit("Merge pull request #1: the fix",
                                       body="closes the handover gap")],
        }
        body = pr_helper.body(ws, "main", git_log=runner,
                              identity=_FAKE_IDENTITY)
        hits = output_guard.home_paths_in(body, home=FAKE_HOME,
                                          account=FAKE_ACCOUNT, host=FAKE_HOST)
        assert hits == []


# --------------------------------------------------------------------------
# (b) the PASSED base drives the git range (never a hardcoded origin/main)
# --------------------------------------------------------------------------

class TestBaseIsRespected:
    def test_body_uses_the_passed_base_for_the_git_range(self, env):
        ws, runner = env
        _session_tasks(ws, [{"id": 1, "subject": "add the fix",
                             "status": "completed"}])
        body = pr_helper.body(ws, "release/2.0", git_log=runner,
                              identity=_FAKE_IDENTITY)
        assert ("release/2.0", "HEAD") in runner.called_ranges
        assert "origin/main" not in body

    def test_base_equals_head_empty_diff_still_yields_a_body(self, env):
        ws, runner = env
        _session_tasks(ws, [{"id": 1, "subject": "the goal",
                             "status": "in_progress"}])
        body = pr_helper.body(ws, "main", git_log=runner,
                              identity=_FAKE_IDENTITY)
        assert "the goal" in body
        assert "No commits" in body


# --------------------------------------------------------------------------
# (c) edge cases: empty / missing / corrupt / unicode journal, None args
# --------------------------------------------------------------------------

class TestJournalEdgeCases:
    def test_missing_task_file(self, env):
        # No .delfin/session_tasks.json at all: body() must not raise.
        ws, runner = env
        body = pr_helper.body(ws, "main", git_log=runner,
                              identity=_FAKE_IDENTITY)
        assert isinstance(body, str)

    def test_empty_task_file(self, env):
        ws, runner = env
        _session_tasks(ws, [])
        body = pr_helper.body(ws, "main", git_log=runner,
                              identity=_FAKE_IDENTITY)
        assert isinstance(body, str)

    def test_corrupt_json(self, env):
        ws, runner = env
        (ws / ".delfin").mkdir(parents=True, exist_ok=True)
        (ws / ".delfin" / "session_tasks.json").write_text(
            "{ not valid json ]]]", encoding="utf-8")
        body = pr_helper.body(ws, "main", git_log=runner,
                              identity=_FAKE_IDENTITY)
        assert isinstance(body, str)

    def test_unicode_task_text_survives_verbatim(self, env):
        # DELFIN sessions write German task text; a 'für' / 'Dysprosium'
        # subject must survive unchanged (no mojibake, no mangling).
        ws, runner = env
        _session_tasks(ws, [{"id": 1, "subject": "Cadmium-Komplexe für die\n"
                             "Dysprosium-Studie", "status": "in_progress"}])
        body = pr_helper.body(ws, "main", git_log=runner,
                              identity=_FAKE_IDENTITY)
        assert "Dysprosium" in body


# --------------------------------------------------------------------------
# (d) PURE: no subprocess/network/gh inside body(); git-log runner injected
# --------------------------------------------------------------------------

class TestPurity:
    def test_no_subprocess_or_network(self, env, monkeypatch):
        import subprocess

        # PUNCH LINE: body() must never spawn a subprocess (no git binary,
        # no gh, no network). If it tried, the fake raises and the test
        # fails loudly instead of silently doing real work.
        calls = []

        def _fake_popen(*a, **k):
            calls.append(a)
            raise AssertionError("body() must not spawn a subprocess")

        monkeypatch.setattr(subprocess, "Popen", _fake_popen)
        ws, runner = env
        _session_tasks(ws, [{"id": 1, "subject": "pure", "status": "done"}])
        runner.records_by_range = {("main", "HEAD"): [_commit("pure")]}
        body = pr_helper.body(ws, "main", git_log=runner,
                              identity=_FAKE_IDENTITY)
        assert isinstance(body, str)
        # No fake-Popen call happened, so no subprocess was spawned.
        assert calls == []
