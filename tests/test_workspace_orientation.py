"""The agent must be oriented at its workspace, not the DELFIN source tree.

Field case 20260729-115601: the session environment block reported the
DELFIN checkout as ``cwd``. Told it was working inside DELFIN, the model
switched to absolute paths under the source tree and built the user's
project (venv, games, notebook) into DELFIN's own checkout instead of
the workspace the dashboard was launched in.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent.prompt_loader import PromptLoader


def _loader(tmp_path: Path) -> PromptLoader:
    pack = tmp_path / "pack_root"
    (pack / "pack" / "shared").mkdir(parents=True)
    (pack / "pack_lite").mkdir(parents=True)
    return PromptLoader(repo_dir=pack)


def test_env_block_reports_the_workspace_not_the_source_tree(tmp_path):
    loader = _loader(tmp_path)
    ws = tmp_path / "user_project"
    ws.mkdir()
    loader.workspace_root = ws
    block = loader._build_session_env_block()
    assert f"cwd: {ws}" in block
    assert str(loader.repo_root) not in block.splitlines()[0]


def test_env_block_warns_when_workspace_differs_from_source_tree(tmp_path):
    loader = _loader(tmp_path)
    ws = tmp_path / "user_project"
    ws.mkdir()
    loader.workspace_root = ws
    block = loader._build_session_env_block()
    assert "do not build the user's project inside it" in block


def test_no_warning_when_working_on_delfin_itself(tmp_path):
    loader = _loader(tmp_path)
    loader.workspace_root = Path(loader.repo_root)
    block = loader._build_session_env_block()
    assert "do not build the user's project inside it" not in block


def test_unset_workspace_falls_back_to_repo_root(tmp_path):
    loader = _loader(tmp_path)
    assert loader.workspace_root is None
    assert f"cwd: {loader.repo_root}" in loader._build_session_env_block()


def test_written_files_count_as_observed_evidence():
    """A file the agent just wrote is grounded — citing it must not be
    flagged as an unverified claim (field case: 5 false flags on files
    the agent had created in the same turn)."""
    from delfin.agent.api_client import _OBSERVATION_TOOLS, _observe_read_files
    for tool in ("write_file", "edit_file", "multi_edit", "notebook_edit"):
        assert tool in _OBSERVATION_TOOLS
    observed: set = set()
    _observe_read_files(observed, "write_file",
                        {"path": "voila_games/snake.py"}, "File created: ...")
    assert observed == {"voila_games/snake.py"}


def test_failed_write_is_not_recorded_as_evidence():
    from delfin.agent.api_client import _observe_read_files
    observed: set = set()
    _observe_read_files(observed, "write_file", {"path": "x.py"},
                        '{"error": "write failed"}')
    assert observed == set()


def test_code_language_rule_ships_in_the_shared_pack():
    """"Ships" means it reaches the model, not that a file contains it.

    This asserted the phrase inside work_cycle_rules.md and was green for
    as long as that file existed — while `load_shared_context` handed
    every role outside the four full-context ones a 16-line slice taken
    entirely from delfin_context.md, so neither the code-English rule nor
    the answer-language rule beside it was ever in the solo, dashboard or
    office prompt. Qwen kept answering English questions in German with
    the fix "shipped". The rule now lives in the honesty addendum, which
    every role does get, and this builds the prompt and looks in it.
    """
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    for role in ("solo_agent", "dashboard_agent", "office_agent",
                 "builder_agent"):
        prompt = loader.build_system_prompt(role, "read calc.py")
        assert "INTO code stays English" in prompt, role
        assert "docstrings" in prompt, role
        assert "language of the user's LATEST message" in prompt, role


def test_delfin_context_suppressed_for_a_user_project(tmp_path):
    """A session in the user's own project must not carry DELFIN's
    product context and module paths — they do not apply, cost tokens,
    and invite the model to drift into the source tree."""
    loader = _loader(tmp_path)
    shared = tmp_path / "pack_root" / "pack" / "shared"
    (shared / "delfin_context.md").write_text(
        "# DELFIN Context\n" + "module guidance\n" * 40, encoding="utf-8")
    loader.is_delfin_workspace = False
    prompt = loader.build_system_prompt("solo_agent", "build a game")
    assert "module guidance" not in prompt
    assert "NOT working on DELFIN's own source" in prompt


def test_delfin_context_kept_when_working_on_delfin(tmp_path):
    loader = _loader(tmp_path)
    shared = tmp_path / "pack_root" / "pack" / "shared"
    (shared / "delfin_context.md").write_text(
        "# DELFIN Context\n" + "module guidance\n" * 40, encoding="utf-8")
    loader.is_delfin_workspace = True
    prompt = loader.build_system_prompt("solo_agent", "fix the pipeline")
    assert "module guidance" in prompt


def test_unknown_workspace_type_keeps_legacy_behaviour(tmp_path):
    loader = _loader(tmp_path)
    shared = tmp_path / "pack_root" / "pack" / "shared"
    (shared / "delfin_context.md").write_text(
        "# DELFIN Context\n" + "module guidance\n" * 40, encoding="utf-8")
    assert loader.is_delfin_workspace is None
    prompt = loader.build_system_prompt("solo_agent", "task")
    assert "module guidance" in prompt


def test_artifact_verification_rule_ships_in_the_contract():
    """Field case 20260729-141616: the agent shipped start.sh as the way
    to run the dashboard but verified voila by calling it directly, so an
    invalid flag in the script reached the user."""
    from pathlib import Path as _P
    text = (_P(__file__).resolve().parent.parent / "delfin" / "agent"
            / "pack" / "shared"
            / "scientific_integrity_addendum.md").read_text(encoding="utf-8")
    assert "Test the artifact you hand over" in text
    assert "not a stand-in" in text
    assert "named as untested" in text


def test_bash_cwd_guidance_matches_reality():
    """The prompt claimed chained `cd` was blocked; it is not. Guidance
    that contradicts the runtime teaches the model to distrust it."""
    from pathlib import Path as _P
    prompt = (_P(__file__).resolve().parent.parent / "delfin" / "agent"
              / "pack" / "agents" / "solo_agent.md").read_text(
                  encoding="utf-8")
    assert "bash(cd <path> && ...)` is blocked" not in prompt
    assert "Prefer `cwd=<path>`" in prompt
