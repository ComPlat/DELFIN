"""The pipeline modes were retired in one place and left standing in four.

quick / reviewed / tdd / cluster / full were dropped from the dashboard
deliberately — the picker is labelled ``Workspace:`` now and answers one
question, which folder the session works in. The code says so in as many
words.

The retirement stopped there. Months later the five were still:

* declared in ``pack_lite/manifest.yaml``, so ``available_modes()``
  returned all nine;
* described in the dashboard's ``_MODE_DESCRIPTIONS``, a dictionary
  nothing looks them up in;
* advertised by the CLI, whose ``--mode`` help read
  "solo / plan / dashboard / quick / …" and whose flag has no
  ``choices=``, so it took them;
* listed in the pack's own README file map;
* present as five .md files, 4 875 bytes.

So the product had dropped a feature the CLI still offered. That is the
shape this codebase keeps producing: a decision applied at one hop and
dropped at the next.

A retired mode still has to be HARMLESS, not fatal — a saved session or
somebody's script may still name one. The engine already falls back to
solo for an unknown mode, and that is asserted below rather than assumed.
"""

from __future__ import annotations

import pathlib

import pytest
import yaml

_ROOT = pathlib.Path(__file__).resolve().parents[1]
_PACK = _ROOT / "delfin" / "agent" / "pack_lite"
_RETIRED = ("quick", "reviewed", "tdd", "cluster", "full")


# ---------------------------------------------------------------------------
# Gone from every place that claimed them
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", _RETIRED)
def test_the_manifest_no_longer_declares_it(mode):
    data = yaml.safe_load((_PACK / "manifest.yaml").read_text(encoding="utf-8"))
    assert mode not in {m.get("id") for m in (data.get("modes") or ())}


@pytest.mark.parametrize("mode", _RETIRED)
def test_the_mode_file_is_gone(mode):
    assert not (_PACK / "modes" / f"{mode}.md").exists()


@pytest.mark.parametrize("mode", _RETIRED)
def test_the_loader_does_not_offer_it(mode):
    from delfin.agent.prompt_loader import PromptLoader
    modes = PromptLoader(repo_dir=_ROOT / "delfin" / "agent").available_modes()
    assert mode not in modes


@pytest.mark.parametrize("mode", _RETIRED)
def test_the_cli_help_does_not_advertise_it(mode):
    """Asserted from the RENDERED help, not from a source-text offset.

    This read `src.index('run.add_argument("--mode"')` and sliced 400
    characters after it — so it broke the moment the flag moved into a
    shared helper, while the property it guards was never in danger. A
    test anchored on where code sits rather than on what it says fails on
    a refactor and stays silent on a real regression, which is the wrong
    way round. It now reads the help every front door actually prints.
    """
    from delfin.agent.cli import build_parser
    import argparse

    parser = build_parser()
    helps = [parser.format_help()]
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            helps.extend(sub.format_help() for sub in action.choices.values())

    for text in helps:
        assert f"/ {mode}" not in text and f"{mode} /" not in text, (
            f"a retired mode is still advertised: {mode}")


@pytest.mark.parametrize("mode", _RETIRED)
def test_the_pack_readme_does_not_list_it(mode):
    readme = (_PACK / "README.md").read_text(encoding="utf-8")
    assert f"modes/{mode}.md" not in readme


# ---------------------------------------------------------------------------
# What remains is what the picker offers
# ---------------------------------------------------------------------------

def test_the_four_that_remain_are_the_four_that_are_declared():
    from delfin.agent.prompt_loader import PromptLoader
    modes = PromptLoader(repo_dir=_ROOT / "delfin" / "agent").available_modes()
    assert sorted(modes) == ["dashboard", "office", "research", "solo"]


@pytest.mark.parametrize("mode", ["dashboard", "solo", "office", "research"])
def test_a_remaining_mode_still_has_its_file(mode):
    assert (_PACK / "modes" / f"{mode}.md").is_file()


# ---------------------------------------------------------------------------
# Retired is harmless, not fatal
# ---------------------------------------------------------------------------

def test_a_session_that_still_names_a_retired_mode_falls_back(tmp_path):
    """Somebody's script or a saved session may still say `quick`. It gets
    solo, not a crash — removing a mode must not strand a caller."""
    from unittest.mock import MagicMock, patch

    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        engine = AgentEngine(repo_dir=tmp_path, backend="api", mode="quick",
                             pack_dir=_ROOT / "delfin" / "agent")
    assert engine.mode == "solo"
    assert engine.route


def test_the_dashboard_describes_only_what_it_offers():
    """The description dictionary and the dropdown are two lists of the
    same thing, and they drifted. Neither may name a retired mode."""
    src = (_ROOT / "delfin" / "dashboard" / "tab_agent.py").read_text(
        encoding="utf-8")
    start = src.index("_MODE_DESCRIPTIONS = {")
    block = src[start:src.index("mode_dropdown = widgets.Dropdown", start)]
    for mode in _RETIRED:
        assert f'"{mode}":' not in block, f"{mode} is still described"
