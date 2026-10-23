"""Where a lazily-loaded block sits decides what a prefix cache can serve.

An OpenAI-compatible endpoint serves a prompt from its cache only up to
the first byte that differs, and the active module set is sticky and
monotonic: it grows the first time a session mentions ORCA or asks for a
permission. A block in the middle of the role file therefore turns every
later byte cold for every session that does not activate it the same way.

Measured on solo_agent.md, 2026-09-07: two sessions differing only in the
KIT-sandbox block (29% into the file) shared 24 343 chars of prefix. With
every block at the tail the worst pairing shares 43 108. On the KIT GLM
deployment a 15k-token prompt the endpoint cannot serve warm costs ~200 s
against ~10 s warm.
"""

import re
from pathlib import Path

import pytest

from delfin.agent.prompt_loader import PromptLoader

_PACK = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack"
_MARK = re.compile(r"^<!--\s*module:([a-zA-Z0-9_-]+)\s*-->\s*$", re.M)

# Every role file that uses the mechanism at all.
_ROLE_FILES = sorted(
    f for f in (_PACK / "agents").glob("*.md") if "<!-- module:" in f.read_text()
)


def test_the_mechanism_is_in_use():
    assert _ROLE_FILES, "no role file uses lazy modules — has it been removed?"


@pytest.mark.parametrize("path", _ROLE_FILES, ids=lambda p: p.name)
def test_every_lazy_block_sits_after_all_eager_text(path):
    """One block in the middle costs the whole tail of the cache."""
    text = path.read_text()
    first = _MARK.search(text)
    assert first is not None
    tail = text[first.start():]
    # From the first marker on, the only H2 sections are module blocks:
    # anything eager down there would sit behind a block that can vanish.
    for m in re.finditer(r"^## .*$", tail, re.M):
        before = tail[:m.start()]
        markers = _MARK.findall(before)
        assert markers, (
            f"{path.name}: '{m.group(0)[:60]}' is eager text after the first "
            f"module marker — move it above them")


def _prefix(task: str, key: str, model: str = "kit.glm-5.3") -> str:
    return PromptLoader().stable_prefix(
        role_id="solo_agent", mode_id="solo", task_text=task,
        session_key=key, model=model)


def _shared(a: str, b: str) -> int:
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


def test_two_sessions_share_their_prompt_up_to_the_blocks():
    """The number that decides a GLM turn: how much of the prompt two
    sessions have in common."""
    cases = {
        "plain": _prefix("Hallo", "P"),
        "permission": _prefix("bitte die berechtigung dauerhaft erlauben", "K"),
        "chemistry": _prefix("rechne mit ORCA die energie", "C"),
        "web": _prefix("recherchier das im internet bitte", "W"),
    }
    names = sorted(cases)
    worst = min(_shared(cases[a], cases[b])
                for i, a in enumerate(names) for b in names[i + 1:])
    # Was 24 343 with the blocks scattered; 43 108 with them at the tail.
    assert worst > 40_000, worst


def test_activating_a_module_only_changes_the_tail():
    """Within one session the set only grows. That growth must cost the
    tail of the cache, not the whole prompt. Measured on a model that
    still pays per module; GLM's profile keeps every module on."""
    plain = _prefix("Hallo", "S1", model="kit.deepseek-v4-flash")
    grown = _prefix("jetzt bitte mit ORCA rechnen und die berechtigung "
                    "dauerhaft erlauben", "S2", model="kit.deepseek-v4-flash")
    shared = _shared(plain, grown)
    assert shared > 40_000, shared
    # ... and the growth is real, or the test proves nothing.
    assert len(grown) > len(plain) + 3_000, (len(plain), len(grown))


def test_the_head_does_not_move_at_all_for_a_model_that_keeps_every_module():
    """GLM pays minutes for a cold head and seconds for a warm one, so its
    profile keeps every module on: the growth above costs it nothing,
    because there is none."""
    plain = _prefix("Hallo", "G1")
    grown = _prefix("jetzt bitte mit ORCA rechnen und die berechtigung "
                    "dauerhaft erlauben", "G2")
    assert plain == grown
