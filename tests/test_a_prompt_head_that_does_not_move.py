"""On the KIT GLM deployment a prompt head the endpoint cannot serve warm
costs 199-266 s against 7-12 s warm (measured 2026-09-07). The lazy
modules sit at the tail of the role file for that reason -- and still, a
module that triggers on turn two is inserted ahead of the ones already
active, and every byte after it goes cold: the second turn of a session
shared 81% of its prompt with the first (measured 2026-09-11). For a
model whose profile says so, every module stays on, and the head is the
same on every turn of every session.
"""

from __future__ import annotations

from delfin.agent.model_profiles import get_profile
from delfin.agent.prompt_loader import PromptLoader


def test_the_profile_says_which_model_keeps_every_module():
    assert get_profile("kit.glm-5.3").all_prompt_modules is True
    assert get_profile("kit.deepseek-v4-flash").all_prompt_modules is False


def _loader():
    return PromptLoader()


def test_everything_means_every_module_in_solo_mode():
    ld = _loader()
    assert ld._detect_active_modules("hallo", "solo") == set()
    assert ld._detect_active_modules("hallo", "solo", everything=True) == set(ld._MODULE_TRIGGERS)


_ROLE = """# Role
Base text that every task gets.

<!-- module:chemistry -->
## ORCA grounding
Read the manual before you cite a keyword.

<!-- module:documents -->
## Spreadsheets, PDFs and Word files
Use the document tools.
"""


def _head(model: str, task: str, session: str) -> str:
    return _loader()._strip_lazy_modules(
        _ROLE, task_text=task, mode_id="solo", model=model,
        session_key=session, role_id="solo_agent")


def test_the_head_is_the_same_whatever_the_task_says_for_that_model():
    a = _head("kit.glm-5.3", "hallo", "s1")
    b = _head("kit.glm-5.3", "read the ORCA output and the xlsx", "s2")
    assert a == b
    assert "ORCA grounding" in a and "Spreadsheets" in a
    assert "<!-- module:" not in a


def test_the_dashboard_router_still_pays_only_for_what_it_uses():
    """Half of its utterances trigger nothing, by measurement, and a
    navigation request must not carry the chemistry rules -- that file
    pins it for this very model."""
    ld = _loader()
    nav = ld._strip_lazy_modules(_ROLE, task_text="wechsel zu Submit", mode_id="dashboard",
                                 model="kit.glm-5.3", session_key="d1", role_id="dashboard")
    assert "ORCA grounding" not in nav


def test_the_other_model_still_pays_only_for_what_it_uses():
    a = _head("kit.deepseek-v4-flash", "hallo", "s3")
    b = _head("kit.deepseek-v4-flash", "read the ORCA output", "s4")
    assert a != b
    assert "ORCA grounding" not in a and "ORCA grounding" in b
