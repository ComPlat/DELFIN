"""Token budget guards for dashboard / solo agent prompts.

These tests fail if a future edit blows the role-prompt size past a
fixed budget. Token estimation uses ``len(text) / 4`` (a tight
upper-bound vs. real BPE tokenization).
"""
from __future__ import annotations

from pathlib import Path

import pytest


_REPO = Path(__file__).resolve().parents[1]
_PROMPT_DIR = _REPO / "delfin" / "agent" / "pack" / "agents"


def _estimate_tokens(text: str) -> int:
    """Fast ~upper-bound token count (1 token ≈ 4 chars in English/German)."""
    return (len(text) + 3) // 4


# Budgets are set just above the CURRENT prompt sizes, so any regrowth
# fails immediately. The prompt diet removed duplicated rules (the same
# contract stated in a role prompt AND a shared addendum), worked examples
# that restated a rule already given, historical incident narrative, and
# blocks that re-listed what a tool schema already declares. No behavioral
# contract was dropped — extending one is fine, but pay for it by trimming
# elsewhere rather than by raising the number.
@pytest.mark.parametrize(
    "filename, max_tokens",
    [
        # 7950 -> 7150: dropped the two worked ACTION examples (they restate
        # the plan-first / verify-after-set rules stated above them), the
        # ORCA counter-example lists, and the duplicated tab-set + command-
        # discovery blocks.
        ("dashboard_agent.md", 7150),
        # 14200 -> 10600: dropped the worked-example dialogs and the
        # "how these compound" walk-through, folded the three separate
        # workspace-location statements into one, compressed the sandbox
        # internals, and removed tool-signature listings that duplicate the
        # tool schemas.
        # 10600 -> 10621, and the 21 are the MEASURED remainder of one rule,
        # not a round number chosen for comfort. Live 2026-08-14: given
        # "Behebe den fehlschlagenden Test" — an instruction that names no
        # test — the agent picked one and EDITED it, fifteen tool calls, an
        # outright forbidden-signal violation. With an empty workspace the
        # same agent asks, in numbered form; with files present it guesses.
        # The autonomy section covered "several valid approaches" and not
        # "no target named, several candidates present", which is the case
        # that has something plausible to do and therefore hides.
        # Paid first: three passages in that same section were compressed,
        # returning 49 of the rule's 82 tokens. 21 is what is left after
        # the sentence was cut twice more, and it buys the narrowest form
        # of the rule — ask before the WRITE, never before the reads.
        # 10621 -> 10672, fifty-one tokens, for the condition under which
        # two routing rows apply. The table sent "Gibbs/SPE many folders"
        # to extract_energy_table with nothing saying the parser is ORCA
        # -specific -- and on xtb it answers status: ok with every value
        # null instead of an error, and reads only the LARGEST .out in the
        # folder, which in the science fixture is the run that did not
        # converge. So a model following the table gets a confident empty
        # answer about the wrong file.
        #
        # Reported by kit.deepseek-v4-flash, asked after a real task what
        # got in its way: it formed the hypothesis that it should have
        # used those parsers, TESTED it against the xtb files, and
        # refuted its own hypothesis with the tool output. Verified here
        # before acting -- and the string form is worse than reported,
        # since properties="scf_converged,single_point" is iterated
        # character by character into keys s, c, f, _, o, n, v ...
        # That part is in delfin/api.py and belongs to another owner.
        ("solo_agent.md", 10672),
        # Written lean from the start: the shared addenda carry the general
        # contracts, so this prompt only states what is specific to working
        # on someone's real records. Raised as the mode's surface grew —
        # series work, record-addressed edits, PDF assembly, remembered
        # folder conventions, the working folder — each replacing something
        # the model would otherwise improvise, so each has to be named. The
        # prose was tightened three times on the way; what is left is
        # contract, not explanation.
        # 1600 -> 1657 on 2026-09-09 for one tool that was in the
        # catalogue and in no routing table. sum_column appeared nowhere in
        # this prompt, while the arithmetic sentence told the model to
        # total a column in bash — so a model that followed the prompt read
        # the CSV with cat and added the amounts in its own text, which is
        # the failure the office module exists to prevent, and it did it
        # while obeying every instruction it had been given. Measured on
        # kit.glm-5.3, office_total_names_what_it_left_out: two bash calls,
        # no document tool, the five amounts written out and summed in the
        # answer.
        #
        # Paid as far as it goes: the rule was cut twice, and the coverage
        # principle above it dropped the list of skipped-row kinds, which
        # the tool now reports itself. 133 tokens became 57, and those 57
        # buy a routing row plus the sentence that stops the shell
        # one-liner. Not a relaxation — the prompt was describing a surface
        # that had changed underneath it.
        #
        # 1657 -> 1679 the same day, for draft_email. Generalising the
        # check above found it: the role could call it and the prompt named
        # no tool for an email, while the egress rule two sections up tells
        # the model that sending anything out is asked about first. So the
        # likely behaviour was to refuse or improvise a text file, when the
        # tool writes a .eml and has no network at all. 22 tokens for the
        # row, and the row says the part that resolves the conflict.
        #
        # 1679 -> 1688, nine tokens, for the words "reads the column's own
        # convention". kit.deepseek-v4-flash gave exactly that as its
        # REASON for going to the shell instead:  "Ich summiere in bash
        # mit Python, um die Konvention korrekt zu beruecksichtigen" --
        # then spent 27 of one task's 50 tool calls writing and debugging
        # the script, for a total sum_column returns in one call and with
        # the coverage attached. The old rule is what it was echoing, so
        # removing that rule was necessary and not sufficient: the model
        # also has to know the tool covers the thing it was worried about.
        ("office_agent.md", 1688),
    ],
)
def test_role_prompt_within_token_budget(filename, max_tokens):
    """Role prompts must stay below their per-role token budget.

    Module markers are removed before counting. They are stripped from
    every composed prompt (``_strip_lazy_modules`` swallows the marker
    line whether or not the module survives, and the disabled path
    substitutes them away), so a marker cannot reach a model and cannot
    cost a token at runtime. Charging the file for them would make the
    mechanism that SHRINKS the prompt read as growth, and would price a
    26-character comment against text the user actually pays for.
    """
    import re as _re

    path = _PROMPT_DIR / filename
    assert path.exists(), f"missing prompt file: {path}"
    text = _re.sub(r"^<!--\s*module:[a-zA-Z0-9_-]+\s*-->\s*$\n?", "",
                   path.read_text(), flags=_re.M)
    actual = _estimate_tokens(text)
    assert actual <= max_tokens, (
        f"{filename}: {actual} tokens (>{max_tokens} budget). "
        f"Trim before extending."
    )


# The file budgets above guard the markdown. This one guards what the model
# actually receives: the CACHEABLE HEAD of the composed prompt (role prompt +
# shared addenda + project context). It is the part that is re-sent verbatim
# on every request of a session, so it is where prompt cost is decided.
@pytest.mark.parametrize(
    "role_id, mode_id, route, max_stable_tokens",
    [
        ("solo_agent", "solo", ["solo_agent"], 11400),
        ("dashboard_agent", "dashboard", ["dashboard_agent"], 10400),
    ],
)
def test_composed_stable_head_within_budget(
        monkeypatch, role_id, mode_id, route, max_stable_tokens):
    from delfin import user_settings
    from delfin.agent.prompt_loader import PromptLoader

    # Pin the lazy-module setting so the budget measures the prompt, not the
    # machine's local configuration.
    monkeypatch.setattr(
        user_settings, "load_settings",
        lambda *a, **k: {"agent": {"slim_prompt": True}})

    report = PromptLoader().prompt_size_report(
        role_id=role_id, mode_id=mode_id, route=route,
        task_text="fix the failing test in foo.py",
        session_key="budget-1")
    actual = report["stable_tokens"]
    assert actual <= max_stable_tokens, (
        f"{role_id}: cacheable head is {actual} tokens "
        f"(>{max_stable_tokens} budget). Trim before extending."
    )


def test_dashboard_prompt_keeps_essential_sections():
    """The prompt must keep the irreducible safety + grounding sections.

    The list pins TODAY's intentional design (dashboard mode is
    guide+UI only since 3a2c802; ORCA claims must be grounded in the
    manual since P4) — update it consciously when the design changes,
    never to silence a failure.
    """
    text = (_PROMPT_DIR / "dashboard_agent.md").read_text()
    must_have = [
        "ACTION:",                       # how commands work
        "Safety rules",                  # safety policy
        "Hard scope limits",             # guide+UI-only contract (3a2c802)
        "Ground every ORCA",             # manual-grounding rule (P4)
        "Tools you may NOT use",         # forbidden tool surface
        "agent_workspace",               # spelled out as NOT available
    ]
    missing = [k for k in must_have if k not in text]
    assert not missing, f"essential sections dropped: {missing}"


def test_a_module_marker_never_reaches_a_model():
    """The budget above excuses markers from the count. That is only
    honest while they really are removed from every composed prompt."""
    import re

    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader()
    for role, mode in (("solo_agent", "solo"),
                       ("dashboard_agent", "dashboard"),
                       ("office_agent", "office")):
        for task in ("Hallo", "rechne mit ORCA die energie und such im netz"):
            built = loader.build_system_prompt(
                role_id=role, mode_id=mode, task_text=task,
                session_key=f"marker-{role}-{len(task)}")
            assert not re.search(r"<!--\s*module:", built), (role, mode, task)
