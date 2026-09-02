"""A stated count is compared to ONE enumeration, not to every bullet added up.

Field report 20260902-080635. An answer about hyperpolarizability said the
sum runs "über alle 27 Tensor-Komponenten" — a fact about a rank-3 tensor
in three dimensions, promising no list at all — and was told:

    ⚠️ Diese Antwort nennt 27, führt aber 11 Einträge auf.

The 11 was two separate Top-5 blocks plus a four-item statistics block,
counted across the whole document with `_LIST_ITEM_RE.findall(text)`. No
enumeration of 11 things existed anywhere in the answer.

It stood directly beside a warning that DID matter — that a count came
from grep output which had been truncated — and that is the real damage.
A reader who is shown one warning that is plainly wrong learns to skip the
next one, and the next one was the true one.

The repair is structural rather than lexical: the enumeration a number
refers to is the longest UNBROKEN run of list items. Nothing here asks
what the number is *about*, which is not decidable, only which list it
could possibly be contradicting.

Also pinned here: the rule that sent the same answer wrong in the first
place. DELFIN computes beta_HRS itself, in parser.py, with the real
orientational average; the agent guessed sqrt(sum(beta^2)/5), disagreed
with DELFIN's own stored value by up to 47%, and called that a
"convention". A convention is a constant factor — its own comparison
(-2.33% to +47.04%) refutes the explanation it was given.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import verify_guard as vg

_PROMPT = Path("delfin/agent/pack/agents/solo_agent.md")


# The answer as it was written, in its actual shape: three short lists,
# each introduced by its own bold lead-in.
FIELD_ANSWER = """## Zusammenfassung

Ich habe den Tensor aus allen **404 Ordnern** extrahiert.

wobei die Summe über alle 27 Tensor-Komponenten läuft.

**Statistik (atomare Einheiten):**
- Minimum: 5.869 au
- Maximum: 393.488 au
- Mittelwert: 102.030 au
- Median: 101.416 au

**Top 5 höchste βHRS-Werte:**
1. `JMC-2021-BTP7_C_245` (TZVP): 393.488 au
2. `JMC-2021-BTP7_C_245` (SVP): 363.184 au
3. `JMC-2021-BTP7_T_325` (TZVP): 350.689 au
4. `JMC-2021-BTP7_T_325` (SVP): 327.802 au
5. `JMC-2021-BTH1_C_1` (TZVP): 324.950 au

**Dateien erstellt:**
- beta_hrs_results.json
- beta_hrs_summary.csv
"""

# What the guard was built for: 31 claimed, 29 actually listed.
TIGHT_LIST = "Ich habe 31 PDF-Dateien verifiziert.\n\n" + "".join(
    f"{i}. Datei_{i}.pdf geprüft\n" for i in range(1, 30))

# The same, written with a blank line between every item.
LOOSE_LIST = "Ich habe 31 Dateien geprüft.\n\n" + "\n".join(
    f"- Datei_{i}.pdf\n" for i in range(1, 30))


# ---------------------------------------------------------------------------
# What the number can be compared against
# ---------------------------------------------------------------------------

def test_separate_lists_are_not_added_together():
    """The reported false positive, gone."""
    assert vg.scan_for_count_vs_enumeration(FIELD_ANSWER) == []


def test_the_longest_run_is_the_longest_single_list():
    """5, not 11 — the Top-5 block, not every bullet in the answer."""
    assert vg._longest_list_block(FIELD_ANSWER) == 5


def test_a_bold_lead_in_ends_the_list_above_it():
    """The detail the first attempt got wrong.

    Treating ``**`` lines as non-breaking merged the three blocks back
    into one run of 11 and the false warning survived the repair — the
    bold line is exactly where one list stops and the next subject
    starts.
    """
    two_lists = "- a\n- b\n- c\n\n**Andere Sache:**\n- d\n- e\n"
    assert vg._longest_list_block(two_lists) == 3


def test_a_paragraph_ends_it_too():
    prose_between = "- a\n- b\n- c\n\nDann ein Satz dazwischen.\n\n- d\n- e\n"
    assert vg._longest_list_block(prose_between) == 3


def test_blank_lines_inside_one_list_do_not_end_it():
    """A loosely spaced list is still one list."""
    assert vg._longest_list_block(LOOSE_LIST) == 29


# ---------------------------------------------------------------------------
# The case the guard exists for still fires
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("answer", [TIGHT_LIST, LOOSE_LIST])
def test_a_count_that_overstates_its_own_list_is_still_caught(answer):
    assert vg.scan_for_count_vs_enumeration(answer) == [(31, 29)]


def test_the_caveat_still_names_both_numbers():
    note = vg.count_vs_enumeration_caveat(
        vg.scan_for_count_vs_enumeration(TIGHT_LIST))
    assert "31" in note and "29" in note


def test_a_short_list_is_not_an_enumeration():
    """A count of 31 followed by two examples contradicts nothing."""
    assert vg.scan_for_count_vs_enumeration(
        "31 Dateien geprüft.\n\n- eine\n- zwei\n") == []


# ---------------------------------------------------------------------------
# The rule that would have prevented the wrong answer
# ---------------------------------------------------------------------------

def test_the_prompt_sends_a_formula_question_to_the_source():
    """search_docs indexes manuals; the derivation lives in the code.

    Both of the agent's queries — the right questions, asked twice — came
    back with config keys at scores 0.26 and 0.12, because parser.py's
    derivation is not in the docs index and cannot be.
    """
    text = _PROMPT.read_text()
    assert "grep -rn" in text and "delfin/" in text
    assert "indexes manuals, not source" in text


def test_the_prompt_says_a_stored_value_outranks_a_recomputation():
    text = _PROMPT.read_text()
    assert "DELFIN_Data.json" in text
    assert "CHECKS DELFIN's value" in text


def test_the_prompt_names_why_a_varying_deviation_is_not_a_convention():
    """The sentence that turns the field case from a judgement call into
    an arithmetic one: the agent's own spread was -2.33% to +47.04%."""
    text = _PROMPT.read_text()
    assert "a convention is a constant factor" in text
