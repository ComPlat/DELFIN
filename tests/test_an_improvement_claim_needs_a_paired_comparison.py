"""A comparative figure in an answer needs a paired comparison behind it.

Pinned here:

* "improved by 12 %", "3x faster", "rose from 6/10 to 9/10", "um 30 %
  schneller" are claims, and get the caveat when no paired comparison
  ran in the session;
* not claims: a hedged sentence, a conditional or plan, a question, a
  figure attributed to a source, fenced code, a quoted line, a table
  row, a figure without a comparative word, an ordinary answer;
* a clean ``delfin-agent compare-arms`` / ``bench compare`` / ``bench
  gate``, ``hyperfine`` with two commands or ``pyperf compare_to`` in
  the ledger silences it; ``hyperfine`` with one command does not;
* the engine appends the caveat at its single exit and keeps it in the
  transcript.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import comparison_claims as C


@pytest.mark.parametrize("text", [
    "The cache change improved the pass rate by 12 %.",
    "With the new parser the suite is now 3x faster.",
    "After the fix the pass rate rose from 6/10 to 9/10.",
    "Latency dropped by 40 percent with the batching switch.",
    "Die Änderung macht den Lauf um 30 % schneller.",
    "Mit dem Patch ist die Fehlerrate um 2 Prozentpunkte gesunken.",
    "Token use went down by half: from 120000 to 60000, a 50% reduction.",
])
def test_a_comparative_figure_is_a_claim(text):
    found = C.scan(text, [])
    assert found, text
    cav = C.caveat(found)
    assert cav.startswith("\n\n[verify] Caveat: comparative figure without")
    assert "delfin-agent compare-arms" in cav


@pytest.mark.parametrize("text", [
    "The change is probably about 10 % faster, but I have not measured it.",
    "This should make it roughly 2x faster once the cache is warm.",
    "We expect the pass rate to rise from 6/10 to 9/10.",
    "Is it 20 % faster now?",
    "According to the paper, the method is 3x faster than CCSD(T).",
    "Laut Handbuch ist RIJCOSX etwa 10x schneller.",
    "```\nimproved by 12 %\n```",
    "> the old report said: 30 % faster",
    "| def2-SVP | ~2-3x faster | coarse geometry |",
    "The run took 12 % of the allocation and 300 s in total.",
    "Fixed the import; all 38 tests pass.",
    "Die Energie liegt bei -1234.5 Eh.",
])
def test_these_are_not_claims(text):
    assert C.scan(text, []) == [], text


@pytest.mark.parametrize("cmd", [
    "bash delfin-agent compare-arms a.jsonl b.jsonl --axis X --stamp-a a --stamp-b b",
    "bash delfin-agent bench compare --run x",
    "bash delfin-agent bench gate --json",
    "bash hyperfine 'python a.py' 'python b.py'",
    "bash hyperfine --warmup 3 'python a.py' 'python b.py'",
    "bash python -m pyperf compare_to a.json b.json",
])
def test_a_paired_comparison_in_the_ledger_silences_it(cmd):
    assert C.scan("The cache change improved the pass rate by 12 %.",
                  [cmd]) == []


def test_one_timing_command_is_no_pair():
    assert C.scan("B is 9 % faster than A.", ["bash hyperfine 'python b.py'"])


def test_the_engine_appends_it_at_its_exit(tmp_path):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        engine = AgentEngine(repo_dir=tmp_path, backend="cli", mode="quick")
    text = "With the switch on, the benchmark is now 25 % faster."
    engine.messages = [{"role": "assistant", "content": text}]
    engine._exec_commands_session = ["bash python bench.py"]
    out = engine._append_answer_caveats(text, functional=[], ambiguous=[])
    assert "[verify] Caveat: comparative figure" in out
    assert "[verify] Caveat: comparative figure" in engine.messages[-1]["content"]
    engine.messages = [{"role": "assistant", "content": text}]
    engine._exec_commands_session = [
        "bash delfin-agent compare-arms a.jsonl b.jsonl --axis X"]
    assert engine._append_answer_caveats(
        text, functional=[], ambiguous=[]) == text
