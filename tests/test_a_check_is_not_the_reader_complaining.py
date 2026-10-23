"""A verification prompt says who is asking and what was read.

Field session, 2026-09-11 (GLM, dashboard): an answer about the CO2
coordinator cited "chain_setup.py" and "CONTROL.txt" bare. Neither
exists at the workspace root, so both were flagged as fabricated and a
correction turn was forced. The prompt said only that the paths did not
exist; the model, whose history holds its answer but not its seven tool
calls, concluded it had read nothing, told the user it had invented the
citations, and re-read the files. The user saw an apology for a
fabrication that never happened.
"""
import inspect
from pathlib import Path

from delfin.agent import verify_guard as vg


def _kinds(text, root):
    return {f.path: f.kind
            for f in vg.scan_for_ungrounded_code_claims(text, repo_root=root)}


def test_a_bare_name_that_exists_somewhere_is_not_a_fabricated_path(tmp_path):
    (tmp_path / "delfin" / "co2").mkdir(parents=True)
    (tmp_path / "delfin" / "co2" / "chain_setup.py").write_text("x = 1\n")
    (tmp_path / "runs" / "job1").mkdir(parents=True)
    (tmp_path / "runs" / "job1" / "CONTROL.txt").write_text("co2=on\n")
    kinds = _kinds("Das Setup steht in chain_setup.py, die Keys in "
                   "CONTROL.txt.", tmp_path)
    assert set(kinds) == {"chain_setup.py", "CONTROL.txt"}
    assert set(kinds.values()) == {"unread"}


def test_a_path_with_a_directory_is_still_checked(tmp_path):
    (tmp_path / "delfin").mkdir()
    kinds = _kinds("Siehe delfin/co2/nope.py:12.", tmp_path)
    assert kinds.get("delfin/co2/nope.py") == "nonexistent"


def test_an_observed_bare_name_is_not_flagged_at_all(tmp_path):
    flags = vg.scan_for_ungrounded_code_claims(
        "Siehe chain_setup.py.", repo_root=tmp_path,
        observed_files={"delfin/co2/chain_setup.py"})
    assert flags == []


def _nonexistent(path):
    return vg.CodeClaimFlag(path=path, line=None, kind="nonexistent")


def test_the_code_claim_prompt_names_who_asks_and_what_was_read():
    text = vg.code_claim_feedback(
        [_nonexistent("chain_setup.py"), _nonexistent("CONTROL.txt")],
        observed={"delfin/co2/CO2_Coordinator6.py", "delfin/cli.py"})
    assert "Automatic check, not a message from the user" in text
    assert "'chain_setup.py'" in text and "'CONTROL.txt'" in text
    assert "delfin/co2/CO2_Coordinator6.py" in text
    assert "delfin/cli.py" in text
    assert "Do not apologise" in text


def test_the_location_claim_prompt_does_the_same():
    flag = vg.LocationClaimFlag(claim="delfin/define.py:119",
                                path="delfin/define.py", line=119,
                                kind="file_line")
    text = vg.location_claim_feedback([flag], observed={"delfin/api.py"})
    assert "Automatic check, not a message from the user" in text
    assert "'delfin/define.py:119'" in text
    assert "delfin/api.py" in text
    assert "Do not apologise" in text


def test_nothing_read_is_said_so():
    text = vg.code_claim_feedback([_nonexistent("foo/bar.py")])
    assert "Nothing was read or grepped this turn." in text


def test_a_long_ledger_is_cut_not_dumped():
    obs = {f"delfin/mod{i}.py" for i in range(20)}
    text = vg.code_claim_feedback([_nonexistent("foo/bar.py")], observed=obs)
    assert "(+14 more)" in text


def test_the_engine_hands_the_ledger_to_the_prompt():
    from delfin.agent import engine
    src = inspect.getsource(engine.AgentEngine._enforce_claim_grounding)
    assert "location_claim_feedback(" in src
    assert 'observed=getattr(self, "_last_observed_files", None)' in src


def test_the_dashboard_says_it_is_an_automatic_check():
    src = Path(__file__).resolve().parents[1].joinpath(
        "delfin", "dashboard", "tab_agent.py").read_text()
    assert "🔎 Automatic check:" in src
    i = src.index("_vg.code_claim_feedback(")
    assert '"_last_observed_files"' in src[i:i + 300]
    assert "the answer is being corrected" not in src


def test_a_bare_name_that_exists_nowhere_is_still_a_fabrication(tmp_path):
    (tmp_path / "delfin").mkdir()
    kinds = _kinds("Laut S1.out liegt S1 bei 2.31 eV.", tmp_path)
    assert kinds.get("S1.out") == "nonexistent"


def test_the_correction_asks_for_the_corrections_only():
    text = vg.code_claim_feedback([_nonexistent("foo/bar.py")])
    assert "do not repeat the rest of the answer" in text
    flag = vg.LocationClaimFlag(claim="a/b.py:1", path="a/b.py", line=1,
                                kind="file_line")
    assert "not the whole answer again" in vg.location_claim_feedback([flag])


def test_the_dashboard_flushes_text_before_thinking_on_a_tool_call():
    """Field case: the answer appeared twice, above the tool calls."""
    src = Path(__file__).resolve().parents[1].joinpath(
        "delfin", "dashboard", "tab_agent.py").read_text()
    i = src.index("def _on_tool_use(tool_name, tool_input):")
    body = src[i:i + 1600]
    text_flush = body.index('_update_last_assistant("".join(chunks), role_label,')
    thinking_flush = body.index("if thinking_chunks:")
    assert text_flush < thinking_flush
