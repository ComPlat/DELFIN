"""Three things the first multi-session run (2026-09-16, three GLM sessions
building `delfin doctor`) showed the chat getting wrong.

1. A message another session sent, and a wake-up the agent had scheduled,
   were shown mid-run as "💬 [you, mid-run]".
2. The glitch-clean note blamed gpt-5.x and advised raising the effort on a
   GLM turn.
3. A file named in the user's task and not yet created was flagged as an
   invented citation, and the forced correction cost a GLM round each time.
"""
import inspect
from pathlib import Path

from delfin.agent import api_client as A
from delfin.agent import verify_guard as vg

TAB = Path(inspect.getfile(__import__("delfin.dashboard.tab_agent", fromlist=["x"]))).read_text()


def test_a_steer_is_labelled_by_its_sender():
    assert A._steer_label("bitte nur mod1").startswith("💬 [you, mid-run]")
    assert A._steer_label('[Message from the session "B" — not from the user.] hi').startswith("✉ [another session")
    assert A._steer_label("[scheduled] Warte auf Session A").startswith("⏰ [wake-up")


def test_the_loop_uses_the_label_everywhere():
    src = inspect.getsource(A.OpenAIClient.stream_message)
    assert "💬 [you, mid-run]" not in src
    assert src.count("_steer_label(") >= 2


def test_the_clean_note_names_the_model_and_blames_gpt_only_for_gpt():
    i = TAB.index("🧹 Cleaned the output of ")
    body = TAB[i - 900:i + 200]
    assert '"gpt" in _mdl.lower()' in body
    assert "_san.leaked_tools and" in body


def test_a_path_named_in_the_task_is_not_an_invented_citation(tmp_path):
    (tmp_path / "delfin").mkdir()
    answer = "Ich lege tests/test_delfin_doctor.py an und nutze delfin/doctor.py."
    task = "Schreibe tests/test_delfin_doctor.py für das Modul delfin/doctor.py, das Session A baut."
    assert vg.scan_for_ungrounded_code_claims(answer, repo_root=tmp_path, named_in=task) == []
    flagged = vg.scan_for_ungrounded_code_claims(answer, repo_root=tmp_path)
    assert "delfin/doctor.py" in {f.path for f in flagged if f.kind == "nonexistent"}


def test_the_dashboard_hands_the_session_text_to_both_scans():
    assert TAB.count("named_in=_session_user_text()") == 2
    assert "def _session_user_text(" in TAB
