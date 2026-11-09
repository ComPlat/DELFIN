"""VR1 live-path hook: verify number+unit claims in a finished answer.

The hook ``engine._append_claim_verification`` runs at the very end of
``engine.stream_response``, before the answer is returned.  It scans the
finished answer for number+unit claims, verifies each against the real
``*.out`` of the turn's calc folder(s) (from the live ledger
``_doc_executor._calc_evidence``), keeps the machine-readable
:class:`VerificationRecord` list on ``engine.last_claim_verification`` and
appends a short "could not be verified" caveat on a mismatch / no-source
claim.

These tests drive the REAL public path (``AgentEngine.stream_response`` with a
fake model backend) and the REAL module ``_doc_executor``, whose
``_calc_evidence`` ledger the hook reads.  The calc output is a synthetic
fixture (no ORCA-generated data in the repo — licence rule).
"""
import time
from unittest.mock import MagicMock, patch

from delfin.agent.api_client import StreamEvent


def _answering_client(text: str):
    """A client that answers with one text_delta ``text`` and a billed delta."""
    def stream(system, messages, max_tokens=4096, session_id="",
               thinking_budget=0, **kw):
        yield StreamEvent(type="text_delta", text=text)
        yield StreamEvent(type="message_delta", input_tokens=100,
                          output_tokens=10, cost_usd=0.0)
    client = MagicMock()
    client.model = "claude"
    client.stream_message = MagicMock(side_effect=stream)
    return client


def _engine(agent_tree, text: str):
    """An AgentEngine on a fresh repo whose backend answers ``text``."""
    from delfin.agent.engine import AgentEngine
    client = _answering_client(text)
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="api", provider="claude")
    eng.client = client
    eng._claim_guard_active = True  # skip verify_guard's own correction turn
    return eng


def _calcdir(tmp_path, energy: str) -> str:
    """A real synthetic calc folder whose largest *.out holds ``energy``."""
    d = tmp_path / "calc_oh"
    d.mkdir()
    (d / "run.out").write_text(
        "ORCA calculation\n"
        f"FINAL SINGLE POINT ENERGY     {energy}\n"
        "\n"
    )
    return str(d)


def _seed_calc_evidence(folder: str) -> None:
    """Write a real-shaped record into the live ledger the hook reads."""
    from delfin.agent.api_client import _doc_executor
    _doc_executor._calc_evidence.append({
        "ts": time.time(),
        "folder": folder,
        "outcome": "run done",
        "worst": "ok",
    })


def _fresh_repo(tmp_path):
    """A throwaway git repo with the minimal agent manifest."""
    import subprocess
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-b", "main"], cwd=repo,
                   check=True, capture_output=True)
    (repo / "manifest.yaml").write_text(
        "agent:\n  mode: solo\n  provider:\n    api: claude\n"
    )
    return repo


def _clear_calc_evidence() -> None:
    from delfin.agent.api_client import _doc_executor
    _doc_executor._calc_evidence.clear()


def test_verified_answer_passes_unchanged_and_record_is_attached(tmp_path):
    """Correct number + unit: no caveat, machine record verdict=verified."""
    folder = _calcdir(tmp_path, "-152.1234567")
    eng = _engine(_fresh_repo(tmp_path),
                  "The energy is -152.1234567 Eh.")
    _seed_calc_evidence(folder)
    try:
        text = eng.stream_response("what is the energy?")
    finally:
        _clear_calc_evidence()
    assert "[verify]" not in text
    recs = eng.last_claim_verification
    assert recs, "a verified claim must leave a machine record"
    assert recs[0]["verdict"] == "verified"
    assert recs[0]["unit"] == "Eh"
    assert recs[0]["value"] == -152.1234567
    assert abs(recs[0]["parsed_value"] - (-152.1234567)) < 1e-6


def test_planted_wrong_number_appends_caveat_and_records_mismatch(tmp_path):
    """A number the real output does NOT hold: caveat + verdict=mismatch."""
    folder = _calcdir(tmp_path, "-152.1234567")  # real value
    eng = _engine(_fresh_repo(tmp_path),
                  "The energy is -152.9 Eh.")  # planted wrong
    _seed_calc_evidence(folder)
    try:
        text = eng.stream_response("what is the energy?")
    finally:
        _clear_calc_evidence()
    assert "[verify]" in text
    assert "could not be verified" in text
    recs = eng.last_claim_verification
    assert recs, "a mismatch must leave a machine record"
    assert recs[0]["verdict"] == "mismatch"


def test_empty_calc_evidence_never_false_alarms(tmp_path):
    """No run this turn: prose numbers produce no caveat and no record."""
    eng = _engine(_fresh_repo(tmp_path),
                  "The notion is 3 years old; index 7 was 2021 Eh input.")
    _clear_calc_evidence()
    text = eng.stream_response("hello")
    assert "[verify]" not in text
    assert eng.last_claim_verification == []
