"""SC1 phase 2: a deviation from workflow defaults needs a stated reason.

Control for result_claims.check_deviations(answer_text, deviations).
All 3 cases are RED on the previous commit (8112f686):
check_deviations does not exist yet (ImportError).
"""

import delfin.agent.result_claims as rcl


def test_deviation_without_reason_is_flagged():
    ans = ("We used B3LYP/def2-TZVP with TightSCF for the final energy. "
           "The electronic energy is -1234.567890 Eh.")
    dev = [{"parameter": "basis", "default": "def2-SVP", "actual": "def2-TZVP"}]
    claims = rcl.check_deviations(ans, dev)
    assert rcl.worst_level(claims) == "error"
    assert any(c.code == "deviation-unexplained" for c in claims)


def test_deviation_with_reason_passes():
    ans = ("We used B3LYP/def2-TZVP instead of the workflow default "
           "def2-SVP because the SVP basis gave a poor description of "
           "the metal centre in the preliminary run.")
    dev = [{"parameter": "basis", "default": "def2-SVP", "actual": "def2-TZVP"}]
    claims = rcl.check_deviations(ans, dev)
    assert not [c for c in claims if c.code == "deviation-unexplained"]


def test_never_raises_on_bad_input():
    assert rcl.check_deviations(None, None) == []
    assert rcl.check_deviations("", [None, "x", {}]) == []
    huge = "basis def2-TZVP " * 500_000
    assert isinstance(rcl.check_deviations(huge, [{"parameter": "basis",
                                                  "default": "def2-SVP",
                                                  "actual": "def2-TZVP"}]),
                      list)
