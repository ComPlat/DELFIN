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


def test_a4_window_collision_excuses_only_near_reason():
    # Red on 4c6668dc: the +-200 window around the grid mention
    # absorbed TightSCF's "because" and silently accepted the
    # unexplained grid deviation. Fixed by sentence-anchoring.
    ans = ("We used TightSCF because the preliminary run had SCF "
           "oscillations. For GridX7 we chose the default GridX5.")
    dev = [{"parameter": "scf", "default": "NormalSCF", "actual": "TightSCF"},
           {"parameter": "grid", "default": "GridX7", "actual": "GridX5"}]
    claims = rcl.check_deviations(ans, dev)
    flagged = {c.message for c in claims if c.code == "deviation-unexplained"}
    assert any("grid=GridX5" in m for m in flagged), flagged
    assert not any("TightSCF" in m for m in flagged), flagged


def test_reason_before_and_after_the_mention_in_the_same_sentence():
    dev = [{"parameter": "basis", "default": "def2-SVP",
            "actual": "def2-TZVP"}]
    before = "Because the SVP basis failed, we used def2-TZVP."
    after = "We used def2-TZVP because the SVP basis failed."
    for ans in (before, after):
        claims = rcl.check_deviations(ans, dev)
        assert not [c for c in claims if c.code == "deviation-unexplained"], ans


def test_record_reason_short_circuits_the_word_list():
    # CB0A<->SC1 contract (QS ruling s29): a non-empty record reason is
    # structure from the workflow itself and passes without any wording.
    ans = "For GridX7 we chose the default GridX5."
    dev = [{"parameter": "grid", "default": "GridX7", "actual": "GridX5",
            "reason": "better convergence in the benchmark"}]
    claims = rcl.check_deviations(ans, dev)
    assert not [c for c in claims if c.code == "deviation-unexplained"]


def test_record_reason_ignored_when_not_a_nonempty_string():
    ans = "For GridX7 we chose the default GridX5."
    for bad in (None, "", "   ", 7, ["x"]):
        dev = [{"parameter": "grid", "default": "GridX7", "actual": "GridX5",
                "reason": bad}]
        claims = rcl.check_deviations(ans, dev)
        assert any(c.code == "deviation-unexplained" for c in claims), bad
