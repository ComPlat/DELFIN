"""Package G, reviewer s6: adversarial tests for phase 4 (reach + additivity).

Attacks the builder's 13 tests do not cover:

  A4.1  reach: duplicate case names in targeted (set() collapses them)
  A4.2  reach: whitespace-only case names are accepted as "cases"
  A4.3  reach: switch_read=True with NO targets at all (empty targeted)
  A4.4  reach: effective subset that is a superset of a REORDERED list
  A4.5  additivity: str instead of bytes is refused (type channel)
  A4.6  additivity: single byte difference (trailing newline) ABORTS
  A4.7  additivity: identical bytes with equal but non-empty label passes
  A4.8  additivity: empty bytes vs empty bytes passes
  A4.9  additivity: label with only whitespace is refused
  A4.10 reach: effective list carries duplicates -- must not mask a miss
"""

import pytest

from delfin.agent.experiment import (
    ExperimentError,
    check_reach,
    require_switch_off_identical,
)


# ── A4.1: duplicate targets collapse in a set() ──────────────────────────

def test_reach_duplicate_target_names(tmp_path=None):
    """targeted = [a, a, b] collapses to {a, b} in the check.  A duplicate
    declaration is a malformed instrument -- document that the check still
    demands both distinct names be hit (no silent passthrough)."""
    with pytest.raises(ExperimentError):
        check_reach(switch_read=True, switch_name="S",
                    targeted=["a", "a", "b"], effective=["a"])


# ── A4.2: whitespace-only case names ─────────────────────────────────────

def test_reach_whitespace_case_name():
    """A whitespace-only case name is not a case.  The check takes
    sequences of names and does not strip them -- 'a ' and 'a' are
    DIFFERENT cases for set membership.  Pin the (surprising) behavior:
    'a ' targeted but 'a' effective is a refusal, because the instrument
    must declare exact names."""
    with pytest.raises(ExperimentError):
        check_reach(switch_read=True, switch_name="S",
                    targeted=["a "], effective=["a"])


# ── A4.3: switch_read=True but zero targeted cases ───────────────────────

def test_reach_empty_targeted_refuses():
    """An empty targeted list means the switch targets nothing -- there is
    no experiment.  Must refuse even with switch_read=True."""
    with pytest.raises(ExperimentError):
        check_reach(switch_read=True, switch_name="S",
                    targeted=[], effective=[])


# ── A4.4: effective beyond targets refuses (superset attack) ─────────────

def test_reach_superset_effective_refuses():
    """effective ⊋ targeted is the 'reaches beyond its declared surface'
    defect.  With partial-attack reversed: every target IS hit, but the
    switch also changed untargeted cases.  Must refuse."""
    with pytest.raises(ExperimentError):
        check_reach(switch_read=True, switch_name="S",
                    targeted=["a"], effective=["a", "b"])


# ── A4.5: additivity with str instead of bytes ───────────────────────────

def test_additivity_str_not_bytes_refuses():
    """'abc' == b'abc' is False in Python, so a str/bytes mix would abort
    with the byte-identity message instead of the clearer type message.
    Either refusal is acceptable; SILENT PASS is not."""
    with pytest.raises(ExperimentError):
        require_switch_off_identical(switch_off="abc", baseline=b"abc",
                                     label="c1")


# ── A4.6: single trailing-newline byte difference ABORTS ─────────────────

def test_additivity_single_byte_diff_aborts():
    """The STRICT byte-identity promise: one byte differs (\\n) -> abort.
    No normalisation, no tolerance."""
    with pytest.raises(ExperimentError):
        require_switch_off_identical(switch_off=b"out\n", baseline=b"out",
                                     label="c1")


# ── A4.7: identical bytes with a real label passes ───────────────────────

def test_additivity_identical_passes():
    require_switch_off_identical(switch_off=b"out", baseline=b"out",
                                 label="case_one")


# ── A4.8: empty bytes vs empty bytes passes ──────────────────────────────

def test_additivity_empty_bytes_pass():
    require_switch_off_identical(switch_off=b"", baseline=b"", label="c0")


# ── A4.9: whitespace-only label ──────────────────────────────────────────

def test_additivity_whitespace_label_refuses():
    with pytest.raises(ExperimentError):
        require_switch_off_identical(switch_off=b"x", baseline=b"x",
                                     label="   ")


# ── A4.10: duplicate effective names must not mask a miss ────────────────

def test_reach_duplicate_effective_no_mask():
    """effective = [a, a] with targeted = [a, b]: the duplicate must not
    make b appear covered.  Set arithmetic refuses here."""
    with pytest.raises(ExperimentError):
        check_reach(switch_read=True, switch_name="S",
                    targeted=["a", "b"], effective=["a", "a"])
