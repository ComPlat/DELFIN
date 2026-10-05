"""Package G, Phase 4: reach and additivity checks.

Two things must be proven before a measurement is believed:

* REACH -- "ran" is not "hit".  Before measuring, prove the switch is read
  and changes something on EVERY case it targets.  A switch that was never
  read, or that only touched a subset of its target cases, is not a hit -- a
  partial reach is exactly the silent self-deception this package stops, so
  `check_reach` refuses it.

* ADDITIVITY -- the switch, off, must reproduce baseline byte-identically
  (a candidate change sits behind ONE switch, default off; additive by
  construction).  Byte-identity is STRICT: no stripping, no normalising.  A
  switch-off output that embeds a timestamp or a counter is NOT identical to
  baseline and the check ABORTS -- the agent must not quietly smooth over a
  real difference.
"""

import pytest

from delfin.agent.experiment import (
    ExperimentError,
    check_reach,
    require_switch_off_identical,
)


# ── Reach: the switch must be read and effective on ALL targets ──────────


def test_reach_refuses_unread_switch():
    with pytest.raises(ExperimentError):
        check_reach(switch_read=False, switch_name="DELFIN_MODE",
                    targeted=["case_a", "case_b"], effective=[])


def test_reach_refuses_zero_effective_cases():
    with pytest.raises(ExperimentError):
        check_reach(switch_read=True, switch_name="DELFIN_MODE",
                    targeted=["case_a", "case_b"], effective=[])


def test_reach_refuses_partial_hit():
    # Switch changed one of two target cases -- that is not a hit on every
    # case it targets; refuse, don't accept the partial.
    with pytest.raises(ExperimentError):
        check_reach(switch_read=True, switch_name="DELFIN_MODE",
                    targeted=["case_a", "case_b"], effective=["case_a"])


def test_reach_passes_when_switch_hits_all_targets():
    # No raise = pass.
    check_reach(switch_read=True, switch_name="DELFIN_MODE",
                targeted=["case_a", "case_b"], effective=["case_a", "case_b"])


def test_reach_refuses_effective_cases_not_in_targets():
    # The switch changed a case it did NOT target -- that reaches beyond the
    # declared surface and is a defect, not a hit.
    with pytest.raises(ExperimentError):
        check_reach(switch_read=True, switch_name="DELFIN_MODE",
                    targeted=["case_a"], effective=["case_a", "unrelated_case"])


# ── Additivity: switch-off must be byte-identical to baseline, else abort ─


def test_identical_output_passes():
    require_switch_off_identical(switch_off=b"abc", baseline=b"abc", label="case_a")


def test_any_byte_difference_aborts():
    with pytest.raises(ExperimentError):
        require_switch_off_identical(switch_off=b"abd", baseline=b"abc", label="case_a")


def test_length_difference_aborts():
    with pytest.raises(ExperimentError):
        require_switch_off_identical(switch_off=b"abcd", baseline=b"abc", label="case_a")


def test_binary_blob_identity_is_byte_exact():
    blob = bytes(range(256))  # every byte value, including NULs
    require_switch_off_identical(switch_off=blob, baseline=blob, label="blob")


def test_embedded_timestamp_field_aborts():
    # A timestamp in the switch-off output (e.g. action_protocol-style) makes
    # it differ from baseline; the check must ABORT, not strip/normalise.
    with pytest.raises(ExperimentError):
        require_switch_off_identical(
            switch_off=b"output timestamp=1234567890 end",
            baseline=b"output timestamp=<STAMP> end",
            label="case_a")


def test_embedded_counter_field_aborts():
    with pytest.raises(ExperimentError):
        require_switch_off_identical(
            switch_off=b"row1\nrow2\nrow3",
            baseline=b"row1\nrow2",
            label="counter")


def test_empty_output_equal_to_empty_baseline_passes():
    # Switch-off producing the same absence as baseline is genuinely
    # additive (empty == empty).
    require_switch_off_identical(switch_off=b"", baseline=b"", label="empty")


def test_additivity_message_names_the_abort():
    try:
        require_switch_off_identical(switch_off=b"A", baseline=b"B", label="case_a")
        assert False, "must raise"
    except ExperimentError as e:
        msg = str(e)
        assert "byte-identical" in msg or "identical" in msg or "abort" in msg.lower()
