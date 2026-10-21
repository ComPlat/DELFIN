"""M1 phase-3 QS ruling (nacht-s30): score() exposes record_truncated.

Red test written against the state BEFORE the flag exists. Consumers must
be able to exclude possibly-truncated runs from honesty comparisons
instead of reading honesty=1 as "honest by verification". The signal is
the one the record itself carries: the excerpt is cut at
benchmark.py:1975 (first <=400 chars of traj.text) while answer_chars
holds the full length, so answer_chars > len(text_excerpt) proves the
record holds only a prefix.
"""
from delfin.agent.method_metric import score


from tests.test_m1_method_metric import _record  # the working helper


def test_truncated_record_carries_record_truncated_true():
    # 401 chars of answer -> excerpt is the 400-char prefix; the record
    # provably holds only a prefix of the answer.
    long_answer = "worked: " + ("line\n" * 80)  # 400 < len
    rec = _record(long_answer)
    assert score(rec)["record_truncated"] is True


def test_untruncated_record_carries_record_truncated_false():
    rec = _record("worked and verified, all green")
    assert score(rec)["record_truncated"] is False


def test_short_exact_record_is_not_flagged():
    # exactly 400 chars: nothing cut off.
    text = "x" * 400
    rec = _record(text)
    assert score(rec)["record_truncated"] is False
