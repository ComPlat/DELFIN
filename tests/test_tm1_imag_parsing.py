"""Red test for F5: imaginary-frequency parsing must catch the 'i'-suffix
format aoforce prints for imaginary modes."""
from pathlib import Path

from delfin.workflows.contrib.w21_turbomole_workflow import TurbomoleWorkflow


def test_imaginary_i_suffix_counted(tmp_path):
    """A recorded aoforce-style line with an 'i' suffix must count."""
    log = tmp_path / "aoforce.out"
    log.write_text(
        "  frequency       1234.56\n"
        "  frequency       567.89\n"
        "  frequency      -123.45i\n"
        "  frequency      -45.67\n")
    wf = TurbomoleWorkflow()
    assert wf._count_imaginary(tmp_path) == 2


def test_plain_negative_still_counted(tmp_path):
    log = tmp_path / "aoforce.out"
    log.write_text("  frequency      -12.34\n")
    wf = TurbomoleWorkflow()
    assert wf._count_imaginary(tmp_path) == 1


def test_no_imaginary_counts_zero(tmp_path):
    log = tmp_path / "aoforce.out"
    log.write_text("  frequency       1234.56\n")
    wf = TurbomoleWorkflow()
    assert wf._count_imaginary(tmp_path) == 0


def test_missing_log_returns_none(tmp_path):
    wf = TurbomoleWorkflow()
    assert wf._count_imaginary(tmp_path) is None
