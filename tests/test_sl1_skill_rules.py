"""SL1 phase 3: the cluster-workflows skill carries the operational rules.

Front matter (domains/name/description) untouched; the sections carry
the throttle rule, the never-do list, the recalc routing and the
sacct artefact — each pinned so a later edit cannot silently drop one.
"""
from __future__ import annotations

from pathlib import Path

SKILL = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "skills" / "cluster-workflows.md"


def _text() -> str:
    return SKILL.read_text(encoding="utf-8")


def test_front_matter_declares_domain():
    text = _text()
    assert text.startswith("---")
    assert "domains: code" in text.split("---")[1]
    assert "name: cluster-workflows" in text.split("---")[1]


def test_read_order_section_names_the_four_entry_points():
    text = _text()
    for needle in ("advise(step, system_size)", "get_pending_start_times",
                   "cluster_efficiency.account", "agent/job_fix.py"):
        assert needle in text, needle


def test_throttle_rule_is_pinned():
    text = _text()
    assert "25 s" in text
    assert "backend_slurm.py:97" in text
    assert "never per id in a loop" in text


def test_never_do_list_is_pinned():
    text = _text()
    assert "never lower pal/maxcore/time" in text.replace(
        "never lower pal/maxcore/time just to", "never lower pal/maxcore/time")
    assert "Never raw" in text
    assert "loosen convergence to fit a queue" in text
    assert "/control path" in text


def test_wait_reason_ruling_is_pinned():
    text = _text()
    assert "ReqNodeNotAvail" in text
    assert "means WAIT" in text


def test_sacct_artefact_is_pinned():
    text = _text()
    assert "TotalCPU" in text
    assert "ARTEFACT" in text
    assert "CPU Efficiency" in text


def test_export_rule_is_pinned():
    text = _text()
    assert "--export=ALL" in text
    assert "never with explicit variables" in text


def test_auto_pal_cap_is_pinned():
    text = _text()
    assert "LOGIN_LOCAL_CORES" in text


def test_no_credential_paths_in_the_skill():
    text = _text().lower()
    for word in ("credential", "password", "token", "api_key", "apikey"):
        assert word not in text, word
