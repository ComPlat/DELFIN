"""EF1 phase-2 red tests: cluster_efficiency.advise().

The suggestion must be built ONLY from DELFIN's own resource logic:
cluster_utils.auto_configure_resources (PAL/maxcore sizing) and
slurm_submit's partition machinery (earliest start = hand SLURM all
matching partitions, never resize to start earlier). Tests inject the
detected cluster info and the partition source, so nothing depends on
the machine the test runs on.
"""
import pytest

from delfin.agent import cluster_efficiency as ce


CLUSTER = {"scheduler": "slurm", "cpus_available": 64,
           "memory_available_mb": 262144, "job_id": None}


def _sizer(config):
    # Stand-in for the real auto_configure_resources, shaped like it.
    from delfin.cluster_utils import auto_configure_resources
    return auto_configure_resources(config)


def test_advise_sizes_from_cluster_utils_logic(monkeypatch):
    calls = []

    def fake_auto(config):
        calls.append(dict(config))
        return {**config, "PAL": 16, "maxcore": 3000}

    monkeypatch.setattr(ce, "_auto_configure", fake_auto)
    adv = ce.advise("opt", 24)
    assert adv["PAL"] == 16
    assert adv["maxcore_mb"] == 3000
    assert calls == [{}]  # the sizer got the empty config it expects


def test_advise_lists_all_matching_partitions_for_earliest_start(monkeypatch):
    monkeypatch.setattr(ce, "_auto_configure",
                        lambda config: {**config, "PAL": 8, "maxcore": 2000})
    monkeypatch.setattr(ce, "_configured_partitions",
                        lambda explicit=None, profile=None: ("cpu", "cpu_il"))
    adv = ce.advise("sp", 10)
    assert adv["partitions"] == "cpu,cpu_il"
    assert "earliest" in adv["partition_rule"]


def test_advise_never_resizes_to_start_earlier(monkeypatch):
    # Even when more partitions are free, PAL/maxcore stay the sizer's
    # values — resources are never changed just to start earlier.
    monkeypatch.setattr(ce, "_auto_configure",
                        lambda config: {**config, "PAL": 32, "maxcore": 4000})
    monkeypatch.setattr(ce, "_configured_partitions",
                        lambda explicit=None, profile=None: ("cpu", "cpu_il"))
    adv = ce.advise("opt", 50)
    assert adv["PAL"] == 32
    assert adv["maxcore_mb"] == 4000
    assert "never" in adv["partition_rule"]


def test_advise_reports_unknown_step_and_size_bounds():
    monkeypatch = pytest.MonkeyPatch()
    try:
        monkeypatch.setattr(ce, "_auto_configure",
                            lambda config: {**config, "PAL": 4,
                                            "maxcore": 1000})
        adv = ce.advise("totally_unknown_step", 3)
        # A workflow DELFIN does not know is said so, not silently sized.
        assert "unknown step" in adv["notes"][0]
    finally:
        monkeypatch.undo()


def test_advise_uses_delfin_step_vocabulary(monkeypatch):
    known = {"opt", "sp", "freq", "ts", "tddft", "gopt"}
    monkeypatch.setattr(ce, "_auto_configure",
                        lambda config: {**config, "PAL": 4, "maxcore": 1000})
    for step in sorted(known):
        adv = ce.advise(step, 12)
        assert not any("unknown step" in n for n in adv["notes"]), step
