"""Analysis-tool version probes must not litter the working directory.

`censo --version` writes an empty `censo.log` into the directory it runs
in, and `c2anmr --version` creates an `anmr/` folder there.  The version
probes in `delfin/analysis_tools` start these CLIs without a cwd, so the
litter lands in whatever directory the caller sits in -- during a suite
run, the checkout itself (seen after tests/test_equatorial_square_kappa4.py
in SLURM 7199892: `<checkout>/anmr` and `<checkout>/censo.log` appeared).
"""

from __future__ import annotations

import os
import shutil

import pytest

from delfin.analysis_tools import _probe_cli_version, collect_analysis_summary


def _litter(directory):
    return sorted(p.name for p in directory.iterdir())


def _stub_tool(bin_dir, name, litter_script):
    """A stand-in binary reproducing what the real tools do: censo
    --version writes censo.log into its working directory, c2anmr
    --version creates an anmr/ folder there (verified by hand on this
    machine against the installed censo/c2anmr). Lives outside the
    directory the litter is checked in."""
    bin_dir.mkdir(exist_ok=True)
    tool = bin_dir / name
    tool.write_text(litter_script, encoding="utf-8")
    tool.chmod(0o755)
    return bin_dir


_CENSO_STUB = "#!/bin/sh\ntouch censo.log\necho 'censo 3.0.6'\n"
_C2ANMR_STUB = "#!/bin/sh\nmkdir -p anmr\necho 'c2anmr 1.0'\n"


def test_censo_version_probe_leaves_no_censo_log(tmp_path, tmp_path_factory, monkeypatch):
    bin_dir = _stub_tool(tmp_path_factory.mktemp("censo-stub"), "censo", _CENSO_STUB)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.chdir(tmp_path)
    version = _probe_cli_version(["censo", "--version"])
    assert version == "censo 3.0.6", "the probe did not read the version"
    assert _litter(tmp_path) == [], (
        f"the probe left { _litter(tmp_path) } in the working directory"
    )


def test_c2anmr_version_probe_leaves_no_anmr_dir(tmp_path, tmp_path_factory, monkeypatch):
    bin_dir = _stub_tool(tmp_path_factory.mktemp("c2anmr-stub"), "c2anmr", _C2ANMR_STUB)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.chdir(tmp_path)
    version = _probe_cli_version(["c2anmr", "--version"])
    assert version == "c2anmr 1.0", "the probe did not read the version"
    assert _litter(tmp_path) == [], (
        f"the probe left { _litter(tmp_path) } in the working directory"
    )


@pytest.mark.skipif(
    not (shutil.which("censo") or shutil.which("c2anmr")),
    reason="no littering analysis tool installed",
)
def test_analysis_summary_leaves_no_litter(tmp_path, monkeypatch):
    """collect_analysis_summary() probes every installed tool's version;
    none of those probes may write into the caller's directory."""
    monkeypatch.chdir(tmp_path)
    summary = collect_analysis_summary()
    assert summary["tools"], "the summary came back empty"
    assert _litter(tmp_path) == [], (
        f"the summary left { _litter(tmp_path) } in the working directory"
    )
