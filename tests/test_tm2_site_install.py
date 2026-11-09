"""Site-install defects from the TM2 real run #608 (operator, 08.11.).

Two findings, each with a red control driven through the public path:

1. The ``define`` dialog can break (non-zero exit) — the old
   ``_setup_control`` never checked its returncode, so the adapter went on
   to run ridft on a missing/partial ``control`` and the failure surfaced
   as a misleading SCF/convergence error. The adapter must report the
   ``define`` failure itself, before any SCF program starts.

2. ``detect_installed_version`` returned None on the site install: the
   TURBODIR path carries no version digits there and no explicit
   module name is passed. The environment module system leaves the
   resolved module name behind in ``LOADEDMODULES``/``_LMFILES_`` —
   detection must read it instead of giving up (never assume a version,
   but also never report unknown when the evidence is present).
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from delfin.tools._registry import get
from delfin.turbomole_version import detect_installed_version


def _make_bin(bin_dir: Path, name: str, body: str) -> Path:
    p = bin_dir / name
    p.write_text(body)
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    return p


def _xyz(work_dir: Path) -> Path:
    p = work_dir / "mol.xyz"
    p.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.8 0.0\n")
    return p


class TestDefineFailureIsSurfaced:
    """A broken define dialog is a define failure, not an SCF failure."""

    def test_define_nonzero_exit_fails_before_ridft(self, tmp_path, monkeypatch):
        work_dir = tmp_path / "work"
        work_dir.mkdir()
        bin_dir = tmp_path / "fakebin"
        bin_dir.mkdir()
        # define breaks (real site: dialog aborts, non-zero exit, stderr).
        _make_bin(bin_dir, "define",
                  "#!/bin/sh\necho 'define: input error' >&2\nexit 1\n")
        # If the adapter wrongly continues, ridft would "converge" quietly:
        _make_bin(bin_dir, "ridft",
                  "#!/bin/sh\nprintf '1  -1.0\\n$end\\n' > energy\n")
        monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
        adapter = get("turbomole_sp")
        assert adapter is not None
        result = adapter.execute(work_dir, geometry=_xyz(work_dir),
                                 cores=1, charge=0)
        assert result.status.name == "FAILED"
        # The error names define, not ridft/SCF (the old misreport).
        assert "define" in (result.error or "").lower(), result.error
        # ridft must never have run: no energy file was produced by it.
        assert not (work_dir / "energy").is_file(), (
            "the adapter ran ridft although define failed")

    def test_define_ok_still_runs_ridft(self, tmp_path, monkeypatch):
        """Nothing new is refused: a working define path still succeeds."""
        work_dir = tmp_path / "work2"
        work_dir.mkdir()
        bin_dir = tmp_path / "fakebin2"
        bin_dir.mkdir()
        _make_bin(bin_dir, "define",
                  f"#!/bin/sh\nprintf '$title\\n$coord\\n  0.0 0.0 0.0 o\\n$end\\n"
                  f"$atoms\\n$basis\\n$basis set\\n$dft\\n  func b3-lyp\\n$rij\\n"
                  f"$closed shells\\n' > {work_dir}/control\n")
        _make_bin(bin_dir, "ridft",
                  "#!/bin/sh\nprintf '1  -1.0\\n$end\\n' > energy\n")
        monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
        adapter = get("turbomole_sp")
        result = adapter.execute(work_dir, geometry=_xyz(work_dir),
                                 cores=1, charge=0)
        assert result.status.name == "SUCCESS", result.error


class TestVersionFromLoadedModules:
    """The module system's own records are a version source."""

    def test_loadedmodules_env_is_read(self, monkeypatch):
        # Site reality: TURBODIR without digits, but the module is loaded.
        monkeypatch.setenv("TURBODIR", "/opt/site/TURBOMOLE")
        monkeypatch.setenv("LOADEDMODULES", "chem/turbomole/7.9:intel/2022")
        ver = detect_installed_version()
        assert ver is not None, (
            "detect_installed_version must read LOADEDMODULES on the site "
            "install instead of returning None")
        assert str(ver) == "7.9"
        assert ver.source == "module"

    def test_lmfiles_env_is_read(self, monkeypatch):
        monkeypatch.delenv("LOADEDMODULES", raising=False)
        monkeypatch.setenv("TURBODIR", "/opt/site/TURBOMOLE")
        monkeypatch.setenv("_LMFILES_", "/opt/modules/chem/turbomole/7.9")
        ver = detect_installed_version()
        assert ver is not None
        assert str(ver) == "7.9"

    def test_explicit_module_name_still_wins_over_env(self, monkeypatch):
        monkeypatch.setenv("TURBODIR", "/opt/site/TURBOMOLE")
        monkeypatch.setenv("LOADEDMODULES", "chem/turbomole/7.9")
        ver = detect_installed_version(module_name="chem/turbomole/8.1")
        assert str(ver) == "8.1"

    def test_no_evidence_still_returns_none(self, monkeypatch):
        for name in ("TURBODIR", "LOADEDMODULES", "_LMFILES_"):
            monkeypatch.delenv(name, raising=False)
        assert detect_installed_version() is None


class TestControlValidatedAfterDefine:
    """(a) A define that 'ends normally' but leaves a partial control
    file is a setup failure, never an SCF failure (real run #608: the
    out-of-step dialog produced a control with only $title / $actual
    step / $end, ridft ran 0 s and the run was reported as 'SCF did not
    converge'). The adapter must validate the control file and name the
    missing sections."""

    def _fake_define(self, bin_dir: Path, work_dir: Path, control_text: str) -> None:
        _make_bin(bin_dir, "define",
                  f"#!/bin/sh\nprintf '{control_text}' > {work_dir}/control\n")

    def _adapter(self, tmp_path, monkeypatch, work_dir, bin_dir):
        monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
        adapter = get("turbomole_sp")
        return adapter.execute(work_dir, geometry=_xyz(work_dir),
                               cores=1, charge=0)

    def test_partial_control_is_refused(self, tmp_path, monkeypatch):
        work_dir = tmp_path / "w_partial"
        work_dir.mkdir()
        bin_dir = tmp_path / "b_partial"
        bin_dir.mkdir()
        # Exactly what the real run produced: no $coord, no $atoms,
        # no $basis, no $dft.
        self._fake_define(bin_dir, work_dir,
                          "$title\\na coord\\n$actual step      ridft\\n$end\\n")
        result = self._adapter(tmp_path, monkeypatch, work_dir, bin_dir)
        assert result.status.name == "FAILED"
        err = (result.error or "").lower()
        assert "define" in err, result.error
        assert "$coord" in err, ("the error must name what is missing: "
                                 + str(result.error))
        assert not (work_dir / "energy").is_file(), "ridft ran on a partial control"

    def test_complete_control_passes(self, tmp_path, monkeypatch):
        work_dir = tmp_path / "w_complete"
        work_dir.mkdir()
        bin_dir = tmp_path / "b_complete"
        bin_dir.mkdir()
        self._fake_define(
            bin_dir, work_dir,
            "$title\\n$coord\\n  0.0 0.0 0.0 o\\n$end\\n$atoms\\n$basis\\n"
            "b def2-SVP\\n$dft\\n  func b3-lyp\\n$rij\\n$closed shells\\n")
        _make_bin(bin_dir, "ridft",
                  "#!/bin/sh\nprintf '1  -1.0\\n$end\\n' > energy\n")
        result = self._adapter(tmp_path, monkeypatch, work_dir, bin_dir)
        assert result.status.name == "SUCCESS", result.error


class TestDefineDialogAnswerSync:
    """(b) define's first two prompts (default-data, title) must receive
    BLANK answers; the old dialog sent them no blank line at all, so
    'a coord' was eaten by the title prompt and every later answer was
    out of step (real run #608)."""

    def test_first_two_answers_are_blank(self, tmp_path, monkeypatch):
        work_dir = tmp_path / "w_sync"
        work_dir.mkdir()
        bin_dir = tmp_path / "b_sync"
        bin_dir.mkdir()
        # The fake define records every stdin line it received.
        _make_bin(bin_dir, "define",
                  "#!/bin/sh\ncat > stdin.log\n"
                  f"printf '$title\\n' > {work_dir}/control\n")
        _make_bin(bin_dir, "ridft",
                  "#!/bin/sh\nprintf '1  -1.0\\n$end\\n' > energy\n")
        # A control define leaves is not enough for validation, so the
        # run fails later — the dialog answers are the point here.
        stdin_log = work_dir / "stdin.log"
        monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
        get("turbomole_sp").execute(work_dir, geometry=_xyz(work_dir),
                                    cores=1, charge=0)
        lines = stdin_log.read_text().splitlines()
        assert len(lines) >= 3, "the dialog sent too few answers: " + str(lines)
        assert lines[0] == "", (
            "answer 1 (default-data prompt) must be BLANK, got %r" % lines[0])
        assert lines[1] == "", (
            "answer 2 (title prompt) must be BLANK, got %r " % lines[1]
            + "-- a non-blank answer here is out of step (real run #608)")
        assert lines[2].strip() == "a coord", lines[2]


class TestMissingWorkDirIsFailed:
    """(c) execute() on a work directory that does not exist must return
    FAILED through the public path, not raise FileNotFoundError."""

    def test_missing_work_dir_returns_failed(self, tmp_path, monkeypatch):
        work_dir = tmp_path / "does_not_exist"
        monkeypatch.setenv("PATH", os.environ.get("PATH", ""))
        result = get("turbomole_sp").execute(work_dir, geometry=None, cores=1)
        assert result.status.name == "FAILED"


class TestDefineSubmenuExitAndOpenShell:
    """(round 3, operator real-run 09.11., transcript
    .gate/rr608_define_transcript_working.txt) — on site TURBOMOLE 7.9
    the dft and ri submenus are LEFT WITH A BLANK line, not '*'; 'on'
    switches each on; the ri prompt does not ask a cores line; and for
    mult >= 3 the eht occupation needs 'n' then 'u <n_unpaired>'."""

    def _stdin_lines(self, tmp_path, monkeypatch, **kw):
        tag = kw.pop("tag")
        work_dir = tmp_path / ("w_" + tag)
        work_dir.mkdir()
        bin_dir = tmp_path / ("b_" + tag)
        bin_dir.mkdir()
        _make_bin(bin_dir, "define",
                  "#!/bin/sh\ncat > stdin.log\n"
                  "printf '$title\\n$coord\\n  0.0 0.0 0.0 o\\n$end\\n"
                  "$atoms\\n$basis\\n$basis set\\n$dft\\n  func b3-lyp\\n$rij\\n"
                  "$closed shells\\n' > control\n")
        _make_bin(bin_dir, "ridft",
                  "#!/bin/sh\nprintf '1  -1.0\\n$end\\n' > energy\n")
        stdin_log = work_dir / "stdin.log"
        monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep
                           + os.environ.get("PATH", ""))
        get("turbomole_sp").execute(work_dir, geometry=_xyz(work_dir),
                                    cores=2, charge=kw.get("charge", 0),
                                    mult=kw.get("mult", 1))
        return stdin_log.read_text().splitlines()

    def test_dft_ri_submenus_left_with_blank_line(self, tmp_path, monkeypatch):
        lines = self._stdin_lines(tmp_path, monkeypatch, tag="subm")
        # 'dft' -> 'on' -> 'func <m>' -> BLANK leaves the submenu.
        i = lines.index("dft")
        assert lines[i + 1] == "on", lines
        assert lines[i + 2].startswith("func "), lines
        assert lines[i + 3] == "", (
            "the dft submenu is left with a BLANK line on 7.9, got %r"
            % lines[i + 3])
        # 'ri' -> 'on' -> BLANK leaves it; no cores line is asked.
        j = lines.index("ri")
        assert lines[j + 1] == "on", lines
        assert lines[j + 2] == "", (
            "the ri submenu is left with a BLANK line and is NOT asked a "
            "cores line on 7.9, got %r" % lines[j + 2])

    def test_mult3_eht_occupation_uses_n_then_u(self, tmp_path, monkeypatch):
        lines = self._stdin_lines(tmp_path, monkeypatch, tag="m3",
                                  mult=3, charge=1)
        # eht -> y -> charge -> 'n' -> 'u <n_unpaired>'.
        k = lines.index("eht")
        tail = lines[k + 1:]
        assert tail[0] == "y", lines
        assert tail[1] == "1", ("the charge from kwargs: " + str(lines))
        assert tail[2] == "n", (
            "for mult >= 3 the eht occupation answer is 'n', got %r" % tail[2])
        assert tail[3] == "u 2", (
            "mult 3 means 2 unpaired electrons: 'u 2', got %r" % tail[3])

    def test_mult2_plain_dialog_sets_uhf_automatically(self, tmp_path,
                                                       monkeypatch):
        """Round 4 (operator real-run, CH3 mult=2): the PLAIN dialog (no
        occupation step) already sets up UHF on 7.9 — the open-shell
        branch must fire only for mult >= 3. For mult 2 the dialog keeps
        the singlet's 'y' (accept the default occupation; define builds
        the alpha/beta shells itself), producing $uhf in the control —
        the old code sent 'n'/'u 1' and fell out of step (its control at
        .gate/rr608_ch3_control_round3.txt held no $dft/$rij)."""
        lines = self._stdin_lines(tmp_path, monkeypatch, tag="m2", mult=2)
        k = lines.index("eht")
        tail = lines[k + 1:]
        assert tail[0] == "y", lines
        assert tail[1] == "0", ("the charge from kwargs: " + str(lines))
        assert tail[2] == "y", (
            "mult 2 keeps the plain 'y' (define sets UHF itself on 7.9), "
            "got %r" % tail[2])
        assert "n" not in tail[:4], (
            "no occupation-decline 'n' for mult 2: " + str(lines))
