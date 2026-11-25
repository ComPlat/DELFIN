"""Install.cmd and Uninstall.cmd lead a person through, not into an error.

Measured on Windows 11 before this change:
- the installed Uninstall.cmd changed into the installed folder; Windows does
  not remove a folder that is a process's current directory, and
  Remove-Item -Recurse then removed nothing at all;
- an update over a running app replaced files under it (and an earlier
  version's DELFIN.exe cannot be removed while it runs);
- the uninstall offered only "everything" or "nothing", and Windows'
  Installed apps ran Uninstall.ps1 in a window that closed at once.

The behaviour itself runs on Windows in Test-InstallRoundTrip.ps1 (CI:
windows-launcher.yml): install, update, a running app, uninstall with each
answer, from inside the installed folder, in a folder of its own. These are
the text-level guards for the Linux test run.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "tools" / "windows-launcher"
INSTALL_CMD = (PKG / "Install.cmd").read_text(encoding="utf-8")
UNINSTALL_CMD = (PKG / "Uninstall.cmd").read_text(encoding="utf-8")
INSTALL = (PKG / "Install.ps1").read_text(encoding="utf-8")
UNINSTALL = (PKG / "Uninstall.ps1").read_text(encoding="utf-8")


def test_the_uninstaller_does_not_stand_in_the_folder_it_removes():
    assert 'cd /d "%~dp0"' not in UNINSTALL_CMD
    assert 'cd /d "%SystemRoot%"' in UNINSTALL_CMD
    assert '-File "%~dp0Uninstall.ps1"' in UNINSTALL_CMD
    i_away = UNINSTALL.index("Set-Location -LiteralPath $env:SystemRoot")
    assert i_away < UNINSTALL.index("Remove-Item -LiteralPath $destination -Recurse")
    assert "[Environment]::CurrentDirectory = $env:SystemRoot" in UNINSTALL


def test_the_uninstall_asks_once_with_keep_as_a_choice():
    assert 'choice /C YKN' in UNINSTALL_CMD
    # In a .cmd file "set" resets ERRORLEVEL: the answer is read in one test.
    assert 'if errorlevel 2 (set "keep=-KeepProfiles") else (set "keep=")' in UNINSTALL_CMD
    assert "%keep%" in UNINSTALL_CMD
    assert "Nothing was changed." in UNINSTALL_CMD


def test_a_removed_uninstall_cmd_ends_cleanly():
    """Without "(goto)" cmd looks for the next line of a file it just
    removed: "batch file cannot be found" and exit code 1 after a complete
    uninstall."""
    assert "(goto) 2>nul & exit /b 0" in UNINSTALL_CMD
    assert "(goto) 2>nul & exit /b 1" in UNINSTALL_CMD


def test_a_running_app_is_asked_about_before_anything_changes():
    for text, first_change in ((INSTALL, "Copy-Item -LiteralPath $source"),
                               (UNINSTALL, "Remove-Item -LiteralPath $destination -Recurse")):
        i_check = text.index("$running = @(Get-RunningApp $destination)")
        assert i_check < text.index(first_change)
        assert "Read-Host 'Close DELFIN now and continue? (y/N)'" in text
        assert "$answer -ne 'y' -and $answer -ne 'Y'" in text, "only an explicit yes closes it"


def test_installed_apps_opens_the_same_uninstall_window():
    assert "(Join-Path $destination 'Uninstall.cmd')" in INSTALL
    assert "System32\\cmd.exe" in INSTALL
    assert "NoModify=1; NoRepair=1" in INSTALL


def test_install_cmd_says_what_to_do_instead_of_failing():
    i_first_powershell = INSTALL_CMD.index("powershell.exe")
    for label, text in (("notextracted", "Extract the ZIP first"),
                        ("nossh", "OpenSSH Client")):
        assert f"goto {label}" in INSTALL_CMD[:i_first_powershell]
        assert f":{label}" in INSTALL_CMD and text in INSTALL_CMD
    assert "start \"\" ms-settings:optionalfeatures" in INSTALL_CMD
    # The developer check runs quietly; its errors still reach the window.
    assert re.search(r'Test-Launcher\.ps1" >nul', INSTALL_CMD)
    assert 'choice /C YN /N /M "Start DELFIN now? [Y/N]: "' in INSTALL_CMD


def test_windows_ci_runs_the_round_trip():
    workflow = (ROOT / ".github" / "workflows" / "windows-launcher.yml").read_text(encoding="utf-8")
    assert "./tools/windows-launcher/Test-InstallRoundTrip.ps1" in workflow
    assert "Test-InstallRoundTrip.ps1" not in (PKG / "build_package.py").read_text(encoding="utf-8"), (
        "a developer test, not part of what a person downloads")
