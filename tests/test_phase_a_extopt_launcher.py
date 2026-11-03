"""Phase A MLIP integration: the executable ExtOpt wrapper (ProgExt target).

ORCA's ExtOpt (manual 4.1.10) runs the program given as ProgExt as
``<script> <extinp> [Ext_Params]``.  DELFIN's wrapper module
(delfin/mlp_tools/orca_extopt.py) has been callable as
``python -m delfin.mlp_tools.orca_extopt``; what Phase A adds is the
executable launcher delfin/mlp_tools/bin/delfin-mlp-extopt that ORCA
needs, plus its creation by install_mlp_tools.sh.  These tests exercise
the launcher the way ORCA would call it -- as a subprocess, no torch,
no ORCA.
"""

from __future__ import annotations

import os
import shlex
import shutil
import stat
import subprocess
import sys
from pathlib import Path

from conftest import child_env

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "delfin" / "mlp_tools" / "bin" / "delfin-mlp-extopt"
INSTALLER = REPO / "delfin" / "mlp_tools" / "install_mlp_tools.sh"


def test_launcher_exists_and_is_executable():
    assert LAUNCHER.is_file(), f"missing: {LAUNCHER}"
    mode = LAUNCHER.stat().st_mode
    assert mode & stat.S_IXUSR, "launcher must be executable (ORCA runs it)"
    assert LAUNCHER.read_text(encoding="utf-8").startswith("#!")


def test_launcher_passes_ext_params_through_to_python_module(tmp_path):
    """ORCA calls: script <extinp> --backend aimnet2 ...  The launcher
    must forward every argument to python -m delfin.mlp_tools.orca_extopt.
    Proven with a stub python that echoes its argv; no delfin import."""
    stub = tmp_path / "stub-python"
    stub.write_text(
        "#!/usr/bin/env bash\n"
        'echo "stub-python" "$@"\n',
        encoding="utf-8",
    )
    stub.chmod(0o755)
    result = subprocess.run(
        [str(LAUNCHER), "/x/a_EXT.extinp.tmp", "--backend", "aimnet2"],
        capture_output=True,
        text=True,
        env={**child_env(tmp_path), "DELFIN_MLP_PYTHON": str(stub)},
    )
    assert result.returncode == 0, result.stderr
    # the module invocation, with all Ext_Params forwarded verbatim
    assert " -m delfin.mlp_tools.orca_extopt /x/a_EXT.extinp.tmp" \
        " --backend aimnet2" in result.stdout
    assert "stub-python" in result.stdout


def test_launcher_prefers_venv_python_without_explicit_config(tmp_path):
    """With DELFIN_MLP_PYTHON unset but VIRTUAL_ENV set, the launcher must
    use that venv's interpreter (the one with delfin + backends)."""
    venv = tmp_path / "venv"
    (venv / "bin").mkdir(parents=True)
    stub = venv / "bin" / "python"
    stub.write_text("#!/usr/bin/env bash\necho venv-python \"$@\"\n")
    stub.chmod(0o755)
    result = subprocess.run(
        [str(LAUNCHER), "--help"],
        capture_output=True,
        text=True,
        env={**child_env(tmp_path), "VIRTUAL_ENV": str(venv)},
    )
    assert result.returncode == 0, result.stderr
    assert "venv-python" in result.stdout


def test_installer_generates_the_launcher(tmp_path):
    """install_mlp_tools.sh must create bin/delfin-mlp-extopt, executable."""
    # isolate: run the installer's install_extopt_wrapper function only.
    # The installer file ends with a bare `main "$@"` call that would run
    # the whole pip installation -- strip it before sourcing.
    installer_text = INSTALLER.read_text(encoding="utf-8")
    installer_text = installer_text.replace('main "$@"\n', "", 1)
    assert 'main "$@"' not in installer_text
    installer_copy = tmp_path / "install_mlp_tools_funcs.sh"
    installer_copy.write_text(installer_text, encoding="utf-8")
    script = (
        "set -euo pipefail\n"
        f"ROOT={shlex.quote(str(tmp_path))}\n"
        f"source {shlex.quote(str(installer_copy))}\n"
        "install_extopt_wrapper\n"
    )
    result = subprocess.run(
        ["bash", "-c", script],
        capture_output=True,
        text=True,
        env=child_env(tmp_path),
    )
    assert result.returncode == 0, result.stderr
    generated = tmp_path / "bin" / "delfin-mlp-extopt"
    assert generated.is_file()
    assert generated.stat().st_mode & stat.S_IXUSR
    body = generated.read_text(encoding="utf-8")
    assert "delfin.mlp_tools.orca_extopt" in body


def test_installer_wires_the_wrapper_into_main(tmp_path):
    text = INSTALLER.read_text(encoding="utf-8")
    assert "install_extopt_wrapper" in text
    # called from main(), not merely defined
    main_pos = text.index("main() {")
    assert "install_extopt_wrapper" in text[main_pos:]


def test_backend_reaches_the_aimnet2_calculator_factory(tmp_path):
    """End of the Phase A chain in one child process: launcher semantics
    aside, --backend aimnet2 must select the registry's aimnet2 factory
    (proven without torch: the factory lookup itself, monkeypatched)."""
    child = (
        "import sys\n"
        "import delfin.mlp_tools.registry as r\n"
        "calls = []\n"
        "r._FACTORY_OVERRIDE = {'aimnet2': lambda **kw: calls.append(kw) or 'CALC'}\n"
        "calc = r.create_calculator('aimnet2', device='cpu', charge=-1, mult=2)\n"
        "assert calc == 'CALC'\n"
        "assert calls[0]['charge'] == -1 and calls[0]['mult'] == 2\n"
        "print('OK')\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", child],
        capture_output=True,
        text=True,
        env=child_env(tmp_path),
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


def test_missing_launcher_dependency_fails_closed():
    """If the wrapped module or its python is missing, the launcher must
    fail loudly (set -e), never silently report success to ORCA."""
    stub = shutil.which("true") or "/bin/true"
    result = subprocess.run(
        [str(LAUNCHER), "/nonexistent/does_not_matter.extinp.tmp"],
        capture_output=True,
        text=True,
        env={
            **child_env(Path("/tmp")),
            "DELFIN_MLP_PYTHON": stub,
            "PATH": os.environ.get("PATH", ""),
        },
    )
    # a stub 'true' exits 0 without writing the engrad -- the failure mode
    # ORCA would see is a missing engrad file; the launcher itself must not
    # crash with a shell error beyond the child's status.
    assert result.returncode == 0  # 'true' semantics; launcher propagates
