"""A job runs the venv as it is now, not as it was when a tar was made.

``delfin_venv.tar`` was packed once, by an installer or a settings button, and
unpacked by every job after it. Months later pip had put fifty-five packages
into the venv that the tar did not have, and nothing said so. The tar is a
cache now, named by what the venv holds, so it changes exactly when the venv
does.
"""

import os
import pathlib
import subprocess
import time

REPO = pathlib.Path(__file__).resolve().parents[1]
TEMPLATE = REPO / "delfin" / "submit_templates" / "submit_delfin.sh"


def _functions() -> str:
    text = TEMPLATE.read_text(encoding="utf-8")
    return text[text.index("venv_cache_key() {"):text.index('VENV_LOCAL=""\nVENV_SRC=""')]


def _bash(tmp_path, body):
    script = tmp_path / "run.sh"
    script.write_text("set -euo pipefail\n" + _functions() + "\n" + body + "\n")
    return subprocess.run(
        ["bash", str(script)], capture_output=True, text=True, timeout=60,
        env={
            "HOME": str(tmp_path / "home"),
            "PATH": os.environ["PATH"],
            "DELFIN_VENV_CACHE_DIR": str(tmp_path / "cache"),
            "SLURM_JOB_ID": "7",
        },
    )


def _venv(tmp_path):
    venv = tmp_path / "anywhere" / "env"
    site = venv / "lib" / "python3.11" / "site-packages"
    (site / "numpy-1.26.4.dist-info").mkdir(parents=True)
    (venv / "bin").mkdir()
    (venv / "pyvenv.cfg").write_text("home = /usr/bin\n")
    return venv, site


def test_the_key_follows_what_is_installed_and_nothing_else(tmp_path):
    venv, site = _venv(tmp_path)

    def key():
        done = _bash(tmp_path, f'venv_cache_key "{venv}"')
        assert done.returncode == 0, done.stderr
        return done.stdout.strip()

    first = key()
    assert len(first) == 16
    assert key() == first

    # Byte-code appearing is not a change of what is installed.
    (site / "__pycache__").mkdir()
    assert key() == first

    (site / "openpyxl-3.1.5.dist-info").mkdir()
    assert key() != first


def test_a_changed_venv_gets_a_new_tar_and_old_ones_do_not_pile_up(tmp_path):
    venv, site = _venv(tmp_path)
    tars = []
    for package in ("a-1.dist-info", "b-1.dist-info", "c-1.dist-info"):
        (site / package).mkdir()
        done = _bash(tmp_path, f'ensure_venv_tar "{venv}"')
        assert done.returncode == 0, done.stderr
        assert "Packing" in done.stderr
        tars.append(pathlib.Path(done.stdout.strip()))
        time.sleep(1.1)

    assert len(set(tars)) == 3
    # The newest and the one before it: a job may still be unpacking that.
    kept = sorted(path.name for path in (tmp_path / "cache").glob("venv-*.tar"))
    assert kept == sorted([tars[1].name, tars[2].name])

    again = _bash(tmp_path, f'ensure_venv_tar "{venv}"')
    assert pathlib.Path(again.stdout.strip()) == tars[2]
    assert "Packing" not in again.stderr, "an unchanged venv was packed again"

    node = tmp_path / "node"
    node.mkdir()
    subprocess.run(["tar", "-xf", str(tars[2]), "--strip-components=1", "-C", str(node)], check=True)
    assert (node / "lib" / "python3.11" / "site-packages" / "c-1.dist-info").is_dir()


def test_the_job_is_not_pointed_at_a_hand_made_tar_in_the_checkout():
    text = TEMPLATE.read_text(encoding="utf-8")
    assert "$DELFIN_DIR/delfin_venv.tar" not in text
    # Where the dashboard runs from is said, not searched for.
    assert 'VENV_SRC="$DELFIN_VENV"' in text
    assert 'DELFIN_DIR="$DELFIN_REPO"' in text
