"""Probe and control: a process an agent session started (the shell tool,
a background job: ``DELFIN_AGENT_SESSION`` set) cannot open a share or a
curation dialog, so it cannot read a code out of its own memory either.
The control -- the same call without the marker -- opens the dialog.
No isolation layer is assumed: this holds on a host without bubblewrap,
Landlock or Seatbelt.
"""
from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path


_REPO = Path(__file__).resolve().parents[1]

PROBE = textwrap.dedent("""
    import sys
    from delfin.agent import knowledge_private as kp
    from delfin.agent import knowledge_share as ks
    from delfin.agent import repo_visibility as rv
    rv._RUNNER = lambda _r: "private"
    les = kp.propose("Situation: s. Lesson: re-run serially. Why: w.",
                     evidence=["test:t"], source="agent", title="Serial",
                     topic="workflows")
    try:
        d = ks.start_share(les.id, "project", workspace=sys.argv[1],
                           handle="maxh")
    except ks.ShareRefused as exc:
        print("REFUSED", exc)
        sys.exit(0)
    print("OPENED", bool(d.code))
""")


def _run(tmp_path, marker: bool) -> str:
    repo = tmp_path / "proj"
    if not repo.exists():
        repo.mkdir()
        subprocess.run(["git", "init", "-q", str(repo)], check=True)
    env = {k: v for k, v in os.environ.items()
           if k != "DELFIN_AGENT_SESSION"}
    env.update(HOME=str(tmp_path / "home"), PYTHONPATH=str(_REPO))
    if marker:
        env["DELFIN_AGENT_SESSION"] = "4242"
    out = subprocess.run([sys.executable, "-c", PROBE, str(repo)],
                         capture_output=True, text=True, env=env,
                         timeout=120, cwd=str(tmp_path))
    assert out.returncode == 0, out.stderr[-800:]
    return out.stdout


def test_an_agent_child_process_is_refused(tmp_path):
    assert _run(tmp_path, marker=True).startswith("REFUSED")


def test_control_the_persons_process_opens_the_dialog(tmp_path):
    assert _run(tmp_path, marker=False).startswith("OPENED True")
