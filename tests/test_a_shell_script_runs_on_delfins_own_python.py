"""A script the agent writes runs on the interpreter DELFIN runs on.

Measured on a Raspberry Pi install: the agent service ran from its venv,
the shell tool's PATH held only /usr/local/bin:/usr/bin:/bin, so `python3`
was the system one without openpyxl. An agent that had written a correct
script to split a workbook asked for a pip install instead of running it.
The shell environment now carries the directory of sys.executable when the
PATH does not already name it; a PATH that does keeps its order.
"""

from __future__ import annotations

import os
import sys

from delfin.agent import api_client as A


def test_the_running_interpreter_is_on_a_bare_path(monkeypatch):
    monkeypatch.setenv("PATH", os.pathsep.join(["/usr/local/bin", "/usr/bin", "/bin"]))
    env = A._scrubbed_bash_env()
    first = env["PATH"].split(os.pathsep)[0]
    assert first == os.path.dirname(sys.executable)
    assert env["PATH"].endswith(os.pathsep.join(["/usr/local/bin", "/usr/bin", "/bin"]))


def test_a_path_that_already_names_it_keeps_its_order(monkeypatch):
    py_bin = os.path.dirname(sys.executable)
    own = os.pathsep.join(["/opt/cluster/env/bin", py_bin, "/usr/bin"])
    monkeypatch.setenv("PATH", own)
    assert A._scrubbed_bash_env()["PATH"] == own
