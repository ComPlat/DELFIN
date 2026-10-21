"""Fault injection for the leave of a ``delfin-agent chat`` child.

Active only for a child whose PYTHONPATH starts with this directory
(``pty_chat.Chat(..., env=leave_shim_env(...))``); the interpreter imports
``sitecustomize`` from sys.path at start-up. ``LEAVE_SHIM`` selects:

  hang-after-save  the first session-end stage (the session report)
                   writes LEAVE_SHIM_DIR/stage-start and sleeps 120 s,
                   after the session is saved
  slow-stage       the same stage writes stage-start, sleeps 15 s,
                   writes stage-done
  slow-save        session_store.save_session writes LEAVE_SHIM_DIR/
                   save-start, sleeps 2 s, saves, writes save-done

Patched in the modules the session end looks up at call time, because
``python -m delfin.agent.cli`` runs the CLI as ``__main__``, whose own
functions a patch of ``delfin.agent.cli`` would not reach.
"""

import os
import time
from pathlib import Path

_MODE = os.environ.get("LEAVE_SHIM", "")
_DIR = Path(os.environ.get("LEAVE_SHIM_DIR", "."))

if _MODE == "hang-after-save":
    from delfin.agent import session_report

    def _hang(*_a, **_k):
        (_DIR / "stage-start").write_text("")
        time.sleep(120)

    session_report.write_session_report = _hang

elif _MODE == "slow-stage":
    from delfin.agent import session_report

    def _slow(*_a, **_k):
        (_DIR / "stage-start").write_text("")
        time.sleep(15)
        (_DIR / "stage-done").write_text("")

    session_report.write_session_report = _slow

elif _MODE == "slow-save":
    from delfin.agent import session_store

    _save = session_store.save_session

    def _slow_save(*a, **k):
        (_DIR / "save-start").write_text("")
        time.sleep(2.0)
        result = _save(*a, **k)
        (_DIR / "save-done").write_text("")
        return result

    session_store.save_session = _slow_save
