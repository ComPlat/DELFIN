"""The line-number lookup `sed -n "$(grep -n PAT F | cut -d: -f1),+Np" F` is a read.

It was the largest group of avoidable dialogs in two audits (2026-09-25,
2026-09-26). Command substitution stays with the confirm gate in general,
because a substitution's OUTPUT becomes code of the outer program. This
one form prints line numbers and nothing else, so it is judged as the
number it yields. Every variation of it still asks.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import KitToolPermissions


@pytest.fixture
def perms(tmp_path):
    return KitToolPermissions(workspace=tmp_path)


FREE = [
    """sed -n "$(grep -n 'def _default_log_path' delfin/agent/audit_log.py | cut -d: -f1),+15p" delfin/agent/audit_log.py""",
    """grep -n "x" a.py | head -3; sed -n "$(grep -n 'def x' a.py | cut -d: -f1),+15p" a.py""",
    """sed -n "$(grep -n '_PREFIXES\\s*=' e.py | head -1 | cut -d: -f1),+6p" e.py""",
]

ASKED = [
    """sed -n "$(cat f)p" f""",
    """sed -n "$(grep -n x f | sh)p" f""",
    """sed -n "$(grep -n x f; rm y)p" f""",
    """sed -n "$(grep -n "$(id)" f | cut -d: -f1)p" f""",
    """sed -n "$(grep -n x f | cut -d: -f1)w out" f""",
    """rm -rf "$(grep -n x f | cut -d: -f1)" """,
    """sed -n "$(grep -n x f | cut -d: -f2)p" f""",
]


@pytest.mark.parametrize("cmd", FREE)
def test_the_lookup_form_runs_free(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd)


@pytest.mark.parametrize("cmd", ASKED)
def test_every_other_substitution_is_asked_about(perms, cmd):
    assert not perms.matches_bash_auto_allow(cmd)
