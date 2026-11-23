"""The suite runs where CI runs: not inside an agent session.

``approval_guard.ENV_MARKER`` (``DELFIN_AGENT_SESSION``) is the mark every
agent child process carries. Six product modules read it to decide that a
person, not an agent, is acting: ``knowledge_share`` refuses to share or
curate, the trust store refuses a write, ``file_confirm`` answers as the
person, ``messenger.pairing`` and the service supervisor read a child as
caged, and ``_classifier_settings`` ignores ``DELFIN_CLASSIFIER=off``.

A DELFIN agent session exports it to everything it starts, so a suite run
from one inherits it -- and then measures "am I inside an agent session"
where it meant to measure the module. Measured 2026-11-20 on a workstation
whose shell carried ``DELFIN_AGENT_SESSION=985466``, on the commit CI runs
green: ``test_zz_u12_second_opinion`` failed on the classifier switch (a
nested agent ignores ``DELFIN_CLASSIFIER=off`` by design) and
``test_knowledge_agent_surface`` failed on sharing ("this runs inside an
agent session"), both passing under ``env -u``.

A host that answers for the code is the class ``_ci_parity_no_secrets``
and ``_no_host_settings_defaults`` already close for keys and the settings
file; this is the third of them, in ``conftest.py``.

The tests below protect two things: a test body sees no AMBIENT marker, and
a test that means to be inside an agent session can still say so -- the
fixture narrows the host, it does not forbid the state. The second half is
what keeps the cage and scratch tests (which set the marker themselves)
measuring anything.
"""

from __future__ import annotations

import os

import conftest


MARKER = conftest.AGENT_SESSION_MARKER


def test_a_test_body_sees_no_ambient_agent_session_marker():
    """The autouse fixture has already run: this process carries no marker.

    This is the property the six product modules depend on. It fails on a
    host whose shell exports the marker and passes on CI, which is the
    whole point: the suite measures the module, not where it is running.
    """
    assert os.environ.get(MARKER) is None


def test_a_test_body_sees_no_ambient_session_scratch_directory():
    """The other variable an agent session exports, same rule.

    ``DELFIN_SESSION_SCRATCH`` names the cage's session scratch directory.
    A test that asserts ``_scrubbed_bash_env`` does not export it when no
    scratch can be made is true on CI and false in a shell that already
    carries one (measured 2026-11-20:
    ``test_sx2_session_scratch_export`` failed here and passed under
    ``env -u DELFIN_SESSION_SCRATCH``).
    """
    assert os.environ.get("DELFIN_SESSION_SCRATCH") is None


def test_the_variables_removed_are_the_ones_an_agent_session_exports():
    """The list here and the product's own names are one set.

    A fixture that removed a misspelt variable would pass the tests above
    and leave every reader measuring the host. Both names are taken from
    the product: ``approval_guard.ENV_MARKER`` for the marker, and the bash
    tool's own export site for the scratch directory.
    """
    from delfin.agent import approval_guard

    assert conftest.AGENT_SESSION_MARKER == approval_guard.ENV_MARKER
    assert conftest.AGENT_SESSION_MARKER in conftest._CAGE_SESSION_VARS
    assert "DELFIN_SESSION_SCRATCH" in conftest._CAGE_SESSION_VARS


def test_a_test_that_means_to_be_inside_an_agent_session_can_still_say_so(
        monkeypatch):
    """The fixture narrows the ambient host, it does not forbid the state.

    The cage tests measure what happens INSIDE an agent session
    (``test_an_agent_session_cannot_act_as_the_operator.py`` and the
    scratch tests set the marker themselves). A fixture that refused the
    marker outright would silence them; one that only removes the host's
    copy leaves them free to set their own -- and this test's own
    ``setenv`` runs after the autouse fixture.
    """
    monkeypatch.setenv(MARKER, "4242")
    assert os.environ.get(MARKER) == "4242"
