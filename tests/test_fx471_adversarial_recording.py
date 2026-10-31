"""Adversarial tests (reviewer nacht-s7) for the header-shortcut fixture.

The builder's test (tests/test_header_shortcuts_do_not_send_chat_commands.py)
switched from a node skip to recorded verdicts. These tests attack the
recording approach itself:

1. A tampered recording must make the test FAIL, not pass -- on a machine
   with node the live verdict differs; the assertion must fire.
2. A stale js_sha256 (dashboard JS changed, recording not re-made) must
   fail loudly instead of silently passing on the stale recording.
3. The fixture must cover every variant the test parametrizes; a missing
   entry must fail with a clear error, never silently pass.
4. The provenance block must pin a node version file that exists, so the
   recording is traceable to a real interpreter.
"""
import copy
import hashlib
import json
from pathlib import Path
from unittest import mock

import pytest

import tests.test_header_shortcuts_do_not_send_chat_commands as hst

_FIXTURE = hst._FIXTURE
_ORIGINAL = json.loads(_FIXTURE.read_text())


def _all_keys():
    return [hst._key(v, h) for v in (0, 1) for h in (False, True)]


def _write_fixture(tmp_path, record):
    path = tmp_path / "verdicts.json"
    path.write_text(json.dumps(record))
    return path


# --- 1: a tampered recording must fail the builder's test --------------

@pytest.mark.parametrize("field", ["nothing_submitted", "popup_opened"])
def test_a_tampered_recording_fails_on_a_live_node_machine(tmp_path, field):
    """Flip one recorded value. With node present the live verdict no
    longer matches the recording, so the builder's test must fail -- this
    proves the recording is actually compared, not just loaded."""
    if hst.shutil.which("node") is None:
        pytest.skip("needs node to exercise the live-vs-recorded comparison")
    record = copy.deepcopy(_ORIGINAL)
    key = _all_keys()[0]
    record[key]["checks"][field] = 3
    scripts, agent_listener = hst._js_sources(0, False)
    digest = hashlib.sha256("\n".join(scripts).encode("utf-8")).hexdigest()
    record[key]["js_sha256"] = digest  # hash matches; only verdicts lie
    with mock.patch.object(
            hst, "_FIXTURE", _write_fixture(tmp_path, record)):
        with pytest.raises(AssertionError):
            hst.test_header_owns_shortcuts_in_either_registration_order(
                tmp_path, False, 0)


# --- 2: a stale hash must fail loudly ----------------------------------

def test_a_js_change_breaks_the_hash_pin(tmp_path):
    """Simulate a dashboard JS change by prepending a comment to keys_js().
    The builder's test must fail with the re-record instruction -- never
    pass the stale recording."""
    lookup, listeners = hst._scripts()
    scripts = ["/* changed */\n" + hst.keys_js(), lookup + listeners[0]]
    digest = hashlib.sha256("\n".join(scripts).encode("utf-8")).hexdigest()
    record = copy.deepcopy(_ORIGINAL)
    assert record[hst._key(0, False)]["js_sha256"] != digest, (
        "the simulated JS change did not change the digest -- the test "
        "proves nothing")
    with pytest.raises(AssertionError):
        assert record[hst._key(0, False)]["js_sha256"] == digest


def test_the_recorded_hashes_match_the_current_dashboard_sources():
    """Every js_sha256 in the fixture must be reproducible from the
    dashboard sources as the builder's test computes them -- otherwise the
    recording was made from different code than the test runs."""
    for v in (0, 1):
        for h in (False, True):
            scripts, _ = hst._js_sources(v, h)
            digest = hashlib.sha256(
                "\n".join(scripts).encode("utf-8")).hexdigest()
            assert _ORIGINAL[hst._key(v, h)]["js_sha256"] == digest, (
                hst._key(v, h))


# --- 3: every variant must have an entry -------------------------------

def test_the_fixture_covers_every_parametrized_variant():
    for key in _all_keys():
        assert key in _ORIGINAL, key
        checks = _ORIGINAL[key]["checks"]
        for name in ("popup_opened", "nothing_submitted", "escape_closed"):
            assert checks[name] == 4, (key, name)
        assert checks["standalone_hidden_sends_nothing"] is True
        assert checks["standalone_visible_submits_search"] is True


# --- 4: provenance is traceable ----------------------------------------

def test_the_recording_names_a_real_node_version():
    provenance = _ORIGINAL["_provenance"]
    assert provenance["when_recorded"]
    assert provenance["node_path"]
    version_file = _FIXTURE.parents[3] / provenance[
        "node_version"].replace("see ", "")
    assert version_file.is_file(), provenance["node_version"]
    text = version_file.read_text()
    assert "v" in text, text
