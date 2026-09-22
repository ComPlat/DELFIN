"""A credential is not written into memory, another session's inbox or a
tool result.

Memory is read back into every later session's prompt (the user-wide
store in every project) and sent to that session's provider. A session
message lands in another session that may run on another provider. A KIT
key has no prefix the shape checks know, so only its exact value finds it.
"""
import json

import pytest

from delfin.agent import api_client as A
from delfin.agent import memory_store as MS
from delfin.agent import output_guard as OG
from delfin.agent import session_messages as SM

KIT_KEY = "a1b2c3d4e5f6a7b8c9d0e1f2a3b4c5d6"          # no recognisable prefix
GH = "ghp_" + "A" * 36


@pytest.fixture
def held_key(monkeypatch):
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", KIT_KEY)


def test_the_exact_value_and_the_shapes_are_scrubbed(held_key):
    out = OG.scrub_secrets(f"key {KIT_KEY} and {GH}, rest stays")
    assert KIT_KEY not in out and GH not in out
    assert "[redacted:credential]" in out and "rest stays" in out


def test_scrubbing_leaves_ordinary_text_alone(held_key):
    text = "B3LYP/def2-SVP converged in 14 cycles"
    assert OG.scrub_secrets(text) == text


def test_a_typed_memory_keeps_no_credential(held_key, tmp_path, monkeypatch):
    monkeypatch.setattr(MS, "_delfin_memory_dir", lambda root: tmp_path / "mem")
    path, slug, mtype = MS.save_typed_memory(
        f"The KIT key is {KIT_KEY}, use it for every run",
        repo_root=tmp_path, title=f"key {GH}")
    written = path.read_text() + (path.parent / "MEMORY.md").read_text()
    assert KIT_KEY not in written and GH not in written


def test_a_legacy_fact_keeps_no_credential(held_key, tmp_path):
    p = tmp_path / "facts.json"
    MS.save_memory(f"token {GH}", path=p)
    assert GH not in p.read_text()


def test_a_session_message_passes_on_no_credential(held_key, tmp_path, monkeypatch):
    monkeypatch.setattr(SM, "_DIR", tmp_path)
    SM.send("abc", f"here is the key: {KIT_KEY}", from_key="me", from_title=GH)
    raw = "".join(f.read_text() for f in tmp_path.iterdir())
    assert KIT_KEY not in raw and GH not in raw
    got = SM.take("abc")
    assert got and "[redacted:credential]" in got[0]["text"]


def test_a_tool_result_carries_no_held_value(held_key):
    out = A._redact_tool_result(json.dumps({"stdout": f"KEY={KIT_KEY}"}))
    assert KIT_KEY not in out
