"""A bug report is the one artifact that deliberately leaves the machine.

``settings_snapshot`` promised "NEVER secrets" in its docstring and then
did ``dict(settings.get("agent") or {})`` -- the whole agent block, no
allow-list, so every field anyone ever adds to it shipped by default. The
settings file on the audited machine already carried an ``agent.api_key``
field: legacy, empty, read by nothing. Empty is luck, not a mechanism.

What happens to that report afterwards is what makes it matter. The writer
sets setgid on the report directory, chgrps every path to the archive's
maintainer group, adds group-read, and rsyncs the result into a shared
archive on a machine with six other accounts. The payload is the system
prompt, the raw traceback, every chat and engine message, and the tool
trace with its raw tool outputs.

``run_output_guards`` existed the whole time and had exactly one caller --
the final answer in the engine. Nothing inspected any of this.

Both halves are tested here: the allow-list, which decides what is
COLLECTED, and the guard, which decides what is WRITTEN. Neither alone is
enough. The allow-list cannot see a credential pasted into a traceback,
and the guard cannot see one in a settings field that is not named like a
secret -- so ``some_future_field`` below carries a real-looking token
under a perfectly innocent name.
"""

from __future__ import annotations

import hashlib
import json
import string

from delfin.agent import bug_report


def _token(n: int, seed: str) -> str:
    """A high-entropy stand-in. The redactor rejects repeated characters on
    purpose, so a lazy probe of "AAAA..." would prove nothing."""
    out: list[str] = []
    h = seed.encode()
    alphabet = string.ascii_letters + string.digits
    while len(out) < n:
        h = hashlib.sha256(h).digest()
        out += [alphabet[b % len(alphabet)] for b in h]
    return "".join(out[:n])


# One credential-shaped value at every depth the settings tree has.
_AT_DEPTH = {
    "agent_api_key": _token(40, "d1"),        # the field really on disk
    "nested_in_allowed": _token(40, "d2"),    # agent.routing.<x>.api_key
    "nested_twice": _token(40, "d3"),         # agent.subagents.<x>.<y>.token
    "unlisted_field": _token(40, "d4"),       # innocent NAME, secret value
    "in_a_list": _token(40, "d5"),            # inside a list value
    "transfer_password": _token(40, "d6"),
    "features_secret": _token(40, "d7"),
    "runtime_secret": _token(40, "d8"),
}

_KEPT = "qwen3-fixture-model"


def _settings() -> dict:
    return {
        "agent": {
            "api_key": _AT_DEPTH["agent_api_key"],
            "model": _KEPT,
            "backend": "openai",
            "routing": {"kit": {"api_key": _AT_DEPTH["nested_in_allowed"]}},
            "subagents": {
                "defaults": {"auth": {"token": _AT_DEPTH["nested_twice"]}},
                "profiles": [{"secret_key": _AT_DEPTH["in_a_list"]}],
            },
            "some_future_field": _AT_DEPTH["unlisted_field"],
        },
        "transfer": {
            "host": "cluster.example",
            "user": "someone",
            "remote_path": "/archive",
            "port": 22,
            "password": _AT_DEPTH["transfer_password"],
        },
        "features": {"beta": {"client_secret": _AT_DEPTH["features_secret"]}},
        "runtime": {"backend": "local",
                    "access_key": _AT_DEPTH["runtime_secret"]},
    }


# ---------------------------------------------------------------------------
# The snapshot collects an allow-list, not "everything under agent"
# ---------------------------------------------------------------------------

def test_the_snapshot_drops_every_credential_at_every_depth():
    snap = bug_report.settings_snapshot(_settings())
    blob = json.dumps(snap, ensure_ascii=False)
    for label, value in _AT_DEPTH.items():
        assert value not in blob, label


def test_the_snapshot_still_carries_what_a_maintainer_needs():
    """An allow-list that empties the report is not a fix, it is a
    different bug: the snapshot exists so a bad turn can be reproduced."""
    snap = bug_report.settings_snapshot(_settings())
    assert snap["agent"]["model"] == _KEPT
    assert snap["agent"]["backend"] == "openai"
    assert snap["transfer"]["host"] == "cluster.example"
    assert snap["runtime_backend"] == "local"


def test_a_field_nobody_listed_is_absent_rather_than_published():
    """The direction of the list is the whole point. This field is named
    like ordinary config and holds a token; no redactor keyed on names can
    see it, so only refusing to collect it works."""
    snap = bug_report.settings_snapshot(_settings())
    assert "some_future_field" not in snap["agent"]
    assert "api_key" not in snap["agent"]


# ---------------------------------------------------------------------------
# ...and the guard runs over what is actually written
# ---------------------------------------------------------------------------

def _write(tmp_path, **kw):
    return bug_report.write_bug_report(
        chat_messages=kw.pop("chat_messages", []),
        settings=_settings(),
        archive_dir=tmp_path / "archive",
        session_id="sess1234",
        **kw,
    )


def test_no_credential_reaches_either_rendered_file(tmp_path):
    """The report.json/report.md pair is what gets rsynced."""
    d = _write(tmp_path)
    rendered = ((d / "report.json").read_text(encoding="utf-8")
                + (d / "report.md").read_text(encoding="utf-8"))
    for label, value in _AT_DEPTH.items():
        assert value not in rendered, label


def test_a_credential_in_the_traceback_is_redacted(tmp_path):
    """The commonest real shape: a key quoted back inside an error."""
    leaked = _token(48, "trace")
    d = _write(
        tmp_path,
        error_text=f"HTTPError 401: OPENAI_API_KEY={leaked} rejected",
    )
    for name in ("report.json", "report.md"):
        assert leaked not in (d / name).read_text(encoding="utf-8"), name


def test_a_credential_in_the_conversation_is_redacted(tmp_path):
    leaked = _token(48, "chat")
    d = _write(tmp_path, chat_messages=[
        {"role": "user", "content": "why does this fail?"},
        {"role": "assistant", "content": f"your token is Bearer {leaked}"},
    ])
    for name in ("report.json", "report.md"):
        assert leaked not in (d / name).read_text(encoding="utf-8"), name


def test_a_credential_in_the_system_prompt_is_redacted(tmp_path):
    leaked = "sk-ant-" + _token(40, "sys")
    d = _write(tmp_path, system_prompt=f"context\napi key: {leaked}\n")
    for name in ("report.json", "report.md"):
        assert leaked not in (d / name).read_text(encoding="utf-8"), name


def test_the_bundled_tool_trace_is_guarded_too(tmp_path, monkeypatch):
    """The tool trace is the file the disk-side finding is about: raw tool
    output, written verbatim, then bundled and given group-read."""
    from delfin.agent import tool_trace

    leaked = _token(48, "tt")
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "traces")
    tool_trace.record("sess1234", tool="bash", tool_input="printenv",
                      output=f"KIT_TOOLBOX_API_KEY={leaked}", ok=True)
    d = _write(tmp_path, trace_session="sess1234")
    shipped = d / "tool_trace.jsonl"
    assert shipped.is_file()
    assert leaked not in shipped.read_text(encoding="utf-8")


def test_the_redacted_report_is_still_machine_readable(tmp_path):
    """report.json is the maintainer-tooling half of the pair. Redacting
    the rendered text rather than the object would be a silent trap if a
    replacement could land across a string boundary."""
    d = _write(
        tmp_path,
        error_text=f"OPENAI_API_KEY={_token(48, 'j')} rejected",
        chat_messages=[{"role": "user", "content": "check this"}],
    )
    payload = json.loads((d / "report.json").read_text(encoding="utf-8"))
    assert payload["schema"] == "delfin-bug-report/1"
    assert payload["session_id"] == "sess1234"
    assert "redacted" in payload["error_text"]


def test_the_tool_trace_stays_one_json_object_per_line(tmp_path,
                                                       monkeypatch):
    from delfin.agent import tool_trace

    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "traces")
    tool_trace.record("sess1234", tool="bash", tool_input="printenv",
                      output=f"api_key={_token(40, 'k')}", ok=True)
    d = _write(tmp_path, trace_session="sess1234")
    for line in (d / "tool_trace.jsonl").read_text(
            encoding="utf-8").splitlines():
        if line.strip():
            json.loads(line)


def test_the_report_says_when_it_redacted_something(tmp_path):
    """A maintainer reading a report needs to know a value was removed,
    rather than concluding the field was empty."""
    d = _write(tmp_path, error_text=f"api_key={_token(40, 'note')}")
    note = d / "REDACTIONS.txt"
    assert note.is_file()
    assert "redacted" in note.read_text(encoding="utf-8").lower()


def test_an_ordinary_report_is_not_littered_with_a_redaction_note(tmp_path):
    d = _write(tmp_path, error_text="ValueError: expected 3 columns, got 2")
    assert not (d / "REDACTIONS.txt").exists()


def test_the_report_survives_a_broken_guard(tmp_path, monkeypatch):
    """Losing the bug report is worse than shipping it unguarded, so the
    guard is best-effort. This pins that it stays best-effort."""
    import delfin.agent.output_guard as og

    def boom(*a, **kw):
        raise RuntimeError("guard exploded")

    monkeypatch.setattr(og, "run_output_guards", boom)
    d = _write(tmp_path, error_text="something went wrong")
    assert (d / "report.json").is_file()
    assert (d / "report.md").is_file()
