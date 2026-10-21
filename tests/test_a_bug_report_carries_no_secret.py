"""A bug report never copies a secret into the shared archive.

Security review 2026-09-16: the dashboard hands the report every path the
agent tried to touch, including reads the deny list refused, and the writer
copied them with shutil.copy2 into a group-readable archive. A denied
read_file("~/.ssh/id_rsa") put the real key there. It also attached the
outcome history, which holds the prompts of every session of the user.
"""
import json

from delfin.agent import bug_report as br


def _ws(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "run.py").write_text("print('ok')\n")
    (ws / ".env").write_text("API_TOKEN=abc123\n")
    (ws / "server.key").write_text("-----BEGIN PRIVATE KEY-----\nMII\n")
    (ws / "credentials.json").write_text(json.dumps({"k": "v"}))
    (ws / "settings.cfg").write_text("OPENAI_API_KEY=sk-proj-" + "A" * 40 + "\n")
    return ws


def test_secret_files_are_skipped_not_copied(tmp_path):
    ws = _ws(tmp_path)
    rep = tmp_path / "rep"; rep.mkdir()
    home_like = tmp_path / "home"
    (home_like / ".ssh").mkdir(parents=True)
    (home_like / ".ssh" / "id_rsa").write_text("PRIVATE")
    (home_like / ".bash_history").write_text("export TOKEN=x")
    recs = br._bundle_files(
        ["run.py", ".env", "server.key", "credentials.json",
         str(home_like / ".ssh" / "id_rsa"), str(home_like / ".bash_history")],
        rep, workspace=ws)
    by = {r["original"]: r["status"] for r in recs}
    assert by["run.py"] == "bundled"
    for secret in (".env", "server.key", "credentials.json", str(home_like / ".ssh" / "id_rsa")):
        assert by[secret] == "skipped-secret", (secret, by[secret])
    assert by[str(home_like / ".bash_history")] == "skipped-private-dotfile"
    copied = {p.name for p in (rep / "workspace").iterdir()}
    assert not any(n.startswith((".env", "server.key", "credentials", "id_rsa", ".bash")) for n in copied)


def test_a_token_inside_a_bundled_text_file_is_redacted(tmp_path):
    ws = _ws(tmp_path)
    rep = tmp_path / "rep"; rep.mkdir()
    recs = br._bundle_files(["settings.cfg"], rep, workspace=ws)
    assert recs[0]["status"] == "bundled"
    text = (rep / recs[0]["bundled"]).read_text()
    assert "sk-proj-" + "A" * 40 not in text


def test_the_report_carries_no_other_sessions_prompts():
    import inspect
    src = inspect.getsource(br)
    assert '"recent_outcomes": recent_outcomes()' not in src


def test_a_binary_file_is_listed_not_copied(tmp_path):
    """Text goes through the redactor; a binary file cannot, so it was
    copied as it was (a keyring blob, an sqlite store, a pickled session,
    a core file with a key in memory). It is not copied: the report names
    it with its size and sha256, enough to ask the reporter for it or to
    check that a copy received later is the same file."""
    import hashlib
    ws = _ws(tmp_path)
    blob = b"\x00\x01SQLite format 3\x00" + b"api_key=sk-proj-" + b"B" * 40
    (ws / "state.db").write_bytes(blob)
    d = br.write_bug_report(chat_messages=[], settings={},
                            archive_dir=tmp_path / "archive",
                            session_id="sess1234",
                            referenced_files=["state.db", "run.py"],
                            workspace=ws)
    meta = json.loads((d / "report.json").read_text(encoding="utf-8"))
    rec = {r["original"]: r for r in meta["referenced_files"]}["state.db"]
    assert rec["status"] == "excluded-binary"
    assert rec["bundled"] == ""
    assert rec["bytes"] == len(blob)
    assert rec["sha256"] == hashlib.sha256(blob).hexdigest()
    copied = [p for p in (d / "workspace").iterdir()]
    assert all(p.read_bytes() != blob for p in copied)
    assert not any(b"SQLite format 3" in p.read_bytes() for p in copied)
    assert rec["sha256"] in (d / "report.md").read_text(encoding="utf-8")
