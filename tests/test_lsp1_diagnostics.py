"""LSP1: diagnostics.check — language-server style findings right after
an edit (red tests written before delfin/agent/diagnostics.py existed).

The fake ruff on PATH is a tiny stand-in linter: it reads the content it
is given (a file argument, or stdin for the old-text pass) and emits a
ruff --output-format json array with one F821 entry per marker name it
finds ("alpha" -> pre-existing marker, "beta" -> newly added marker).
All tests use throwaway tmp_path workspaces; nothing touches a real
home or the live web.
"""
import json
import os
import stat
import textwrap
import time

import pytest

import delfin.agent.diagnostics as diag

FAKE_RUFF = r"""#!/usr/bin/python3
import json, re, sys

args = sys.argv[1:]
if "--stdin-filename" in args:
    content = sys.stdin.read()
else:
    path = args[-1]
    with open(path, encoding="utf-8", errors="replace") as fh:
        content = fh.read()

out = []
for i, line in enumerate(content.splitlines(), start=1):
    for name in ("alpha", "beta"):
        if re.search(rf"\b{name}\(\)", line):
            col = line.index(name)
            out.append({
                "code": "F821",
                "message": f"Undefined name `{name}`",
                "location": {"row": i, "column": col},
            })
print(json.dumps(out))
"""


@pytest.fixture()
def fake_ruff(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    exe = bin_dir / "ruff"
    exe.write_text(FAKE_RUFF, encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
    return bin_dir


def _make_file(tmp_path, text, name="mod.py"):
    p = tmp_path / name
    p.write_text(textwrap.dedent(text), encoding="utf-8")
    return p


def test_f821_reported_at_right_line(tmp_path, fake_ruff):
    p = _make_file(tmp_path, """\
        def go():
            beta()
    """)
    found = diag.check(p, tmp_path)
    f821 = [d for d in found if d.code == "F821"]
    assert len(f821) == 1
    assert f821[0].line == 2
    assert "beta" in f821[0].message


def test_preexisting_finding_not_repeated(tmp_path, fake_ruff):
    # The edit ADDS beta(); alpha() was already in the pre-edit text.
    old = "def go():\n    alpha()\n"
    new = "def go():\n    alpha()\n    beta()\n"
    p = _make_file(tmp_path, new)
    found = diag.check(p, tmp_path, old_text=old)
    names = [d.message for d in found if d.code == "F821"]
    assert names == ["Undefined name `beta`"]


def test_clean_edit_reports_nothing(tmp_path, fake_ruff):
    text = "def go():\n    return 1\n"
    p = _make_file(tmp_path, text)
    assert diag.check(p, tmp_path, old_text=text) == []


# --- degraded and hostile paths (reviewer s15 requests, brief phase 6) ---


def test_no_ruff_on_path_is_one_line_note(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    p = _make_file(tmp_path, "def go():\n    beta()\n")
    found = diag.check(p, tmp_path)
    assert found == []


def test_ruff_is_a_directory_not_a_tool(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    (bin_dir / "ruff").mkdir(parents=True)
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
    p = _make_file(tmp_path, "def go():\n    beta()\n")
    found = diag.check(p, tmp_path)
    assert found == []


def test_slow_tool_cut_at_5s(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    exe = bin_dir / "ruff"
    exe.write_text(
        "#!/usr/bin/python3\n"
        "import time\n"
        "time.sleep(30)\n"
        "print('[]')\n",
        encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
    p = _make_file(tmp_path, "def go():\n    beta()\n")
    start = time.monotonic()
    found = diag.check(p, tmp_path)
    elapsed = time.monotonic() - start
    assert found == []
    assert elapsed < 10.0, f"slow tool not cut: {elapsed:.1f}s"


def test_ansi_and_huge_output_bounded_and_sanitised(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    exe = bin_dir / "ruff"
    # ANSI-escape-laden message plus 10 MB of padding after the JSON.
    exe.write_text(
        "#!/usr/bin/python3\n"
        "import json, sys\n"
        "msg = '\\x1b[31mUndefined name \\x1b[1m`beta`\\x1b[0m'\n"
        "row = json.dumps([{"
        '"code": "F821", "message": msg,'
        ' "location": {"row": 2, "column": 4}}])\n'
        "sys.stdout.write(row + ' ' * (10 * 1024 * 1024))\n",
        encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
    p = _make_file(tmp_path, "def go():\n    beta()\n")
    found = diag.check(p, tmp_path)
    # Bounded: the run must come back at all (capture_output handles 10 MB,
    # but json.loads on trailing padding must not crash the module).
    for d in found:
        assert "\x1b" not in d.message, "ANSI escapes leaked into diagnostics"
    assert found, "the 10 MB padded output must still yield its finding"


def test_fake_json_output_rejected_cleanly(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    exe = bin_dir / "ruff"
    exe.write_text("#!/usr/bin/python3\nprint('this is not json')\n",
                   encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
    p = _make_file(tmp_path, "def go():\n    beta()\n")
    found = diag.check(p, tmp_path)
    assert found == []


def test_file_outside_workspace_never_checked(tmp_path, fake_ruff):
    outside = tmp_path.parent / "outside_mod.py"
    outside.write_text("def go():\n    beta()\n", encoding="utf-8")
    p = _make_file(tmp_path, "def go():\n    beta()\n")
    found = diag.check(p, tmp_path)  # sanity: inside is checked
    assert isinstance(found, list)
    assert diag.check(outside, tmp_path) == []


def test_non_python_file_returns_empty_list(tmp_path, fake_ruff):
    p = tmp_path / "notes.md"
    p.write_text("# beta()\n", encoding="utf-8")
    assert diag.check(p, tmp_path) == []


def test_cap_at_20_diagnostics(tmp_path, fake_ruff):
    lines = "\n".join(f"    beta{i}()" for i in range(30))
    p = _make_file(tmp_path, "def go():\n" + lines + "\n")
    found = diag.check(p, tmp_path)
    assert len(found) <= diag.MAX_DIAGNOSTICS
