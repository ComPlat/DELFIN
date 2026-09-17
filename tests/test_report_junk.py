"""Tests for delfin/agent/report_junk.py against a fabricated tree."""

import os
from pathlib import Path

import pytest

from delfin.agent.report_junk import collect, format_text, main


@pytest.fixture()
def fake_workspace(tmp_path):
    """Fabricate a workspace with junk git status would (mostly) hide."""
    # normal tracked-looking files
    (tmp_path / "src" / "code.py").parent.mkdir(parents=True)
    (tmp_path / "src" / "code.py").write_text("print('hi')")
    (tmp_path / "README.md").write_text("hello")

    # __pycache__ + bytecode
    pc = tmp_path / "src" / "__pycache__"
    pc.mkdir()
    (pc / "code.cpython-311.pyc").write_bytes(b"\x00" * 10)
    # stray bytecode OUTSIDE __pycache__ (the dir itself is not descended)
    (tmp_path / "stray.pyc").write_bytes(b"\x00" * 4)

    # QM leftovers: files
    qm = tmp_path / "run1"
    qm.mkdir()
    (qm / "job.out").write_text("ORCA 5")
    (qm / "job.engrad").write_text("gradient")
    (qm / "NORMAL_TERMINATION").write_text("")
    (qm / "xtbrestart").write_text("restart")
    # QM scratch dir
    (tmp_path / "run1" / "orca.tmp").mkdir()
    (tmp_path / "run1" / "orca.tmp" / "part.prop").write_text("x")

    # lock file
    (tmp_path / "run1" / "lease.lock").write_text("")

    # large file (threshold is 100 MiB, so patch the module-level constant
    # is avoided by writing a file just over it)
    big = tmp_path / "run1" / "big.trj"
    big.write_bytes(b"\x00" * (100 * 1024 * 1024 + 1))

    # secret-looking files: must NOT appear anywhere in the report
    (tmp_path / "run1" / ".env").write_text("SECRET")
    (tmp_path / "run1" / "my.key").write_text("SECRET")

    return tmp_path


def test_collect_finds_all_junk_categories(fake_workspace):
    data = collect(fake_workspace)
    rels = lambda lst: [p if isinstance(p, str) else p["path"]
                        for p in lst]
    assert "src/__pycache__" in data["pycache_dirs"]
    assert "stray.pyc" in data["pyc_files"]
    assert "run1/job.out" in data["qm_leftovers"]
    assert "run1/job.engrad" in data["qm_leftovers"]
    assert "run1/NORMAL_TERMINATION" in data["qm_leftovers"]
    assert "run1/xtbrestart" in data["qm_leftovers"]
    assert "run1/orca.tmp/" in data["qm_leftovers"]
    assert "run1/lease.lock" in data["lock_files"]
    assert "run1/big.trj" in rels(data["large_files"])


def test_collect_skips_secrets(fake_workspace):
    data = collect(fake_workspace)
    blob = repr(data)
    assert ".env" not in blob
    assert "my.key" not in blob


def test_collect_folder_sizes_and_total(fake_workspace):
    data = collect(fake_workspace)
    assert set(data["folder_sizes"]) >= {"src", "run1"}
    assert data["total_bytes"] >= 100 * 1024 * 1024
    # README.md at top level counts under its own name
    assert data["folder_sizes"]["README.md"] == 5


def test_collect_empty_dir(tmp_path):
    data = collect(tmp_path)
    assert data["large_files"] == []
    assert data["qm_leftovers"] == []
    assert data["total_bytes"] == 0


def test_collect_default_root_is_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    data = collect()
    assert Path(data["root"]) == tmp_path.resolve()


def test_collect_is_read_only(fake_workspace):
    before = sorted(
        (str(p), p.stat().st_size) for p in fake_workspace.rglob("*")
    )
    collect(fake_workspace)
    after = sorted(
        (str(p), p.stat().st_size) for p in fake_workspace.rglob("*")
    )
    assert before == after


def test_format_text_worst_first(fake_workspace):
    data = collect(fake_workspace)
    text = format_text(data)
    assert "Junk report for" in text
    assert "QM leftovers" in text
    assert "big.trj" in text
    assert "100.0 MiB" in text
    # large-file section is sorted worst-first: big.trj before anything else
    large = [ln for ln in text.splitlines()
             if ln.startswith("  ") and "MiB" in ln]
    assert large and "big.trj" in large[0]
    assert "nothing was modified or deleted" in text.lower()
    # secrets never rendered
    assert ".env" not in text
    assert "my.key" not in text


def test_format_text_empty(tmp_path):
    text = format_text(collect(tmp_path))
    assert "Nothing suspicious found." in text


def test_main_cli(capsys, fake_workspace):
    rc = main([str(fake_workspace)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "run1/job.out" in out
    assert "Size per top-level entry" in out


def test_module_runnable_via_python_m(fake_workspace):
    import subprocess, sys
    repo_root = os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))
    env = dict(os.environ, PYTHONPATH=repo_root)
    proc = subprocess.run(
        [sys.executable, "-m", "delfin.agent.report_junk", str(fake_workspace)],
        capture_output=True, text=True, cwd=repo_root, env=env,
    )
    assert proc.returncode == 0, proc.stderr
    assert "QM leftovers" in proc.stdout


def test_symlinks_not_followed(fake_workspace):
    target = fake_workspace / "src"
    (fake_workspace / "link_to_src").symlink_to(target, target_is_directory=True)
    data = collect(fake_workspace)
    # only one __pycache__ counted (via src), link itself ignored
    assert data["pycache_dirs"].count("src/__pycache__") == 1
    assert "link_to_src" not in str(data["pyc_files"])
