#!/usr/bin/env python3
"""Setup for generalist_merge_conflict: a merge stopped on a conflict.

``configlib/`` -- base commit B: ``load(path)`` reads ``key=value`` lines.
main (M, on B): comment and blank lines are skipped; adds a test.
feature/env-override (F, on B): ``load(path, env=None)`` lets an
environment variable ``CONFIGLIB_<KEY>`` override a value; adds a test.
Both rewrite the same loop. The setup runs ``git merge
feature/env-override`` on main and leaves it stopped on the conflict in
configlib/loader.py (tests/ merges cleanly).

M and F are recorded in .git/delfin-bench-base and
.git/delfin-bench-feature for the acceptance script.
Preconditions: the merge is in progress (MERGE_HEAD exists) and
loader.py is the one conflicted path.
argv[1] is the workspace.
"""
from __future__ import annotations

from generalist_repo import (GITIGNORE_PY, commit, expect, git, main,
                             mark_base, new_repo, write)

LOADER_BASE = '''"""Read key=value configuration files."""


def load(path):
    """The settings in ``path`` as a dict of strings."""
    settings = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            key, _, value = line.strip().partition("=")
            settings[key.strip()] = value.strip()
    return settings
'''

LOADER_MAIN = '''"""Read key=value configuration files."""


def load(path):
    """The settings in ``path`` as a dict of strings.

    Blank lines and lines starting with ``#`` are skipped.
    """
    settings = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            key, _, value = line.partition("=")
            settings[key.strip()] = value.strip()
    return settings
'''

LOADER_FEATURE = '''"""Read key=value configuration files."""

import os


def load(path, env=None):
    """The settings in ``path`` as a dict of strings.

    A variable ``CONFIGLIB_<KEY>`` in ``env`` (default: os.environ)
    overrides the value of ``key``.
    """
    env = os.environ if env is None else env
    settings = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            key, _, value = line.strip().partition("=")
            key = key.strip()
            settings[key] = env.get("CONFIGLIB_" + key.upper(), value.strip())
    return settings
'''

TEST_BASE = '''from configlib.loader import load


def test_reads_pairs(tmp_path):
    p = tmp_path / "a.conf"
    p.write_text("host = example.org\\nport=8080\\n", encoding="utf-8")
    assert load(p) == {"host": "example.org", "port": "8080"}
'''

TEST_MAIN = '''from configlib.loader import load


def test_skips_comments_and_blank_lines(tmp_path):
    p = tmp_path / "a.conf"
    p.write_text("# settings\\n\\nhost=example.org\\n", encoding="utf-8")
    assert load(p) == {"host": "example.org"}
'''

TEST_FEATURE = '''from configlib.loader import load


def test_environment_overrides_a_value(tmp_path):
    p = tmp_path / "a.conf"
    p.write_text("port=8080\\n", encoding="utf-8")
    assert load(p, env={"CONFIGLIB_PORT": "9090"}) == {"port": "9090"}
'''


def build(workspace: str):
    repo = new_repo(workspace, "configlib")
    write(repo, {
        ".gitignore": GITIGNORE_PY,
        "README.md": "# configlib\n\nRead key=value configuration files.\n",
        "pyproject.toml": ("[project]\nname = \"configlib\"\n"
                           "version = \"0.2.0\"\n\n"
                           "[tool.pytest.ini_options]\n"
                           "pythonpath = [\".\"]\ntestpaths = [\"tests\"]\n"),
        "configlib/__init__.py": "\"\"\"configlib.\"\"\"\n",
        "configlib/loader.py": LOADER_BASE,
        "tests/test_loader.py": TEST_BASE,
    })
    commit(repo, "Add the key=value loader")
    git(repo, "checkout", "-q", "-b", "feature/env-override")
    write(repo, {"configlib/loader.py": LOADER_FEATURE,
                 "tests/test_env.py": TEST_FEATURE})
    feature = commit(repo, "Let CONFIGLIB_<KEY> override a value")
    git(repo, "checkout", "-q", "main")
    write(repo, {"configlib/loader.py": LOADER_MAIN,
                 "tests/test_comments.py": TEST_MAIN})
    mark_base(repo, commit(repo, "Skip comments and blank lines"))
    (repo / ".git" / "delfin-bench-feature").write_text(feature + "\n",
                                                        encoding="utf-8")
    import subprocess
    proc = subprocess.run(["git", "merge", "--no-edit",
                           "feature/env-override"], cwd=str(repo),
                          capture_output=True, text=True)
    expect("the merge stops on a conflict", proc.returncode != 0,
           proc.stdout + proc.stderr)
    expect("MERGE_HEAD exists", (repo / ".git" / "MERGE_HEAD").is_file())
    conflicted = git(repo, "diff", "--name-only", "--diff-filter=U")
    expect("loader.py is the one conflicted path",
           conflicted.split() == ["configlib/loader.py"], conflicted)
    return repo


if __name__ == "__main__":
    main(build)
