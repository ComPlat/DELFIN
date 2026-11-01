#!/usr/bin/env python3
"""Task setup entry: resolve a count-like merge conflict by summing
(gh1_count_conflict).

Both branches grew the shared ``metrics/counts.py`` module: main added
a feature and bumped its own counter, the partner branch added its
feature and bumped the SAME counter key to a different value. The
agent must merge the partner branch and resolve the conflict by
SUMMING the two sides (count semantics), not by picking one, and the
suite must be green afterwards.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gh1_repo import counts_conflict, main, mark_base, new_repo, write  # noqa: E402
from generalist_repo import commit  # noqa: E402

COUNTS = {
    "conftest.py": (
        'import sys\nfrom pathlib import Path\n\n'
        'sys.path.insert(0, str(Path(__file__).parent / "src"))\n'),
    "metrics/counts.py": (
        '"""Measured counts per task class; summed at merge time."""\n'
        'COUNTS = {"count": 0}\n'),
    "metrics/__init__.py": "",
}

MAIN_FEATURE = {
    "src/widgetkit/rename.py": (
        '"""Rename widgets in place."""\n\n'
        'def rename(widgets, mapping):\n'
        '    """Apply a name mapping, unknown names kept."""\n'
        '    return [mapping.get(w, w) for w in widgets]\n'),
}

PARTNER_FEATURE = {
    "src/widgetkit/lookup.py": (
        '"""Look up widget settings."""\n\n'
        'def load_config(path):\n'
        '    """Read config.json; missing file -> empty dict."""\n'
        '    try:\n'
        '        import json\n'
        '        with open(path, encoding="utf-8") as fh:\n'
        '            return json.load(fh)\n'
        '    except FileNotFoundError:\n'
        '        return {}\n'),
    "tests/test_features.py": (
        'from widgetkit.lookup import load_config\n'
        'from widgetkit.rename import rename\n\n\n'
        'def test_rename_maps_known_names():\n'
        '    assert rename(["a", "b"], {"a": "A"}) == ["A", "b"]\n\n\n'
        'def test_load_config_missing_file_is_empty():\n'
        '    assert load_config("no-such-file.json") == {}\n'),
}


def build(workspace: str) -> Path:
    repo = new_repo(workspace, "widgetkit")
    write(repo, {**COUNTS, **MAIN_FEATURE})
    base = commit(repo, "initial widgetkit package")
    mark_base(repo, base)
    # main: bump count 0 -> 10 in its own commit
    counts = repo / "metrics/counts.py"
    counts.write_text(
        '"""Measured counts per task class; summed at merge time."""\n'
        'COUNTS = {"count": 10}\n', encoding="utf-8")
    commit(repo, "main: rename feature, counts 10")
    # partner branch from base: bump count 0 -> 15 + its feature
    counts_conflict(repo, "metrics/counts.py", 0, 15,
                    "src/widgetkit/lookup.py",
                    PARTNER_FEATURE["src/widgetkit/lookup.py"],
                    partner_files={"tests/test_features.py":
                                   PARTNER_FEATURE["tests/test_features.py"]})
    return repo


if __name__ == "__main__":
    main(build)
