#!/usr/bin/env python3
"""Task setup entry: cross-review another branch on a scratch branch
(gh1_cross_review).

The fixture is a repo with the agent's own small ``main`` (with a
local commit) and a ``partner/feature`` branch that adds a module with
a real defect the review should catch (a bare ``except: pass`` that
swallows a real failure) and a module that is fine. The agent must
read the partner's diff without merging it into its branch, and write
the review to answer.txt (the graded artifact).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gh1_repo import feature_branch, main, mark_base, new_repo, write  # noqa: E402
from generalist_repo import commit  # noqa: E402

OWN = {
    "conftest.py": (
        'import sys\nfrom pathlib import Path\n\n'
        'sys.path.insert(0, str(Path(__file__).parent / "src"))\n'),
    "src/widgetkit/__init__.py": 'from .core import normalise\n\n'
                                 '__all__ = ["normalise"]\n',
    "src/widgetkit/core.py": (
        '"""Normalise widget labels."""\n'
        'import re\n\n'
        '_SEP = re.compile(r"[-_\\s]+")\n\n'
        'def normalise(label: str) -> str:\n'
        '    """A label becomes lower-case, separators single spaces."""\n'
        '    return _SEP.sub(" ", label.strip().lower())\n'),
}

PARTNER = {
    "src/widgetkit/lookup.py": (
        '"""Look up widget settings."""\n'
        'import json\n'
        'from pathlib import Path\n\n'
        'def load_config(path):\n'
        '    """Read config.json; missing file -> defaults."""\n'
        '    try:\n'
        '        with Path(path).open(encoding="utf-8") as fh:\n'
        '            return json.load(fh)\n'
        '    except:  # noqa: E722\n'
        '        pass\n'
        '    return {}\n'),
    "src/widgetkit/pretty.py": (
        '"""Format a widget table."""\n\n'
        'def rows(items):\n'
        '    """Two columns, aligned."""\n'
        '    width = max((len(k) for k in items), default=0)\n'
        '    for k, v in items:\n'
        '        yield f"{k.ljust(width)}  {v}"\n'),
}


def build(workspace: str) -> Path:
    repo = new_repo(workspace, "widgetkit")
    write(repo, OWN)
    base = commit(repo, "initial widgetkit package")
    mark_base(repo, base)
    # The partner's work sits on partner/feature; own branch is main.
    feature_branch(repo, "partner/feature", PARTNER,
                   "partner: config loader and table formatter")
    from generalist_repo import git
    git(repo, "checkout", "-q", "main")
    return repo


if __name__ == "__main__":
    main(build)
