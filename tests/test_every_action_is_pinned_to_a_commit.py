"""Every action a workflow runs is named by a full commit SHA.

Input: every ``uses:`` in .github/workflows/*.yml (step actions and job-level
reusable workflows), read as text so the trailing version comment is seen,
and cross-checked against the parsed YAML so no ``uses:`` is missed by the
line pattern. Also .github/dependabot.yml.

Output: a failure that names the file, the line and the reference for:
- a reference by tag or branch (``actions/checkout@v7``): a tag is a
  mutable pointer, and whoever can move it decides what this repository's
  CI token runs (OpenSSF Scorecard, Pinned-Dependencies);
- a SHA shorter than 40 hex digits: ambiguous, and resolved by GitHub as a
  ref name first;
- a pin without a ``# vX.Y.Z`` comment: Dependabot rewrites SHA pins only
  when the comment names the version, and a reader cannot tell what the
  SHA is without it.
Local actions (``./path``) and ``docker://...@sha256:`` images are allowed.

Rejected alternative: the repository setting ``sha_pinning_required``
alone. It refuses an unpinned workflow at run time, after the change has
merged; this test refuses it in the pull request, and also covers the
version comment the setting does not look at.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = ROOT / ".github" / "workflows"

_USES = re.compile(r"^\s*(?:-\s+)?uses:\s*(?P<ref>\S+)(?P<rest>.*)$")
_PINNED = re.compile(r"^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+(?:/[^@\s]+)?@[0-9a-f]{40}$")
_VERSION_COMMENT = re.compile(r"^\s+#\s*v\d+(?:\.\d+)*\s*$")


def unpinned(text: str) -> list[tuple[int, str]]:
    """(line number, reason) for every ``uses:`` line that is not pinned."""
    found = []
    for number, line in enumerate(text.splitlines(), 1):
        if line.lstrip().startswith("#"):
            continue
        match = _USES.match(line)
        if not match:
            continue
        ref = match.group("ref").strip("'\"")
        rest = match.group("rest")
        if ref.startswith("./"):
            continue
        if ref.startswith("docker://"):
            if "@sha256:" not in ref:
                found.append((number, f"{ref}: image not pinned by digest"))
            continue
        if not _PINNED.match(ref):
            found.append((number, f"{ref}: not a full commit SHA"))
        elif not _VERSION_COMMENT.match(rest):
            found.append((number, f"{ref}: no '# vX.Y.Z' comment"))
    return found


def _parsed_uses(node) -> list[str]:
    if isinstance(node, dict):
        out = [str(v) for k, v in node.items() if k == "uses"]
        for value in node.values():
            out += _parsed_uses(value)
        return out
    if isinstance(node, list):
        return [u for item in node for u in _parsed_uses(item)]
    return []


_FILES = sorted(WORKFLOWS.glob("*.y*ml"))


@pytest.mark.parametrize("path", _FILES, ids=[p.name for p in _FILES])
def test_every_uses_is_a_commit_sha_with_its_version(path):
    text = path.read_text(encoding="utf-8")
    bad = unpinned(text)
    assert not bad, "\n".join(f"{path.name}:{n}: {why}" for n, why in bad)


@pytest.mark.parametrize("path", _FILES, ids=[p.name for p in _FILES])
def test_the_line_pattern_sees_every_uses_the_parser_sees(path):
    """A ``uses:`` written in a form the line pattern does not match (a flow
    mapping, a quoted key) would pass the test above by not being seen."""
    text = path.read_text(encoding="utf-8")
    parsed = sorted(_parsed_uses(yaml.safe_load(text)))
    by_line = sorted(m.group("ref").strip("'\"") for m in map(_USES.match, text.splitlines())
                     if m and not m.group(0).lstrip().startswith("#"))
    assert parsed == by_line, (parsed, by_line)


@pytest.mark.parametrize("line, ok", [
    ("      - uses: actions/checkout@v7", False),
    ("        uses: actions/checkout@main", False),
    ("        uses: actions/checkout@3d3c42e", False),
    ("        uses: actions/checkout@3d3c42e5aac5ba805825da76410c181273ba90b1", False),
    ("        uses: actions/checkout@3d3c42e5aac5ba805825da76410c181273ba90b1 # v7.0.1", True),
    ("        uses: github/codeql-action/init@6f2ac81f5cd81fe351a197ba0b00d91d03cd2aee # v4.39.0", True),
    ("    uses: org/repo/.github/workflows/x.yml@v1", False),
    ("        uses: ./.github/actions/local", True),
    ("        uses: docker://alpine:3", False),
    ("        uses: docker://alpine@sha256:" + "0" * 64, True),
])
def test_the_check_refuses_what_it_must_and_passes_a_pin(line, ok):
    assert (not unpinned(line)) is ok, unpinned(line)


def _updates() -> dict[str, dict]:
    config = yaml.safe_load((ROOT / ".github" / "dependabot.yml").read_text(encoding="utf-8"))
    return {f"{u['package-ecosystem']}:{u.get('directory', '/')}": u for u in config["updates"]}


def test_dependabot_moves_the_pins_and_the_window_dependencies():
    """A pin nobody moves is a known advisory left in place. Dependabot's
    github-actions ecosystem updates SHA pins together with their comment;
    the bun ecosystem reads tui/bun.lock."""
    updates = _updates()
    assert "github-actions:/" in updates, sorted(updates)
    assert "bun:/tui" in updates, sorted(updates)
    assert (ROOT / "tui" / "bun.lock").is_file(), "the bun ecosystem reads bun.lock"
    assert "pip:/" in updates, sorted(updates)
