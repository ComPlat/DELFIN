"""Every specialist preset the pack ships passes ``validate_preset``.

This is the CI half of the pack/user split. A specialist is one markdown
file with YAML frontmatter, discovered from ``~/.delfin/subagents/`` as
well as from the packaged ``delfin/agent/pack/agents/``. The two halves
are refused in different places on purpose:

* a USER preset that fails the check is refused at load, with its reason,
  because it is the user's file and the user can fix it;
* a PACK preset that fails must fail HERE, in CI, because at runtime it
  would simply not appear. A specialist the product ships and then drops
  silently is worse than one that never existed: the feature is
  documented and absent, and nothing says why.

Input: the ``*_subagent.md`` files found in the pack's preset directory.
Output: for each one, ``subagents.validate_preset``'s list of problems --
required fields, a unique kebab-case id, tools from the known catalogue,
``write_globs`` inside the repo and never a protected path,
``may_delegate_to`` only known presets, a valid ``output_schema``, no
override of a built-in or pack preset. Empty means valid.

Semantics pinned here: the set of files is DISCOVERED, so a preset added
later is checked without touching this file, and every problem is named
in the failure message so a broken preset is actionable rather than
merely red.

The directory comes from ``subagents._md_preset_search_dirs()`` -- the
product's own list -- and not from the repository root. The reason is
that the product reads presets relative to the imported
``delfin.agent`` package, which on a CI host need not be the checkout
the tests were collected from; deriving the path from the repo root
would check files the running product never opens. The pack entry is
picked by comparing against the imported package's own
``pack/agents``, so the test does not depend on the order of that list.

The rejected alternative: a hard-coded list of preset names, or a
hand-written parametrisation. Both pass while checking nothing once a
file is renamed or added, which is the defect this file exists to
prevent -- hence also the explicit non-empty assertion below.

First run of this file, against commit 93c498537: 13 of the 13 live
pack presets were refused with ``name <id> overrides an existing pack
preset <id>`` -- each one reported as overriding ITSELF. The uniqueness
block exempted the definition being introduced by discarding its own
stem from ``pack_ids``, and then re-added every live registry id, the
preset's own among them, so the exemption was undone for exactly the
presets that were loaded. Two variants of the same block are worth
knowing when a failure here is read: the exemption keys on the file
STEM while the collision keys on the frontmatter NAME, so a preset
whose id differs from its filename (``smiles-builder_subagent.md`` is
named ``metal-smiles-builder``) is not covered by it; and the exemption
is skipped outright when a definition arrives as TEXT rather than as a
path. A failure naming only the preset's own id is this shape -- the
check, not the preset.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import subagents as S
from delfin.agent.skills import _parse_frontmatter

# The glob the product uses in ``_load_md_presets`` / ``_build_preset_registry``.
_PRESET_GLOB = "*_subagent.md"


def _pack_preset_dir() -> Path:
    """The pack directory out of the product's own search list."""
    want = (Path(S.__file__).resolve().parent / "pack" / "agents").resolve()
    for d in S._md_preset_search_dirs():
        try:
            if Path(d).resolve() == want:
                return Path(d)
        except OSError:  # pragma: no cover - unreadable path
            continue
    return want


def _preset_id(path: Path) -> str:
    """The id the product gives this file.

    Mirrors ``_build_preset_registry``: the frontmatter ``name``, falling
    back to the file stem without the ``_subagent`` suffix. Read with the
    product's own frontmatter parser, so the id is taken from the file's
    content rather than from its name.
    """
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:  # pragma: no cover - unreadable file
        text = ""
    meta, _body = _parse_frontmatter(text)
    return (meta.get("name") or path.stem.replace("_subagent", "")).strip()


def _is_skeleton(preset_id: str) -> bool:
    """The product's underscore rule, in one place.

    ``_build_preset_registry`` skips a preset whose id begins with an
    underscore: that marks a copy-and-fill skeleton
    (``_template_subagent.md``), not a live type. Such a file is not a
    shipped specialist, so it is not checked as one -- and it could not
    pass, because an underscore-prefixed id is not kebab-case. The rule
    is applied here rather than a filename being excluded by hand, so a
    second skeleton added later is excluded for the same reason.
    """
    return preset_id.startswith("_")


PACK_DIR = _pack_preset_dir()
DISCOVERED = sorted(PACK_DIR.glob(_PRESET_GLOB)) if PACK_DIR.is_dir() else []
SHIPPED = [p for p in DISCOVERED if not _is_skeleton(_preset_id(p))]
SKELETONS = [p for p in DISCOVERED if _is_skeleton(_preset_id(p))]


def test_the_discovery_itself_is_not_empty():
    """A glob that matches nothing makes every check below vacuous."""
    assert PACK_DIR.is_dir(), (
        f"the pack preset directory does not exist: {PACK_DIR} -- the "
        "product reads presets from here, so nothing was checked")
    assert DISCOVERED, (
        f"no {_PRESET_GLOB} file found under {PACK_DIR}: this file would "
        "pass while checking no preset at all")
    assert SHIPPED, (
        f"every {_PRESET_GLOB} file under {PACK_DIR} was excluded as a "
        "skeleton; no shipped preset was checked")


def test_a_skeleton_is_excluded_by_the_products_underscore_rule():
    """Excluded because the product excludes it, not because of its name.

    Each file left out above is left out because its id begins with an
    underscore, and the product's registry agrees: the id is absent from
    ``_build_preset_registry()``. Each file kept has an id the product
    would accept as a live type.
    """
    registry = S._build_preset_registry()
    assert not [k for k in registry if k.startswith("_")], (
        "the registry offers an underscore-prefixed id as a delegateable "
        f"type: {sorted(k for k in registry if k.startswith('_'))}")
    for path in SKELETONS:
        preset_id = _preset_id(path)
        assert _is_skeleton(preset_id)
        assert preset_id not in registry, (
            f"{path.name} has the skeleton id {preset_id!r} and the "
            "product still loaded it as a live preset")
    for path in SHIPPED:
        assert not _is_skeleton(_preset_id(path))


@pytest.mark.parametrize(
    "path", SHIPPED, ids=[p.name for p in SHIPPED])
def test_every_shipped_preset_is_valid(path: Path):
    """One shipped preset, one check; the problems are the message."""
    problems = S.validate_preset(str(path))
    assert problems == [], (
        f"shipped preset {path.name} fails validate_preset "
        f"({len(problems)} problem(s)):\n  - "
        + "\n  - ".join(problems)
        + f"\n  file: {path}")
