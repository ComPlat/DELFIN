"""Tests for the ORCA keyword extractor (anti-hallucination ground-truth).

The extractor walks the indexed ORCA manual and produces a structured
keyword namespace per computation block.  Tests verify:

- The four canonical CASSCF keywords (nel, norb, mult, nroots) are
  identified — proves the extractor catches the right things.
- Common hallucination patterns (Nactel, Nactorb, Multiplicity) are
  NOT in the extracted namespace — proves the manual is the actual
  source of truth.
- ``is_real_keyword`` API is case-insensitive and rejects unknown
  blocks gracefully.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.orca_keyword_extractor import (
    extract_keyword_namespace,
    is_real_keyword,
)


# ---------------------------------------------------------------------------
# Live extraction (uses the actual ~/.delfin/doc_index.json on the box).
# These tests are integration tests in spirit — they rely on a real
# index existing on the developer machine.  Skip if absent.
# ---------------------------------------------------------------------------


_INDEX = Path.home() / ".delfin" / "doc_index.json"
_has_index = _INDEX.exists()


@pytest.fixture(scope="module")
def namespace():
    if not _has_index:
        pytest.skip(f"No doc_index at {_INDEX}")
    return extract_keyword_namespace()


@pytest.mark.skipif(not _has_index, reason="No indexed ORCA manual")
def test_extractor_finds_casscf_block(namespace):
    assert "casscf" in namespace
    info = namespace["casscf"]
    assert info["block"] == "%casscf"
    assert info["n_chapter_refs"] > 0
    assert info["n_keywords"] > 10  # ~200 in practice


@pytest.mark.skipif(not _has_index, reason="No indexed ORCA manual")
def test_extractor_includes_canonical_casscf_keywords(namespace):
    """The four keywords every CASSCF input needs MUST be extracted —
    if these go missing, our anti-hallucination ground-truth is broken."""
    kws = {k.lower() for k in namespace["casscf"]["keywords"]}
    for required in ("nel", "norb", "mult", "nroots"):
        assert required in kws, (
            f"Canonical keyword '{required}' missing from extracted "
            f"CASSCF namespace — extraction is broken"
        )


@pytest.mark.skipif(not _has_index, reason="No indexed ORCA manual")
def test_extractor_excludes_known_hallucinations(namespace):
    """Common hallucination patterns observed in real production
    sessions MUST NOT appear in the extracted namespace.  If any of
    these match, the model is being given a false ground-truth."""
    kws = {k.lower() for k in namespace["casscf"]["keywords"]}
    for hallu in ("nactel", "nactorb", "multiplicity"):
        assert hallu not in kws, (
            f"Hallucinated keyword '{hallu}' leaked into CASSCF "
            f"namespace — extraction is too permissive"
        )


@pytest.mark.skipif(not _has_index, reason="No indexed ORCA manual")
def test_extractor_covers_major_method_families(namespace):
    """All these should have at least a few sections in the manual."""
    for block in ("casscf", "mp2", "scf", "method", "freq", "tddft"):
        assert block in namespace, f"Missing block: {block}"
        assert namespace[block]["n_chapter_refs"] > 0


# ---------------------------------------------------------------------------
# is_real_keyword API
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _has_index, reason="No indexed ORCA manual")
def test_is_real_keyword_case_insensitive(namespace):
    """API must be case-insensitive — agents emit casing variants."""
    assert is_real_keyword("casscf", "nel", namespace=namespace) is True
    assert is_real_keyword("casscf", "NEL", namespace=namespace) is True
    assert is_real_keyword("casscf", "Nel", namespace=namespace) is True
    assert is_real_keyword("CASSCF", "nel", namespace=namespace) is True


@pytest.mark.skipif(not _has_index, reason="No indexed ORCA manual")
def test_is_real_keyword_rejects_hallucinations(namespace):
    """The pre-set hallucination synonyms must be rejected."""
    assert is_real_keyword("casscf", "Nactel", namespace=namespace) is False
    assert is_real_keyword("casscf", "Nactorb", namespace=namespace) is False


def test_is_real_keyword_unknown_block_returns_false():
    """Querying a block that doesn't exist in the manual must not
    crash — just return False."""
    ns = {"casscf": {"keywords": ["nel"]}}
    assert is_real_keyword("nonexistent", "anykey", namespace=ns) is False


def test_is_real_keyword_works_with_passed_namespace():
    """Caller can supply a pre-extracted namespace for speed."""
    ns = {
        "casscf": {"keywords": ["nel", "norb"]},
        "tddft":  {"keywords": ["nroots"]},
    }
    assert is_real_keyword("casscf", "nel", namespace=ns) is True
    assert is_real_keyword("casscf", "nroots", namespace=ns) is False
    assert is_real_keyword("tddft",  "nroots", namespace=ns) is True


# ---------------------------------------------------------------------------
# Committed snapshot (in pack/benchmark/) — must stay in sync
# ---------------------------------------------------------------------------


def test_committed_snapshot_exists():
    """The repo-committed snapshot JSON is the artefact downstream
    consumers (benchmark, dashboard hints) actually load.  Verify it
    exists + is well-formed."""
    p = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
         / "pack" / "benchmark" / "keywords_groundtruth_orca.json")
    assert p.exists(), "ORCA keyword snapshot missing"
    data = json.loads(p.read_text(encoding="utf-8"))
    assert data.get("version") or data.get("orca_version")
    assert "blocks" in data
    assert "casscf" in data["blocks"]
    casscf = data["blocks"]["casscf"]
    assert "nel" in {k.lower() for k in casscf["keywords"]}
    assert "norb" in {k.lower() for k in casscf["keywords"]}
