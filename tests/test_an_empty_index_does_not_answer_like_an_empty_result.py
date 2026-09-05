""""You have no calculations" and "I scanned nothing" were one sentence.

``build_calc_index`` skips a root that does not exist, so a home without a
``calc`` directory produced a SUCCESS with zero records — which made the
explicit "Calc index could not be built" error unreachable, ``search_calcs``
return ``[]`` and ``calc_summary`` return all zeros. Byte-identical to a
populated home whose query simply missed. Two more silences sat behind it:
the engine is cached for the life of the process, so a summary claiming to
cover ALL indexed calculations is a snapshot from the first call; and
``calc`` is scanned one level deep while the archives get three, so a
calculation in ``calc/group/name/`` is outside every count.

Every answer now carries where it came from.

The doc index had the matching pair. A PDF that yields no text was chunked
into ONE section with an empty body, so the document reported
``section_count: 1``, the "is the ORCA manual indexed" check passed on a
title match alone, and every search against it returned nothing forever —
which is what happens with no pypdf installed and for any scanned or
DRM-protected manual. And ``built_at`` was written by both indexers and
compared to nothing, so an answer from a section deleted from the source
last month was indistinguishable from a fresh hit.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.doc_server import indexer as I
from delfin.doc_server.calc_indexer import build_calc_index
from delfin.doc_server.calc_search import CalcSearchEngine
from delfin.doc_server.search import DocSearchEngine


# ---------------------------------------------------------------------------
# 1. Calc scan provenance
# ---------------------------------------------------------------------------

def test_an_index_over_nothing_says_it_scanned_nothing():
    idx = build_calc_index(calc_dir=None, archive_dir=None, quiet=True)
    prov = idx["provenance"]
    assert prov["any_root_scanned"] is False
    assert prov["roots_scanned"] == 0

    eng = CalcSearchEngine(idx)
    answer = eng.search(query="PBE0")
    assert answer["results"] == []
    assert answer["index"]["any_root_scanned"] is False
    assert "empty index, not an empty result" in answer["index"]["note"]


def test_a_missing_root_is_named_as_missing(tmp_path):
    idx = build_calc_index(calc_dir=tmp_path / "nope",
                           archive_dir=tmp_path / "also_nope", quiet=True)
    roots = {r["source"]: r for r in idx["provenance"]["roots"]}
    assert roots["calc"]["status"] == "missing"
    assert roots["archive"]["status"] == "missing"
    assert roots["remote_archive"]["status"] == "not_configured"


def _calc(root: Path, name: str) -> Path:
    d = root / name
    d.mkdir(parents=True)
    (d / "CONTROL.txt").write_text("NAME = x\nmethod = classic\n",
                                   encoding="utf-8")
    return d


def test_a_populated_index_says_which_roots_it_read(tmp_path):
    calc = tmp_path / "calc"
    _calc(calc, "flat_one")
    idx = build_calc_index(calc_dir=calc, quiet=True)
    prov = idx["provenance"]
    assert prov["any_root_scanned"] is True
    roots = {r["source"]: r for r in prov["roots"]}
    assert roots["calc"]["status"] == "scanned"
    assert roots["calc"]["records"] == 1
    assert roots["calc"]["path"] == str(calc)
    assert "note" not in CalcSearchEngine(idx).provenance()


def test_an_empty_result_over_a_real_root_is_not_an_empty_index(tmp_path):
    calc = tmp_path / "calc"
    _calc(calc, "flat_one")
    answer = CalcSearchEngine(build_calc_index(calc_dir=calc, quiet=True)
                              ).search(query="definitely_not_present")
    assert answer["results"] == []
    assert answer["index"]["any_root_scanned"] is True
    assert answer["index"]["record_count"] == 1


def test_a_calculation_below_the_depth_limit_is_reported_as_unscanned(tmp_path):
    """calc/ is scanned one level deep against the archives' three, so a
    calculation nested deeper is outside every count. It stays outside —
    but the count now says where it stopped looking."""
    calc = tmp_path / "calc"
    _calc(calc, "flat_one")
    _calc(calc, "group/sub/too_deep")
    idx = build_calc_index(calc_dir=calc, quiet=True)
    found = {r["calc_id"] for r in idx["records"]}
    assert "flat_one" in found
    assert "too_deep" not in found
    roots = {r["source"]: r for r in idx["provenance"]["roots"]}
    assert roots["calc"]["max_depth"] == 1
    assert idx["provenance"]["truncated"], (
        "a calculation was left out of the scan without a word about it")
    assert any("group/sub" in t for t in idx["provenance"]["truncated"])


def test_a_fully_scanned_tree_reports_no_truncation(tmp_path):
    """A calculation folder that simply has no subfolders is not a place
    the scan gave up on."""
    calc = tmp_path / "calc"
    d = _calc(calc, "flat_one")
    (d / "scratch").mkdir()      # a calc's own subfolder, not a grouping one
    idx = build_calc_index(calc_dir=calc, quiet=True)
    assert {r["calc_id"] for r in idx["records"]} == {"flat_one"}
    assert idx["provenance"]["truncated"] == []


def test_the_summary_carries_the_same_provenance(tmp_path):
    calc = tmp_path / "calc"
    _calc(calc, "flat_one")
    summary = CalcSearchEngine(build_calc_index(calc_dir=calc, quiet=True)
                               ).summary()
    assert summary["total_calculations"] == 1
    assert summary["index"]["built_at"]
    assert summary["index"]["any_root_scanned"] is True


def test_an_empty_summary_says_it_covers_nothing():
    summary = CalcSearchEngine(
        build_calc_index(calc_dir=None, quiet=True)).summary()
    assert summary["total_calculations"] == 0
    assert summary["index"]["any_root_scanned"] is False
    assert summary["index"]["note"]


# ---------------------------------------------------------------------------
# 2. A file that yields no text is a failure, not a document
# ---------------------------------------------------------------------------

def test_no_extracted_text_is_not_a_one_section_document():
    assert I._chunk_pdf_into_sections([]) == []


def test_a_pdf_that_yields_nothing_is_not_indexed(tmp_path, monkeypatch):
    lit = tmp_path / "literature"
    lit.mkdir()
    manual = lit / "ORCA_6_1_1.pdf"
    manual.write_bytes(b"%PDF-1.4 not really a pdf")

    def _nothing(path, quiet=False, report=None):
        if report is not None:
            report["reason"] = "pypdf is not installed"
        return []

    monkeypatch.setattr(I, "_extract_pdf_text", _nothing)
    idx = I.build_index(lit, quiet=True)

    assert idx["documents"] == {}, (
        "a PDF with no extractable text was indexed as present-and-empty")
    assert idx["document_count"] == 0
    failed = idx["failed_documents"]
    assert len(failed) == 1
    assert failed[0]["source_path"] == str(manual)
    assert "pypdf" in failed[0]["reason"]


def test_the_manual_check_no_longer_passes_on_a_title(tmp_path, monkeypatch):
    """The "is the ORCA manual indexed" check matches a title against
    ``documents``. With the empty entry gone it answers honestly."""
    lit = tmp_path / "literature"
    lit.mkdir()
    (lit / "ORCA_6_1_1.pdf").write_bytes(b"%PDF-1.4")
    monkeypatch.setattr(I, "_extract_pdf_text",
                        lambda p, quiet=False, report=None: [])
    idx = I.build_index(lit, quiet=True)
    titles = [d.get("title", "") for d in idx["documents"].values()]
    assert not any("orca" in t.lower() for t in titles)
    assert any("orca" in f["source_path"].lower()
               for f in idx["failed_documents"])


def test_a_real_document_is_still_indexed(tmp_path):
    lit = tmp_path / "literature"
    lit.mkdir()
    (lit / "notes.md").write_text("# Heading\nreal content here\n",
                                  encoding="utf-8")
    idx = I.build_index(lit, quiet=True)
    assert idx["document_count"] == 1
    assert idx["failed_documents"] == []
    doc = next(iter(idx["documents"].values()))
    assert doc["total_chars"] > 0


def test_a_pdf_with_text_is_still_indexed(tmp_path, monkeypatch):
    lit = tmp_path / "literature"
    lit.mkdir()
    (lit / "paper.pdf").write_bytes(b"%PDF-1.4")
    monkeypatch.setattr(
        I, "_extract_pdf_text",
        lambda p, quiet=False, report=None: [{"page": 1, "text": "RIJCOSX "
                                              "approximation details"}])
    idx = I.build_index(lit, quiet=True)
    assert idx["document_count"] == 1
    assert idx["failed_documents"] == []


# ---------------------------------------------------------------------------
# 3. A one-section index must not break search outright
# ---------------------------------------------------------------------------

def _index_with(sections: dict) -> dict:
    return {
        "version": 1, "built_at": "2026-01-01T00:00:00+00:00",
        "document_count": 1,
        "documents": {"d": {"title": "T", "source_path": "",
                            "sections": sections}},
        "failed_documents": [],
    }


def test_a_single_section_corpus_does_not_raise():
    """max_df=0.95 of one document floors to zero, below min_df=1, and the
    vectoriser raised ValueError — so the FIRST PDF a user indexed broke
    search_docs on every call."""
    eng = DocSearchEngine(_index_with(
        {"s": {"title": "Only", "text": "hello world def2-TZVP"}}))
    answer = eng.search("def2-TZVP")
    assert [r["section_id"] for r in answer["results"]] == ["s"]
    assert answer["index"]["corpus_sections"] == 1


def test_a_two_section_corpus_still_ranks():
    eng = DocSearchEngine(_index_with({
        "s1": {"title": "A", "text": "hello world"},
        "s2": {"title": "B", "text": "goodbye moon"},
    }))
    assert [r["section_id"]
            for r in eng.search("hello")["results"]] == ["s1"]


def test_an_empty_index_answers_with_its_own_emptiness():
    answer = DocSearchEngine(_index_with({})).search("anything")
    assert answer["results"] == []
    assert answer["index"]["document_count"] == 1
    assert answer["index"]["corpus_sections"] == 0


# ---------------------------------------------------------------------------
# 4. Staleness is representable
# ---------------------------------------------------------------------------

def test_a_changed_source_is_reported_as_stale(tmp_path):
    lit = tmp_path / "literature"
    lit.mkdir()
    src = lit / "notes.md"
    src.write_text("# Heading\noriginal content\n", encoding="utf-8")
    idx = I.build_index(lit, quiet=True)
    assert I.stale_documents(idx) == []

    import os
    src.write_text("# Heading\nrewritten content\n", encoding="utf-8")
    future = os.path.getmtime(src) + 600
    os.utime(src, (future, future))

    stale = I.stale_documents(idx)
    assert [s["doc_id"] for s in stale] == ["notes"]
    assert "changed" in stale[0]["reason"]


def test_a_deleted_source_is_reported_as_stale(tmp_path):
    lit = tmp_path / "literature"
    lit.mkdir()
    src = lit / "notes.md"
    src.write_text("# Heading\ncontent\n", encoding="utf-8")
    idx = I.build_index(lit, quiet=True)
    src.unlink()
    stale = I.stale_documents(idx)
    assert [s["doc_id"] for s in stale] == ["notes"]
    assert "no longer exists" in stale[0]["reason"]


def test_every_search_payload_carries_the_build_time(tmp_path):
    lit = tmp_path / "literature"
    lit.mkdir()
    (lit / "notes.md").write_text("# Heading\nRIJCOSX approximation\n",
                                  encoding="utf-8")
    idx = I.build_index(lit, quiet=True)
    answer = DocSearchEngine(idx).search("RIJCOSX")
    assert answer["results"]
    assert answer["index"]["built_at"] == idx["built_at"]
    assert answer["index"]["stale_documents"] == []

    # ...and a miss carries it too, so "no hit" can be told from "nothing
    # to hit" and from "built before the source changed".
    miss = DocSearchEngine(idx).search("nothing_matches_this_token")
    assert miss["results"] == []
    assert miss["index"]["built_at"] == idx["built_at"]


def test_the_doctor_warns_about_drift(tmp_path, monkeypatch):
    import os
    from delfin.agent import doctor

    lit = tmp_path / "literature"
    lit.mkdir()
    src = lit / "notes.md"
    src.write_text("# Heading\ncontent\n", encoding="utf-8")
    idx = I.build_index(lit, quiet=True)
    idx_path = tmp_path / "doc_index.json"
    idx_path.write_text(json.dumps(idx), encoding="utf-8")
    monkeypatch.setattr("delfin.doc_server.indexer.get_default_index_path",
                        lambda: idx_path)

    rows = doctor._check_doc_index({})
    assert all(r["status"] == doctor.PASS for r in rows), rows

    future = os.path.getmtime(src) + 600
    src.write_text("# Heading\nchanged\n", encoding="utf-8")
    os.utime(src, (future, future))
    rows = doctor._check_doc_index({})
    assert any(r["status"] == doctor.WARN and "changed" in r["detail"]
               for r in rows), rows


def test_the_doctor_names_a_document_that_could_not_be_extracted(
        tmp_path, monkeypatch):
    from delfin.agent import doctor

    idx = {
        "version": 1, "built_at": "2026-01-01T00:00:00+00:00",
        "document_count": 1,
        "documents": {"ok": {"title": "OK", "source_path": "",
                             "section_count": 1, "total_chars": 10,
                             "sections": {}}},
        "failed_documents": [{"doc_id": "orca_manual",
                              "source_path": "/lit/ORCA.pdf",
                              "reason": "pypdf is not installed"}],
    }
    idx_path = tmp_path / "doc_index.json"
    idx_path.write_text(json.dumps(idx), encoding="utf-8")
    monkeypatch.setattr("delfin.doc_server.indexer.get_default_index_path",
                        lambda: idx_path)

    rows = doctor._check_doc_index({})
    warn = [r for r in rows if r["status"] == doctor.WARN]
    assert warn and "orca_manual" in warn[0]["detail"]
    assert "pypdf" in warn[0]["detail"]
