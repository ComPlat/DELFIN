"""Literature tab: online search -> results -> "Add to my literature".

Widget level, no browser and no network: the sources answer from the
recorded responses of test_openalex_sources, and the tab's background
runner is replaced by an inline call so each step is deterministic.
"""

from __future__ import annotations

import json
import types
import urllib.error

import pytest

from delfin.agent import scholarly_sources as ss
from delfin.agent import source_apis, web_tools
from delfin.dashboard import tab_literature

from tests.test_openalex_sources import FakeNet


@pytest.fixture
def tab(monkeypatch, tmp_path):
    net = FakeNet()
    monkeypatch.setattr(ss, "_get", net)
    monkeypatch.setattr(ss, "contact_email", lambda: "")
    monkeypatch.setattr(ss, "_credential", lambda name: "")
    monkeypatch.setattr(web_tools, "outbound_secret_reason", lambda text: None)
    monkeypatch.setattr(tab_literature, "_spawn", lambda fn: fn())
    from delfin.doc_server import indexer
    monkeypatch.setattr(indexer, "build_index",
                        lambda d, quiet=True: {"document_count": 0, "documents": {}})
    monkeypatch.setattr(indexer, "get_default_index_path",
                        lambda: tmp_path / "index" / "index.json")
    source_apis.clear_cache()
    ss._COOLDOWN.clear()
    ctx = types.SimpleNamespace(repo_dir=tmp_path / "repo", run_js=lambda s: None,
                                add_init_js=lambda s: None)
    (tmp_path / "repo").mkdir()
    _widget, refs = tab_literature.create_tab(ctx)
    for key, box in refs["lit_search_sources"].items():
        box.value = key == "openalex"
    yield refs, net
    source_apis.clear_cache()


def _rows(refs):
    return list(refs["lit_search_results"].children)


def _row_title(row):
    return row.children[0].value


def test_search_panel_opens_and_searches_on_enter_only(tab):
    refs, net = tab
    refs["lit_mode"].value = "search"
    assert refs["lit_search_panel"].layout.display == "flex"
    # continuous_update=False: the value (and so a search) changes on
    # Enter, never on each keystroke.
    assert refs["lit_search_query"].continuous_update is False
    refs["lit_search_query"].value = "exact exchange"
    rows = _rows(refs)
    assert len(rows) == 2
    html = _row_title(rows[0])
    assert "Density-functional thermochemistry" in html
    assert 'href="https://doi.org/10.1063/1.464913"' in html
    assert "Cited by 100000" in html and "Found in: OpenAlex" in html
    assert 'target="_blank" rel="noopener noreferrer"' in html
    assert "2 result(s)" in refs["lit_search_status"].value
    assert len(net.calls) == 1


def test_sort_reorders_without_a_new_request(tab):
    refs, net = tab
    refs["lit_search_query"].value = "spin"
    refs["lit_search_sort"].value = "year"
    assert "2022" in _row_title(_rows(refs)[0])
    refs["lit_search_sort"].value = "citations"
    assert "Cited by 100000" in _row_title(_rows(refs)[0])
    assert len(net.calls) == 1


def test_same_search_again_comes_from_the_cache(tab):
    refs, net = tab
    refs["lit_search_query"].value = "spin"
    refs["lit_search_button"].click()
    assert len(net.calls) == 1 and len(_rows(refs)) == 2


def test_add_stores_record_and_open_pdf(tab, monkeypatch):
    refs, net = tab
    monkeypatch.setattr(ss, "download_open_pdf", lambda rec, **k: b"%PDF-1.7 test")
    refs["lit_search_query"].value = "spin"
    oa_row = next(r for r in _rows(refs) if "[PDF]" in _row_title(r))
    oa_row.children[1].click()
    lit = refs["lit_dir"]
    md = lit / "Doe2022_An_open_benchmark_of_spin_state.md"
    pdf = lit / "Doe2022_An_open_benchmark_of_spin_state.pdf"
    assert md.is_file() and pdf.read_bytes().startswith(b"%PDF-")
    text = md.read_text()
    assert "DOI: [10.1039/test.2]" in text and "Local copy: " + pdf.name in text
    assert "@article{Doe2022," in text
    assert oa_row.children[1].description == "Added"
    assert "Added" in refs["lit_search_status"].value
    # adding the same work again creates nothing new and deletes nothing
    refs["lit_search_add"](next(r for r in refs["lit_search_state"]["results"]
                                if r["doi"] == "10.1039/test.2"))
    assert "Already in your literature" in refs["lit_search_status"].value
    assert sorted(p.name for p in lit.iterdir()) == [md.name, pdf.name]


def test_add_without_open_pdf_keeps_only_the_record(tab, monkeypatch):
    refs, net = tab
    called = []
    monkeypatch.setattr(ss, "download_open_pdf", lambda rec, **k: called.append(rec))
    refs["lit_search_query"].value = "exact exchange"
    becke = next(r for r in _rows(refs) if "thermochemistry" in _row_title(r))
    becke.children[1].click()
    files = sorted(p.name for p in refs["lit_dir"].iterdir())
    assert files == ["Becke1993_Density_functional_thermochemistry_III_The_role.md"]
    assert "no open-access PDF found" in (refs["lit_dir"] / files[0]).read_text()
    assert called == []


def test_pdf_with_closed_licence_is_linked_not_stored(tab):
    refs, _net = tab
    rec = {"title": "Closed", "authors": ["A B"], "year": 2020, "doi": "10.1/x",
           "pdf_url": "https://p.example.org/x.pdf", "pdf_license": "publisher-specific-oa",
           "sources": ["openalex"]}
    refs["lit_search_add"](rec)
    files = [p.name for p in refs["lit_dir"].iterdir()]
    assert files == ["B2020_Closed.md"]
    assert "does not allow keeping a copy" in (refs["lit_dir"] / files[0]).read_text()


def test_network_error_is_one_line_and_tab_stays_usable(tab):
    refs, net = tab
    net.fail["api.openalex.org"] = urllib.error.URLError("Temporary failure in name resolution")
    refs["lit_search_query"].value = "spin"
    status = refs["lit_search_status"].value
    assert "OpenAlex: network error" in status and "<br" not in status
    assert "0 result(s)" in status and _rows(refs) == []
    del net.fail["api.openalex.org"]
    refs["lit_search_button"].click()
    assert len(_rows(refs)) == 2


def test_empty_query_and_no_source_send_nothing(tab):
    refs, net = tab
    refs["lit_search_button"].click()
    assert "Type a query" in refs["lit_search_status"].value
    for box in refs["lit_search_sources"].values():
        box.value = False
    refs["lit_search_query"].value = "spin"
    assert "at least one source" in refs["lit_search_status"].value
    assert net.calls == []


def test_mode_switch_shows_search_or_preview(tab):
    refs, _net = tab
    refs["lit_mode"].value = "search"
    assert refs["lit_search_panel"].layout.display == "flex"
    assert refs["lit_content"].layout.display == "none"
    refs["lit_mode"].value = "preview"
    assert refs["lit_search_panel"].layout.display == "none"
    assert refs["lit_content"].layout.display == "block"


def test_result_html_escapes_record_text():
    rec = {"title": "<script>x</script>", "authors": ["<b>A</b>"], "year": 2020,
           "venue": "V", "doi": "10.1/a\"b", "url": "javascript:alert(1)",
           "pdf_url": "javascript:x", "sources": ["openalex"]}
    html = tab_literature._format_hit_html(rec)
    assert "<script>" not in html and "<b>A</b>" not in html
    assert "javascript:" not in html
