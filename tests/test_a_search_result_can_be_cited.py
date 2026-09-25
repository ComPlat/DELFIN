"""A literature hit names the work it came from.

Input: a document's own metadata and its first page. Output: ``title``,
``authors``, ``year``, ``journal``, ``doi`` and a one-line ``reference``,
attached to the document record and carried into every search result.

Why. The index already held papers and named them by filename -- a search
returned ``1 s2.0 S0021979724009044 main``. The agent could read a paper
and repeat what it said, and could not say where it read it, which is the
difference between an opinion and a claim somebody can check.

Where the fields come from is measured over the real corpus, not assumed:
``/Subject`` is the richest source because publishers put the citation
itself there, ``/Title`` held the real title in every file measured, and
the first page is consulted for a DOI only, and only when the metadata has
none.

Offline by construction. A bibliographic service would resolve more and
would also send a list of what the user reads to a third party and fail on
a machine with no route out. A field left empty is better than either.
"""

from __future__ import annotations

import pytest

from delfin.doc_server import bibliography as B


# -- the fields ------------------------------------------------------------

@pytest.mark.parametrize("text, want", [
    ("see doi:10.1016/j.jcis.2024.04.176 for details",
     "10.1016/j.jcis.2024.04.176"),
    ("https://doi.org/10.1186/s13321-025-01008-1",
     "10.1186/s13321-025-01008-1"),
    ("ends a sentence 10.1021/acs.jctc.4c00619.", "10.1021/acs.jctc.4c00619"),
    ("10.1021/ACS.JCTC.4C00619", "10.1021/acs.jctc.4c00619"),
    ("no identifier here", ""),
    ("", ""),
])
def test_a_doi_is_found_and_normalised(text, want):
    assert B.doi_in(text) == want


def test_a_doi_may_contain_a_dot_and_still_end_cleanly():
    """The trailing strip cannot simply forbid a dot: plenty of DOIs have
    one inside, and only the last one is punctuation."""
    assert B.doi_in("10.1021/acs.jctc.4c00619.") == "10.1021/acs.jctc.4c00619"


@pytest.mark.parametrize("text, want", [
    ("Journal of Colloid And Interface Science, 668 (2024) 366-374", "2024"),
    ("J. Chem. Theory Comput. 2024.20:8367-8377", "2024"),
    ("no year", ""),
])
def test_the_year_is_read_when_it_is_there(text, want):
    assert B.year_in(text) == want


def test_the_journal_is_the_head_of_a_publisher_subject_line():
    assert B.journal_in(
        "Journal of Colloid And Interface Science, 668 (2024) 366-374. "
        "doi:10.1016/j.jcis.2024.04.176"
    ) == "Journal of Colloid And Interface Science"


def test_a_subject_that_is_an_abstract_is_not_a_journal():
    """Some producers put the abstract in /Subject. A sentence is not a
    journal name, and printing one as though it were makes the reference
    look like data it is not."""
    assert B.journal_in(
        "In this work we present a systematic study of the electronic "
        "structure of a series of transition metal complexes and their "
        "reduction potentials under aqueous conditions."
    ) == ""


# -- putting them together -------------------------------------------------

def test_a_paper_describes_itself_from_its_metadata():
    fields = B.describe(
        {"/Title": "SMILES all around",
         "/Author": "Maria H. Rasmussen",
         "/Subject": "Journal of Cheminformatics, "
                     "https://doi.org/10.1186/s13321-025-01008-1"},
        first_page="", fallback_title="s13321 025 01008 1 1")
    assert fields["title"] == "SMILES all around"
    assert fields["doi"] == "10.1186/s13321-025-01008-1"
    assert fields["journal"] == "Journal of Cheminformatics"


def test_the_first_page_supplies_a_doi_the_metadata_lacks():
    fields = B.describe({"/Title": "A paper"},
                        first_page="... https://doi.org/10.1000/xyz123 ...",
                        fallback_title="file")
    assert fields["doi"] == "10.1000/xyz123"


def test_a_file_that_says_nothing_keeps_the_title_it_had():
    """The floor is the behaviour that existed before any of this."""
    fields = B.describe({}, first_page="", fallback_title="some file name")
    assert fields["title"] == "some file name"
    assert B.reference(fields) == "some file name"


@pytest.mark.parametrize("junk", [
    "Untitled", "Microsoft Word - manuscript.docx", "document1", "Print Job",
])
def test_a_producer_placeholder_is_not_a_title(junk):
    fields = B.describe({"/Title": junk}, "", "the real file name")
    assert fields["title"] == "the real file name"


# -- the line a reader acts on --------------------------------------------

def test_several_authors_become_the_first_surname_and_the_rest():
    assert B.reference({"authors": "Hongni Jin and Kenneth M. Merz Jr.",
                        "year": "2024", "title": "T"}).startswith(
        "Jin et al. (2024)")


def test_one_author_is_left_exactly_as_printed():
    """Taking the last word of a lone author turned an institute into
    "Kohlenforschung" -- measured on the corpus. With one name there is no
    list to abbreviate, and guessing a family name guesses about people
    and organisations at once."""
    line = B.reference({"authors": "Max-Planck-Institut für Kohlenforschung",
                        "title": "ORCA Manual"})
    assert line.startswith("Max-Planck-Institut für Kohlenforschung")
    assert line.split(".")[0] != "Kohlenforschung", "the institute was mangled"


def test_a_missing_field_is_left_out_rather_than_filled():
    """"n.d." and "Anon." look like data and are not."""
    line = B.reference({"title": "A paper", "doi": "10.1/x"})
    assert "n.d." not in line and "Anon" not in line
    assert line == "A paper. doi:10.1/x"


def test_nothing_known_is_an_empty_reference():
    assert B.reference({}) == ""
    assert B.reference(None) == ""


# -- offline by construction ----------------------------------------------

def test_the_module_reaches_no_network():
    """Not a policy in a docstring: a lookup would send a list of what the
    user reads to a third party and fail where there is no route out."""
    import inspect

    src = inspect.getsource(B)
    for forbidden in ("requests", "urllib", "http://", "https://", "socket"):
        assert forbidden not in src, forbidden


# -- the wiring ------------------------------------------------------------

def test_a_search_result_carries_the_reference():
    from delfin.doc_server.search import DocSearchEngine

    index = {
        "documents": {
            "p1": {
                "title": "SMILES all around",
                "reference": "Rasmussen. SMILES all around. doi:10.1186/x",
                "sections": {"s1": {"title": "Methods",
                                    "text": "conversion of transition metal "
                                            "complexes to SMILES strings"}},
            },
            "p2": {
                "title": "Something else",
                "reference": "",
                "sections": {"s1": {"title": "Other", "text": "unrelated"}},
            },
        }
    }
    hits = DocSearchEngine(index).search("SMILES conversion", max_results=5)
    results = hits.get("results", [])
    assert results, "the engine returned nothing to cite"
    assert any(r.get("reference", "").startswith("Rasmussen")
               for r in results), "the citation did not travel with the hit"


def test_a_document_is_findable_by_its_own_title():
    """Measured before this: the corpus held per-section headings only, so
    a document's real title was searchable text nowhere. Ranks of the
    document a query names went >20->16, 2->1, 11->4, 1->1."""
    from delfin.doc_server.search import DocSearchEngine

    index = {
        "documents": {
            "paper": {
                "title": "Partial to Total Generation of 3D "
                         "Transition-Metal Complexes",
                "sections": {"s1": {"title": "2 Results",
                                    "text": "tables of computed values"}},
            },
            "manual": {
                "title": "A manual",
                "sections": {"s1": {"title": "Chapter",
                                    "text": "generation and complexes and "
                                            "transition and metal " * 20}},
            },
        }
    }
    hits = DocSearchEngine(index).search(
        "Partial to Total Generation of 3D Transition-Metal Complexes",
        max_results=5)
    titles = [r.get("doc_title", "") for r in hits.get("results", [])]
    assert titles and titles[0].startswith("Partial to Total"), titles


def test_the_indexer_attaches_a_reference_to_what_it_indexes():
    import inspect

    from delfin.doc_server import indexer

    src = inspect.getsource(indexer)
    assert '"reference": reference' in src
    assert '"bibliography": bib' in src
    assert "_bibliography.describe(" in src
