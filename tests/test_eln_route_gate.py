"""The ELN route gate refuses by name and allows reads.

Foundations for a later Chemotion ELN integration: nothing here (and nothing
in delfin/lab/eln_routes.py) sends a request to any ELN instance -- these
tests only run the decision function offline.

Why named assertions instead of one loop over the table: a predicate passing
on a looked-up row is not the gate passing. Every test here drives
:func:`delfin.lab.eln_routes.decide_route` with the request's real route
string (the strings the server's own routes spell) and asserts the refusal
the gate must give, by name.
"""

from __future__ import annotations

import ast
from pathlib import Path


def _decide(method, path, body=None):
    from delfin.lab.eln_routes import decide_route

    return decide_route(method, path, body)


def _module_source() -> str:
    import delfin.lab.eln_routes as eln_routes

    return Path(eln_routes.__file__).read_text(encoding="utf-8")


# -- refusals BY NAME: the three POST-destructive routes -----------------------

def test_eln_gate_refuses_versions_revert_as_a_write():
    """POST /api/v1/versions/revert replaces an element with an earlier
    revision and elements are acts_as_paranoid (no undelete via the API). A
    verb-keyed gate ("POST is a write") would arrive at "refuse" here for the
    wrong reason; the search test below is the other half of that."""
    decision = _decide("POST", "/api/v1/versions/revert",
                       {"element_type": "Sample", "id": 7, "version_id": 2})

    assert decision.allowed is False
    assert decision.writes is True
    assert decision.reason.split(":", 1)[0] == "refused_write"
    assert "/api/v1/versions/revert" in decision.reason


def test_eln_gate_refuses_generic_elements_delete_revision_as_a_write():
    """Deletes a revision. acts_as_paranoid again means the delete stands."""
    decision = _decide("POST", "/api/v1/generic_elements/delete_revision",
                       {"element_type": "GenericElement", "id": 7,
                        "revision_id": 3})

    assert decision.allowed is False
    assert decision.writes is True
    assert decision.reason.split(":", 1)[0] == "refused_write"
    assert "/api/v1/generic_elements/delete_revision" in decision.reason


def test_eln_gate_refuses_molecules_delete_name_as_a_write():
    """Deletes a name of a molecule; there is no undo in the API."""
    decision = _decide("POST", "/api/v1/molecules/delete_name",
                       {"molecule_id": 7, "name": "old label"})

    assert decision.allowed is False
    assert decision.writes is True
    assert decision.reason.split(":", 1)[0] == "refused_write"
    assert "/api/v1/molecules/delete_name" in decision.reason


# -- the other half of the verdict: a POST that is a read ----------------------

def test_eln_gate_allows_a_search_post_as_a_read():
    """Searches are POSTs on this server (family /api/v1/search/*): the same
    verb as the destructive operations above, opposite verdict."""
    decision = _decide("POST", "/api/v1/search/samples",
                       {"search_query": {"sample": {"name": "aspirin"}}})

    assert decision.allowed is True
    assert decision.writes is False
    assert decision.reason.split(":", 1)[0] == "read"
    assert "/api/v1/search/samples" in decision.reason


# -- fail closed ---------------------------------------------------------------

def test_eln_gate_refuses_an_unknown_route():
    """A route the table does not hold is refused. /api/v1/reactions is a
    plausible route for this server and is deliberately NOT in the table:
    that is what refusing unlisted routes must feel like."""
    for path in ("/api/v1/reactions", "/api/v1/definitely_not_a_route",
                 "/api/v1/search",            # the family with no scope
                 "/api/v1/search/samples/",   # trailing slash: not its string
                 "/api/v1/search/samples/extra"):
        decision = _decide("POST", path, {"id": 1})
        assert decision.allowed is False, path
        assert decision.reason.split(":", 1)[0] == "refused_not_allowlisted", (
            path, decision.reason)
        assert path in decision.reason


def test_eln_gate_refuses_an_unparseable_request():
    """An unparseable request is refused: unknown verb, a path the table
    cannot hold, a body that was supplied but is not a JSON object."""
    bad = [
        ("FETCH", "/api/v1/search/samples", {"q": "*"}),
        (None, "/api/v1/search/samples", None),
        ("POST", "/api/v1/search/samples", "not json at all"),
        ("POST", "/api/v1/search/samples", "[1, 2]"),
        ("DELETE", "/api/v1/../etc/passwd", None),
        ("DELETE", "/api/v1/ui_state?all=1", None),
    ]
    for method, path, body in bad:
        decision = _decide(method, path, body)
        assert decision.allowed is False, (method, path, body)
        assert decision.reason.split(":", 1)[0] == "refused_not_parsed", (
            method, path, body, decision.reason)


def test_eln_route_table_has_one_entry_per_verb_and_path():
    """Data-level guard: the table itself must be a clean mapping, so a row
    cannot shadow another and silently change a decision."""
    from delfin.lab.eln_routes import ELN_ROUTES

    keys = [(row.verb, row.path) for row in ELN_ROUTES]
    assert len(keys) == len(set(keys)), keys


def test_eln_route_module_sends_nothing_and_pulls_no_client_in():
    """The gate decides offline. No HTTP surface in the module at all, and
    none of the AGPL client (LICENSE_GUARD.md) may be imported or hinted at
    by name -- this repository cannot carry that licence."""
    tree = ast.parse(_module_source())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported |= {alias.name.split(".")[0] for alias in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert not (imported & {"requests", "urllib", "http", "httpx", "socket",
                            "aiohttp", "urllib3"}), sorted(imported)
    source = _module_source()
    for forbidden in ("chemotion_api", "chemotion-api", "chemotionapi",
                      "ChemAssistent", "chemassistent"):
        assert forbidden not in source, forbidden


# -- the body-dependent blast radius -------------------------------------------

def test_eln_gate_needs_the_body_to_scope_ui_state():
    """DELETE /api/v1/ui_state without currentCollection.id in the body
    targets everything the user owns (server fact, 3657b15). The decision
    needs the request body here -- without one it cannot classify, and it
    refuses rather than guessing."""
    from delfin.lab.eln_routes import ELN_ROUTES

    row = next(r for r in ELN_ROUTES if r.path == "/api/v1/ui_state")
    assert row.verb == "DELETE"
    assert row.writes is True
    assert row.needs_body is True
    assert row.scope_fields == ("currentCollection.id",)

    no_body = _decide("DELETE", "/api/v1/ui_state")
    assert no_body.allowed is False
    assert no_body.reason.split(":", 1)[0] == "refused_needs_body"

    no_scope = _decide("DELETE", "/api/v1/ui_state", {"currentCollection": {}})
    assert no_scope.allowed is False
    assert no_scope.writes is True
    assert no_scope.reason.split(":", 1)[0] == "refused_write_unscoped"
    assert "currentCollection.id" in no_scope.reason

    scoped = _decide("DELETE", "/api/v1/ui_state",
                     {"currentCollection": {"id": 7}})
    assert scoped.allowed is False
    assert scoped.writes is True
    assert scoped.reason.split(":", 1)[0] == "refused_write_scoped"
    assert scoped.body_consulted is True
    assert scoped.scope == {"currentCollection.id": 7}


# -- the whole table, driven with the real route strings -----------------------

def test_eln_gate_allows_reads_and_refuses_writes_across_the_whole_table():
    """Every entry, driven through :func:`decide_route` with the route string
    as the table spells it (the trailing search wildcard turned into one
    concrete scope: /api/v1/search/samples). Reads are allowed, every write
    is refused as a write, and the reason names the entry."""
    from delfin.lab.eln_routes import ELN_ROUTES, decide_route

    scope_body = {"currentCollection": {"id": 7}}
    for row in ELN_ROUTES:
        path = (row.path[:-1] + "samples" if row.path.endswith("/*")
                else row.path)
        body = scope_body if row.scope_fields else None
        decision = decide_route(row.verb, path, body)
        token = decision.reason.split(":", 1)[0]
        if row.writes:
            if row.scope_fields:
                token_expected = "refused_write_scoped"
            elif row.needs_body:
                token_expected = "refused_needs_body"
            else:
                token_expected = "refused_write"
            assert decision.allowed is False, (row.path, decision)
            assert token == token_expected, (row.path, decision.reason)
        else:
            assert decision.allowed is True, (row.path, decision)
            assert token == "read", (row.path, decision.reason)
        assert decision.writes is row.writes, (row.path, decision.writes)
        assert row.reason, f"{row.path} carries no reason"
        assert row.source in ("@3657b15", "exporter"), row.path