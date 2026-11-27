"""The ELN client and its approval: reads answer, writes cannot be sent.

What each test binds to:

* the credential never leaves the store: a token that is set appears in NO
  message, no ``repr``, no URL and no error, while the request still carries
  it as a bearer header;
* a read answers from the server's real shapes, pages through more than one
  page, and reports when a row has no analysis;
* a write WITHOUT a person's approval sends nothing -- including in
  ``bypassPermissions``, which is the mode a standing permission comes from;
* a write WITH an approval sends once, and the same approval cannot be used
  twice, cannot carry a different body, and cannot be granted by an agent
  session.

The far end is ``tests/_eln_fake.py``: an in-process function, no socket,
no real server. Every shape it serves is transcribed from the API's source
and annotated there.
"""

from __future__ import annotations

import json

import pytest

from _eln_fake import FAKE_URL, FakeELN


@pytest.fixture
def eln(tmp_path, monkeypatch):
    """A scratch HOME, lab dir, person and credential store."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("DELFIN_LAB_DIR", str(tmp_path / "lab"))
    monkeypatch.setenv("DELFIN_PERSON", "alice")
    monkeypatch.setenv("DELFIN_ELN_URL", FAKE_URL)
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    from delfin.agent import trust_store
    trust_store.forget_cached_keys()
    from delfin.agent import credentials
    # The path conftest's autouse guard installed for THIS test. Working
    # out a path of my own would write a file the product never reads,
    # which is how "the client found the token" passes while the real
    # lookup finds nothing.
    assert credentials._DEFAULT_PATH.parent == tmp_path
    assert not credentials._DEFAULT_PATH.exists()
    credentials.set_credential("CHEMOTION_API_TOKEN", "chemtoken_secret_XYZ",
                               store="file")
    return tmp_path


def _client(fake: FakeELN):
    from delfin.lab import eln_api
    from delfin.lab.eln_tools import _credential
    return eln_api.Client(FAKE_URL, _credential, transport=fake.send)


# ---------------------------------------------------------------------------
# The credential
# ---------------------------------------------------------------------------

def test_the_token_is_sent_as_a_bearer_and_never_appears_in_a_message(eln):
    from delfin.lab import eln_api

    fake = FakeELN()
    client = _client(fake)
    listing = client.collections()
    assert listing["own"][0]["label"] == "Synthesis 2026"
    assert fake.calls[-1]["auth"] == "Bearer chemtoken_secret_XYZ"
    # The secret appears in no string this object can produce.
    secret = "chemtoken_secret_XYZ"
    assert secret not in repr(client)
    assert "set" in repr(client)
    # ... and not in the display of the address either.
    assert secret not in eln_api.display_url(
        "https://user:chemtoken_secret_XYZ@eln.example.org/api")


def test_a_missing_credential_names_the_key_and_not_the_person(eln, tmp_path):
    from delfin.lab import eln_api

    fake = FakeELN()
    client = eln_api.Client(FAKE_URL, lambda name: "", transport=fake.send)
    with pytest.raises(eln_api.ClientError) as excinfo:
        client.collections()
    text = str(excinfo.value)
    assert "CHEMOTION_API_TOKEN" in text          # what to store
    assert "alice" not in text                    # not who
    assert not fake.calls                         # nothing was sent


def test_a_wrong_token_quotes_the_servers_401(eln):
    from delfin.lab import eln_api

    fake = FakeELN()
    client = eln_api.Client(FAKE_URL, lambda name: "not-a-chemtoken",
                            transport=fake.send)
    with pytest.raises(eln_api.ClientError) as excinfo:
        client.collections()
    assert excinfo.value.status == 401
    assert "401 Unauthorized" in str(excinfo.value)


def test_a_server_error_is_quoted_not_paraphrased(eln):
    from delfin.lab import eln_api

    fake = FakeELN()
    fake.refuse_with = {"status": 403, "error": "403 Forbidden"}
    client = _client(fake)
    with pytest.raises(eln_api.ClientError) as excinfo:
        client.collections()
    assert excinfo.value.status == 403
    assert "403 Forbidden" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------

def test_a_search_reads_one_page_and_reports_the_total(eln):
    fake = FakeELN()
    page = _client(fake).search("samples", {}, per_page=2, page=1)
    assert len(page.rows) == 2
    assert page.total == 5
    assert page.total_pages == 3
    assert fake.calls[-1]["path"] == "/api/v1/search/samples"
    assert fake.calls[-1]["body"]["page"] == 1


def test_a_search_follows_the_pages_rather_than_answering_from_page_one(eln):
    """The population of "which samples have no analysis" spans pages.

    With four per page five samples need two, and page one alone would
    have reported S-003 and S-004 while missing S-005 -- and would have
    missed any sample with no analysis on page two.
    """
    fake = FakeELN()
    rows = _client(fake).search_all("samples", {}, per_page=4)
    assert [r["id"] for r in rows] == [11, 12, 13, 14, 15]
    assert len(fake.calls) == 2


def test_an_element_read_carries_its_analyses(eln):
    fake = FakeELN()
    row = _client(fake).element("samples", 11)
    assert row["name"] == "S-001"
    kinds = [c["extended_metadata"]["kind"] for c in
             row["container"]["children"]
             if c["container_type"] == "analysis"]
    assert kinds == ["nmr"]


# ---------------------------------------------------------------------------
# The route gate: what may be sent at all
# ---------------------------------------------------------------------------

def test_a_collection_read_is_allowed_by_the_gate(eln):
    """The gate allowed only searches on the base; a collection name could
    not be resolved to an id, so no search could be scoped."""
    from delfin.lab import eln_routes
    decision = eln_routes.decide_route("GET", "/api/v1/collections")
    assert decision.allowed is True
    assert decision.writes is False


def test_a_create_is_refused_by_the_gate_before_any_transport(eln):
    """The client cannot write at all: only the approval path sends a POST
    to /reactions, and it does so through ``_send`` after the grant."""
    from delfin.lab import eln_api

    fake = FakeELN()
    client = _client(fake)
    with pytest.raises(eln_api.ClientError) as excinfo:
        client._send("POST", "/api/v1/reactions", body={"name": "x"})
    assert "refused before sending" in str(excinfo.value)
    assert "refused_write" in str(excinfo.value)
    assert not fake.calls, "nothing may reach the server"


def test_a_delete_is_refused_by_the_gate(eln):
    from delfin.lab import eln_api

    fake = FakeELN()
    with pytest.raises(eln_api.ClientError) as excinfo:
        _client(fake)._send("DELETE", "/api/v1/samples/11",
                            body={"id": 11})
    assert excinfo.value.status is None        # never sent
    assert not fake.calls


# ---------------------------------------------------------------------------
# The write: a draft is filed, and only a person's grant sends it
# ---------------------------------------------------------------------------

def _draft(fake: FakeELN, **kwargs):
    from delfin.lab import eln_write
    client = _client(fake)
    return eln_write.create_draft(
        client, "reaction",
        {"name": "Oxidation of 4a", "collection_id": 3,
         "materials": {"reactants": [1]}}, **kwargs)


def test_a_draft_without_approval_sends_nothing_and_names_the_draft(eln):
    from delfin.lab import eln_write

    fake = FakeELN()
    with pytest.raises(eln_write.WriteRefused) as excinfo:
        _draft(fake)
    exc = excinfo.value
    assert exc.code == "needs_approval"
    assert not fake.calls, "a filed request must not reach the server"
    # What it says: the element, its name, the collection, the command.
    text = str(exc)
    assert "reaction draft 'Oxidation of 4a'" in text
    assert "collection 3" in text
    assert "1 reactants" in text
    assert f"delfin-agent eln approve {exc.request_id}" in text
    assert exc.digest[:12] in text


def test_bypass_permissions_cannot_approve_a_draft(eln):
    """The mode is not consulted. Two calls, two modes, one refusal."""
    from delfin.lab import eln_tools

    fake = FakeELN()
    outs = []
    # A real project directory, not the pytest tmp root: the workspace
    # guard refuses a home directory or a system root, and a fixture that
    # trips it would be measuring the guard instead of the tool.
    ws = eln / "ws"
    ws.mkdir(exist_ok=True)
    for mode in ("bypassPermissions", "plan", "default"):
        from delfin.agent import api_client as A
        perms = A.KitToolPermissions(workspace=ws)
        perms.mode = mode
        perms.confirm_callback = lambda *a, **k: True   # answers every dialog
        raw = eln_tools.execute(
            "eln_draft",
            {"element": "reaction", "name": "Oxidation of 4a",
             "collection_id": 3, "_transport": fake.send},
            perms)
        outs.append(json.loads(raw))
    for out in outs:
        assert out["code"] == "needs_approval", out
        assert "not decided by the permission mode" in out["note"]
    assert not fake.calls
    # All three filed the SAME request (one digest), so a person approves
    # one thing rather than three spellings of it.
    assert len({o["request_id"] for o in outs}) == 1


def test_an_approved_draft_is_sent_once_and_the_grant_cannot_be_reused(eln):
    from delfin.lab import eln_write

    fake = FakeELN()
    with pytest.raises(eln_write.WriteRefused) as excinfo:
        _draft(fake)
    rid = excinfo.value.request_id
    eln_write.grant(rid)                        # the person, in their terminal

    out = _draft(fake, approval=rid)
    assert out["ok"] is True
    assert out["sent"]["name"] == "Oxidation of 4a"
    posts = [c for c in fake.calls if c["method"] == "POST"]
    assert len(posts) == 1
    assert posts[0]["path"] == "/api/v1/reactions"
    assert posts[0]["body"]["collection_id"] == 3

    with pytest.raises(eln_write.WriteRefused) as again:
        _draft(fake, approval=rid)
    assert again.value.code == "refused_approval"
    assert "already used" in str(again.value)
    assert len([c for c in fake.calls if c["method"] == "POST"]) == 1


def test_an_approval_cannot_carry_a_different_body(eln):
    from delfin.lab import eln_write

    fake = FakeELN()
    with pytest.raises(eln_write.WriteRefused) as excinfo:
        _draft(fake)
    rid = excinfo.value.request_id
    eln_write.grant(rid)
    # The same request id, a DIFFERENT name: the digest no longer matches.
    with pytest.raises(eln_write.WriteRefused) as excinfo2:
        eln_write.create_draft(
            _client(fake), "reaction",
            {"name": "Boc deprotection", "collection_id": 3,
             "materials": {"reactants": [1]}}, approval=rid)
    assert excinfo2.value.code == "refused_approval"
    assert "different write" in str(excinfo2.value)
    assert not [c for c in fake.calls if c["method"] == "POST"]


def test_an_agent_session_cannot_grant_its_own_write(eln, monkeypatch):
    from delfin.lab import eln_write

    fake = FakeELN()
    with pytest.raises(eln_write.WriteRefused) as excinfo:
        _draft(fake)
    rid = excinfo.value.request_id
    monkeypatch.setenv("DELFIN_AGENT_SESSION", "4242")
    from delfin.agent import trust_store
    trust_store.forget_cached_keys()
    with pytest.raises(eln_write.WriteRefused) as refused:
        eln_write.grant(rid)
    assert refused.value.code == "refused_agent_session"
    assert "person" in str(refused.value)


def test_the_server_refusal_of_an_approved_write_is_quoted(eln):
    """An approval gets the request SENT; it does not promise the server
    will accept it. A 403 comes back as the server's own sentence."""
    from delfin.lab import eln_api, eln_write

    fake = FakeELN()
    fake.own = []                               # nothing writable now
    with pytest.raises(eln_write.WriteRefused) as excinfo:
        _draft(fake)
    eln_write.grant(excinfo.value.request_id)
    with pytest.raises(eln_api.ClientError) as sent:
        _draft(fake, approval=excinfo.value.request_id)
    assert sent.value.status == 403
    assert "403 Forbidden" in str(sent.value)
    # The grant was consumed even though the send failed: a retry needs a
    # new approval, so a failed send cannot silently become two elements.
    with pytest.raises(eln_write.WriteRefused) as again:
        _draft(fake, approval=excinfo.value.request_id)
    assert again.value.code == "refused_approval"


def test_an_unnamed_or_unknown_element_is_refused_before_filing(eln):
    from delfin.lab import eln_write

    fake = FakeELN()
    with pytest.raises(eln_write.WriteRefused) as excinfo:
        eln_write.create_draft(_client(fake), "reaction", {"collection_id": 3})
    assert excinfo.value.code == "refused_unnamed"
    with pytest.raises(eln_write.WriteRefused) as other:
        eln_write.create_draft(_client(fake), "wellplate", {"name": "x"})
    assert other.value.code == "refused_unknown_element"
    assert not fake.calls
