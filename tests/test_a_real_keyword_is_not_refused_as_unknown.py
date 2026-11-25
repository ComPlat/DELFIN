"""A keyword ORCA defines must not be refused as "not a keyword".

``orca_check_input`` reports every described keyword the shipped dictionary
does not define with ``severity: "error"`` (``delfin/api.py``,
``_unknown_keywords``), and an error-severity finding makes ``ok`` false. The
dictionary covers 1,481 records, and ORCA knows more keywords than that, so a
real spelling the dictionary lacks was refused with the statement *"X is not a
keyword ... the input cannot run as written"* -- about an input that runs.

Measured on the branch, same call on both trees::

    keywords ["B3LYP","def2-SVP","STABILITY"], job sp
      main   -> ok=True   ids=[]
      branch -> ok=False  ids=['unknown-keyword.stability']

``STABILITY`` is not a typo. The shipped ORCA manual, section 2.16 "SCF
Stability Analysis", says verbatim: *"An SCF stability analysis with default
settings can be requested via STABILITY, SCFSTABILITY, SCFSTAB or STAB on the
simple input line."* The dictionary resolves ``STAB`` and ``SCFSTAB`` but not
``STABILITY``/``SCFSTABILITY``.

The fix keeps the finding for a name nothing knows and drops it when ORCA's own
keyword catalog (``delfin.orca_knowledge.known_keyword``, which merges the
pinned OPI list with the words ORCA's own parser accepted) does define the
word. A misspelt or invented keyword must still be an error -- the whole point
of the rule -- so both directions are pinned here.

Calls go through the product's own dispatch (``_DocToolExecutor.execute``),
the path the model's call takes.
"""
from __future__ import annotations

import json

import pytest


@pytest.fixture()
def call(tmp_path, monkeypatch):
    """One call through the executor's dispatch, in a temporary HOME+workspace."""
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    ws = tmp_path / "ws"
    ws.mkdir()
    from delfin.agent import api_client as A

    perms = A.KitToolPermissions(workspace=str(ws))
    perms.mode = "bypassPermissions"
    perms.confirm_callback = lambda *a, **k: True
    executor = A._DocToolExecutor()

    def _call(description) -> dict:
        raw = executor.execute(
            "orca_check_input", {"calc_description": description}, perms)
        out = json.loads(raw) if isinstance(raw, str) else raw
        assert isinstance(out, dict), f"dispatch answered {raw!r}"
        return out

    return _call


def _desc(**kw) -> str:
    base = {"keywords": [], "blocks": {}, "job": "sp", "charge": 0, "mult": 1}
    base.update(kw)
    return json.dumps(base)


@pytest.mark.parametrize("word", ["STABILITY", "SCFSTABILITY"])
def test_a_stability_keyword_the_manual_names_is_accepted(call, word):
    """STABILITY/SCFSTABILITY are on the simple input line per manual 2.16."""
    out = call(_desc(keywords=["B3LYP", "def2-SVP", word]))

    refused = [f for f in out["findings"]
               if str(f.get("id", "")).startswith("unknown-keyword.")]
    assert not refused, (
        f"{word} is a real ORCA keyword (manual 2.16) but was refused as "
        f"unknown: {refused}")
    assert out["ok"] is True, out


def test_a_misspelt_keyword_is_still_refused(call):
    """The control: widening the witness must not turn a typo into a pass."""
    out = call(_desc(keywords=["B3LYP", "def2-SVP", "B3LYPP"]))

    refused = [f for f in out["findings"]
               if str(f.get("id", "")).startswith("unknown-keyword.")]
    assert refused, f"a misspelt keyword must still be a finding: {out}"
    assert out["ok"] is False, (
        "ok must stay false when a genuinely unknown keyword fired an error")


def test_an_invented_basis_is_still_refused(call):
    """The reproduction from the real run that created the rule, kept green."""
    out = call(_desc(keywords=["nonexistent-basis-xyz"]))

    refused = [f for f in out["findings"]
               if str(f.get("id", "")).startswith("unknown-keyword.")]
    assert refused, out
    assert out["ok"] is False, out
