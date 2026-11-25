"""``orca_check_input`` must not answer ``ok`` when it verified nothing.

Measured in a real run (package f5, wave 42): a model sent
``{"keywords": ["nonexistent-basis-xyz"], ...}`` and a bare ``{}``, got
``{"ok": true, "n_findings": 0, ...}`` for both, and wrote in its own
answer that *"a checker that passes deliberate garbage validates nothing"*.
It then checked the tool's controls itself and concluded the result meant
nothing. The tool was not wrong about the card store -- every shipped card
is a v1 statement card (``test_w04_orca_knowledge_tools.py:63`` says so) --
but ``ok: true`` said "checked" about a description nothing had checked.

``ok`` is true only when a check actually ran and produced its verdict, so
an unknown keyword is reported instead of silently passing. The keyword
dictionary (1,481 records) is the checker that exists: every keyword
DELFIN's own templates write resolves through it, and names no manual
defines do not.

The calls go through the product's own dispatch
(``_DocToolExecutor.execute``), the path the model's call takes.
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


#: A description shaped like the tool's documented argument.
def _desc(**kw) -> str:
    base = {"keywords": [], "blocks": {}, "job": "opt", "charge": 0, "mult": 1}
    base.update(kw)
    return json.dumps(base)


def test_an_unknown_keyword_is_a_finding_not_a_pass(call):
    """The input that must fail: a keyword no manual defines.

    Red on origin/main: ``ok`` was True and ``findings`` empty for this
    description, so a model could report a nonexistent basis set as
    "nothing wrong".
    """
    out = call(_desc(keywords=["def2-SVP", "nonexistent-basis-xyz"]))
    unknown = [f for f in out.get("findings", [])
               if f.get("kind") == "unknown_keyword"]
    assert unknown, (
        "a name the dictionary does not know was not reported: "
        f"{json.dumps(out)[:300]}")
    assert any("nonexistent-basis-xyz" in str(f.get("statement", "")).lower()
               for f in unknown), "the finding does not name the unknown word"
    assert out["ok"] is False, (
        "ok must be false when a check fired with severity error")
    assert out["n_findings"] == len(out["findings"])


def test_every_keyword_of_a_template_input_resolves(call):
    """The control: a real template route is not accused of an unknown keyword.

    Passes on BOTH trees -- on origin/main because nothing is ever reported,
    and after the change because every one of these words is in the
    dictionary. Without this case the first test could be satisfied by a
    check that accuses every input.
    """
    out = call(_desc(
        keywords=["B3LYP", "def2-TZVPP", "DLPNO-CCSD(T)", "TightPNO",
                  "TightSCF", "RIJCOSX", "def2/J", "D4", "CPCM"],
        job="opt"))
    unknown = [f for f in out.get("findings", [])
               if f.get("kind") == "unknown_keyword"]
    assert not unknown, (
        "a keyword DELFIN's own templates write was called unknown: "
        f"{[f.get('statement') for f in unknown]}")
    assert out["ok"] is True
    assert out["n_findings"] == 0


def test_an_empty_description_verifies_nothing_and_says_so(call):
    """``{}`` was the other input the real run measured as a false ``ok``."""
    out = call("{}")
    assert out["ok"] is False
    assert "error" in out and "nothing to check" in out["error"].lower()
