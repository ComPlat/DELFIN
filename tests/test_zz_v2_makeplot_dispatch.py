"""make_plot through the executor's real dispatch: the write gate, the read
gate and plan mode refuse before anything is written, in both directions
(a refused call writes nothing, an admitted one writes the figure).
"""

from __future__ import annotations

import json


from delfin.agent.api_client import KitToolPermissions


def _perms(ws):
    return KitToolPermissions(workspace=ws, mode="default")


def _perms_plan(ws):
    return KitToolPermissions(workspace=ws, mode="plan")


def _perms_glob(ws, globs):
    return KitToolPermissions(workspace=ws, mode="default",
                              write_allow_globs=globs)


def _spec():
    return {"kind": "line", "title": "E per method", "x": ["a", "b"],
            "y": [1.0, 2.0]}


def test_tool_is_registered_in_the_advertised_schema():
    from delfin.agent import api_client as A
    names = [t.get("function", {}).get("name") for t in A._DOC_TOOLS_OPENAI]
    assert "make_plot" in names


def test_dispatch_returns_a_delfin_card_result(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    ex = _DocToolExecutor()
    out = ex._dispatch("make_plot", {"spec": _spec()}, _perms(str(tmp_path)))
    assert out.startswith("DELFIN_CARD:")
    obj = json.loads(out[len("DELFIN_CARD:"):])
    assert obj["card"] == "plot"
    assert obj["path"] == str((tmp_path / "plot_line.svg").resolve())
    assert "<" not in out, "the result carries no markup"


def test_dispatch_writes_only_inside_the_workspace(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    ws = tmp_path / "ws"
    ws.mkdir()
    ex = _DocToolExecutor()
    spec = _spec()
    spec["filename"] = "../escape.svg"
    out = ex._dispatch("make_plot", {"spec": spec}, _perms(str(ws)))
    assert not out.startswith("DELFIN_CARD:")  # error, not a card
    assert "error" in out.lower()
    assert not (tmp_path / "escape.svg").exists()


def test_bad_spec_is_refused_clearly(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    ex = _DocToolExecutor()
    out = ex._dispatch("make_plot", {"spec": {"kind": "pie"}},
                       _perms(str(tmp_path)))
    assert not out.startswith("DELFIN_CARD:")
    assert "kind" in out.lower()


def test_requires_permissions():
    from delfin.agent.api_client import _DocToolExecutor
    ex = _DocToolExecutor()
    out = ex._dispatch("make_plot", {"spec": _spec()}, None)
    assert "workspace" in out.lower() or "permission" in out.lower()


def test_plan_mode_is_refused(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    out = _DocToolExecutor()._dispatch("make_plot", {"spec": _spec()},
                                       _perms_plan(str(tmp_path)))
    assert not out.startswith("DELFIN_CARD:")
    assert "plan" in out.lower()
    assert not list(tmp_path.glob("*.svg")), "plan mode must not write"


def test_data_file_secret_path_is_refused_by_read_gate(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    ws = tmp_path / "ws"
    ws.mkdir()
    dd = ws / ".ssh"
    dd.mkdir()
    (dd / "data.csv").write_text("x,y\n1,2", encoding="utf-8")
    spec = {"kind": "line", "data": {"file": ".ssh/data.csv"},
            "title": "t", "filename": "from_secret.svg"}
    out = _DocToolExecutor()._dispatch("make_plot", {"spec": spec},
                                       _perms(str(ws)))
    assert not out.startswith("DELFIN_CARD:")
    assert "read denied" in out.lower() or "denied" in out.lower()
    assert not (ws / "from_secret.svg").exists()


def test_filename_outside_write_globs_is_refused(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    ws = tmp_path / "ws"
    ws.mkdir()
    spec = _spec()
    spec["filename"] = "x.svg"  # a .svg, but globs only allow allowed/*
    out = _DocToolExecutor()._dispatch("make_plot", {"spec": spec},
                                       _perms_glob(str(ws), ("allowed/*",)))
    assert not out.startswith("DELFIN_CARD:")
    assert not (ws / "x.svg").exists(), "refused by the write gate"


def test_existing_allowed_file_written_only_with_gate_ok(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    ws = tmp_path / "ws"
    ws.mkdir()
    spec = _spec()
    spec["filename"] = "plot.svg"
    out = _DocToolExecutor()._dispatch("make_plot", {"spec": spec},
                                       _perms(str(ws)))
    assert out.startswith("DELFIN_CARD:")  # gate allowed it
    assert (ws / "plot.svg").exists()


def test_a_huge_bin_count_is_refused_through_dispatch(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    spec = {"kind": "histogram", "data": [0, 1], "bins": 3_000_000}
    out = _DocToolExecutor()._dispatch("make_plot", {"spec": spec},
                                       _perms(str(tmp_path)))
    assert "bins" in json.loads(out)["error"]
    assert not list(tmp_path.iterdir())


def test_the_advertised_kinds_are_the_kinds_it_draws(tmp_path):
    """The schema named bar_by_method and hist, which the renderer refused."""
    import re
    from delfin.agent import api_client as A
    from delfin.agent import chat_plots
    desc = next(t["function"]["description"] for t in A._DOC_TOOLS_OPENAI
                if t["function"]["name"] == "make_plot")
    kinds = set(re.search(r"kind: ([a-z_|]+)", desc).group(1).split("|"))
    assert kinds == set(chat_plots._KINDS)
