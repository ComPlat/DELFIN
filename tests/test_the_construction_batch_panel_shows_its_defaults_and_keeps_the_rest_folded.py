"""The construction batch panel: real defaults on show, an automatic run name, the rarely
needed options in a folded Advanced section, and readable messages.

Synthetic input only (invented IDs, textbook SMILES).
"""
from __future__ import annotations

import datetime
import json
from types import SimpleNamespace

from delfin.cluster_bench.provenance import TOOLS
from delfin.dashboard import construction_batch as cb

LIST = ("cisplatin;[Cl][Pt-2]([Cl])([NH3+])[NH3+]\n"
        "hexaaqua_fe;[OH2+][Fe-3]([OH2+])([OH2+])([OH2+])([OH2+])[OH2+]\n")


def _cb_panel(tmp_path, batch_text=""):
    ctx = SimpleNamespace(calc_dir=tmp_path / "dash", backend=None)
    acc = cb.create_construction_batch_panel(ctx, SimpleNamespace(value=batch_text))
    return acc.construction_batch_widgets


def _cb_descendants(w):
    yield w
    for c in getattr(w, "children", ()):
        yield from _cb_descendants(c)


def test_only_the_basic_fields_are_open_the_rest_is_folded(tmp_path):
    w = _cb_panel(tmp_path)
    basic = set(map(id, _cb_descendants(w["basic"])))
    advanced = set(map(id, _cb_descendants(w["advanced"])))
    for key in ("tool", "mode", "run_name", "list_path"):
        assert id(w[key]) in basic and id(w[key]) not in advanced, key
    for key in ("select_path", "specs_path", "tool_python", "shard_size", "timeout", "speed",
                "throttle", "repeat", "wall"):
        assert id(w[key]) in advanced and id(w[key]) not in basic, key
    assert w["advanced"].selected_index is None                 # folded
    assert w["advanced"].get_title(0) == "Advanced"


def test_the_fields_show_real_defaults_and_follow_the_builder(tmp_path):
    w = _cb_panel(tmp_path, batch_text=LIST)
    assert w["shard_size"].value == TOOLS["manta"]["shard_size"]
    assert w["mode"].value == "champion" and "shipped MANTA" in w["mode_help"].value
    assert "= 6.0 h" in w["timeout_help"].value and "21600 s" in w["timeout_help"].value
    assert "auto = 72 h (computed: 6.0 h" in w["wall_help"].value
    assert "2 molecules, 1 shard(s)" in w["plan"].value
    w["tool"].value = "mace"
    assert w["shard_size"].value == TOOLS["mace"]["shard_size"]
    assert [v for _, v in w["mode"].options] == ["paper", "extended"]
    assert w["mode"].value == "paper" and "square-planar" in w["mode_help"].value
    w["mode"].value = "extended"
    assert "TBP" in w["mode_help"].value
    w["shard_size"].value = 7                                   # a user's value stays
    w["tool"].value = "architector"
    assert w["shard_size"].value == 7
    w["speed"].value = "1.5"
    assert "9.0 h (32400 s)" in w["timeout_help"].value


def test_a_default_is_left_out_of_the_command():
    f = {"tool": "manta", "mode": "champion", "input": "in.txt", "run_dir": "RUN",
         "shard_size": TOOLS["manta"]["shard_size"], "timeout": 21600, "speed_factor": "1.0",
         "repeat": 0}
    assert cb.construction_batch_argv("prepare", f) == [
        "prepare", "--tool", "manta", "--input", "in.txt", "--run-dir", "RUN"]
    assert cb.construction_batch_argv("slurm", {"run_dir": "RUN", "throttle": 40,
                                                "time_limit": cb.CB_DEFAULT_WALL}) == [
        "slurm", "RUN", "--submit"]


def test_an_empty_run_name_becomes_builder_date_and_count(tmp_path):
    root = tmp_path / "runs"
    day = datetime.date(2026, 10, 5)
    assert cb.construction_batch_default_run_name("mace", 3, root, day) == "mace_20261005_3mol"
    (root / "mace_20261005_3mol").mkdir(parents=True)
    assert cb.construction_batch_default_run_name("mace", 3, root, day) == "mace_20261005_3mol_2"

    w = _cb_panel(tmp_path, batch_text=LIST)
    today = datetime.date.today().strftime("%Y%m%d")
    assert w["run_name"].placeholder == f"auto: manta_{today}_2mol"
    w["buttons"]["Prepare"].click()
    assert w["run_name"].value == f"manta_{today}_2mol"
    run = tmp_path / "dash" / cb.CB_RUNS_SUBDIR / w["run_name"].value
    man = json.loads((run / "manifest.json").read_text())
    assert man["n_systems"] == 2


def test_an_existing_run_is_named_before_prepare_and_never_overwritten(tmp_path):
    w = _cb_panel(tmp_path, batch_text=LIST)
    w["run_name"].value = "r1"
    w["buttons"]["Prepare"].click()
    run = tmp_path / "dash" / cb.CB_RUNS_SUBDIR / "r1"
    before = (run / "manifest.json").read_bytes()
    assert "exists" in w["run_help"].value
    w["run_name"].value = "r2"
    assert "exists" not in w["run_help"].value
    w["run_name"].value = "r1"
    w["buttons"]["Prepare"].click()
    assert (run / "manifest.json").read_bytes() == before


def test_bad_lines_of_the_field_are_named_by_their_line_number():
    rows, errors = cb.construction_batch_field_rows(
        "a;[Cl][Pt-2]([Cl])([NH3+])[NH3+]\n\nno id here\na;CCO\nb;C C\n")
    assert rows == [("a", "[Cl][Pt-2]([Cl])([NH3+])[NH3+]")]
    assert [e.split(":")[0] for e in errors] == ["line 3", "line 4", "line 5"]
    assert "twice" in errors[1]


def test_the_prepare_summary_and_the_status_line_read_like_sentences(tmp_path):
    man = {"n_systems": 1200, "settings": {"timeout_s": 3600, "workers": 48},
           "sets": {"main": {"n_shards": 5, "shards": [{"n": 250}] * 4 + [{"n": 200}]},
                    "repeat": {"n_shards": 1, "shards": [{"n": 10}]}}}
    s = cb.construction_batch_summary(man, tmp_path / "run", throttle=2)
    assert s.startswith("1200 molecules, 5 shard(s) + 1 repeat shard(s);")
    assert "~18.0 h wall with 2 parallel job(s)" in s and "1,728 core-hours" in s
    assert str(tmp_path / "run") in s
    st = {"set": "main", "n_done": 7, "n_systems": 10, "shards_done": 1, "n_shards": 3,
          "by_class": {"ok": 5, "timeout": 1, "fail": 1},
          "shards": [{"state": "done"}, {"state": "running/partial"}, {"state": "pending"}]}
    assert cb.construction_batch_status_line(st) == (
        "main: 7/10 done -- ok 5, timeout 1, failed 1 -- shards 1/3 complete, 1 in progress")
