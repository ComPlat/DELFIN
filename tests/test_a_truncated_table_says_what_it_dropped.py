"""A typed tool's table, cut in the middle of a row, told the model
nothing about what was missing.

The tool-result window is a hard cap, and the head-and-tail rule that
keeps tracebacks readable is the wrong knife for a JSON table: both
models read a ranking with its middle cut out and re-fetched what they
took for missing folders, and one asked for "a hint like request fewer
folders" (2026-09-11). A table is now cut at a row boundary, from the
end of the longest list of rows, and carries how many rows of how many
went and what to do about it. What is not a table keeps the old rule.
"""

from __future__ import annotations

import json

from delfin.agent import api_client as A


def _rows(n, width=60):
    return [{"folder": f"calc/run_{i:03d}", "method": "PBE0/def2-SVP",
             "single_point": -113.0 - i / 1000, "note": "x" * width} for i in range(n)]


def test_a_grouped_table_is_cut_at_a_row_and_says_so():
    obj = {"note": "compare within a method", "root": "/archive",
           "groups": [{"method": "PBE0/def2-SVP", "rows": _rows(30)},
                      {"method": "B3LYP/def2-SVP", "rows": _rows(30)}],
           "skipped": [{"folder": "calc/run_x", "reason": "no output"}]}
    text = json.dumps(obj)
    cut = A._truncate_table_result(text, cap=2500)
    assert cut and len(cut) <= 2500
    out = json.loads(cut)                                   # still JSON, whole rows only
    assert out["truncated"] is True
    kept = sum(len(g["rows"]) for g in out["groups"])
    assert out["rows_omitted"] == 60 - kept > 0
    assert out["rows_total"] == 61                          # skipped rows count too
    assert "fewer folders" in out["hint"]
    assert out["note"] == obj["note"] and out["root"] == obj["root"]
    for g in out["groups"]:
        for r in g["rows"]:
            assert set(r) == {"folder", "method", "single_point", "note"}


def test_the_cut_spreads_across_groups():
    obj = {"groups": [{"method": "A", "rows": _rows(40)}, {"method": "B", "rows": _rows(40)}]}
    out = json.loads(A._truncate_table_result(json.dumps(obj), cap=3000))
    sizes = [len(g["rows"]) for g in out["groups"]]
    assert all(n > 0 for n in sizes), sizes
    assert abs(sizes[0] - sizes[1]) <= 1


def test_a_plain_list_is_cut_and_ends_with_the_note():
    text = json.dumps(_rows(50))
    out = json.loads(A._truncate_table_result(text, cap=2000))
    assert out[-1]["truncated"] is True and out[-1]["rows_omitted"] > 0
    assert all("folder" in r for r in out[:-1])


def test_a_result_that_fits_is_left_alone():
    assert A._truncate_table_result(json.dumps(_rows(3)), cap=5000) == ""


def test_what_is_not_a_table_keeps_the_head_and_tail_rule():
    trace = "Traceback (most recent call last):\n" + "  frame\n" * 400 + "ValueError: boom\n"
    assert A._truncate_table_result(trace, cap=2000) == ""
    text = A._smart_truncate(trace, cap=2000, label="tool_result")
    assert "ValueError: boom" in text and "truncated" in text
    assert A._truncate_table_result(json.dumps({"error": "x" * 9000}), cap=500) == ""


def test_the_flag_the_event_reads_is_set():
    """The engine reports whether the model saw the whole result by
    looking for '"truncated": true'; a table cut this way says it in the
    JSON, and the compact dump spells it the way the check expects."""
    text = json.dumps(_rows(50))
    cut = A._truncate_table_result(text, cap=2000)
    assert '"truncated":true' in cut
