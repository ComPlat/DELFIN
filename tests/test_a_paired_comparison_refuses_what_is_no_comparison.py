"""A paired comparison refuses what is no comparison and lists every
never-worse case.

Input: two arm result lists, their provenance stamps and the declared
axis, through ``paired_compare.compare`` and the ``delfin-agent
provenance`` / ``compare-arms`` commands.

Pinned here:

* a different population (an item only in one arm, an item twice in an
  arm) is ``skewed_pair`` and no numbers are compared;
* the code or the platform differing outside the declared axis is
  ``skewed_pair``; a behaviour variable of the project that differs and
  was not declared is ``skewed_pair`` even when its name shares nothing
  with the axis;
* the axis at the same effective value in both arms (a switch at its
  default, resolved through ``--default``) is ``null_test_suspected``;
  null -> value without a declared default is ``axis_unknown``;
* a behaviour switch set identically in both arms is refused unless the
  caller opts in; a path value and an empty value are no switch, a number
  is;
* identical leading items flag a null test;
* one ok -> not ok item, one item whose defect count rose, one quality
  regression or one nondeterministic rerun item blocks, each listed;
* fewer paired items than the minimum is ``insufficient_n``;
* a clean pair is ``needs_human_check`` with a stratified sample, never
  a landing;
* the runner's stamp names code that changed after the runner started;
* a saved stamp is written once under one plain name and reads back as
  sealed; a stamp read from any other file is reported as unsealed;
* both commands are read-only to the shell gate in an ordinary
  workspace, and the CLI's exit code says whether a comparison stands.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.agent import paired_compare as P


def _stamp(**over):
    base = {
        "schema": P.SCHEMA, "taken_at": 1000.0, "platform": "Linux-x",
        "machine": "x86_64", "python": "3.11.9",
        "git": {"toplevel": "/r", "head": "abc", "dirty_hash": ""},
        "code": {"hash": "h1", "files": 3, "missing": []},
        "prefixes": ["DELFIN_", "PROJ_"], "env_source": "command",
        "env": {},
    }
    for k, v in over.items():
        base[k] = v
    return base


def _arm(n=120, *, ok=lambda i: True, extra=lambda i: {}):
    return [{"id": f"i{i:03d}", "ok": ok(i), "value": i, **extra(i)}
            for i in range(n)]


def _changed(n=120):
    return _arm(n, extra=lambda i: {"value": i + (1 if i % 2 else 0)})


def _pair(**kw):
    a = kw.pop("a", None) or _arm()
    b = kw.pop("b", None) or _changed()
    sa = kw.pop("stamp_a", None) or _stamp(env={"PROJ_FAST": "0"})
    sb = kw.pop("stamp_b", None) or _stamp(env={"PROJ_FAST": "1"})
    kw.setdefault("axis", ["PROJ_FAST"])
    kw.setdefault("ok", "ok")
    return P.compare(a, b, stamp_a=sa, stamp_b=sb, env={}, **kw)


# ---------------------------------------------------------------------------
# the population
# ---------------------------------------------------------------------------

def test_a_clean_pair_asks_a_person_and_never_lands():
    v = _pair()
    assert v["state"] == "needs_human_check", v["reasons"]
    assert v["totals"]["n"] == 120 and v["totals"]["ok_a"] == 120
    assert v["never_worse"] == []
    assert len(v["sample"]) == 100
    # every changed item comes first in the sample
    assert set(v["sample"][:60]) == {f"i{i:03d}" for i in range(1, 120, 2)}
    text = P.format_report(v)
    assert "Not a landing" in text and "user's own message" in text
    assert P.exit_code(v) == 0


def test_an_item_in_one_arm_only_is_a_different_population():
    v = _pair(b=_changed()[:-1])
    assert v["state"] == "skewed_pair"
    assert "only in A" in " ".join(v["reasons"])
    assert "totals" not in v
    assert P.exit_code(v) == 3


def test_an_item_listed_twice_moves_the_denominator():
    b = _changed()
    b.append(dict(b[0]))
    v = _pair(b=b)
    assert v["state"] == "skewed_pair"
    assert "more than once" in " ".join(v["reasons"])


# ---------------------------------------------------------------------------
# one instrument, the declared axis
# ---------------------------------------------------------------------------

def test_code_that_differs_outside_the_axis_is_no_pair():
    v = _pair(stamp_b=_stamp(env={"PROJ_FAST": "1"},
                             code={"hash": "h2", "missing": []}))
    assert v["state"] == "skewed_pair"
    assert "code differs" in " ".join(v["reasons"])


def test_code_as_the_axis_with_equal_code_is_a_null_test():
    v = _pair(axis=["code"], stamp_a=_stamp(), stamp_b=_stamp())
    assert v["state"] == "null_test_suspected"


def test_the_platform_differing_is_no_pair_unless_declared():
    sb = _stamp(env={"PROJ_FAST": "1"}, python="3.13.1")
    assert _pair(stamp_b=sb)["state"] == "skewed_pair"
    v = _pair(stamp_b=sb, axis=["PROJ_FAST", "platform"])
    assert v["state"] == "needs_human_check"


def test_an_undeclared_variable_of_the_project_under_another_name():
    # The field case: the declared axis was one variable, the arms also
    # differed in another variable of the same project with a different
    # name. The checkout prefix (PROJ_) catches it, not the axis prefix.
    sa = _stamp(env={"PROJ_FAST": "0", "PROJ_CACHE_MODE": "warm"})
    sb = _stamp(env={"PROJ_FAST": "1", "PROJ_CACHE_MODE": "cold"})
    v = _pair(stamp_a=sa, stamp_b=sb)
    assert v["state"] == "skewed_pair"
    assert "PROJ_CACHE_MODE" in " ".join(v["reasons"])


def test_a_determinism_variable_differing_is_no_pair():
    sa = _stamp(env={"PROJ_FAST": "0", "OMP_NUM_THREADS": "1"})
    sb = _stamp(env={"PROJ_FAST": "1", "OMP_NUM_THREADS": "8"})
    assert _pair(stamp_a=sa, stamp_b=sb)["state"] == "skewed_pair"


def test_a_switch_at_its_default_is_no_ab():
    # Arm A leaves the switch unset, arm B sets it to the code's default:
    # both arms ran the same behaviour.
    v = _pair(stamp_a=_stamp(env={}), stamp_b=_stamp(env={"PROJ_FAST": "1"}),
              defaults={"PROJ_FAST": "1"})
    assert v["state"] == "null_test_suspected"
    assert "same effective value" in " ".join(v["reasons"])


def test_null_to_value_without_a_default_is_a_new_measurement():
    v = _pair(stamp_a=_stamp(env={}), stamp_b=_stamp(env={"PROJ_FAST": "1"}))
    assert v["state"] == "axis_unknown"
    assert "new measurement" in " ".join(v["reasons"])


def test_zero_is_a_value_not_off():
    # Whether 0 means off depends on the code's default; with default 1 an
    # explicit 0 is a real axis.
    v = _pair(stamp_a=_stamp(env={"PROJ_FAST": "0"}),
              stamp_b=_stamp(env={}), defaults={"PROJ_FAST": "1"})
    assert v["state"] == "needs_human_check", v["reasons"]


def test_a_shared_switch_is_refused_unless_opted_in():
    sa = _stamp(env={"PROJ_FAST": "0", "PROJ_THRESHOLD": "0.5"})
    sb = _stamp(env={"PROJ_FAST": "1", "PROJ_THRESHOLD": "0.5"})
    v = _pair(stamp_a=sa, stamp_b=sb)
    assert v["state"] == "undeclared_switch"
    assert "PROJ_THRESHOLD" in " ".join(v["reasons"])
    ok = P.compare(_arm(), _changed(), stamp_a=sa, stamp_b=sb,
                   axis=["PROJ_FAST"], ok="ok",
                   env={P.ALLOW_SHARED_ENV: "1"})
    assert ok["state"] == "needs_human_check"


def test_a_data_path_an_empty_value_and_the_harness_marker_are_no_switch():
    shared = {"PROJ_DATA": "/data/pool.csv", "PROJ_EMPTY": "",
              "DELFIN_AGENT_SESSION": "4711"}
    v = _pair(stamp_a=_stamp(env={"PROJ_FAST": "0", **shared}),
              stamp_b=_stamp(env={"PROJ_FAST": "1", **shared}))
    assert v["state"] == "needs_human_check", v["reasons"]


def test_an_arm_whose_environment_could_not_be_read_is_checked_by_hand():
    v = _pair(stamp_b=_stamp(env=None, env_source="unavailable"))
    assert v["state"] == "axis_unknown"


def test_no_stamp_or_no_axis_is_axis_unknown():
    v = P.compare(_arm(), _changed(), ok="ok", env={})
    assert v["state"] == "axis_unknown"
    assert P.exit_code(v) == 3


# ---------------------------------------------------------------------------
# reach, never-worse, minimum n
# ---------------------------------------------------------------------------

def test_identical_arms_are_a_null_test():
    v = _pair(b=_arm())
    assert v["state"] == "null_test_suspected"
    assert "none of the 120 items differs" in " ".join(v["reasons"])


def test_identical_leading_items_flag_a_null_test():
    b = _arm(extra=lambda i: {"value": i + (1 if i >= 30 else 0)})
    v = _pair(b=b)
    assert v["state"] == "null_test_suspected"
    assert "first 24 items" in " ".join(v["reasons"])


def test_volatile_fields_do_not_count_as_a_change():
    a = _arm(extra=lambda i: {"duration_s": 1.0, "ts": 5})
    b = _arm(extra=lambda i: {"duration_s": 2.5, "ts": 9})
    assert _pair(a=a, b=b)["state"] == "null_test_suspected"


def test_one_lost_item_blocks_and_is_named():
    b = _changed()
    b[7]["ok"] = False
    v = _pair(b=b)
    assert v["state"] == "blocked"
    assert v["never_worse"] == [{"item": "i007", "term": "ok: ok -> not ok"}]
    assert "never by loosening a threshold" in " ".join(v["reasons"])
    assert "i007" in P.format_report(v)


def test_a_rising_defect_count_blocks_even_when_the_rate_improves():
    a = _arm(ok=lambda i: i % 10 != 0, extra=lambda i: {"crashes": 0})
    b = _changed()
    for r in b:
        r["crashes"] = 0
    b[3]["crashes"] = 2
    v = _pair(a=a, b=b, defects=["crashes"])
    assert v["totals"]["ok_b"] > v["totals"]["ok_a"]
    assert v["state"] == "blocked"
    assert {"item": "i003", "term": "crashes: 0 -> 2"} in v["never_worse"]
    assert "absolute crashes rose: 0 -> 2" in v["reasons"]


def test_a_quality_regression_on_one_item_blocks():
    a = _arm(extra=lambda i: {"score": 10})
    b = _arm(extra=lambda i: {"score": 11 if i != 5 else 9})
    v = _pair(a=a, b=b, quality="score")
    assert v["state"] == "blocked"
    assert v["never_worse"] == [{"item": "i005", "term": "score: 10 -> 9"}]
    low = _pair(a=a, b=b, quality="score", lower_is_better=True)
    assert low["never_worse_count"] == 119


def test_a_nondeterministic_rerun_blocks():
    rerun = _arm()
    rerun[11]["value"] = -1
    v = _pair(rerun_a=rerun)
    assert v["state"] == "blocked"
    assert v["totals"]["nondeterministic_a"] == 1


def test_below_the_minimum_n_the_user_decides():
    v = _pair(a=_arm(40), b=_changed(40))
    assert v["state"] == "insufficient_n"
    assert "the user decides" in " ".join(v["reasons"])
    assert _pair(a=_arm(40), b=_changed(40), min_n=20)["state"] == (
        "needs_human_check")


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------

def test_the_command_assignments_are_the_arms_environment():
    env = P.command_env("cd run && PROJ_FAST=1 env -u PROJ_OLD PROJ_N=3 "
                        "python bench.py", {"PROJ_OLD": "x", "HOME": "/h"})
    assert env["PROJ_FAST"] == "1" and env["PROJ_N"] == "3"
    assert "PROJ_OLD" not in env
    exported = P.command_env("export PROJ_FAST=0; python bench.py", {})
    assert exported["PROJ_FAST"] == "0"
    # assignments of an earlier program do not leak into the next one
    assert "PROJ_A" not in P.command_env("PROJ_A=1 true && python b.py", {})


def test_a_stamp_takes_the_project_prefixes_and_redacts_secrets(tmp_path):
    code = tmp_path / "proj"
    code.mkdir()
    (code / "m.py").write_text("x = 1\n")
    s = P.provenance(code=[str(code)], command="PROJ_FAST=1 python m.py",
                     prefixes=["PROJ_"], cwd=str(tmp_path),
                     env={"PROJ_API_TOKEN": "s3cr3t", "UNRELATED": "1",
                          "OMP_NUM_THREADS": "4"})
    assert s["env"] == {"OMP_NUM_THREADS": "4", "PROJ_API_TOKEN": "<redacted>",
                        "PROJ_FAST": "1"}
    assert s["code"]["files"] == 1 and s["code"]["missing"] == []
    h1 = s["code"]["hash"]
    (code / "m.py").write_text("x = 2\n")
    assert P.code_hash([str(code)])["hash"] != h1
    assert P.code_hash([str(tmp_path / "nope")])["missing"]


def test_a_runner_started_before_its_code_changed_is_named(tmp_path):
    code = tmp_path / "tool.py"
    code.write_text("x = 1\n")
    os.utime(code, (time.time() - 30, time.time() - 30))
    proc = subprocess.Popen([sys.executable, "-c",
                             "import sys; sys.stdin.read()"],
                            stdin=subprocess.PIPE,
                            env={**os.environ, "PROJ_FAST": "1"})
    try:
        deadline = time.time() + 10
        stamp = None
        while time.time() < deadline:
            stamp = P.provenance(code=[str(code)], pid=proc.pid,
                                 prefixes=["PROJ_"], cwd=str(tmp_path))
            if stamp.get("process", {}).get("start_epoch"):
                break
            time.sleep(0.05)
        start = stamp["process"]["start_epoch"]
        assert start and abs(start - time.time()) < 120
        assert stamp["code_changed_after_start"] == []
        os.utime(code, (start + 5, start + 5))
        later = P.provenance(code=[str(code)], pid=proc.pid,
                             prefixes=["PROJ_"], cwd=str(tmp_path))
        assert later["code_changed_after_start"] == [str(code)]
        if later["env_source"] == "pid":
            assert later["env"].get("PROJ_FAST") == "1"
        else:
            assert later["env"] is None
        ax = P.axis_check(later, later, ["code"])
        assert ax["state"] == "skewed_pair"
        assert "runner started before" in " ".join(ax["reasons"])
    finally:
        proc.stdin.close()
        proc.wait(timeout=10)


def test_records_load_from_jsonl_json_list_and_wrapped_list(tmp_path):
    rows = [{"id": "a"}, {"id": "b"}]
    (tmp_path / "x.jsonl").write_text("\n".join(json.dumps(r) for r in rows))
    (tmp_path / "x.json").write_text(json.dumps(rows))
    (tmp_path / "w.json").write_text(json.dumps({"items": rows, "n": 2}))
    for name in ("x.jsonl", "x.json", "w.json"):
        assert P.load_records(tmp_path / name) == rows


def test_a_saved_stamp_is_written_once_and_reads_back_sealed(tmp_path):
    path = P.save_stamp("arm-a.1", _stamp())
    assert path.parent == P.stamp_dir() and path.stat().st_mode & 0o077 == 0
    assert P.load_stamp("arm-a.1")["sealed"] is True
    with pytest.raises(ValueError, match="written once"):
        P.save_stamp("arm-a.1", _stamp(python="3.13"))
    for bad in ("../x", "a/b", "", "-x", "x" * 65):
        with pytest.raises(ValueError):
            P.save_stamp(bad, _stamp())
    loose = tmp_path / "s.json"
    loose.write_text(json.dumps(_stamp()))
    assert P.load_stamp(str(loose))["sealed"] is False
    v = _pair(stamp_a=P.load_stamp(str(loose)), stamp_b=P.load_stamp(
        str(loose)), axis=["platform"])
    assert v["stamps_sealed"] is False
    assert "not from the saved store" in P.format_report(v)


# ---------------------------------------------------------------------------
# the commands
# ---------------------------------------------------------------------------

def _cli(*argv, cwd, env):
    return subprocess.run([sys.executable, "-m", "delfin.agent.cli", *argv],
                          cwd=cwd, env=env, capture_output=True, text=True,
                          timeout=120)


def test_the_commands_print_and_exit_by_the_verdict(tmp_path):
    from conftest import child_env
    env = child_env(tmp_path)
    env["PYTHONPATH"] = str(Path(P.__file__).resolve().parents[2])
    a, b = tmp_path / "a.jsonl", tmp_path / "b.jsonl"
    a.write_text("\n".join(json.dumps(r) for r in _arm()))
    b.write_text("\n".join(json.dumps(r) for r in _changed()))
    for arm, val in (("a", "0"), ("b", "1")):
        r = _cli("provenance", "--code", str(a), "--prefix", "PROJ_",
                 "--command", f"PROJ_FAST={val} python run.py",
                 "--save", f"t-{arm}", cwd=str(tmp_path), env=env)
        assert r.returncode == 0, r.stderr
        assert json.loads(r.stdout)["env"]["PROJ_FAST"] == val
        saved = [p for p in tmp_path.rglob(f"provenance/t-{arm}.json")]
        assert len(saved) == 1, saved
    again = _cli("provenance", "--save", "t-a", cwd=str(tmp_path), env=env)
    assert again.returncode == 2 and "written once" in again.stderr
    r = _cli("compare-arms", str(a), str(b), "--stamp-a", "t-a",
             "--stamp-b", "t-b", "--axis", "PROJ_FAST", "--ok", "ok",
             cwd=str(tmp_path), env=env)
    assert r.returncode == 0, r.stdout + r.stderr
    assert r.stdout.startswith("state: needs_human_check")
    assert "saved store" not in r.stdout
    r = _cli("compare-arms", str(a), str(a), "--stamp-a", "t-a",
             "--stamp-b", "t-a", "--axis", "PROJ_FAST", "--json",
             cwd=str(tmp_path), env=env)
    assert r.returncode == 3
    assert json.loads(r.stdout)["state"] == "null_test_suspected"


@pytest.mark.parametrize("cmd", [
    "delfin-agent provenance --code src --command 'PROJ_FAST=1 python run.py'",
    "delfin-agent provenance --pid 4242 --prefix PROJ_ --save run-b.2",
    "delfin-agent compare-arms a.jsonl b.jsonl --stamp-a a.json --stamp-b "
    "b.json --axis PROJ_FAST --key id --ok ok --defect crashes --min-n 50",
    "delfin-agent compare-arms a.jsonl b.jsonl --json",
])
def test_both_commands_read_only_in_any_workspace(tmp_path, cmd):
    import re
    from delfin.agent import api_client as A
    assert any(re.search(p, cmd) for p in A._DEFAULT_BASH_AUTO_ALLOW), cmd


@pytest.mark.parametrize("cmd", [
    "delfin-agent provenance --code src; rm -rf x",
    "delfin-agent compare-arms a b --json && curl http://x",
    "delfin-agent compare-arms a b --stamp-a $(cat key)",
    "delfin-agent provenance --code src > /etc/passwd",
])
def test_nothing_rides_along(tmp_path, cmd):
    import re
    from delfin.agent import api_client as A
    assert not any(re.search(p, cmd) for p in A._DEFAULT_BASH_AUTO_ALLOW), cmd
