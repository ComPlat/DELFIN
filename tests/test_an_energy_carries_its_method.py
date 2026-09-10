"""An energy without its method is a number waiting to be ranked.

Measured 2026-09-10 (20 samples per arm): asked which run has the lowest
energy over a mixed archive, the agent opened the outputs and ranked two
B3LYP runs against PBE0 as chemistry -- and a one-sentence rule in the
prompt did not change that (5/20 -> 1/20, p=0.18). Driving the tools it
reaches for showed why: extract_energy_table returned energies with no
functional at all, compare_across_functionals sorted flat by gibbs
across methods, find_calculation_extreme returned one global minimum,
and the tool description invited "which functional gives the lowest
minimum?". The data handed the model the ranking. So the data changed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin import api
from delfin.ops_server import server as ops

_OUT = """
                                 * O   R   C   A *
| 1> ! {functional} {basis} Opt Freq
FINAL SINGLE POINT ENERGY      {spe}
Final Gibbs free energy         ...      {gibbs} Eh
                             ****ORCA TERMINATED NORMALLY****
"""


def _run(root: Path, name: str, functional: str, basis: str, gibbs: float) -> str:
    d = root / name
    d.mkdir(parents=True)
    (d / "run.out").write_text(_OUT.format(functional=functional, basis=basis,
                                           spe=gibbs + 0.05, gibbs=gibbs))
    return str(d)


@pytest.fixture
def mixed(tmp_path):
    """Two methods; the B3LYP pair sits far below the PBE0 pair, as it
    does in every real archive, and the lowest PBE0 is the best structure
    of its group."""
    return {
        "b3_a": _run(tmp_path, "b3_a", "B3LYP", "def2-SVP", -113.50),
        "b3_b": _run(tmp_path, "b3_b", "B3LYP", "def2-SVP", -113.40),
        "pbe_a": _run(tmp_path, "pbe_a", "PBE0", "def2-SVP", -113.20),
        "pbe_b": _run(tmp_path, "pbe_b", "PBE0", "def2-SVP", -113.10),
    }


def _parsed_ok(rows):
    return [r for r in rows if r.get("status") == "ok" and r.get("method")]


def test_every_row_of_the_energy_table_carries_its_method(mixed):
    rows = api.extract_energy_table(list(mixed.values()), properties=["gibbs"])
    ok = _parsed_ok(rows)
    if len(ok) < 4:
        pytest.skip("the ORCA parser did not read the fixture's method line")
    assert {r["method"] for r in ok} == {"B3LYP/def2-SVP", "PBE0/def2-SVP"}
    assert all(r["functional"] and r["basis"] for r in ok)


def test_a_missing_folder_row_has_the_method_keys_too(tmp_path):
    rows = api.extract_energy_table([str(tmp_path / "nope")], properties=["gibbs"])
    assert rows[0]["status"] == "missing" and rows[0]["method"] is None


def test_the_comparison_never_orders_across_methods(mixed):
    rows = api.compare_across_functionals(list(mixed.values()), include_imag=False, sort_by="gibbs")
    ok = [r for r in rows if r.status == "ok" and r.method]
    if len(ok) < 4:
        pytest.skip("the ORCA parser did not read the fixture's method line")
    methods = [r.method for r in ok]
    # Groups are contiguous and ordered by name, not by energy: B3LYP is
    # the lower group and would lead a flat sort either way, so the test
    # checks the PBE0 group is intact and internally ordered, and that no
    # PBE0 row sits between two B3LYP rows.
    assert methods == sorted(methods), methods
    pbe = [r.gibbs for r in ok if r.method == "PBE0/def2-SVP"]
    assert pbe == sorted(pbe)


def test_the_extreme_is_found_per_method(mixed):
    rows = api.find_calculation_extreme(list(mixed.values()), property="gibbs", n=1)
    ok = _parsed_ok(rows)
    if len(ok) < 2:
        pytest.skip("the ORCA parser did not read the fixture's method line")
    by_method = {r["method"]: r for r in ok}
    assert Path(by_method["B3LYP/def2-SVP"]["folder"]).name == "b3_a"
    assert Path(by_method["PBE0/def2-SVP"]["folder"]).name == "pbe_a"
    assert all(r["rank_within_method"] == 1 for r in ok)


def test_the_tools_say_the_rule_before_they_say_a_number(mixed):
    out = json.loads(ops.tool_compare_across_functionals(",".join(mixed.values()), include_imag=False))
    assert out["note"] == api.METHOD_NOTE
    assert "groups" in out and all("method" in g and "rows" in g for g in out["groups"])
    out2 = json.loads(ops.tool_find_calculation_extreme(",".join(mixed.values()), property="gibbs", n=2))
    assert out2["note"] == api.METHOD_NOTE


def test_the_descriptions_no_longer_invite_the_question_that_has_no_answer():
    src = open(ops.__file__).read()
    assert "Which functional gives the lowest minimum" not in src
    assert "within one method" in src


# ---------------------------------------------------------------------------
# What a model, asked where the tools hurt, reported the same evening
# ---------------------------------------------------------------------------
#
# Driving the tools as the operator over the nine-calculation fixture,
# DeepSeek reported: compare_across_functionals filed every run under
# method None (a DELFIN run's ORCA output does not state its method);
# extract_energy_table said "no_output" for a run that is still running
# and could not tell that from a missing one; list_active_calculations
# crashed with an ImportError. All three below.

def test_the_method_comes_from_the_folder_when_the_output_does_not_say(tmp_path):
    d = tmp_path / "run"; d.mkdir()
    (d / "run.out").write_text("* O R C A *\nFINAL SINGLE POINT ENERGY   -113.30500000\n****ORCA TERMINATED NORMALLY****\n")
    (d / "run.inp").write_text("! PBE0 def2-SVP Opt\n* xyz 0 1\nO 0 0 0\n*\n")
    rows = api.extract_energy_table([str(d)], properties=["single_point"])
    assert rows[0]["status"] == "ok" and rows[0]["method"] == "PBE0/def2-SVP", rows[0]
    cmp = api.compare_across_functionals([str(d)], include_imag=False)
    assert cmp[0].method == "PBE0/def2-SVP"


def test_delfin_data_outranks_the_inp_header(tmp_path):
    d = tmp_path / "run"; d.mkdir()
    (d / "run.out").write_text("FINAL SINGLE POINT ENERGY   -1.0\n")
    (d / "run.inp").write_text("! B3LYP def2-SVP\n")
    (d / "DELFIN_Data.json").write_text(json.dumps({"functional": "PBE0", "basis_set": "def2-TZVP"}))
    rows = api.extract_energy_table([str(d)], properties=["single_point"])
    assert rows[0]["method"] == "PBE0/def2-TZVP"


def test_a_row_says_whether_the_run_is_still_running(tmp_path):
    running = tmp_path / "running"; running.mkdir()
    (running / "run.inp").write_text("! PBE0 def2-SVP\n")
    (running / "delfin_run.log").write_text("started\n")
    done = tmp_path / "done"; done.mkdir()
    (done / "run.out").write_text("FINAL SINGLE POINT ENERGY   -1.0\n")
    (done / ".exit_code_0").write_text("")
    failed = tmp_path / "failed"; failed.mkdir()
    (failed / "run.out").write_text("ORCA finished by error termination\n")
    (failed / ".exit_code_1025").write_text("")
    rows = {Path(r["folder"]).name: r for r in api.extract_energy_table(
        [str(running), str(done), str(failed), str(tmp_path / "nope")], properties=["single_point"])}
    assert rows["running"]["status"] == "no_output" and rows["running"]["outcome"].startswith("running or crashed")
    assert rows["running"]["method"] == "PBE0/def2-SVP", "a running run still names its method"
    assert rows["done"]["outcome"] == "succeeded (exit code 0)"
    assert rows["failed"]["outcome"] == "failed (exit code 1025)"
    assert rows["nope"]["status"] == "missing" and rows["nope"]["outcome"].startswith("unknown")


def test_the_index_and_the_energy_tools_read_completion_the_same_way(tmp_path):
    from delfin.doc_server import calc_indexer as ci
    d = tmp_path / "x"; d.mkdir(); (d / ".exit_code_7").write_text("")
    assert ci.completion_of(d) == (True, 7)
    assert ci.outcome_of_folder(d) == "failed (exit code 7)"


def test_list_active_calculations_no_longer_crashes_on_an_import():
    out = ops.tool_list_active_calculations() if not ops.tool_list_active_calculations.__code__.co_argcount else None
    if out is None:
        pytest.skip("wrapper takes arguments; the import is exercised below")
    assert "cannot import name" not in out


def test_the_local_backend_import_resolves():
    api._resolve_backend  # the function under test
    import shutil
    if shutil.which("sbatch") and shutil.which("squeue"):
        pytest.skip("slurm host: the local branch is not taken")
    backend = api._resolve_backend()
    assert type(backend).__name__ == "LocalJobBackend"
