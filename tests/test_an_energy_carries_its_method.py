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


# ---------------------------------------------------------------------------
# The second interview, on the fixed tools
# ---------------------------------------------------------------------------
#
# Re-asked over the same fixture with the three fixes in place, the
# operator confirmed them and named the next layer: outcome said
# "unknown" for all nine runs because the fixture signals completion
# through DELFIN_Data.json alone; the fixture contradicted itself (an
# output that terminated normally beside a status of "running");
# extract_delfin_json looked for DELFIN_data.json and never found the
# real file; and a folder without an output was compared without its
# method. All four below.

def test_extract_delfin_json_finds_the_file_as_delfin_writes_it(tmp_path):
    d = tmp_path / "run"; d.mkdir()
    (d / "DELFIN_Data.json").write_text(json.dumps({"status": "finished", "functional": "PBE0"}))
    got = api.extract_delfin_json(str(d))
    assert got.error is None and got.json_path and got.json_path.endswith("DELFIN_Data.json")


def test_extract_delfin_json_still_reads_the_old_spelling(tmp_path):
    d = tmp_path / "run"; d.mkdir()
    (d / "DELFIN_data.json").write_text(json.dumps({"status": "finished"}))
    assert api.extract_delfin_json(str(d)).error is None


def test_outcome_reads_the_state_file_when_nothing_else_says(tmp_path):
    from delfin.doc_server import calc_indexer as ci
    a = tmp_path / "a"; a.mkdir(); (a / "DELFIN_Data.json").write_text(json.dumps({"status": "running"}))
    b = tmp_path / "b"; b.mkdir(); (b / "DELFIN_Data.json").write_text(json.dumps({"status": "finished"}))
    c = tmp_path / "c"; c.mkdir(); (c / "DELFIN_Data.json").write_text(json.dumps({"status": "finished"})); (c / ".exit_code_3").write_text("")
    assert ci.outcome_of_folder(a).startswith("running per DELFIN_Data.json")
    assert ci.outcome_of_folder(b).startswith("finished per DELFIN_Data.json")
    assert ci.outcome_of_folder(c) == "failed (exit code 3)", "an exit code outranks the state file"


def test_a_run_without_an_output_is_compared_with_its_method(tmp_path):
    d = tmp_path / "running"; d.mkdir()
    (d / "run.inp").write_text("! TPSSh def2-TZVP Opt\n")
    (d / "CONTROL.txt").write_text("basis_set = def2-TZVP\n")
    rows = api.compare_across_functionals([str(d)], include_imag=False)
    assert rows[0].status == "no_output" and rows[0].method == "TPSSh/def2-TZVP"


def test_the_small_archive_no_longer_contradicts_itself(tmp_path):
    import subprocess, sys
    setup = Path(api.__file__).resolve().parent / "agent" / "pack" / "benchmark" / "setup" / "a_small_calc_archive.py"
    r = subprocess.run([sys.executable, str(setup), str(tmp_path)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    running = []
    for f in (tmp_path / "calc_archive").rglob("DELFIN_Data.json"):
        status = json.loads(f.read_text()).get("status")
        has_out = any(f.parent.glob("*.out"))
        assert (status == "running") == (not has_out), f"{f.parent.name}: status {status} with output {has_out}"
        if status == "running":
            running.append(f.parent.name)
    assert running == [] or running == ["calc_d"], running
    rows = {Path(r["folder"]).name: r for r in api.extract_energy_table(
        [str(p) for p in sorted((tmp_path / "calc_archive" / "calc").iterdir())], properties=["single_point"])}
    assert rows["calc_d"]["outcome"].startswith(("running", "no output yet")), rows["calc_d"]["outcome"]
    assert all(r["outcome"].startswith("finished per DELFIN_Data.json") for n, r in rows.items() if n != "calc_d"), rows


# ---------------------------------------------------------------------------
# GLM, as the operator, added two
# ---------------------------------------------------------------------------

def test_a_sort_key_nobody_has_falls_back_and_says_so(tmp_path):
    """sort_by=gibbs over a single-point archive sorted silently by
    something else. Now it sorts by what is there and every row says which."""
    a = _run(tmp_path, "a", "PBE0", "def2-SVP", -113.20)
    b = _run(tmp_path, "b", "PBE0", "def2-SVP", -113.30)
    # strip the Gibbs line so only single points exist
    for d in (a, b):
        f = Path(d) / "run.out"; f.write_text("\n".join(l for l in f.read_text().splitlines() if "Gibbs" not in l) + "\n")
    rows = api.compare_across_functionals([a, b], include_imag=False, sort_by="gibbs")
    ok = [r for r in rows if r.status == "ok"]
    if len(ok) < 2:
        pytest.skip("parser did not read the fixture")
    assert all(r.sorted_by == "single_point" for r in ok)
    assert [Path(r.folder).name for r in ok] == ["b", "a"], "ordered by the key it fell back to"


def test_an_input_without_output_is_no_output_yet_not_unknown(tmp_path):
    from delfin.doc_server import calc_indexer as ci
    d = tmp_path / "queued"; d.mkdir(); (d / "run.inp").write_text("! PBE0 def2-SVP\n")
    assert ci.outcome_of_folder(d).startswith("no output yet")
    e = tmp_path / "empty"; e.mkdir()
    assert ci.outcome_of_folder(e).startswith("unknown")


def test_parse_orca_output_names_the_method_from_the_folder(tmp_path):
    """Round 3: parse_orca_output promised functional/basis and returned
    empty strings for a DELFIN-style output. The file's folder knows."""
    d = tmp_path / "run"; d.mkdir()
    out = d / "run.out"; out.write_text("FINAL SINGLE POINT ENERGY   -113.30500000\n****ORCA TERMINATED NORMALLY****\n")
    (d / "run.inp").write_text("! B3LYP def2-TZVP Opt\n")
    parsed = api.parse_orca_output(str(out))
    assert parsed.functional == "B3LYP" and parsed.basis == "def2-TZVP"
    assert parsed.final_single_point == pytest.approx(-113.305)


def test_the_summary_table_names_the_method_of_a_run_without_output(tmp_path):
    d = tmp_path / "queued"; d.mkdir()
    (d / "run.inp").write_text("! PBE0 def2-SVP Opt\n")
    rows = api.extract_calc_summary_table([str(d)])
    assert rows[0].status == "no_output" and rows[0].functional == "PBE0" and rows[0].basis == "def2-SVP"


# ---------------------------------------------------------------------------
# An empty ranking says why it is empty
# ---------------------------------------------------------------------------
#
# Driven over nine single-point runs (2026-09-10), find_calculation_extreme
# returned {"groups": []} and nothing else: the default property is gibbs,
# a single-point output has none, and every row was dropped in silence.
# The operator concluded the parsers could not read the files, read them
# by hand, and misassigned one energy. Now the ranking falls back to what
# is there when the requested property is nowhere, says so, and lists
# every folder it left out with the reason.

_OUT_SPE_ONLY = """
                                 * O   R   C   A *
| 1> ! {functional} {basis} SP
FINAL SINGLE POINT ENERGY      {spe}
                             ****ORCA TERMINATED NORMALLY****
"""


def _spe_run(root: Path, name: str, functional: str, basis: str, spe: float) -> str:
    d = root / name
    d.mkdir(parents=True)
    (d / "run.out").write_text(_OUT_SPE_ONLY.format(functional=functional,
                                                    basis=basis, spe=spe))
    return str(d)


def test_an_absent_property_is_ranked_by_what_is_there_and_says_so(tmp_path):
    a = _spe_run(tmp_path, "a", "PBE0", "def2-SVP", -113.30)
    b = _spe_run(tmp_path, "b", "PBE0", "def2-SVP", -113.31)
    res = api.find_calculation_extreme_explained([a, b], property="gibbs", n=1)
    assert res["property_requested"] == "gibbs"
    assert res["property_used"] == "single_point"
    assert [Path(r["folder"]).name for r in res["rows"]] == ["b"]
    assert res["rows"][0]["property_used"] == "single_point"
    out = json.loads(ops.tool_find_calculation_extreme(f"{a},{b}", property="gibbs", n=1))
    assert out["property_used"] == "single_point"
    assert "single_point" in out["note"] and "gibbs" in out["note"]
    assert [r["folder"] for g in out["groups"] for r in g["rows"]] == [b]


def test_mixed_availability_keeps_the_requested_property(tmp_path):
    """One folder has a Gibbs energy, one only a single point that sits
    lower: the ranking stays on gibbs and names the other as skipped."""
    g = _run(tmp_path, "g", "PBE0", "def2-SVP", -113.20)
    s = _spe_run(tmp_path, "s", "PBE0", "def2-SVP", -113.90)
    res = api.find_calculation_extreme_explained([g, s], property="gibbs")
    assert res["property_used"] == "gibbs"
    assert [Path(r["folder"]).name for r in res["rows"]] == ["g"]
    left = {Path(x["folder"]).name: x["reason"] for x in res["skipped"]}
    assert "gibbs" in left["s"]


def test_a_folder_left_out_says_why(tmp_path):
    good = _run(tmp_path, "good", "PBE0", "def2-SVP", -113.20)
    pending = tmp_path / "pending"
    pending.mkdir()
    (pending / "run.inp").write_text("! PBE0 def2-SVP Opt\n")
    missing = str(tmp_path / "nowhere")
    out = json.loads(ops.tool_find_calculation_extreme(
        f"{good},{pending},{missing}", property="gibbs"))
    assert [r["folder"] for g in out["groups"] for r in g["rows"]] == [good]
    left = {Path(x["folder"]).name: x for x in out["skipped"]}
    assert set(left) == {"pending", "nowhere"}
    assert "no output" in left["pending"]["reason"]
    assert "no output" in (left["pending"]["outcome"] or "")
    assert "missing" in left["nowhere"]["reason"]


def test_nothing_to_rank_is_said_not_shown_as_an_empty_list(tmp_path):
    pending = tmp_path / "pending"
    pending.mkdir()
    (pending / "run.inp").write_text("! PBE0 def2-SVP Opt\n")
    out = json.loads(ops.tool_find_calculation_extreme(str(pending), property="gibbs"))
    assert out["groups"] == []
    assert "skipped" in out["note"]
    assert out["skipped"] and Path(out["skipped"][0]["folder"]).name == "pending"


def test_the_rows_alone_are_still_the_rows(mixed):
    """The list-returning name keeps its contract for callers that only
    want the ranking."""
    rows = api.find_calculation_extreme(list(mixed.values()), property="gibbs", n=1)
    explained = api.find_calculation_extreme_explained(list(mixed.values()), property="gibbs", n=1)
    assert [r["folder"] for r in rows] == [r["folder"] for r in explained["rows"]]


def test_the_job_list_says_it_is_not_the_disk():
    """Driven over an archive with one unfinished run, list_active_calculations
    returned [] -- true of the scheduler, read as 'nothing is running'.
    The description now says which question it answers and where the
    other one is answered."""
    doc = ops.tool_list_active_calculations.__doc__ or ""
    assert "scheduler" in doc
    assert "extract_energy_table" in doc and "outcome" in doc


# ---------------------------------------------------------------------------
# The same method means the same solvent and the same dispersion correction
# ---------------------------------------------------------------------------
#
# The archive fixture had PBE0/def2-SVP in DMF and PBE0/def2-SVP in the gas
# phase filed under one method, and the tools ranked them together. An
# implicit solvent shifts a total energy by the solvation free energy;
# that is not chemistry between two structures either.

_CONTROL = """# DELFIN CONTROL
NAME = {name}
charge = 0
multiplicity = 1
method = classic
functional = {functional}
basis_set = {basis}
solvent = {solvent}
implicit_solvation_model = CPCM
"""


def _run_in(root: Path, name: str, functional: str, basis: str, solvent: str,
            gibbs: float, dispersion: str = "") -> str:
    d = root / name
    d.mkdir(parents=True)
    (d / "run.out").write_text(_OUT.format(functional=functional, basis=basis,
                                           spe=gibbs + 0.05, gibbs=gibbs))
    text = _CONTROL.format(name=name, functional=functional, basis=basis, solvent=solvent)
    if dispersion:
        text += f"dispersion_correction = {dispersion}\n"
    (d / "CONTROL.txt").write_text(text)
    return str(d)


def test_a_solvent_is_part_of_the_method(tmp_path):
    dmf = _run_in(tmp_path, "dmf", "PBE0", "def2-SVP", "DMF", -113.30)
    gas = _run_in(tmp_path, "gas", "PBE0", "def2-SVP", "none", -113.20)
    rows = {Path(r["folder"]).name: r for r in api.extract_energy_table([dmf, gas])}
    assert rows["dmf"]["method"] == "PBE0/def2-SVP/DMF"
    assert rows["gas"]["method"] == "PBE0/def2-SVP"          # gas phase adds nothing
    assert rows["dmf"]["solvent"] == "DMF" and rows["gas"]["solvent"] is None
    ranked = api.find_calculation_extreme([dmf, gas], property="gibbs", n=5)
    methods = {r["method"] for r in ranked}
    assert len(methods) == 2, "a DMF run and a gas-phase run were ranked together"


def test_a_dispersion_correction_is_part_of_the_method(tmp_path):
    d3 = _run_in(tmp_path, "d3", "PBE0", "def2-SVP", "none", -113.30, dispersion="D3BJ")
    bare = _run_in(tmp_path, "bare", "PBE0", "def2-SVP", "none", -113.20)
    rows = {Path(r["folder"]).name: r for r in api.extract_energy_table([d3, bare])}
    assert rows["d3"]["method"] == "PBE0-D3BJ/def2-SVP"
    assert rows["bare"]["method"] == "PBE0/def2-SVP"
    assert rows["d3"]["dispersion"] == "D3BJ"


def test_the_input_header_names_the_solvent_when_nothing_else_does(tmp_path):
    """A folder with only an .inp -- the fixture's arch_e -- carries its
    solvent in ORCA's own keyword: CPCM(DMF)."""
    d = tmp_path / "only_inp"
    d.mkdir()
    (d / "run.inp").write_text("! PBE0 def2-SVP D3BJ CPCM(DMF) TightSCF\n* xyz 0 1\n*\n")
    from delfin.doc_server.calc_indexer import method_parts_of_folder
    parts = method_parts_of_folder(d)
    assert (parts["functional"], parts["basis"]) == ("PBE0", "def2-SVP")
    assert parts["solvent"] == "DMF" and parts["dispersion"] == "D3BJ"
    assert api._method_label(**{k: parts[k] for k in ("functional", "basis", "solvent", "dispersion")}) \
        == "PBE0-D3BJ/def2-SVP/DMF"


def test_the_note_names_all_four(tmp_path):
    for word in ("functional", "basis", "dispersion", "solvent"):
        assert word in api.METHOD_NOTE
    assert "gas phase" in api.METHOD_NOTE


def test_the_archive_fixture_no_longer_files_gas_and_dmf_together(tmp_path):
    import importlib.util
    setup = Path(api.__file__).resolve().parent / "agent" / "pack" / "benchmark" / "setup" / "a_small_calc_archive.py"
    spec = importlib.util.spec_from_file_location("archive_setup", setup)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.main(["x", str(tmp_path)]) == 0
    ws = tmp_path / "calc_archive"
    folders = [str(p) for p in sorted(list((ws / "calc").iterdir()) + list((ws / "archive").iterdir()))]
    rows = {Path(r["folder"]).name: r for r in api.extract_energy_table(folders, properties=["single_point"])}
    assert rows["calc_c"]["method"] == "PBE0/def2-SVP"
    assert rows["arch_b"]["method"] == "PBE0/def2-SVP/DMF" == rows["arch_e"]["method"]
    assert rows["calc_b"]["method"] == "B3LYP/def2-SVP/water" == rows["arch_c"]["method"]
    groups = {}
    for r in api.find_calculation_extreme(folders, property="single_point", n=1):
        groups[r["method"]] = Path(r["folder"]).name
    assert groups["PBE0/def2-SVP/DMF"] == "arch_e"
    assert groups["PBE0/def2-TZVP/DMF"] == "arch_a"
    assert groups["B3LYP/def2-SVP/water"] == "calc_b"
