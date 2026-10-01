"""delfin cluster: sharding, the Slurm script, the tool adapters, resume and the merged classes.

Everything here is synthetic: invented IDs, textbook SMILES, fake builders.  No structure
database content is needed or allowed in the repository.
"""
from __future__ import annotations

import json
import sys
import textwrap
from pathlib import Path

import pytest

from delfin.cluster_bench import prepare as cp
from delfin.cluster_bench import provenance as cprov
from delfin.cluster_bench import report as crep
from delfin.cluster_bench import runner as crun
from delfin.cluster_bench.slurm_script import cbatch_render_sbatch

WORKERS = Path(cprov.WORKERS_DIR)

CISPLATIN = "[Cl][Pt-2]([Cl])([NH3+])[NH3+]"
HEXAAQUA = "[OH2+][Fe-3]([OH2+])([OH2+])([OH2+])([OH2+])[OH2+]"
COCL4N2 = "[Cl][Co-4]([Cl])([Cl])([Cl])([NH3+])[NH3+]"
CUBIPY = "C1=CC=[N+]2C(=C1)C1=CC=CC=[N+]1[Cu-2]2([Cl])[Cl]"
NO_METAL = "CC(C)=O"


def _write(path: Path, text: str) -> Path:
    path.write_text(textwrap.dedent(text).lstrip())
    return path


# ---------------------------------------------------------------- input and sharding
def test_a_bad_input_line_writes_nothing(tmp_path):
    p = _write(tmp_path / "in.txt", "A1;CCO\nA2 CCO\nA1;CC\nA3;C C\n")
    with pytest.raises(SystemExit) as e:
        cp.cbatch_parse_list(p)
    assert "3 invalid" in str(e.value)


def test_identical_smiles_stay_in_one_shard_and_the_cut_is_deterministic():
    rows = [(f"S{i}", f"C{'C' * (i % 7)}O") for i in range(40)]
    a = cp.cbatch_split_grouped(rows, 5, "salt")
    b = cp.cbatch_split_grouped(list(rows), 5, "salt")
    assert a == b
    assert sorted(r for sh in a for r, _ in sh) == sorted(r for r, _ in rows)
    where = {}
    for k, sh in enumerate(a):
        for _, smi in sh:
            where.setdefault(smi, set()).add(k)
    assert all(len(v) == 1 for v in where.values())
    assert all(len(sh) >= 5 for sh in a[:-1])


def test_the_selection_keeps_its_order_and_refuses_unknown_ids(tmp_path):
    rows = [("A", "C"), ("B", "CC"), ("C", "CCC")]
    sel = _write(tmp_path / "sel.txt", "C\nA;whatever\n")
    assert cp.cbatch_apply_selection(rows, cp.cbatch_read_selection(sel)) == [("C", "CCC"), ("A", "C")]
    with pytest.raises(SystemExit):
        cp.cbatch_apply_selection(rows, ["Z"])
    assert cp.cbatch_split_ordered(rows, 2) == [[("A", "C"), ("B", "CC")], [("C", "CCC")]]


def test_the_limit_is_scaled_in_exact_decimal():
    assert cprov.cbatch_effective_timeout(21600, "1.3") == 28080
    assert cprov.cbatch_effective_timeout(21600, 1.0) == 21600
    assert cprov.cbatch_effective_timeout(21600, "1.00001") == 21601


def test_the_repeat_subset_takes_one_id_per_smiles():
    rows = [("A", "C"), ("B", "C"), ("C", "CC"), ("D", "CCC")]
    sub = cp.cbatch_repeat_subset(rows, 3, "s")
    assert len(sub) == 3 and len({s for _, s in sub}) == 3


# ---------------------------------------------------------------- specs and adapters
def test_specs_and_the_tool_adapters_on_textbook_complexes():
    sys.path.insert(0, str(WORKERS))
    try:
        import bench_convert as bc
        import mace_convert as mc
    finally:
        sys.path.remove(str(WORKERS))
    specs = {s: cp.cbatch_spec_of(("X", s)) for s in (CISPLATIN, HEXAAQUA, COCL4N2, CUBIPY, NO_METAL)}
    assert specs[CISPLATIN]["status"] == "ok" and specs[CISPLATIN]["cn"] == 4
    assert specs[CISPLATIN]["metal_ox"] == 2
    assert specs[HEXAAQUA]["metal_ox"] == 3 and specs[HEXAAQUA]["cn"] == 6
    assert specs[NO_METAL]["status"] == "split_failed" and specs[NO_METAL]["bail"] == "no_metal"
    assert json.dumps(specs[CUBIPY]).startswith('{"refcode": "X"')  # the workers find specs by this prefix
    assert bc.bench_common_check(specs[NO_METAL]) == "spec_not_ok"
    assert bc.bench_common_check(specs[COCL4N2]) is None
    job, err = mc.mace_convert(specs[HEXAAQUA], "paper")
    assert err is None and job["geoms"] == ["OH"] and job["CA"] == "[Fe+3]"
    job, err = mc.mace_convert(specs[CUBIPY], "extended")
    assert err is None and job["geoms"] == ["SP", "TET"] and job["info"]["n_sites"] == 4
    assert mc.mace_convert(specs[NO_METAL])[1] == "not_expressible:spec_no_metal"


# ---------------------------------------------------------------- slurm script
def test_the_array_script_asks_for_a_full_node_for_72_hours(tmp_path):
    man = {"tool": "architector", "label": "arch_x",
           "settings": {"timeout_base_s": 21600, "speed_factor": "1.0", "workers": 48, "threads": 1},
           "sets": {"main": {"n_shards": 85}}}
    text = cbatch_render_sbatch(tmp_path, man, speed_factor="1.3", throttle=20,
                                setup=["module load chem/xyz"], python="/ws/env/bin/python")
    assert "#SBATCH --cpus-per-task=48" in text
    assert "#SBATCH --time=72:00:00" in text
    assert "#SBATCH --array=0-84%20" in text
    assert "limit 28080 s" in text
    assert "-m delfin.cluster_bench run-shard" in text and '--shard "$SLURM_ARRAY_TASK_ID"' in text
    assert "--speed-factor 1.3" in text and "module load chem/xyz" in text
    assert crep.cbatch_array_spec([0, 1, 2, 5, 7, 8]) == "0-2,5,7-8"


# ---------------------------------------------------------------- MANTA: resume, duplicates, collect
def test_a_manta_shard_resumes_serves_duplicates_and_collects(tmp_path, monkeypatch):
    inp = _write(tmp_path / "in.txt", f"""
        M1;{CISPLATIN}
        M2;{CISPLATIN}
        M3;{HEXAAQUA}
        T1;{COCL4N2}
        E1;{CUBIPY}
        """)
    man = cp.cbatch_prepare(tool="manta", input_list=inp, run_dir=tmp_path / "run", size=50,
                            repeat=2, salt="t")
    assert man["sets"]["main"]["n_shards"] == 1 and man["settings"]["threads"] == 7
    calls = []

    def fake_child(cmd, env, timeout):
        rid, out = cmd[3], Path(cmd[5])
        calls.append(rid)
        assert env["DELFIN_DETERMINISTIC"] == "1" and env["PYTHONHASHSEED"] == "0"
        if rid.startswith("T"):
            return None, "", "", True
        if rid.startswith("E"):
            return 0, json.dumps({"rid": rid, "status": "empty"}) + "\n", "", False
        out.mkdir(parents=True, exist_ok=True)
        (out / f"{rid}.xyz").write_text(f"1\n{rid} frame0 SP-4\nPt 0.0 0.0 0.0\n")
        return 0, json.dumps({"rid": rid, "status": "ok", "niso": 1}) + "\n", "", False

    monkeypatch.setattr(crun, "cbatch_run_child", fake_child)
    chunk = cp.cbatch_chunk_dir(tmp_path / "run", "main", "main", 0)
    chunk.mkdir(parents=True)
    (chunk / "status.jsonl").write_text(json.dumps({"rid": "M3", "status": "ok", "time_s": 1}) + "\n")
    (chunk / "archive").mkdir()
    (chunk / "archive" / "M3.xyz").write_text("1\nM3 frame0 OC-6\nFe 0 0 0\n")
    rc = crun.cbatch_run_shard(tmp_path / "run", 0, log=lambda m: None)
    assert rc == 0
    assert sorted(calls) == ["E1", "M1", "T1"]          # M3 resumed, M2 served from M1
    assert (chunk / "archive" / "M2.xyz").read_text() == "1\nM2 frame0 SP-4\nPt 0.0 0.0 0.0\n"
    st = crep.cbatch_status(tmp_path / "run")
    assert st["n_done"] == 5 and st["shards_done"] == 1
    s = crep.cbatch_collect(tmp_path / "run")
    assert s["n_problems"] == 0, s["problems"]
    assert s["by_class"] == {"ok": 3, "timeout": 1, "empty": 1}
    assert s["n_xyz"] == 3
    with pytest.raises(SystemExit):
        crep.cbatch_collect(tmp_path / "run")              # never overwritten
    assert crun.cbatch_run_shard(tmp_path / "run", 0, set_name="repeat", log=lambda m: None) == 0
    rs = crep.cbatch_repeat_stats(tmp_path / "run")
    assert rs["identical"] + rs["different"] + rs["neither"] + rs["only_first"] + rs["only_second"] == 2


def test_a_changed_shard_is_refused(tmp_path):
    inp = _write(tmp_path / "in.txt", f"M1;{CISPLATIN}\n")
    cp.cbatch_prepare(tool="manta", input_list=inp, run_dir=tmp_path / "run")
    shard = tmp_path / "run" / "shards" / "main" / "shard_0000.txt"
    shard.write_text(shard.read_text() + f"M9;{HEXAAQUA}\n")
    assert crun.cbatch_run_shard(tmp_path / "run", 0, log=lambda m: None) == 2


# ---------------------------------------------------------------- external builder: classes
FAKE_WORKER = '''
import json, os, sys, time
tool, ref, specs, archive, work = sys.argv[1:6]
os.makedirs(os.path.join(archive, "_meta"), exist_ok=True)
if ref.startswith("SLOW"):
    time.sleep(30)
st, n = {"OK": ("ok", 1), "NE": ("convert:spec_not_ok", 0), "FAIL": ("convert:ms_cn_unsupported", 0),
         "EMPTY": ("no_structure", 0)}[ref.rstrip("0123456789")]
if n:
    open(os.path.join(archive, ref + ".xyz"), "w").write("1\\n%s frame0 oct tool=%s E=na\\nFe 0 0 0\\n" % (ref, tool))
json.dump({"refcode": ref, "tool": tool, "status": st, "n_frames": n},
          open(os.path.join(archive, "_meta", ref + ".json"), "w"))
'''


def test_an_external_builder_shard_gives_every_class(tmp_path, monkeypatch):
    wd = tmp_path / "workers"
    wd.mkdir()
    (wd / "bench_worker.py").write_text(FAKE_WORKER)
    (wd / "env_probe.py").write_text((WORKERS / "env_probe.py").read_text())
    monkeypatch.setattr(cprov, "WORKERS_DIR", wd)
    monkeypatch.setattr(crun, "WORKERS_DIR", wd)
    ids = ["OK1", "OK2", "NE1", "FAIL1", "EMPTY1", "SLOW1"]
    inp = _write(tmp_path / "in.txt", "".join(f"{i};{CISPLATIN}\n" for i in ids))
    specs = tmp_path / "specs.jsonl"
    specs.write_text("".join(json.dumps({"refcode": i, "status": "split_failed" if i == "NE1" else "ok",
                                         "bail": "radical_in_ligand", "cn": 4}) + "\n" for i in ids))
    cp.cbatch_prepare(tool="molsimplify", input_list=inp, run_dir=tmp_path / "run", specs_file=specs,
                      tool_python=sys.executable, timeout_base=2, size=4, check_tool=False)
    assert crun.cbatch_run_shard(tmp_path / "run", 0, log=lambda m: None) == 0
    assert crun.cbatch_run_shard(tmp_path / "run", 1, log=lambda m: None) == 0
    s = crep.cbatch_collect(tmp_path / "run")
    assert s["n_problems"] == 0, s["problems"]
    assert s["by_class"] == {"ok": 2, "not_expressible": 1, "fail": 1, "empty": 1, "timeout": 1}
    assert s["n_expressible"] == 5 and s["coverage_of_expressible"] == 0.4
    assert s["not_expressible_reasons"] == {"radical_in_ligand": 1}
    meta = json.loads((Path(tmp_path / "run" / "collected" / s["label"] / f"archive_{s['label']}"
                            / "_meta" / "SLOW1.json")).read_text())
    assert meta["class"] == "timeout" and meta["timeout_s"] == 2


def test_the_mace_status_marks_not_expressible_itself():
    assert crep.cbatch_class({"status": "convert:not_expressible:sites_5_no_paper_geometry",
                              "n_frames": 0}, {"status": "ok"}) == "not_expressible"
    assert crep.cbatch_class({"status": "convert:arch_cn_unsupported", "n_frames": 0},
                             {"status": "ok"}) == "fail"
    assert crep.cbatch_class({"status": "ok", "n_frames": 3}) == "ok"
    assert crep.cbatch_class({"status": "crash(rc=-9)", "n_frames": 0}) == "empty"


def test_the_code_is_compared_by_content_and_by_commit_only_where_both_know_it():
    a = {"code": {"commit": "abc", "code_sha256": "1", "dirty": False}, "tool_env": {"python": "3.7"}}
    assert cprov.cbatch_provenance_mismatch(a, a) == []
    no_git = {"code": {"commit": None, "code_sha256": "1", "dirty": None}, "tool_env": {"python": "3.7"}}
    assert cprov.cbatch_provenance_mismatch(a, no_git) == []
    other = {"code": {"commit": "abc", "code_sha256": "2"}, "tool_env": {"python": "3.11"}}
    assert cprov.cbatch_provenance_mismatch(a, other) == ["code", "tool_env"]


# ---------------------------------------------------------------- the tool's interpreter


def test_the_tool_interpreter_defaults_to_the_one_delfins_builders_use(tmp_path, monkeypatch):
    fake = tmp_path / "py"
    fake.write_text("#!/bin/sh\n")
    fake.chmod(0o755)
    monkeypatch.setenv("DELFIN_MOLSIMPLIFY_PYTHON", str(fake))
    assert cp.cbatch_default_tool_python("molsimplify") == str(fake)
    monkeypatch.delenv("DELFIN_MACE_PYTHON", raising=False)
    monkeypatch.setenv("DELFIN_AI_TOOLS_ROOT", str(tmp_path / "no_tools_here"))
    with pytest.raises(SystemExit, match="epic-mace"):
        cp.cbatch_default_tool_python("mace")


def test_an_environment_without_the_tool_is_refused_before_anything_is_written():
    assert cp.cbatch_tool_missing("architector", {"packages": {"architector": "0.0.10"}}) is None
    assert "molSimplify" in cp.cbatch_tool_missing("molsimplify", {"packages": {}, "executable": "x"})
    assert cp.cbatch_tool_missing("mace", {"mace_files_sha256": "ab"}) is None
    assert "epic-MACE" in cp.cbatch_tool_missing("mace", {"packages": {}})


def test_the_batch_code_and_its_guide_name_no_private_place():
    root = Path(cp.__file__).resolve().parent
    texts = {p: p.read_text() for p in root.rglob("*.py")}
    guide = root.parent.parent / "docs" / "CONSTRUCTION_BATCH.md"
    texts[guide] = guide.read_text()
    for path, text in texts.items():
        low = text.lower()
        for word in ("agent_workspace", "/home/", "weddell", "delfin-backup", "heldout",
                     "batch_v2", "ccdc", "csd "):
            assert word not in low, f"{path}: {word}"
