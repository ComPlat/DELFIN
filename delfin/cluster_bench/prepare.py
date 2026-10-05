"""Cut an ``ID;SMILES`` list into shards and write the run directory with its manifest.

Run directory (never overwritten; ``prepare`` refuses an existing one)::

    RUN/manifest.json            tool, settings, input/selection sha256, provenance, shard table
    RUN/shards/<set>/shard_NNNN.txt   "ID;SMILES" lines          (set = main | repeat)
    RUN/shards/<set>/specs_NNNN.jsonl tool-neutral specs (external builders only)
    RUN/shards/<set>/SHA256SUMS
    RUN/out/<set>_<run>/chunk_NNNN/   written by run-shard
    RUN/slurm/, RUN/logs/             written by the slurm step

Sharding:
* ``manta``: identical SMILES stay in ONE shard (the runner builds a SMILES once and serves the
  other IDs from that build, relabelled); groups are ordered by sha256(salt + SMILES), because an
  input sorted by metal or family would otherwise put all heavy systems into the same shards;
  a shard is closed as soon as it holds >= size IDs.
* external builders: the order of the list (or of the selection file) is kept -- a stratified,
  interleaved selection stays stratified in every prefix -- and cut into blocks of ``size``.

Specs (external builders): metal, oxidation state, CN and free ligands with donor indices, cut by
``delfin.common.external_builders.split_complex_smiles`` -- the same cut DELFIN's own external
builder module uses.  The tools' input limits are applied later, by the tool adapters in their
own environment (``tool_workers/``).
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from delfin.cluster_bench.provenance import (MODES, TOOLS, cbatch_effective_timeout,
                                             cbatch_provenance, cbatch_sha256_file)

ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+-]*$")
FORMAT = "delfin-cluster-batch/1"

#: split_complex_smiles messages -> category (the reason a system is not expressible as
#: metal + free ligands + oxidation state).
SPEC_BAIL = (("could not be parsed", "unparseable"), ("contains no metal", "no_metal"),
             ("metal atoms; only mononuclear", "multinuclear"),
             ("not bound to the metal", "spectator_fragment"),
             ("could not be sanitised", "ligand_sanitize_failed"),
             ("hydrogen is bound", "donor_is_hydrogen"), ("round-trip", "ligand_roundtrip_failed"),
             ("left as a radical", "radical_in_ligand"), ("outside 0..8", "ox_out_of_range"),
             ("has no ligands", "no_ligands"))


def cbatch_check_line(ln: str):
    """One non-blank list line -> (id, smiles), or the reason it is invalid (a str)."""
    sep = ";" if ";" in ln else ("|" if "|" in ln else None)
    if sep is None or ln.count(sep) != 1:
        return "need exactly one ';' (or '|') between ID and SMILES"
    rid, smi = (x.strip() for x in ln.split(sep))
    if not ID_RE.match(rid):
        return f"ID {rid!r} (letters, digits, _ . + - only)"
    if not smi or any(c.isspace() for c in smi):
        return "empty SMILES or whitespace inside it"
    return rid, smi


def cbatch_parse_list(path) -> list:
    """``ID;SMILES`` (or ``ID|SMILES``) lines -> [(id, smiles)].  Any bad line aborts."""
    rows, errors, seen = [], [], set()
    if not Path(path).is_file():
        raise SystemExit(f"{path}: no such file")
    text = Path(path).read_text(encoding="utf-8")
    for i, ln in enumerate(text.split("\n"), 1):
        ln = ln.rstrip("\r")
        if not ln.strip():
            continue
        row = cbatch_check_line(ln)
        if isinstance(row, str):
            errors.append(f"line {i}: {row}")
            continue
        rid, smi = row
        if rid in seen:
            errors.append(f"line {i}: ID {rid} twice")
            continue
        seen.add(rid)
        rows.append((rid, smi))
    if errors:
        head = "\n  ".join(errors[:20])
        raise SystemExit(f"{len(errors)} invalid line(s) in {path} -- nothing written:\n  {head}")
    if not rows:
        raise SystemExit(f"{path}: no systems")
    return rows


def cbatch_read_selection(path) -> list:
    """IDs, one per line (anything after ';', '|' or whitespace is ignored), order kept."""
    out, seen = [], set()
    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        rid = re.split(r"[;|\s]", ln.strip(), maxsplit=1)[0] if ln.strip() else ""
        if rid and rid not in seen:
            seen.add(rid)
            out.append(rid)
    return out


def cbatch_apply_selection(rows, selection) -> list:
    by_id = dict(rows)
    missing = [r for r in selection if r not in by_id]
    if missing:
        raise SystemExit(f"{len(missing)} selected ID(s) not in the input list, e.g. {missing[:5]}")
    return [(r, by_id[r]) for r in selection]


def cbatch_split_grouped(rows, size, salt) -> list:
    """MANTA shards: SMILES groups never split, groups in sha256(salt + SMILES) order."""
    groups: dict = {}
    for rid, smi in rows:
        groups.setdefault(smi, []).append(rid)
    order = sorted(groups, key=lambda s: hashlib.sha256((salt + s).encode()).hexdigest())
    shards, cur = [], []
    for smi in order:
        cur.extend((r, smi) for r in groups[smi])
        if len(cur) >= size:
            shards.append(cur)
            cur = []
    if cur:
        shards.append(cur)
    return shards


def cbatch_split_ordered(rows, size) -> list:
    return [rows[k:k + size] for k in range(0, len(rows), size)]


def cbatch_repeat_subset(rows, n, salt) -> list:
    """n systems for the second build, deterministic: smallest sha256(salt + 'repeat' + ID),
    one per distinct SMILES (a duplicate would be served from the first build, not rebuilt)."""
    if n <= 0:
        return []
    firsts, seen = [], set()
    for rid, smi in rows:
        if smi not in seen:
            seen.add(smi)
            firsts.append((rid, smi))
    firsts.sort(key=lambda t: hashlib.sha256((salt + "repeat" + t[0]).encode()).hexdigest())
    chosen = {r for r, _ in firsts[:n]}
    return [(r, s) for r, s in rows if r in chosen]


def _git_blob_sha1(path) -> str:
    data = Path(path).read_bytes()
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def cbatch_spec_of(item) -> dict:
    """Tool-neutral spec of one system (the schema the tool adapters read)."""
    rid, smi = item
    from rdkit import Chem, RDLogger

    from delfin.common import external_builders as eb

    RDLogger.DisableLog("rdApp.*")
    rec = {"refcode": rid, "smiles": smi,
           "converter": "external_builders@" + _git_blob_sha1(eb.__file__)[:12]}
    try:
        sp = eb.split_complex_smiles(smi)
    except eb.BuildError as e:
        msg = str(e)
        rec.update(status="split_failed", bail_msg=msg[:200],
                   bail=next((c for k, c in SPEC_BAIL if k in msg), "other"))
        return rec
    except Exception as e:  # noqa: BLE001
        rec.update(status="split_failed", bail="exception:" + type(e).__name__,
                   bail_msg=str(e)[:200])
        return rec
    ligs, hapto = [], False
    for lig in sp["ligands"]:
        chk = Chem.MolFromSmiles(lig["smiles"])
        if chk is None:
            n_heavy, de = 0, []
        else:
            n_heavy = chk.GetNumAtoms()
            de = [chk.GetAtomWithIdx(c).GetSymbol() for c in lig["coordList"]]
            cl = set(lig["coordList"])
            for c in cl:
                for nb in chk.GetAtomWithIdx(c).GetNeighbors():
                    if nb.GetIdx() in cl:
                        hapto = True
        if lig["smiles"] == "[H-]":
            n_heavy, de = 0, ["H"]
        ligs.append(dict(lig, n_heavy=n_heavy, donor_elems=de))
    rec.update(status="ok", metal=sp["metal"], cn=sp["cn"], metal_ox=sp["metal_ox"],
               total_charge=sp["total_charge"], hapto=hapto, ligands=ligs,
               denticities=sorted((lig["denticity"] for lig in ligs), reverse=True))
    rec["max_dent"] = rec["denticities"][0]
    return rec


def cbatch_specs_for(rows, workers, specs_file=None) -> dict:
    """{id: spec JSON line}.  From a given specs file (its lines kept verbatim) or computed."""
    want = {r for r, _ in rows}
    if specs_file:
        out = {}
        for ln in open(specs_file, encoding="utf-8"):
            if ln.strip():
                rid = json.loads(ln)["refcode"]
                if rid in want:
                    out[rid] = ln.rstrip("\n") + "\n"
        missing = [r for r, _ in rows if r not in out]
        if missing:
            raise SystemExit(f"{len(missing)} ID(s) without a spec in {specs_file}, e.g. {missing[:5]}")
        return out
    out = {}
    if workers <= 1:
        recs = map(cbatch_spec_of, rows)
        for rec in recs:
            out[rec["refcode"]] = json.dumps(rec) + "\n"
        return out
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for rec in ex.map(cbatch_spec_of, rows, chunksize=32):
            out[rec["refcode"]] = json.dumps(rec) + "\n"
    return out


def _write_set(root: Path, shards, specs) -> dict:
    root.mkdir(parents=True)
    table, sums = [], []
    for k, sh in enumerate(shards):
        body = "".join(f"{r};{s}\n" for r, s in sh).encode()
        sp = root / f"shard_{k:04d}.txt"
        sp.write_bytes(body)
        ent = {"shard": k, "file": sp.name, "n": len(sh), "n_unique_smiles": len({s for _, s in sh}),
               "sha256": hashlib.sha256(body).hexdigest()}
        sums.append(f"{ent['sha256']}  {sp.name}")
        if specs is not None:
            jb = "".join(specs[r] for r, _ in sh).encode()
            jp = root / f"specs_{k:04d}.jsonl"
            jp.write_bytes(jb)
            ent["specs_sha256"] = hashlib.sha256(jb).hexdigest()
            sums.append(f"{ent['specs_sha256']}  {jp.name}")
        table.append(ent)
    (root / "SHA256SUMS").write_text("\n".join(sums) + "\n")
    return {"n_shards": len(shards), "n_systems": sum(len(s) for s in shards), "shards": table}


#: Distribution (or, for epic-MACE, the probed file hash) that shows the tool is installed.
_TOOL_PACKAGE = {"architector": "architector", "molsimplify": "molSimplify"}


def cbatch_default_tool_python(tool) -> str:
    """The interpreter an external builder runs in when ``--tool-python`` is not given: the
    one DELFIN's own builders use (``DELFIN_<TOOL>_PYTHON``, else the environment DELFIN's
    installer built for the tool, else this interpreter)."""
    from delfin.common import external_builders as eb

    py = eb.tool_python(tool)
    if tool == "mace" and py == sys.executable:
        raise SystemExit("epic-MACE runs in an environment of its own (Python 3.7): install it "
                         "with  python -m delfin.installer --install epic-mace  or pass "
                         "--tool-python (or set DELFIN_MACE_PYTHON)")
    return py


def cbatch_tool_missing(tool, tool_env: dict):
    """A message when the probed tool environment does not have the tool, else None."""
    if tool == "mace":
        if tool_env.get("mace_files_sha256"):
            return None
    elif _TOOL_PACKAGE[tool] in tool_env.get("packages", {}):
        return None
    from delfin.common import external_builders as eb

    return eb.not_installed_message(tool, tool_env.get("executable", "?"))


def cbatch_prepare(*, tool, input_list, run_dir, selection=None, specs_file=None, size=None,
                   salt="delfin-cluster-v1", mode=None, tool_python=None, timeout_base=21600,
                   speed_factor=1.0, workers=None, threads=None, repeat=0, spec_workers=8,
                   label=None, check_tool=True) -> dict:
    if tool not in TOOLS:
        raise SystemExit(f"unknown tool {tool!r}; one of {sorted(TOOLS)}")
    t = TOOLS[tool]
    mode = mode or t["mode"]
    if mode not in MODES[tool]:
        raise SystemExit(f"mode {mode!r} not valid for {tool}: {MODES[tool]}")
    if t["needs_specs"] and not tool_python:
        tool_python = cbatch_default_tool_python(tool)
    tool_python = str(Path(tool_python).absolute()) if tool_python else sys.executable
    run = Path(run_dir)
    if run.exists():
        raise SystemExit(f"{run} exists -- a run directory is never overwritten; choose a new one")
    rows = cbatch_parse_list(input_list)
    n_input = len(rows)
    if selection:
        rows = cbatch_apply_selection(rows, cbatch_read_selection(selection))
    size = int(size or t["shard_size"])
    if tool == "manta":
        shards = cbatch_split_grouped(rows, size, salt)
    else:
        shards = cbatch_split_ordered(rows, size)
    rep_rows = cbatch_repeat_subset(rows, int(repeat), salt)
    specs = None
    if t["needs_specs"]:
        specs = cbatch_specs_for(rows, spec_workers, specs_file)
    prov = cbatch_provenance(tool, tool_python, mode)
    for key in ("delfin_env", "tool_env"):
        if "error" in prov.get(key, {}):
            raise SystemExit(f"{key}: the interpreter could not be probed: {prov[key]['error']}")
    if t["needs_specs"] and check_tool:
        missing = cbatch_tool_missing(tool, prov.get("tool_env", {}))
        if missing:
            raise SystemExit(missing)
    run.mkdir(parents=True)
    sets = {"main": _write_set(run / "shards" / "main", shards, specs)}
    if rep_rows:
        sets["repeat"] = _write_set(run / "shards" / "repeat", cbatch_split_ordered(rep_rows, size)
                                    if tool != "manta" else cbatch_split_grouped(rep_rows, size, salt),
                                    specs)
    for d in ("out", "slurm", "logs"):
        (run / d).mkdir()
    man = {
        "format": FORMAT, "tool": tool, "label": label or f"{tool}_{run.name}",
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "input": {"path": str(Path(input_list).resolve()), "sha256": cbatch_sha256_file(input_list),
                  "n_lines": n_input},
        "selection": ({"path": str(Path(selection).resolve()), "sha256": cbatch_sha256_file(selection)}
                      if selection else None),
        "specs_file": ({"path": str(Path(specs_file).resolve()), "sha256": cbatch_sha256_file(specs_file)}
                       if specs_file else None),
        "n_systems": len(rows), "n_unique_smiles": len({s for _, s in rows}),
        "settings": {"mode": mode, "tool_python": tool_python, "shard_size": size, "salt": salt,
                     "workers": int(workers or t["workers"]), "threads": int(threads or t["threads"]),
                     "timeout_base_s": int(timeout_base), "speed_factor": str(speed_factor),
                     "timeout_s": cbatch_effective_timeout(timeout_base, speed_factor),
                     "repeat": int(repeat)},
        "provenance": prov,
        "sets": sets,
    }
    (run / "manifest.json").write_text(json.dumps(man, indent=1))
    return man


def cbatch_load_manifest(run_dir) -> dict:
    p = Path(run_dir) / "manifest.json"
    if not p.exists():
        raise SystemExit(f"{p} not found -- not a run directory (delfin cluster prepare)")
    man = json.loads(p.read_text())
    if man.get("format") != FORMAT:
        raise SystemExit(f"{p}: unknown format {man.get('format')!r}")
    return man


def cbatch_chunk_dir(run_dir, set_name, run_name, shard) -> Path:
    return Path(run_dir) / "out" / f"{set_name}_{run_name}" / f"chunk_{int(shard):04d}"


def cbatch_shard_rows(run_dir, set_name, shard) -> list:
    p = Path(run_dir) / "shards" / set_name / f"shard_{int(shard):04d}.txt"
    return [tuple(ln.split(";", 1)) for ln in p.read_text().splitlines() if ln.strip()]

