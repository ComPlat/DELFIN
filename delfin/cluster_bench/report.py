"""Progress of a run (``status``), the merged result (``collect``) and the repeat comparison.

Classes in ``build_<label>.json`` and in every ``_meta`` record:

    ok               at least one frame
    timeout          killed at the per-system wall limit
    empty            the tool ran and returned no frame (no structure, exception, crash)
    fail             the tool's own input limits refuse the system (CN or denticity without a
                     geometry, metal unsupported, ...) -- a failure of the tool
    not_expressible  the input cannot be written in the tool's input language (the system cannot
                     be cut into metal + free ligands + oxidation state 0..8: radical ligand, H as
                     donor, ...; for epic-MACE also a site count without a MACE geometry).
                     NOT a failure of the tool and not in its denominator.

MANTA reads SMILES directly, so its classes are ok / empty / timeout / fail.
Coverage in the summary: ok / (n - not_expressible), and ok / n for reference.
"""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

from delfin.cluster_bench.prepare import cbatch_chunk_dir, cbatch_load_manifest, cbatch_shard_rows
from delfin.cluster_bench.runner import cbatch_read_final_metas, cbatch_read_status

NOT_EXPRESSIBLE = "not_expressible"
CLASSES = ("ok", "empty", "timeout", "fail", NOT_EXPRESSIBLE)


def cbatch_class(meta: dict, spec: dict | None = None) -> str:
    """Class of one external-builder record (see the module doc)."""
    st = str(meta.get("status", ""))
    if meta.get("n_frames", 0) > 0:
        return "ok"
    if st == "timeout":
        return "timeout"
    if st.startswith("convert:not_expressible:"):
        return NOT_EXPRESSIBLE
    if st.startswith("convert:"):
        if spec is not None and spec.get("status") != "ok":
            return NOT_EXPRESSIBLE
        return "fail"
    return "empty"


def _load_specs(path: Path) -> dict:
    out = {}
    if path.exists():
        for ln in path.read_text().splitlines():
            if ln.strip():
                r = json.loads(ln)
                out[r["refcode"]] = r
    return out


def cbatch_chunk_progress(run_dir, man, set_name, run_name, shard) -> dict:
    rows = cbatch_shard_rows(run_dir, set_name, shard)
    ids = [r for r, _ in rows]
    chunk = cbatch_chunk_dir(run_dir, set_name, run_name, shard)
    by = {}
    if man["tool"] == "manta":
        st = cbatch_read_status(chunk / "status.jsonl")
        for r in ids:
            if r in st:
                by[st[r]["status"]] = by.get(st[r]["status"], 0) + 1
    else:
        specs = _load_specs(Path(run_dir) / "shards" / set_name / f"specs_{int(shard):04d}.jsonl")
        for r, m in cbatch_read_final_metas(chunk / "archive", ids).items():
            c = cbatch_class(m, specs.get(r))
            by[c] = by.get(c, 0) + 1
    n_done = sum(by.values())
    state = ("done" if (chunk / "DONE.json").exists() and n_done == len(ids)
             else "running/partial" if n_done else ("started" if chunk.exists() else "pending"))
    return {"shard": int(shard), "n": len(ids), "done": n_done, "by_class": by, "state": state}


def cbatch_status(run_dir, set_name="main", run_name="main") -> dict:
    man = cbatch_load_manifest(run_dir)
    if set_name not in man["sets"]:
        raise SystemExit(f"shard set {set_name!r} not in this run ({sorted(man['sets'])})")
    rows = [cbatch_chunk_progress(run_dir, man, set_name, run_name, e["shard"])
            for e in man["sets"][set_name]["shards"]]
    tot: dict = {}
    for r in rows:
        for k, v in r["by_class"].items():
            tot[k] = tot.get(k, 0) + v
    return {"tool": man["tool"], "label": man["label"], "set": set_name, "run": run_name,
            "timeout_s": man["settings"]["timeout_s"], "n_systems": sum(r["n"] for r in rows),
            "n_done": sum(r["done"] for r in rows), "by_class": tot,
            "shards_done": sum(1 for r in rows if r["state"] == "done"), "n_shards": len(rows),
            "incomplete": [r["shard"] for r in rows if r["state"] != "done"], "shards": rows}


def cbatch_format_status(st: dict, verbose=False) -> str:
    lines = [f"{st['label']} ({st['tool']}), set {st['set']}, run {st['run']}, "
             f"limit {st['timeout_s']} s",
             f"systems {st['n_done']}/{st['n_systems']} final, shards {st['shards_done']}/"
             f"{st['n_shards']} complete",
             "classes: " + ", ".join(f"{k} {st['by_class'].get(k, 0)}" for k in CLASSES
                                     if st["by_class"].get(k) or k in ("ok", "timeout"))]
    if st["incomplete"]:
        inc = st["incomplete"]
        lines.append(f"incomplete shards ({len(inc)}): {cbatch_array_spec(inc)}")
    if verbose:
        for r in st["shards"]:
            lines.append(f"  shard {r['shard']:4d}  {r['state']:<15} {r['done']:4d}/{r['n']:<4d} "
                         + " ".join(f"{k}={v}" for k, v in sorted(r["by_class"].items())))
    return "\n".join(lines)


def cbatch_array_spec(indices) -> str:
    """[0,1,2,5,7,8] -> '0-2,5,7-8' (sbatch --array syntax)."""
    idx = sorted(set(int(i) for i in indices))
    out, i = [], 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and idx[j + 1] == idx[j] + 1:
            j += 1
        out.append(str(idx[i]) if i == j else f"{idx[i]}-{idx[j]}")
        i = j + 1
    return ",".join(out)


def _link_or_copy(src: Path, dst: Path) -> None:
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def cbatch_collect(run_dir, set_name="main", run_name="main", dest=None, label=None) -> dict:
    """Merge the complete chunks of one set/run into ``archive_<label>`` (hard links where the
    file system allows, the chunks stay untouched).  ``dest`` must not exist."""
    run_dir = Path(run_dir).resolve()
    man = cbatch_load_manifest(run_dir)
    if set_name not in man["sets"]:
        raise SystemExit(f"shard set {set_name!r} not in this run ({sorted(man['sets'])})")
    label = label or (man["label"] + ("" if (set_name, run_name) == ("main", "main")
                                      else f"_{set_name}_{run_name}"))
    dest = Path(dest) if dest else run_dir / "collected" / label
    if dest.exists():
        raise SystemExit(f"{dest} exists -- choose a new directory (never overwritten)")
    arch = dest / f"archive_{label}"
    tool = man["tool"]
    arch.mkdir(parents=True)
    if tool != "manta":
        (arch / "_meta").mkdir()
    classes, times, mem, reasons = {}, {}, {}, {}
    problems, missing, provs, tmos = [], [], {}, {}
    for ent in man["sets"][set_name]["shards"]:
        k = ent["shard"]
        chunk = cbatch_chunk_dir(run_dir, set_name, run_name, k)
        if not (chunk / "DONE.json").exists():
            missing.append(k)
            continue
        ids = [r for r, _ in cbatch_shard_rows(run_dir, set_name, k)]
        tj = chunk / "timeout.json"
        tmos[k] = json.loads(tj.read_text()) if tj.exists() else None
        metas = sorted(chunk.glob("run_meta_*.json"))
        if metas:
            pv = json.loads(metas[-1].read_text()).get("provenance", {})
            provs[k] = json.dumps(pv, sort_keys=True)
        if tool == "manta":
            st = cbatch_read_status(chunk / "status.jsonl")
            for r in ids:
                if r not in st:
                    problems.append(f"chunk {k}: {r} without status")
                    continue
                if r in classes:
                    problems.append(f"chunk {k}: {r} seen twice")
                    continue
                c = st[r]["status"]
                classes[r] = c
                times[r] = st[r].get("time_s")
                if "peak_gb" in st[r]:
                    mem[r] = st[r]["peak_gb"]
                x = chunk / "archive" / f"{r}.xyz"
                if (c == "ok") != x.exists():
                    problems.append(f"chunk {k}: {r} status {c} xyz={'yes' if x.exists() else 'no'}")
                if x.exists():
                    _link_or_copy(x, arch / f"{r}.xyz")
        else:
            specs = _load_specs(run_dir / "shards" / set_name / f"specs_{k:04d}.jsonl")
            fin = cbatch_read_final_metas(chunk / "archive", ids)
            for r in ids:
                m = fin.get(r)
                if m is None:
                    problems.append(f"chunk {k}: {r} without final _meta")
                    continue
                if r in classes:
                    problems.append(f"chunk {k}: {r} seen twice")
                    continue
                sp = specs.get(r)
                if sp is None:
                    problems.append(f"chunk {k}: {r} without spec")
                c = cbatch_class(m, sp)
                m["class"] = c
                if c == NOT_EXPRESSIBLE:
                    m["not_expressible_reason"] = (m.get("not_expressible_reason")
                                                   or (sp or {}).get("bail") or (sp or {}).get("status"))
                    reasons[m["not_expressible_reason"]] = reasons.get(m["not_expressible_reason"], 0) + 1
                classes[r] = c
                times[r] = m.get("wall_total_s", m.get("wall_s"))
                mem[r] = m.get("rss_gb")
                x = chunk / "archive" / f"{r}.xyz"
                if (m.get("n_frames", 0) > 0) != x.exists():
                    problems.append(f"chunk {k}: {r} n_frames={m.get('n_frames')} "
                                    f"xyz={'yes' if x.exists() else 'no'}")
                if x.exists():
                    _link_or_copy(x, arch / f"{r}.xyz")
                (arch / "_meta" / f"{r}.json").write_text(json.dumps(m, indent=1))
    pset = sorted(set(provs.values()))
    if len(pset) > 1:
        problems.append(f"chunks ran with {len(pset)} different code/environment states")
    tset = sorted({json.dumps(t, sort_keys=True) for t in tmos.values()})
    if len(tset) > 1:
        problems.append(f"chunks ran with {len(tset)} different limits: {tset}")
    counts: dict = {}
    for c in classes.values():
        counts[c] = counts.get(c, 0) + 1
    n = len(classes)
    nexpr = n - counts.get(NOT_EXPRESSIBLE, 0)
    (dest / f"build_{label}.json").write_text(json.dumps(classes, indent=1, sort_keys=True))
    (dest / f"buildtime_{label}.json").write_text(json.dumps(times, indent=0, sort_keys=True))
    (dest / f"buildmem_{label}.json").write_text(json.dumps(mem, indent=0, sort_keys=True))
    (dest / "resubmit.txt").write_text(cbatch_array_spec(missing) + "\n")
    summary = {
        "label": label, "tool": tool, "set": set_name, "run": run_name,
        "n_systems_merged": n, "n_systems_in_set": man["sets"][set_name]["n_systems"],
        "n_expressible": nexpr, "by_class": counts, "not_expressible_reasons": reasons,
        "coverage_of_expressible": round(counts.get("ok", 0) / nexpr, 4) if nexpr else None,
        "coverage_of_all": round(counts.get("ok", 0) / n, 4) if n else None,
        "n_xyz": sum(1 for _ in arch.glob("*.xyz")),
        "chunks_merged": len(tmos), "chunks_missing": missing,
        "resubmit_array": cbatch_array_spec(missing),
        "timeout": [json.loads(t) for t in tset],
        "provenance": [json.loads(p) for p in pset][:3],
        "manifest_settings": man["settings"],
        "problems": problems[:500], "n_problems": len(problems)}
    (dest / f"summary_{label}.json").write_text(json.dumps(summary, indent=1))
    return summary


def cbatch_repeat_stats(run_dir, run_a="main", run_b="main") -> dict:
    """Second build of the ``repeat`` subset against the first build of the same IDs in the
    ``main`` set: byte-identical / different / only in one run / in neither."""
    run_dir = Path(run_dir).resolve()
    man = cbatch_load_manifest(run_dir)
    if "repeat" not in man["sets"]:
        raise SystemExit("this run has no repeat set (prepare --repeat N)")
    where = {}
    for ent in man["sets"]["main"]["shards"]:
        for r, _ in cbatch_shard_rows(run_dir, "main", ent["shard"]):
            where[r] = ent["shard"]
    res = {"identical": [], "different": [], "only_first": [], "only_second": [], "neither": []}
    detail = {}
    for ent in man["sets"]["repeat"]["shards"]:
        for r, _ in cbatch_shard_rows(run_dir, "repeat", ent["shard"]):
            a = cbatch_chunk_dir(run_dir, "main", run_a, where[r]) / "archive" / f"{r}.xyz"
            b = cbatch_chunk_dir(run_dir, "repeat", run_b, ent["shard"]) / "archive" / f"{r}.xyz"
            ea, eb = a.exists(), b.exists()
            if ea and eb:
                ba, bb = a.read_bytes(), b.read_bytes()
                if ba == bb:
                    res["identical"].append(r)
                else:
                    res["different"].append(r)
                    la = [ln for ln in ba.decode().splitlines() if " frame" in ln]
                    lb = [ln for ln in bb.decode().splitlines() if " frame" in ln]
                    detail[r] = {"frames": [len(la), len(lb)], "same_labels": la == lb}
            elif ea:
                res["only_first"].append(r)
            elif eb:
                res["only_second"].append(r)
            else:
                res["neither"].append(r)
    out = {k: len(v) for k, v in res.items()}
    both = out["identical"] + out["different"]
    out["identical_fraction_of_both"] = round(out["identical"] / both, 4) if both else None
    out["ids"] = res
    out["different_detail"] = detail
    (run_dir / f"repeat_stats_{run_a}_{run_b}.json").write_text(json.dumps(out, indent=1))
    return out
