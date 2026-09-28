"""Architector and molSimplify as SMILES constructors, one definition per tool.

The dashboard (Submit Job and ORCA Builder, through the shared structure
editor) and the CONTROL pipeline (``smiles_converter=ARCHITECTOR`` /
``MOLSIMPLIFY``) all build through :func:`build_frames`, the same way the MANTA
entry points all build through :mod:`delfin.common.manta_build`.

How a build runs
----------------
1. The SMILES is split into metal + free ligands (:func:`split_complex_smiles`):
   every metal-donor bond is cut in the SMILES *as written*; a dative bond
   leaves the donor neutral (L-type), a covalent bond of order *n* leaves it
   with charge -*n* (X-type), a donor carbon with two sigma bonds left becomes
   a singlet carbene, ``[MHn]`` becomes *n* hydrides.  The coordination number
   is the sum of the denticities and the oxidation state is the total charge
   minus the ligand charges.  Neither tool is left to guess either of them:
   without ``coreCN`` Architector fills the core up with water, and without
   ``metal_ox`` it reads the SMILES charge as the metal charge.
2. The spec is turned into the tool's own input (:func:`architector_input`,
   :func:`molsimplify_inputs`).
3. The tool runs in a subprocess (``python -m delfin.common.external_builders``)
   so that its crashes, its ``quit()`` calls and its stdout never reach the
   caller (a Voila kernel, a pipeline), and so that it can live in another
   Python environment than DELFIN: ``DELFIN_ARCHITECTOR_PYTHON`` /
   ``DELFIN_MOLSIMPLIFY_PYTHON`` name that interpreter; default is the running
   one.
4. Every frame the tool returns comes back, lowest energy first when the tool
   gives energies.  Nothing falls back to another builder: a missing tool or a
   failed build is an error message.

Known pitfalls handled here (Architector 0.0.10, molSimplify 2.0.0):

* Architector computes a default spin through an API that mendeleev >= 1.0 no
  longer has, so the spin is always computed here and passed explicitly.
* rdkit is imported before the tools, otherwise a foreign libstdc++ from
  ``LD_LIBRARY_PATH`` (e.g. ORCA's) can be loaded first and rdkit fails.
* molSimplify's ``-isomers`` only works with its own ligand dictionary, so it
  builds one structure per geometry of the coordination number.
* molSimplify fills sites in ligand order: ligands go by decreasing denticity.
* ``-keepHs yes``: the protonation written in the SMILES is final.
* molSimplify drops implicit H on very small SMILES ligands and treats any
  ligand string similar to one of its dictionary names as that ligand
  (``CO`` -> carbonyl); ligands are re-spelled until neither happens.

This module has no import-time dependency beyond the standard library, so the
worker runs in a tool environment that only has rdkit, openbabel and the tool.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

#: Tools this module knows, keyed by the name used everywhere in DELFIN.
TOOLS = {
    "architector": {
        "display": "Architector",
        "module": "architector",
        "pip": "architector",
        "python_env": "DELFIN_ARCHITECTOR_PYTHON",
    },
    "molsimplify": {
        "display": "molSimplify",
        "module": "molSimplify",
        "pip": "molSimplify",
        "python_env": "DELFIN_MOLSIMPLIFY_PYTHON",
    },
}

#: molSimplify geometries (its coordinations.dict names) per coordination number.
MS_GEOMS = {2: ["li"], 3: ["tpl"], 4: ["sqp", "thd"], 5: ["spy", "tbp"],
            6: ["oct", "tpr"], 7: ["pbp"], 8: ["sqap", "tdhd"]}

#: Wall-clock limit of one build subprocess, seconds.
DEFAULT_TIMEOUT = float(os.environ.get("DELFIN_EXTERNAL_BUILDER_TIMEOUT", "1800"))

# Everything that is not a metal for the purpose of "which atom is the core".
_NON_METALS = {
    "H", "He", "B", "C", "N", "O", "F", "Ne", "Si", "P", "S", "Cl", "Ar",
    "Ge", "As", "Se", "Br", "Kr", "Sb", "Te", "I", "Xe", "At", "Rn",
}


class BuildError(ValueError):
    """A SMILES that cannot be handed to the tool, with the reason."""


def tool_key(tool: str) -> str:
    key = str(tool or "").strip().lower()
    if key not in TOOLS:
        raise ValueError(f"unknown external builder {tool!r} "
                         f"(known: {', '.join(sorted(TOOLS))})")
    return key


def tool_python(tool: str, environ=None) -> str:
    """The interpreter the tool runs in: its env variable, else this one."""
    env = os.environ if environ is None else environ
    return (env.get(TOOLS[tool_key(tool)]["python_env"]) or "").strip() or sys.executable


def not_installed_message(tool: str, python: str) -> str:
    info = TOOLS[tool_key(tool)]
    return (f"{info['display']} is not installed in {python}. Install it with "
            f"pip install 'delfin-complat[ai-complex]' (or pip install {info['pip']}), "
            f"or set {info['python_env']} to a Python interpreter that has it.")


# ---------------------------------------------------------------------------
# SMILES -> tool-neutral spec
# ---------------------------------------------------------------------------

def split_complex_smiles(smiles: str) -> dict:
    """Metal, oxidation state, coordination number and free ligands of a SMILES.

    Returns ``{'metal', 'metal_ox', 'cn', 'total_charge', 'ligands': [...]}``;
    each ligand is ``{'smiles', 'coordList' (0-based), 'charge', 'denticity',
    'flags'}``.  Raises :class:`BuildError` for what neither tool can build: no
    metal, more than one metal, a fragment not bound to the metal, a ligand
    left as a radical.
    """
    from rdkit import Chem

    smi = str(smiles or "").strip()
    raw = Chem.MolFromSmiles(smi, sanitize=False)
    if raw is None or raw.GetNumAtoms() == 0:
        raise BuildError("the SMILES could not be parsed")
    metals = [a.GetIdx() for a in raw.GetAtoms() if a.GetSymbol() not in _NON_METALS]
    if not metals:
        raise BuildError("the SMILES contains no metal atom (Architector and "
                         "molSimplify build metal complexes)")
    if len(metals) != 1:
        raise BuildError(f"the SMILES contains {len(metals)} metal atoms; only "
                         "mononuclear complexes can be built with this tool")
    total_charge = int(sum(a.GetFormalCharge() for a in raw.GetAtoms()))
    # Freeze the hydrogen count every atom has as written.  An unbracketed
    # donor (the Cl of [Pt](Cl)...) would otherwise be given a fresh implicit H
    # once its metal bond is cut and come back as HCl instead of chloride.
    raw.UpdatePropertyCache(strict=False)
    for a in raw.GetAtoms():
        a.SetNumExplicitHs(a.GetTotalNumHs())
        a.SetNoImplicit(True)
    m = metals[0]
    metal_sym = raw.GetAtomWithIdx(m).GetSymbol()
    rw = Chem.RWMol(raw)
    for a in rw.GetAtoms():
        a.SetIntProp("_ridx", a.GetIdx())
    n_hydride = rw.GetAtomWithIdx(m).GetNumExplicitHs()
    donors, flags = [], []
    for b in list(rw.GetAtomWithIdx(m).GetBonds()):
        donors.append(b.GetOtherAtomIdx(m))
    for d in donors:
        rw.RemoveBond(m, d)
    for b in rw.GetBonds():  # dative bonds inside a ligand -> charge separated
        if b.GetBondType() == Chem.BondType.DATIVE:
            b.SetBondType(Chem.BondType.SINGLE)
            s, e = b.GetBeginAtom(), b.GetEndAtom()
            s.SetFormalCharge(s.GetFormalCharge() + 1)
            e.SetFormalCharge(e.GetFormalCharge() - 1)
    rw.RemoveAtom(m)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.AssignRadicals(mol)
    dset = set(donors)
    for a in mol.GetAtoms():
        if a.GetIntProp("_ridx") not in dset:
            continue
        r = a.GetNumRadicalElectrons()
        if r == 0:
            continue                      # closed shell without the metal: L-type
        val = sum(bd.GetBondTypeAsDouble() for bd in a.GetBonds()) + a.GetTotalNumHs()
        if a.GetSymbol() == "C" and int(round(val)) == 2:
            a.SetFormalCharge(0)          # two sigma bonds only: singlet carbene
            a.SetNumRadicalElectrons(2)
            a.SetNoImplicit(True)
            continue
        a.SetFormalCharge(a.GetFormalCharge() - r)   # X-type: fill with electrons
        a.SetNumRadicalElectrons(0)
        a.SetNoImplicit(True)
    ligands = []
    for fi in Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False):
        fdon = [i for i in fi if mol.GetAtomWithIdx(i).GetIntProp("_ridx") in dset]
        if not fdon:
            raise BuildError("the SMILES has a fragment that is not bound to the "
                             "metal (counter-ion or solvent); give the complex alone")
        sub = Chem.RWMol(mol)
        for i in sorted(set(range(mol.GetNumAtoms())) - set(fi), reverse=True):
            sub.RemoveAtom(i)
        fm = sub.GetMol()
        try:
            Chem.SanitizeMol(fm)
        except Exception as exc:
            raise BuildError(f"a ligand could not be sanitised after cutting the "
                             f"metal off: {exc}") from exc
        rid = {a.GetIntProp("_ridx"): a.GetIdx() for a in fm.GetAtoms()}
        dl = [rid[mol.GetAtomWithIdx(i).GetIntProp("_ridx")] for i in fdon]
        if any(fm.GetAtomWithIdx(i).GetAtomicNum() == 1 for i in dl):
            raise BuildError("a hydrogen is bound to the metal as a donor atom")
        fh = Chem.RemoveHs(fm)
        rid2 = {a.GetIntProp("_ridx"): a.GetIdx() for a in fh.GetAtoms()}
        dl = [rid2[fm.GetAtomWithIdx(i).GetIntProp("_ridx")] for i in dl]
        mk = Chem.Mol(fh)
        Chem.Kekulize(mk, clearAromaticFlags=True)
        s = Chem.MolToSmiles(mk, kekuleSmiles=True, canonical=True)
        order = list(mk.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"])
        pos = {old: i for i, old in enumerate(order)}
        coord = sorted(pos[i] for i in dl)
        chk = Chem.MolFromSmiles(s)
        if chk is None:
            raise BuildError(f"ligand {s} does not round-trip through RDKit")
        lflags = []
        for c in coord:
            ca = chk.GetAtomWithIdx(c)
            if ca.GetSymbol() == "C" and ca.GetNumRadicalElectrons() == 2:
                lflags.append("carbene")
        if sum(a.GetNumRadicalElectrons() for a in chk.GetAtoms()) and "carbene" not in lflags:
            raise BuildError(f"ligand {s} is left as a radical after cutting the metal off")
        ligands.append({"smiles": s, "coordList": coord,
                        "charge": int(Chem.GetFormalCharge(chk)),
                        "denticity": len(coord), "flags": lflags + flags})
    for _ in range(n_hydride):
        ligands.append({"smiles": "[H-]", "coordList": [0], "charge": -1,
                        "denticity": 1, "flags": ["hydride"]})
    if not ligands:
        raise BuildError("the metal has no ligands in the SMILES")
    metal_ox = total_charge - sum(l["charge"] for l in ligands)
    if not 0 <= metal_ox <= 8:
        raise BuildError(f"oxidation state {metal_ox} of {metal_sym} read from the "
                         "SMILES is outside 0..8")
    return {"metal": metal_sym, "metal_ox": int(metal_ox),
            "cn": int(sum(l["denticity"] for l in ligands)),
            "total_charge": total_charge, "ligands": ligands}


# ---------------------------------------------------------------------------
# spec -> tool input (needs the tool itself; runs in the worker)
# ---------------------------------------------------------------------------

def default_unpaired(metal: str, ox: int) -> int:
    """Architector's own default number of unpaired electrons.

    Its reference table when the oxidation state is the reference one, else
    the aufbau rule of mendeleev -- through the current mendeleev API, because
    the one Architector 0.0.10 calls is gone in mendeleev >= 1.0.
    """
    from architector import io_ptable
    import mendeleev

    if ox == io_ptable.metal_charge_dict.get(metal, 100):
        return int(io_ptable.metal_spin_dict[metal])
    return int(mendeleev.element(metal).ec.ionize(ox).unpaired_electrons())


def _ob_parses(smi: str) -> bool:
    from openbabel import openbabel as ob

    conv = ob.OBConversion()
    conv.SetInFormat("smi")
    m = ob.OBMol()
    return bool(conv.ReadString(m, smi)) and m.NumAtoms() > 0


def architector_input(spec: dict, params=None) -> dict:
    from architector import io_ptable
    from architector.io_core import Geometries

    if spec["metal"] not in io_ptable.all_metals:
        raise BuildError(f"Architector does not support the metal {spec['metal']}")
    if spec["cn"] not in Geometries().cn_geo_dict:
        raise BuildError(f"Architector has no geometry for coordination number {spec['cn']}")
    ligs = []
    for l in spec["ligands"]:
        if not _ob_parses(l["smiles"]):
            raise BuildError(f"OpenBabel cannot read ligand {l['smiles']}")
        ligs.append({"smiles": l["smiles"], "coordList": list(l["coordList"])})
    p = {"metal_ox": int(spec["metal_ox"]),
         "metal_spin": default_unpaired(spec["metal"], int(spec["metal_ox"]))}
    p.update(params or {})
    return {"core": {"metal": spec["metal"], "coreCN": int(spec["cn"])},
            "ligands": ligs, "parameters": p}


_LICORES = None


def _ms_collides(s: str) -> bool:
    import difflib

    global _LICORES
    if _LICORES is None:
        from molSimplify.Scripts.io import getlicores
        _LICORES = list(getlicores().keys())
    if s in _LICORES:
        return True
    return max(difflib.SequenceMatcher(None, s, k).ratio() for k in _LICORES) > 0.6


def _smiles_variants(smi, coord):
    """Equivalent spellings of a ligand SMILES, donor positions re-mapped."""
    from rdkit import Chem

    m = Chem.MolFromSmiles(smi)
    if m is None:
        return
    Chem.Kekulize(m, clearAromaticFlags=True)
    for allh in ((True, False) if m.GetNumAtoms() <= 3 else (False, True)):
        for root in range(-1, m.GetNumAtoms()):
            s = Chem.MolToSmiles(m, kekuleSmiles=True, canonical=True,
                                 rootedAtAtom=root, allHsExplicit=allh)
            order = list(m.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"])
            pos = {old: i for i, old in enumerate(order)}
            yield s, sorted(pos[c] for c in coord)


def molsimplify_inputs(spec: dict) -> list:
    """One molSimplify input per geometry of the coordination number."""
    if spec["cn"] not in MS_GEOMS:
        raise BuildError(f"molSimplify has no geometry for coordination number {spec['cn']}")
    if max(l["denticity"] for l in spec["ligands"]) > 6:
        raise BuildError("molSimplify cannot place a ligand with more than six donors")
    smis, cats = [], []
    for l in sorted(spec["ligands"], key=lambda l: -l["denticity"]):
        if not _ob_parses(l["smiles"]):
            raise BuildError(f"OpenBabel cannot read ligand {l['smiles']}")
        chosen = next(((s, c) for s, c in _smiles_variants(l["smiles"], l["coordList"])
                       if not _ms_collides(s)), None)
        if chosen is None:
            raise BuildError(f"ligand {l['smiles']} cannot be spelled so that molSimplify "
                             "does not mistake it for one of its dictionary ligands")
        smis.append(chosen[0])
        cats.append(chosen[1])
    spin = default_unpaired(spec["metal"], int(spec["metal_ox"]))
    base = {"-core": spec["metal"],
            "-lig": ",".join(smis),
            "-ligocc": ",".join("1" for _ in smis),
            "-smicat": "[" + ",".join("[" + ",".join(str(c + 1) for c in cc) + "]"
                                      for cc in cats) + "]",
            "-coord": str(spec["cn"]),
            "-oxstate": str(int(spec["metal_ox"])),
            "-spinmultiplicity": str(spin + 1),
            "-ff": "uff", "-ffoption": "BA",
            "-keepHs": ",".join("yes" for _ in smis)}
    return [(g, dict(base, **{"-geometry": g})) for g in MS_GEOMS[spec["cn"]]]


# ---------------------------------------------------------------------------
# the worker (runs in the tool's interpreter)
# ---------------------------------------------------------------------------

def _run_architector(spec: dict, work: str, options: dict) -> tuple:
    from architector import build_complex
    from architector.io_process_input import inparse
    from architector import io_ptable

    params = {"temp_prefix": work.rstrip("/") + "/"}
    if options.get("mode", "full") == "full":
        # every distinct symmetry (isomer) per core geometry is relaxed and returned
        params.update({"n_symmetries": 10, "n_conformers": 10})
    inp = architector_input(spec, params)
    out = build_complex(inp)
    if not out and max(l["denticity"] for l in spec["ligands"]) >= 2:
        # Architector's own rescue for chelates: scaled metal radii.
        from architector.complex_construction import build_complex_driver
        for larger in (True, False):
            out = build_complex_driver(io_ptable.map_metal_radii(inparse(inp), larger=larger))
            out = {k: v for k, v in out.items() if "_init_only" not in k}
            if out:
                break
    frames = []
    for key, v in out.items():
        at = v.get("ase_atoms")
        if at is None:
            continue
        e = v.get("energy")
        frames.append({"label": str(key), "symbols": at.get_chemical_symbols(),
                       "coords": at.get_positions().tolist(),
                       "energy": float(e) if e is not None else None})
    return frames, {}


def _run_molsimplify(spec: dict, work: str, options: dict) -> tuple:
    import glob
    from molSimplify.Scripts.generator import startgen_pythonic

    frames, per_geo = [], {}
    for geo, d in molsimplify_inputs(spec):
        rd = os.path.join(work, geo)
        name = f"delfin_{geo}"
        d = dict(d, **{"-rundir": rd + "/", "-name": name})
        try:
            startgen_pythonic(d, write=True)
        except (Exception, SystemExit) as exc:   # molSimplify also quit()s
            # Its bond-length ANN (octahedral only) is a Keras 2.2 model that
            # current Keras refuses to load; the build itself does not need it.
            try:
                import shutil
                shutil.rmtree(rd, ignore_errors=True)
                startgen_pythonic(dict(d, **{"-skipANN": "True"}), write=True)
                per_geo[geo + "_ann"] = f"ANN skipped after {type(exc).__name__}"
            except (Exception, SystemExit) as exc2:
                per_geo[geo] = f"{type(exc2).__name__}: {str(exc2)[:160]}"
                continue
        xyzs = sorted(glob.glob(os.path.join(rd, "**", name + ".xyz"), recursive=True))
        if not xyzs:
            per_geo[geo] = "no structure written"
            continue
        per_geo[geo] = "ok"
        for f in xyzs:
            ln = Path(f).read_text().splitlines()
            n = int(ln[0].split()[0])
            syms, xyz = [], []
            for row in ln[2:2 + n]:
                p = row.split()
                syms.append(p[0])
                xyz.append([float(p[1]), float(p[2]), float(p[3])])
            frames.append({"label": geo, "symbols": syms, "coords": xyz, "energy": None})
    return frames, {"per_geometry": per_geo}


def _worker(tool: str, request_path: str, result_path: str) -> int:
    # rdkit before any tool: see the module docstring (libstdc++).
    result = {"frames": [], "error": None}
    try:
        from rdkit import RDLogger  # noqa: F401
        RDLogger.DisableLog("rdApp.*")
        req = json.loads(Path(request_path).read_text())
        key = tool_key(tool)
        import importlib.util
        if importlib.util.find_spec(TOOLS[key]["module"]) is None:
            result["error"] = not_installed_message(key, sys.executable)
            result["not_installed"] = True
        else:
            spec = split_complex_smiles(req["smiles"])
            result["spec"] = spec
            run = _run_architector if key == "architector" else _run_molsimplify
            frames, extra = run(spec, req["work"], req.get("options") or {})
            result.update(extra)
            result["frames"] = frames
            if not frames:
                detail = "; ".join(f"{k}: {v}" for k, v in extra.get("per_geometry", {}).items())
                result["error"] = (f"{TOOLS[key]['display']} returned no structure"
                                   + (f" ({detail})" if detail else ""))
    except BuildError as exc:
        result["error"] = str(exc)
    except (Exception, SystemExit) as exc:
        import traceback
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["traceback"] = traceback.format_exc()[-3000:]
    Path(result_path).write_text(json.dumps(result))
    return 0


# ---------------------------------------------------------------------------
# the one entry point
# ---------------------------------------------------------------------------

def _xyz_body(symbols, coords) -> str:
    return "\n".join(f"{s}  {x:.6f}  {y:.6f}  {z:.6f}" for s, (x, y, z) in zip(symbols, coords))


def build_frames(tool: str, smiles: str, *, options=None, python=None,
                 timeout=None, workdir=None, environ=None):
    """Build *smiles* with *tool*; ``([(xyz_body, n_atoms, label), ...], error)``.

    Frames come back lowest energy first when the tool gives energies (it is
    Architector's order anyway), else in the tool's order.  *error* is ``None``
    on success, else a message for the user.  *options* for Architector:
    ``{'mode': 'full'}`` (default: every isomer it finds) or ``'default'``
    (its own defaults, one structure per core geometry).
    """
    key = tool_key(tool)
    info = TOOLS[key]
    python = python or tool_python(key, environ)
    if python == sys.executable:
        import importlib.util
        if importlib.util.find_spec(info["module"]) is None:
            return [], not_installed_message(key, python)
    import delfin

    repo_root = str(Path(delfin.__file__).resolve().parent.parent)
    env = dict(os.environ if environ is None else environ)
    env["PYTHONPATH"] = repo_root + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env["PYTHONHASHSEED"] = "0"
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        env.setdefault(var, "1")
    env.setdefault("OMP_STACKSIZE", "1G")
    timeout = DEFAULT_TIMEOUT if timeout is None else float(timeout)

    tmp = tempfile.mkdtemp(prefix=f"delfin_{key}_", dir=workdir)
    work = os.path.join(tmp, "work")
    os.makedirs(work, exist_ok=True)
    request_path = os.path.join(tmp, "request.json")
    result_path = os.path.join(tmp, "result.json")
    log_path = os.path.join(tmp, "tool.log")
    Path(request_path).write_text(json.dumps(
        {"smiles": str(smiles).strip(), "work": work, "options": dict(options or {})}))
    cmd = [python, "-m", "delfin.common.external_builders", key, request_path, result_path]
    keep = os.environ.get("DELFIN_EXTERNAL_BUILDER_KEEP", "0") == "1"
    try:
        with open(log_path, "w") as log:
            try:
                proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                                        cwd=work, env=env, start_new_session=True)
            except OSError as exc:
                return [], f"{info['display']} could not be started with {python}: {exc}"
            try:
                proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                from delfin.common.manta_build import kill_process_group
                kill_process_group(proc)
                return [], f"{info['display']} did not finish within {int(timeout)} s"
        if not os.path.exists(result_path):
            tail = Path(log_path).read_text(errors="replace")[-600:]
            keep = True
            return [], (f"{info['display']} ended without a result (exit code "
                        f"{proc.returncode}); log kept in {log_path}:\n{tail}")
        res = json.loads(Path(result_path).read_text())
        frames = res.get("frames") or []
        if res.get("traceback"):
            keep = True
        if res.get("error") and not frames:
            return [], f"{info['display']}: {res['error']}"
        if frames and all(f.get("energy") is not None for f in frames):
            frames = sorted(frames, key=lambda f: f["energy"])
        out = []
        for i, f in enumerate(frames, start=1):
            e = f.get("energy")
            label = f"{info['display']} {f.get('label') or i}"
            if e is not None:
                label += f" (E = {e:.4f} eV)"   # Architector's xTB energy, ASE units
            out.append((_xyz_body(f["symbols"], f["coords"]), len(f["symbols"]), label))
        return out, None
    finally:
        if not keep:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)


def build_first_xyz(tool: str, smiles: str, **kwargs):
    """The first (lowest-energy) frame as an XYZ body, for the CONTROL pipeline."""
    frames, error = build_frames(tool, smiles, **kwargs)
    if error:
        return None, error
    if not frames:
        return None, f"{TOOLS[tool_key(tool)]['display']} returned no structure"
    return frames[0][0], None


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: python -m delfin.common.external_builders TOOL REQUEST.json RESULT.json")
        sys.exit(2)
    sys.exit(_worker(*sys.argv[1:]))
