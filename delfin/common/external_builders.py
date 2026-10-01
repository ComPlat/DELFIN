"""Architector, molSimplify and epic-MACE as SMILES constructors, one definition per tool.

The dashboard (Submit Job and ORCA Builder, through the shared structure
editor) and the CONTROL pipeline (``smiles_converter=ARCHITECTOR`` /
``MOLSIMPLIFY`` / ``MACE``) all build through :func:`build_frames`, the same
way the MANTA entry points all build through :mod:`delfin.common.manta_build`.

epic-MACE (Chernyshov & Pidko, J. Chem. Theory Comput. 2024, 20, 2313;
github.com/EPiCs-group/epic-mace, GPL-3.0) is called as an external program in
an environment of its own -- Python 3.7 with RDKit 2020.09 -- and is never
imported into DELFIN.  ``python -m delfin.installer --install epic-mace``
builds that environment; ``DELFIN_MACE_PYTHON`` overrides where it is.

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
   :func:`molsimplify_inputs`, :func:`mace_job`).  For epic-MACE the spec is
   cut in the caller's interpreter (current RDKit) and handed to the worker,
   which only writes it into MACE's input.
3. The tool runs in a subprocess (``python -m delfin.common.external_builders``)
   so that its crashes, its ``quit()`` calls and its stdout never reach the
   caller (a Voila kernel, a pipeline), and so that it can live in another
   Python environment than DELFIN: ``DELFIN_ARCHITECTOR_PYTHON`` /
   ``DELFIN_MOLSIMPLIFY_PYTHON`` / ``DELFIN_MACE_PYTHON`` name that
   interpreter; default is the tool's managed environment when DELFIN built
   one (epic-MACE), else the running one.
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

epic-MACE (GitHub commit efb5778e, the hapto centroids and the TET / SPY /
TBP / SAN geometries of its unreleased 0.6.0; PyPI 0.5.0 has OH and SP only):

* Each ligand of the spec becomes a ligand SMILES with atom-mapped donors, and
  ``mace.ComplexFromLigands`` puts them on the central atom ``[M+ox]``.
* Donors of one ligand bonded to each other form one hapto group, written as
  one centroid dummy atom: bonded to one carbon of a 5- or 6-membered carbon
  ring (MACE expands it to the ring), else bonded to every group atom with the
  bonds inside the group removed.  The centroids (element X in MACE's output)
  are dropped from the frames.
* The geometry follows from the number of sites, not from the CN (a hapto
  group is one site): ``paper`` = OH for 6, SP for 4 (the geometries of the
  paper); ``extended`` (default) adds SPY + TBP for 5, TET for 4 and SAN for
  two hapto centroids.  Every geometry that fits is built.
* Stereomers: ``GetStereomers(regime='all', dropEnantiomers=False,
  minTransCycle=None, merRule=False)`` -- every arrangement at the metal and at
  unassigned ligand stereocentres, enantiomers kept; ``merRule=True`` (the
  library default) returns no stereomer at all for fac-only tripods.
* 3D: ``AddConformers(numConfs=10, maxAttempts=10, rmsThresh=-1)`` per
  stereomer, ordered by MACE's force-field energy (UFF, kcal/mol).  MACE has
  no random seed, so two builds can differ.

This module has no import-time dependency beyond the standard library and runs
under Python 3.7, so the worker runs in a tool environment that only has
rdkit and the tool (and openbabel for Architector and molSimplify).
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
    "mace": {
        "display": "epic-MACE",
        "module": "mace",
        "python_env": "DELFIN_MACE_PYTHON",
        # Built by install_ai_tools.sh (INSTALL_EPIC_MACE) under the AI tools
        # root: Python 3.7 + RDKit 2020.09, which DELFIN itself cannot run in.
        "managed_env": ".mamba_env/epic_mace",
        "install": "python -m delfin.installer --install epic-mace",
        "spec_in_caller": True,
        "energy_unit": "kcal/mol, UFF",
        # Geometry by geometry as listed (the paper's OH / SP first), lowest
        # energy first within each: force-field energies of a square and a
        # tetrahedron are no ranking of the two.
        "keep_order": True,
    },
}

#: epic-MACE settings, the same for the dashboard and for CONTROL.
MACE_DEFAULTS = {"geometries": "extended", "num_confs": 10, "max_attempts": 10}

#: epic-MACE geometries per number of donor sites (a hapto group is one site).
MACE_GEOMS = {"paper": {6: ["OH"], 4: ["SP"]},
              "extended": {6: ["OH"], 5: ["SPY", "TBP"], 4: ["SP", "TET"]}}

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


def managed_python(tool: str, environ=None):
    """The interpreter of the environment DELFIN's installer builds for *tool*.

    ``<AI tools root>/<managed_env>/bin/python`` (the root is
    ``DELFIN_AI_TOOLS_ROOT``, else ``~/.delfin/ai_tools``, where the installer
    puts it), or ``None`` when the tool has no environment of its own.
    """
    rel = TOOLS[tool_key(tool)].get("managed_env")
    if not rel:
        return None
    env = os.environ if environ is None else environ
    root = (env.get("DELFIN_AI_TOOLS_ROOT") or "").strip()
    base = Path(root).expanduser() if root else Path.home() / ".delfin" / "ai_tools"
    return str(base / rel / "bin" / "python")


def tool_python(tool: str, environ=None) -> str:
    """The interpreter the tool runs in.

    Its env variable, else its managed environment when that exists, else this
    one.
    """
    env = os.environ if environ is None else environ
    named = (env.get(TOOLS[tool_key(tool)]["python_env"]) or "").strip()
    if named:
        return named
    managed = managed_python(tool, env)
    if managed and os.access(managed, os.X_OK):
        return managed
    return sys.executable


def not_installed_message(tool: str, python: str) -> str:
    info = TOOLS[tool_key(tool)]
    if info.get("install"):
        return (f"{info['display']} is not installed in {python}. Install it with "
                f"{info['install']} (it builds an environment of its own), "
                f"or set {info['python_env']} to a Python interpreter that has it.")
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


def _architector_sites(lig: dict) -> int:
    """Core sites Architector gives one ligand: 3 for a 'sandwich' (every donor in one aromatic
    ring, more than two donors -- Architector's own test), else one per donor.  Architector fills
    coreCN minus its site count with water, so coreCN is counted the same way."""
    from architector import io_obabel

    c = lig["coordList"]
    if len(c) > 2:
        obmol = io_obabel.get_obmol_smiles(lig["smiles"])  # keep alive: rings point into it
        for ring in obmol.GetSSSR():
            if all(ring.IsInRing(x + 1) for x in c) and ring.IsAromatic():
                return 3
    return len(c)


def architector_input(spec: dict, params=None) -> dict:
    from architector import io_ptable
    from architector.io_core import Geometries

    if spec["metal"] not in io_ptable.all_metals:
        raise BuildError(f"Architector does not support the metal {spec['metal']}")
    ligs = []
    for l in spec["ligands"]:
        if not _ob_parses(l["smiles"]):
            raise BuildError(f"OpenBabel cannot read ligand {l['smiles']}")
        ligs.append({"smiles": l["smiles"], "coordList": list(l["coordList"])})
    core_cn = sum(_architector_sites(l) for l in ligs)
    if core_cn not in Geometries().cn_geo_dict:
        raise BuildError(f"Architector has no geometry for coordination number {core_cn}")
    p = {"metal_ox": int(spec["metal_ox"]),
         "metal_spin": default_unpaired(spec["metal"], int(spec["metal_ox"]))}
    p.update(params or {})
    return {"core": {"metal": spec["metal"], "coreCN": int(core_cn)},
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


def _mace_donor_groups(mol, coord):
    """Connected components of the donor atoms in the ligand graph (sorted)."""
    cs = set(coord)
    seen, groups = set(), []
    for d in sorted(cs):
        if d in seen:
            continue
        comp, stack = [], [d]
        seen.add(d)
        while stack:
            a = stack.pop()
            comp.append(a)
            for n in mol.GetAtomWithIdx(a).GetNeighbors():
                j = n.GetIdx()
                if j in cs and j not in seen:
                    seen.add(j)
                    stack.append(j)
        groups.append(sorted(comp))
    return groups


def _mace_carbon_ring(mol, group):
    """The ring a hapto group is, when it is a whole 5- or 6-membered carbon ring."""
    if len(group) not in (5, 6):
        return None
    gs = set(group)
    if not all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in group):
        return None
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == len(gs) and set(ring) == gs:
            return list(ring)
    return None


def _mace_ring_anchor(rw, ring):
    """Give one ring carbon a free valence for the centroid; its index.

    Only the Lewis structure changes, never the atoms or their hydrogens:
    adjacent radical carbons are joined into a double bond; a radical carbon
    left without a double bond, or a charged / radical carbon, is the anchor;
    else a ring double bond at the anchor becomes single and its partner a
    carbocation (an sp2 cation keeps the ring planar for RDKit 2020.09, where a
    radical would be typed sp3 and MACE would read eta1 instead of the ring).
    """
    from rdkit import Chem

    rs = set(ring)
    n = len(ring)
    for k in range(n):
        i, j = ring[k], ring[(k + 1) % n]
        ai, aj = rw.GetAtomWithIdx(i), rw.GetAtomWithIdx(j)
        b = rw.GetBondBetweenAtoms(i, j)
        if (b is not None and b.GetBondType() == Chem.BondType.SINGLE
                and ai.GetNumRadicalElectrons() > 0 and aj.GetNumRadicalElectrons() > 0
                and not any(x.GetBondType() == Chem.BondType.DOUBLE
                            for x in list(ai.GetBonds()) + list(aj.GetBonds()))):
            b.SetBondType(Chem.BondType.DOUBLE)
            ai.SetNumRadicalElectrons(ai.GetNumRadicalElectrons() - 1)
            aj.SetNumRadicalElectrons(aj.GetNumRadicalElectrons() - 1)
    for i in ring:
        a = rw.GetAtomWithIdx(i)
        if a.GetNumRadicalElectrons() > 0 and not any(
                x.GetBondType() == Chem.BondType.DOUBLE for x in a.GetBonds()):
            a.SetNumRadicalElectrons(a.GetNumRadicalElectrons() - 1)
            return i
    for i in ring:
        a = rw.GetAtomWithIdx(i)
        if a.GetFormalCharge() != 0:
            a.SetFormalCharge(0)
            return i
        if a.GetNumRadicalElectrons() > 0:
            a.SetNumRadicalElectrons(a.GetNumRadicalElectrons() - 1)
            return i
    for i in ring:
        for b in rw.GetAtomWithIdx(i).GetBonds():
            j = b.GetOtherAtomIdx(i)
            if j in rs and b.GetBondType() == Chem.BondType.DOUBLE:
                b.SetBondType(Chem.BondType.SINGLE)
                p = rw.GetAtomWithIdx(j)
                p.SetFormalCharge(p.GetFormalCharge() + 1)
                return i
    return ring[0]


def mace_ligand(lig: dict) -> tuple:
    """A spec ligand as epic-MACE's ligand SMILES; ``(smiles, sites)``.

    Donor atoms carry atom-map number 1; each hapto group becomes one centroid
    dummy ``[*:1]`` (see the module docstring).  Raises :class:`BuildError`.
    """
    from rdkit import Chem

    m = Chem.MolFromSmiles(lig["smiles"])
    if m is None:
        raise BuildError(f"ligand {lig['smiles']} cannot be read by the RDKit of epic-MACE")
    try:
        Chem.Kekulize(m, clearAromaticFlags=True)
    except Exception as exc:
        raise BuildError(f"ligand {lig['smiles']} cannot be kekulized: {exc}") from exc
    for a in m.GetAtoms():          # the hydrogens are final
        a.SetNumExplicitHs(a.GetTotalNumHs())
        a.SetNoImplicit(True)
    rw = Chem.RWMol(m)
    sites = []
    for g in _mace_donor_groups(m, lig["coordList"]):
        if len(g) == 1:
            rw.GetAtomWithIdx(g[0]).SetAtomMapNum(1)
            sites.append({"kind": "atom", "elem": m.GetAtomWithIdx(g[0]).GetSymbol()})
            continue
        ring = _mace_carbon_ring(m, g)
        star = rw.AddAtom(Chem.Atom(0))
        rw.GetAtomWithIdx(star).SetAtomMapNum(1)
        if ring is not None:
            rw.AddBond(star, _mace_ring_anchor(rw, ring), Chem.BondType.SINGLE)
            sites.append({"kind": "hapto", "eta": len(g), "enc": "anchor"})
        else:
            gs = set(g)
            for b in list(m.GetBonds()):
                i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                if i in gs and j in gs:
                    rw.RemoveBond(i, j)
            for i in g:
                rw.AddBond(star, i, Chem.BondType.SINGLE)
            sites.append({"kind": "hapto", "eta": len(g), "enc": "star"})
    mol = rw.GetMol()
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol)
    except Exception as exc:
        raise BuildError(f"ligand {lig['smiles']} cannot be written for epic-MACE: "
                         f"{str(exc)[:120]}") from exc
    return Chem.MolToSmiles(mol), sites


def mace_job(spec: dict, geometries: str = "extended") -> dict:
    """The epic-MACE input of a spec: ligands, central atom, geometries.

    ``{'geoms', 'ligands', 'CA', 'info'}``.  Raises :class:`BuildError` when no
    epic-MACE geometry has the number of sites.
    """
    if geometries not in MACE_GEOMS:
        raise ValueError(f"unknown epic-MACE geometry set {geometries!r} "
                         f"(known: {', '.join(sorted(MACE_GEOMS))})")
    ligs, sites = [], []
    for lig in spec["ligands"]:
        s, info = mace_ligand(lig)
        ligs.append(s)
        sites += info
    n_sites = len(sites)
    hapto = [(s["eta"], s["enc"]) for s in sites if s["kind"] == "hapto"]
    ox = spec.get("metal_ox")
    ca = (f"[{spec['metal']}+{ox}]" if isinstance(ox, int) and 0 < ox <= 8
          else f"[{spec['metal']}]")
    geoms = list(MACE_GEOMS[geometries].get(n_sites, []))
    if geometries == "extended" and n_sites == 2 and len(hapto) == 2:
        geoms = ["SAN"]
    if not geoms:
        have = ", ".join(f"{k} ({'/'.join(v)})" for k, v in sorted(MACE_GEOMS[geometries].items()))
        raise BuildError(f"epic-MACE has no geometry for {n_sites} donor sites "
                         f"(a hapto ligand is one site; it has: {have}"
                         + ("; two hapto sites: SAN" if geometries == "extended" else "") + ")")
    return {"geoms": geoms, "ligands": ligs, "CA": ca,
            "info": {"n_sites": n_sites, "hapto": hapto}}


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


def _run_mace(spec: dict, work: str, options: dict) -> tuple:
    import mace

    opts = dict(MACE_DEFAULTS, **(options or {}))
    job = mace_job(spec, opts["geometries"])
    frames, per_geo = [], {}
    for geom in job["geoms"]:
        try:
            X = mace.ComplexFromLigands(job["ligands"], job["CA"], geom)
            # As MACE's own command line does before a stereomer search: every
            # donor mapped (and marked) alike, the complex rebuilt from that.
            for idx in X._DAs:
                X.mol.GetAtomWithIdx(idx).SetAtomMapNum(1)
                X.mol.GetAtomWithIdx(idx).SetIsotope(1)
            X = mace.ComplexFromMol(X.mol, X.geom)
            stereomers = X.GetStereomers("all", False, None, False)
        except Exception as exc:
            text = str(exc).strip().splitlines()
            per_geo[geom] = f"{type(exc).__name__}: {text[0][:160] if text else ''}"
            continue
        n_before = len(frames)
        for k, x in enumerate(stereomers):
            try:
                x.AddConformers(numConfs=int(opts["num_confs"]),
                                maxAttempts=int(opts["max_attempts"]), rmsThresh=-1)
            except Exception:
                continue
            if not x.GetNumConformers():
                continue
            x.OrderConfsByEnergy()
            for j in range(x.GetNumConformers()):
                block = x.ToXYZBlock(j).splitlines()
                n = int(block[0])
                try:
                    e = json.loads(block[1]).get("E")
                except ValueError:
                    e = None
                syms, xyz = [], []
                for row in block[2:2 + n]:
                    p = row.split()
                    if p[0] == "X":      # a hapto centroid, not an atom
                        continue
                    syms.append(p[0])
                    xyz.append([float(p[1]), float(p[2]), float(p[3])])
                frames.append({"label": f"{geom}-iso{k}-conf{j}", "symbols": syms,
                               "coords": xyz, "energy": float(e) if e is not None else None})
        mine = frames[n_before:]
        mine.sort(key=lambda f: (f["energy"] is None, f["energy"] or 0.0))
        frames[n_before:] = mine
        per_geo[geom] = (f"{len(stereomers)} stereomers, {len(frames) - n_before} frames"
                         if len(frames) > n_before else
                         f"{len(stereomers)} stereomers, no conformer embedded")
    return frames, {"per_geometry": per_geo, "mace_input": {k: job[k] for k in ("ligands", "CA")}}


def _is_the_tool(key: str) -> bool:
    """Whether the importable module is the tool and not a namesake.

    ``import mace`` is also the MACE machine-learning potential (mace-torch),
    which DELFIN installs into its own environment; epic-MACE is the one with
    ``ComplexFromLigands``.
    """
    if key != "mace":
        return True
    try:
        import mace
    except Exception:
        return False
    return hasattr(mace, "ComplexFromLigands")


_RUNNERS = {"architector": _run_architector, "molsimplify": _run_molsimplify,
            "mace": _run_mace}


def _worker(tool: str, request_path: str, result_path: str) -> int:
    # rdkit before any tool: see the module docstring (libstdc++).
    result = {"frames": [], "error": None}
    try:
        from rdkit import RDLogger  # noqa: F401
        RDLogger.DisableLog("rdApp.*")
        req = json.loads(Path(request_path).read_text())
        key = tool_key(tool)
        import importlib.util
        if importlib.util.find_spec(TOOLS[key]["module"]) is None or not _is_the_tool(key):
            result["error"] = not_installed_message(key, sys.executable)
            result["not_installed"] = True
        else:
            # A spec cut by the caller (spec_in_caller) is taken as it is.
            spec = req.get("spec") or split_complex_smiles(req["smiles"])
            result["spec"] = spec
            frames, extra = _RUNNERS[key](spec, req["work"], req.get("options") or {})
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
    Architector's order anyway), else in the tool's order (epic-MACE: geometry
    by geometry, lowest energy first within each).  *error* is ``None``
    on success, else a message for the user.  *options* for Architector:
    ``{'mode': 'full'}`` (default: every isomer it finds) or ``'default'``
    (its own defaults, one structure per core geometry).  For epic-MACE:
    :data:`MACE_DEFAULTS` (``geometries`` ``'extended'`` or ``'paper'``,
    ``num_confs``, ``max_attempts``).
    """
    key = tool_key(tool)
    info = TOOLS[key]
    python = python or tool_python(key, environ)
    if python == sys.executable:
        import importlib.util
        if importlib.util.find_spec(info["module"]) is None:
            return [], not_installed_message(key, python)
    spec = None
    if info.get("spec_in_caller"):
        # Cut here, with DELFIN's RDKit; the tool's environment may carry an
        # RDKit years older (epic-MACE: 2020.09) that only has to write it out.
        try:
            spec = split_complex_smiles(smiles)
        except BuildError as exc:
            return [], f"{info['display']}: {exc}"
    import delfin

    repo_root = str(Path(delfin.__file__).resolve().parent.parent)
    env = dict(os.environ if environ is None else environ)
    if info.get("managed_env"):
        # Another Python version: nothing of this one's paths may reach it.
        env["PYTHONPATH"] = repo_root
        env["PYTHONNOUSERSITE"] = "1"
    else:
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
        {"smiles": str(smiles).strip(), "work": work, "options": dict(options or {}),
         "spec": spec}))
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
        if (frames and not info.get("keep_order")
                and all(f.get("energy") is not None for f in frames)):
            frames = sorted(frames, key=lambda f: f["energy"])
        out = []
        for i, f in enumerate(frames, start=1):
            e = f.get("energy")
            label = f"{info['display']} {f.get('label') or i}"
            if e is not None:
                # Architector: its xTB energy in ASE units; epic-MACE: UFF.
                label += f" (E = {e:.4f} {info.get('energy_unit', 'eV')})"
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
