#!/usr/bin/env python3
"""Generate the two task pools for a specialist-versus-generalist comparison.

Input: a source tree (default: the repository this file sits in).
Output: JSON on ``--out`` plus a printed breakdown; exit status 0 even when a
pool is underpowered, because an underpowered pool is a result about the
corpus and not an error of this script.

Semantics
---------
The corpus is produced by a query, never by a person choosing tasks: whoever
chooses knows which domain the specialist covers, which confounds arm with
domain. One work type is held fixed and only the domain varies.

The work type: *a module-level callable in the product that nothing in the
product reaches.* Two tiers, as the reachability test defines them:

* ``nowhere``   -- nothing under ``delfin/`` or ``tests/`` names it;
* ``test-only`` -- only the suite names it.

The task sentence attached to every candidate is identical across pools and
arms: "connect this capability to the product, or delete it together with
its tests."

Two pools, so that the domain varies and the work type does not:

* ``infra`` -- defined under ``delfin/agent/`` and not in a chemistry file.
  Outside every chemistry specialist's domain.
* ``chem``  -- defined under ``delfin/tools/`` (all of it, by the brief), or
  in a chemistry file anywhere else under ``delfin/``.

Re-use of the scanning logic
----------------------------
The scanning functions are **copied**, not imported, from
``tests/test_every_capability_has_a_caller.py`` as of blob
``f290b9753151c5c9b323865bf37b96e1c09bcf66`` on branch
``commit/a-commit-records-the-model-that-wrote-it``. Copied for three
reasons, each of which blocks an import:

1. that file does not exist on ``origin/main``, which is the tree this
   experiment runs against, so an import would fail on the experiment tree
   and succeed only on one branch;
2. its ``_definitions`` hard-codes ``delfin/agent`` as the definition root,
   and this script needs the same question asked of three other roots;
3. it is a pytest module with a module-scoped fixture, so importing it
   drags pytest collection into a measurement script.

The copies are faithful: ``_definitions`` takes its root as an argument and
accepts a file filter, and ``_mentions``/``_count_dotted_dispatch`` are
unchanged apart from the baseline-file skip, which is expressed as a set of
paths because the baseline lives on another branch. Any behavioural
divergence between copy and original is a defect of this script;
``tests/test_ab_apparatus.py`` pins the pair whenever both are present in
one tree.

Neither that test nor its fixture is modified or written by this script.

The chemistry rule
------------------
Classification is per **file**, not per package, because a pool entry is a
definition at a file and a line, and because a package verdict mislabels the
whole of ``delfin/agent/`` (609 files) on the strength of the 5 of them that
import a chemistry library.

A file under ``delfin/`` is a *chemistry file* when either signal holds:

* **S1, library**: it imports a third-party library whose own subject is
  chemistry (``CHEM_LIBRARIES``). Objective: the list names packages, not
  modules of this repository.
* **S2, external vocabulary**: it names third-party chemistry *software*,
  standard chemistry *file formats*, chemistry *notations* or published
  *method* names (``CHEM_EXTERNAL``) at least ``S2_MIN_HITS`` times and at a
  rate of at least ``S2_MIN_PER_1K`` per thousand words. The absolute count
  travels with the rate on purpose: a rate alone promoted a 20-word banner
  module.

Every term in ``CHEM_EXTERNAL`` is the name of something that exists outside
this repository and can be checked against its own documentation. Generic
words ("molecule", "solvent", "functional") are deliberately absent: they
occur in user-facing strings all over the product and made the rule select
modules with no chemistry in them.

A derived vocabulary was measured as an alternative -- words enriched in the
S1 files against the rest of the tree -- and rejected: the S1 set contains
the dashboard widgets, so the top of the derived list is
``input_textarea``, ``model_dropdown``, ``canvas``, ``mouse``.

Exclusions
----------
Applied in the order below, each counted and reported:

1. ``documented_only`` -- tier ``nowhere`` and the name occurs under
   ``docs/`` or ``test-evidence/``. A documented entry point is a reference,
   not a caller, and it is also the shape that makes deletion wrong.
2. ``defined_in_test`` -- the definition itself sits in a test file.
3. ``module_absent`` -- the defining file does not exist in the comparison
   tree (``--main-ref``).
4. ``short_name`` -- fewer than 4 characters, so a grep for it is ambiguous.

A fifth exclusion was written and then dropped: a name defined in both
pools, which no name-based scan could attribute to one domain. It is
structurally impossible. The mention count skips only the definition line of
the pool's own definition, so a second module-level definition of the same
name anywhere under ``delfin/`` counts as a mention and the name is never a
candidate. Every candidate therefore has exactly one module-level
definition under ``delfin/``. The field ``candidates_with_a_twin`` reports
that as a number instead of asserting it, because it is the kind of
assumption that stops holding when the scan changes.

The same property is a limitation of the measure rather than of this
script: two unreachable twins hide each other, in the original instrument
as well as here.
"""
from __future__ import annotations

import argparse
import ast
import datetime as dt
import io
import json
import re
import subprocess
import tokenize
from collections import Counter
from pathlib import Path

# --------------------------------------------------------------------------
# Copied from tests/test_every_capability_has_a_caller.py (see the docstring)
# --------------------------------------------------------------------------

#: Not source text. Reading these as text would spray random identifiers
#: into the reference count and hide real orphans.
_BINARY = frozenset({
    ".png", ".jpg", ".jpeg", ".gif", ".pdf", ".ico", ".woff", ".woff2",
    ".zip", ".gz", ".so", ".pyc", ".xlsx", ".docx", ".pptx", ".odt", ".ods",
})

_WORD = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def _definitions(root: Path, tree_root: Path,
                 keep=None) -> dict[str, tuple[str, int]]:
    """``{name: (file, line)}`` for every module-level callable under ``root``.

    Module level only: a method is reached through its class and asking the
    same question of it would need the class's own reachability first.
    ``keep`` filters the files, and paths are reported relative to
    ``tree_root``.
    """
    out: dict[str, tuple[str, int]] = {}
    for path in sorted(root.rglob("*.py")):
        if keep is not None and not keep(path):
            continue
        try:
            tree = ast.parse(path.read_text(errors="replace"))
        except SyntaxError:
            continue
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) \
                    and not node.name.startswith("__"):
                out.setdefault(node.name,
                               (str(path.relative_to(tree_root)), node.lineno))
    return out


def _mentions(tree_root: Path, defs: dict, skip_own_definition: bool,
              root: Path, skip_paths=frozenset()) -> dict:
    """How often each name is mentioned under ``tree_root``.

    Python files are tokenised, so a name in a comment does not count while
    a name in a string literal does -- that is name-based dispatch, which is
    a real reference. Everything else readable as text is scanned for words,
    which is how a module named from YAML is found. With
    ``skip_own_definition`` the ``def`` line itself does not count as a
    mention of the thing it defines. ``skip_paths`` holds files that list
    the names being looked for; scanning one makes the instrument measure
    itself.
    """
    hits = {name: 0 for name in defs}
    for path in sorted(tree_root.rglob("*")):
        if not path.is_file() or path.suffix.lower() in _BINARY:
            continue
        if path in skip_paths:
            continue
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        rel = str(path.relative_to(root))
        if path.suffix == ".py":
            try:
                for tok in tokenize.generate_tokens(io.StringIO(text).readline):
                    if tok.type == tokenize.NAME and tok.string in hits:
                        name = tok.string
                        if (skip_own_definition and rel == defs[name][0]
                                and tok.start[0] == defs[name][1]):
                            continue
                        hits[name] += 1
                    elif tok.type == tokenize.STRING:
                        for word in _WORD.findall(tok.string):
                            if word in hits:
                                hits[word] += 1
                continue
            except (tokenize.TokenError, IndentationError, SyntaxError):
                pass        # unparseable: fall through to the word scan
        for word in _WORD.findall(text):
            if word in hits:
                hits[word] += 1
    _count_dotted_dispatch(tree_root, hits, skip_paths)
    return hits


def _count_dotted_dispatch(tree_root: Path, hits: dict,
                           skip_paths=frozenset()) -> None:
    """Also count a name spelled with dots instead of underscores.

    A dispatcher that does ``getattr(obj, method.replace(".", "_"))`` means
    the protocol declares ``"swarm.observe.feed"`` while the underscored
    function name appears nowhere. This compares literal text, because a
    grep that leaves the dot unescaped matches every name against itself.
    """
    dotted = {n.replace("_", "."): n for n in hits if "_" in n}
    if not dotted:
        return
    runs = re.compile(r"[A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)+")
    for path in sorted(tree_root.rglob("*")):
        if not path.is_file() or path.suffix.lower() in _BINARY:
            continue
        if path in skip_paths:
            continue
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        for run in set(runs.findall(text)):
            name = dotted.get(run)
            if name:
                hits[name] += 1


# --------------------------------------------------------------------------
# The chemistry rule
# --------------------------------------------------------------------------

#: Third-party packages whose own subject is chemistry. Signal S1.
CHEM_LIBRARIES = (
    "rdkit", "ase", "opi", "openbabel", "pybel", "pyscf", "cclib",
    "pymatgen", "spglib", "openmm", "MDAnalysis", "xtb", "morfeus",
    "geometric", "py3Dmol", "nglview", "dscribe", "mace", "nequip",
)

#: Names of chemistry software, file formats, notations and published
#: methods -- each checkable outside this repository. Signal S2.
CHEM_EXTERNAL = (
    # programs
    "orca", "xtb", "crest", "turbomole", "gaussian", "molpro", "psi4",
    "nwchem", "gamess", "cp2k", "vasp", "gromacs", "openmm", "goat",
    # notations and identifiers
    "smiles", "smarts", "inchi",
    # methods, functionals, basis sets
    "b3lyp", "pbe0", "wb97", "tpss", "cam_b3lyp", "def2", "cc_pvdz",
    "ccsd", "mp2", "nevpt2", "casscf", "caspt2", "dlpno", "gfn2", "gfn_ff",
    "tddft", "dft", "nbo", "nics",
    # implicit solvation models
    "cpcm", "smd", "cosmo",
    # file formats
    "xyz", "mol2", "sdf", "cif", "pdb", "gbw", "hess", "molden",
    # orbital and spectroscopy terms with no non-chemical reading
    "homo", "lumo",
)

S2_MIN_HITS = 25
S2_MIN_PER_1K = 2.0

_CHEM_IMPORT = re.compile(
    r"^\s*(?:from|import)\s+(?:" + "|".join(CHEM_LIBRARIES) + r")\b",
    re.MULTILINE)
_CHEM_EXTERNAL = re.compile(
    r"\b(?:" + "|".join(CHEM_EXTERNAL) + r")\b", re.IGNORECASE)


def chemistry_signals(text: str) -> dict:
    """The two signals for one file's text, counts beside the rate.

    Output: ``{"s1": bool, "hits": int, "words": int, "per_1k": float,
    "s2": bool, "chemistry": bool}``.
    """
    hits = len(_CHEM_EXTERNAL.findall(text))
    words = max(len(text.split()), 1)
    per_1k = hits / words * 1000
    s1 = bool(_CHEM_IMPORT.search(text))
    s2 = hits >= S2_MIN_HITS and per_1k >= S2_MIN_PER_1K
    return {"s1": s1, "hits": hits, "words": words,
            "per_1k": round(per_1k, 2), "s2": s2, "chemistry": s1 or s2}


def is_chemistry_file(path: Path) -> bool:
    """True when ``path`` is a chemistry file by the rule in the docstring."""
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return False
    return chemistry_signals(text)["chemistry"]


# --------------------------------------------------------------------------
# Pools
# --------------------------------------------------------------------------

def _is_test_file(rel: str) -> bool:
    name = Path(rel).name
    return (name.startswith("test_") or name.endswith("_test.py")
            or "/tests/" in "/" + rel)


def _module_of(rel: str) -> str:
    return rel[:-3].replace("/", ".") if rel.endswith(".py") else rel


def _package_of(rel: str) -> str:
    return str(Path(rel).parent)


def _chem_file_cache(root: Path) -> dict[Path, bool]:
    """``{path: is_chemistry}`` for every ``.py`` file under ``delfin/``."""
    out: dict[Path, bool] = {}
    for path in sorted((root / "delfin").rglob("*.py")):
        out[path] = is_chemistry_file(path)
    return out


def pool_roots(root: Path, chem: dict[Path, bool]):
    """The file filters for the two pools.

    ``infra``: under ``delfin/agent/`` and not a chemistry file.
    ``chem``: all of ``delfin/tools/``, plus any chemistry file elsewhere
    under ``delfin/`` (``delfin/agent/`` included, which is how a chemistry
    file in the agent tree leaves ``infra`` and is not simply dropped).
    """
    agent = (root / "delfin" / "agent").resolve()
    tools = (root / "delfin" / "tools").resolve()

    def infra(path: Path) -> bool:
        return not chem.get(path, False)

    def chemistry(path: Path) -> bool:
        p = path.resolve()
        if tools in p.parents:
            return True
        return chem.get(path, False)

    return {
        "infra": (agent, infra),
        "chem": ((root / "delfin").resolve(), chemistry),
    }


def _twin_definitions(root: Path) -> set[str]:
    """Names with more than one module-level definition under ``delfin/``.

    Output: the set of such names. A candidate may not be one of them, and
    the count of those that are is reported rather than assumed away.
    """
    seen = Counter()
    for path in sorted((root / "delfin").rglob("*.py")):
        try:
            tree = ast.parse(path.read_text(errors="replace"))
        except SyntaxError:
            continue
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) \
                    and not node.name.startswith("__"):
                seen[node.name] += 1
    return {n for n, c in seen.items() if c > 1}


def _head_sha(root: Path) -> str | None:
    """HEAD of ``root``, or None when it is not a readable git tree.

    Recorded in the document so that the draw reports the sha the scan ran
    on, not the sha the tree happens to be at when it is drawn from.
    """
    try:
        out = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() or None if out.returncode == 0 else None


def _tracked_paths(root: Path, ref: str) -> set[str] | None:
    """Paths in ``ref``, or None when the ref cannot be read."""
    try:
        out = subprocess.run(["git", "-C", str(root), "ls-tree", "-r",
                              "--name-only", ref],
                             capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0:
        return None
    return set(out.stdout.split())


def build_pools(root: Path, main_ref: str = "origin/main") -> dict:
    """Both pools with their candidates, tiers, exclusions and counts.

    Input: a source tree and a git ref to check module existence against.
    Output: the JSON document described in ``--out``.
    """
    root = root.resolve()
    chem_cache = _chem_file_cache(root)
    roots = pool_roots(root, chem_cache)
    # The baseline file lists every name the scan looks for. Scanning it
    # makes every listed name read as reached and empties the sharp tier.
    skip = {(root / "tests" / "fixtures" / "unmounted_capabilities.txt")}

    raw: dict[str, dict[str, dict]] = {}
    for pool, (def_root, keep) in roots.items():
        defs = _definitions(def_root, root, keep)
        in_product = _mentions(root / "delfin", defs, True, root, skip)
        unmounted = sorted(n for n in defs if in_product[n] == 0)
        in_suite = _mentions(root / "tests", defs, False, root, skip)
        raw[pool] = {
            n: {
                "name": n,
                "tier": "test-only" if in_suite.get(n) else "nowhere",
                "file": defs[n][0],
                "line": defs[n][1],
                "module": _module_of(defs[n][0]),
                "package": _package_of(defs[n][0]),
            }
            for n in unmounted
        }
        raw[pool + "__ndefs"] = len(defs)      # type: ignore[assignment]

    # Names mentioned in prose. Needed for exclusion 1, and the scan does
    # not read these trees, so a documented name is in tier `nowhere`.
    prose_names: set[str] = set()
    wanted = set(raw["infra"]) | set(raw["chem"])
    for sub in ("docs", "test-evidence"):
        tree = root / sub
        if not tree.is_dir():
            continue
        for path in sorted(tree.rglob("*")):
            if not path.is_file() or path.suffix.lower() in _BINARY:
                continue
            try:
                text = path.read_text(errors="replace")
            except OSError:
                continue
            for word in _WORD.findall(text):
                if word in wanted:
                    prose_names.add(word)

    tracked = _tracked_paths(root, main_ref)
    twins = _twin_definitions(root)

    removed = Counter()
    removed_names: dict[str, list[str]] = {k: [] for k in (
        "documented_only", "defined_in_test", "module_absent",
        "short_name")}
    pools: dict[str, list[dict]] = {}
    for pool in ("infra", "chem"):
        kept = []
        for name, row in sorted(raw[pool].items()):
            reason = None
            if row["tier"] == "nowhere" and name in prose_names:
                reason = "documented_only"
            elif _is_test_file(row["file"]):
                reason = "defined_in_test"
            elif tracked is not None and row["file"] not in tracked:
                reason = "module_absent"
            elif len(name) < 4:
                reason = "short_name"
            if reason:
                removed[f"{pool}:{reason}"] += 1
                removed[reason] += 1
                removed_names[reason].append(f"{pool}:{name}")
                continue
            kept.append(row)
        pools[pool] = kept

    return {
        "generated_for": str(root),
        "tree_sha": _head_sha(root),
        "scanned_at": dt.datetime.now(dt.timezone.utc).isoformat(
            timespec="seconds"),
        "main_ref": main_ref,
        "main_ref_readable": tracked is not None,
        "chemistry_rule": {
            "granularity": "file",
            "libraries": list(CHEM_LIBRARIES),
            "external_terms": list(CHEM_EXTERNAL),
            "s2_min_hits": S2_MIN_HITS,
            "s2_min_per_1k": S2_MIN_PER_1K,
            "chemistry_files_under_delfin":
                sum(1 for v in chem_cache.values() if v),
            "py_files_under_delfin": len(chem_cache),
            "chemistry_files_under_delfin_agent": sum(
                1 for p, v in chem_cache.items()
                if v and "delfin/agent/" in str(p).replace("\\", "/")),
        },
        "task_sentence": ("connect this capability to the product, or delete "
                          "it together with its tests"),
        "callables_scanned": {p: raw[p + "__ndefs"] for p in ("infra", "chem")},
        "candidates_with_a_twin": {
            p: sorted(r["name"] for r in pools[p] if r["name"] in twins)
            for p in pools},
        "exclusions_removed": dict(removed),
        "exclusions_removed_names": removed_names,
        "pools": pools,
        "counts": {p: len(pools[p]) for p in pools},
        "counts_by_tier": {
            p: dict(Counter(r["tier"] for r in pools[p])) for p in pools},
        "counts_by_package": {
            p: dict(sorted(Counter(r["package"] for r in pools[p]).items()))
            for p in pools},
        "underpowered": {p: len(pools[p]) < 10 for p in pools},
    }


def report(doc: dict) -> None:
    """Print the counts, the per-package breakdown and the exclusions."""
    rule = doc["chemistry_rule"]
    print(f"chemistry rule: {rule['chemistry_files_under_delfin']} of "
          f"{rule['py_files_under_delfin']} .py files under delfin/ are "
          f"chemistry files "
          f"({rule['chemistry_files_under_delfin_agent']} of them under "
          "delfin/agent/, which therefore leave the infra pool)")
    print(f"callables scanned: {doc['callables_scanned']}")
    print()
    for pool in ("infra", "chem"):
        rows = doc["pools"][pool]
        print(f"pool {pool}: {len(rows)} candidates "
              f"{doc['counts_by_tier'][pool]}")
        for pkg, n in doc["counts_by_package"][pool].items():
            print(f"    {n:4}  {pkg}")
    print()
    print("exclusions removed:")
    for key in ("documented_only", "defined_in_test", "module_absent",
                "short_name"):
        print(f"    {doc['exclusions_removed'].get(key, 0):4}  {key}"
              f"   (infra {doc['exclusions_removed'].get(f'infra:{key}', 0)},"
              f" chem {doc['exclusions_removed'].get(f'chem:{key}', 0)})")
    twins = {p: v for p, v in doc["candidates_with_a_twin"].items() if v}
    print(f"    candidates with a second module-level definition: "
          f"{twins or 'none'}")
    if not doc["main_ref_readable"]:
        print(f"    note: {doc['main_ref']} could not be read, so "
              "module_absent removed nothing by construction")
    print()
    for pool, under in doc["underpowered"].items():
        if under:
            print(f"UNDERPOWERED: pool {pool} has "
                  f"{doc['counts'][pool]} candidates, fewer than 10")
    if not any(doc["underpowered"].values()):
        print("neither pool is underpowered (both at 10 or more)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[1]),
                    help="source tree to scan")
    ap.add_argument("--main-ref", default="origin/main",
                    help="ref whose file list decides exclusion module_absent")
    ap.add_argument("--out", default=None, help="write the JSON document here")
    args = ap.parse_args(argv)

    doc = build_pools(Path(args.root), args.main_ref)
    if args.out:
        Path(args.out).write_text(json.dumps(doc, indent=2, sort_keys=True))
        print(f"wrote {args.out}")
    report(doc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
