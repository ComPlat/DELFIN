"""Grade a workflow SPEC an agent built, before anything runs (package WB0).

The wave-21 goal is that agents build their OWN physically realistic
workflows. What was missing is the measurement: a function that takes the
spec the agent produced, the chemist's request it answers, and the answer
text the agent wrote, and returns a verdict — without any calculation.

The grader owns NO chemistry rules of its own. Every check routes through
the real platform call sites:

- :func:`delfin.tools.platform.resolve_spec` (:256) — the minimal spec must
  resolve to a complete one;
- :func:`delfin.tools.platform.validate_spec` (:368) — structural
  validation (unknown steps, missing required params, port wiring);
- the science lint via :func:`delfin.tools.platform.scientific_lint` (:239,
  rules in :mod:`delfin.tools._science`, extended by other packages) —
  every ``warning`` must be justified in the answer text; only findings the
  answer justifies are tolerated, a ``error``-level finding never is;
- :meth:`delfin.tools.pipeline.Pipeline.dry_run` (:1648) — the resolved
  pipeline must build and walk without executing anything.

On top of the platform checks the grader adds the two physics rules that
the structural validator and the current lint rule set cannot see, because
they live ACROSS steps or in the request, not inside one step's kwargs:

- one energy difference, one level of theory: two (or more) steps compared
  against each other must agree on method AND basis;
- request coverage: every capability the request asks for (optimise,
  frequencies, conformers, spectrum, spin states, ...) must be covered by
  a step.

And the two resource rules the wave treats as correctness:

- a transition state must be verified: a frequency check (exactly one
  imaginary mode) or an IRC must follow the TS optimisation;
- charged open-shell species need the basis matched to the physics: an
  anion needs diffuse functions, a 5d metal needs a relativistic
  Hamiltonian, a flexible molecule asked for its conformer landscape needs
  a conformer step before the expensive level.

The cost block is an ESTIMATE computed from the resolved spec (documented
as such in the verdict): the grader has no access to the scheduler, and a
run's true cost is EF1's accounting (:mod:`delfin.agent.cluster_efficiency`)
which exists only after the run.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

__all__ = ["grade"]

# Rough per-step core-hour weights for the budget ESTIMATE. These are
# order-of-magnitude anchors for <= 20-atom systems at def2-SVP..TZVP on
# 4 cores, not measured costs; the verdict marks them ``estimate`` and the
# real-run tasks check the actuals with EF1's accounting.
_CORE_HOUR_ESTIMATE = {
    "smiles_to_xyz": 0.0,
    "crest_conformers": 0.05,
    "xtb_opt": 0.02,
    "xtb_goat": 0.1,
    "orca_sp": 0.1,
    "orca_opt": 0.3,
    "orca_freq": 0.3,
    "orca_opt_freq": 0.5,
    "orca_ts_opt": 1.0,
    "orca_ts_freq_check": 0.4,
    "orca_ts_irc": 0.8,
    "orca_tddft": 0.5,
    "cclib_parse": 0.0,
    "cclib_energies": 0.0,
    "cclib_vibrations": 0.0,
    "imag_fix": 0.2,
    "uv_vis_parse": 0.0,
    "ir_parse": 0.0,
    "xtb_solvator": 0.05,
    "esd": 0.3,
}

#: Request-coverage vocabulary. Each entry maps a keyword family in the
#: chemist's request to the step names that satisfy it. A request keyword
#: nothing in the spec matches fails the coverage check.
_COVERAGE = (
    (("optimise", "optimize", "optimization", "optimisation"), ("orca_opt", "orca_opt_freq", "xtb_opt", "orca_ts_opt")),
    (("frequency", "frequencies", "thermochem", "thermochemical", "minimum", "saddle"), ("orca_freq", "orca_opt_freq", "cclib_vibrations", "ir_parse")),
    (("conformer", "conformers", "flexible", "rotamer"), ("crest_conformers", "xtb_goat")),
    # ir_parse satisfies a spectrum request only when the request names
    # IR; a UV/visible spectrum needs the TDDFT chain (s25 A2 review run).
    (("ir spectrum", "infrared"), ("ir_parse", "orca_opt_freq")),
    (("spectrum", "excitation", "absorption", "optically", "bright"), ("orca_tddft", "uv_vis_parse", "esd", "ir_parse")),
    # a UV/visible request needs the TDDFT/ESD chain; ir_parse alone
    # cannot deliver it (the IR entry above covers the infrared side).
    (("uv", "uv-vis"), ("orca_tddft", "uv_vis_parse", "esd")),
    (("spin",), ("orca_sp", "orca_opt", "orca_opt_freq")),
    (("potential", "redox", "ionisation", "ionization", "oxidation", "reduction"), ("orca_sp", "cclib_energies")),
    (("solvation", "solvent", "solvated", "solvate", "medium"), ("xtb_solvator", "esd")),
    (("transition state", "ts search", "reaction path", "irc", "barrier"), ("orca_ts_opt", "orca_ts_irc", "orca_ts_freq_check", "imag_fix")),
)

# Level-of-theory normalisation for the "one difference, one level" rule.
_LEVEL_CANON = {"pbe0": "pbe0", "b3lyp": "b3lyp", "b3lyp-g": "b3lyp",
                "r2scan-3c": "r2scan-3c", "pbeh-3c": "pbeh-3c",
                "wb97x-3c": "wb97x-3c", "hf-3c": "hf-3c",
                "b97-3c": "b97-3c", "tpss": "tpss", "pbe": "pbe"}
_BASIS_CANON = {"def2-svp": "def2-svp", "def2-tzvp": "def2-tzvp",
                "def2-tzvpp": "def2-tzvpp", "def2-qzvp": "def2-qzvp",
                "cc-pvtz": "cc-pvtz", "cc-pvdz": "cc-pvdz",
                "ma-def2-tzvp": "ma-def2-tzvp", "6-31g*": "6-31g*"}

# Diffuse-function markers an anion's basis must carry (wave-21 target 5).
_DIFFUSE_MARKERS = ("aug-", "ma-", "def2-tzvpd", "def2-tzvpdd", "6-31+g",
                    "cc-pvtz-d", "sapporo", "cdnz", "aug-cc")

# 5d metals whose chemistry is meaningless without scalar relativistic
# effects (ZORA/DKH/X2C) — a plain non-relativistic Hamiltonian is a
# physics defect, not a style choice. (symbol, element name) pairs so a
# chemist's request ("a gold-containing complex") is caught like a
# formula ("Au").
_5D_METALS = (
    ("au", "gold"), ("pt", "platinum"), ("ir", "iridium"),
    ("os", "osmium"), ("re", "rhenium"), ("w", "tungsten"),
    ("ta", "tantalum"), ("hf", "hafnium"), ("hg", "mercury"),
    ("tl", "thallium"), ("pb", "lead"), ("bi", "bismuth"),
)
#: (case-sensitive symbol, case-insensitive element name) — the symbol
#: match is case-sensitive so 'IR' (infrared) never reads as 'Ir'
#: (iridium) or 'HF' (hydrogen fluoride) as 'Hf' (hafnium). Symbols are
#: stored capitalised; the element name in lowercase.
_5D_SYMBOLS = tuple(
    (sym.capitalize(), name) for sym, name in _5D_METALS)

_RELATIVISTIC = ("zora", "dkh", "x2c", "sarc")


def _norm(value: Any) -> str:
    return str(value or "").strip().lower()


def _level_of(step: dict) -> Optional[tuple]:
    """(method, basis) of an orca step, else None."""
    name = _norm(step.get("step"))
    if not name.startswith("orca"):
        return None
    kw = step.get("kwargs", {}) or {}
    method, basis = _norm(kw.get("method")), _norm(kw.get("basis"))
    if not method and not basis:
        return None
    return (_LEVEL_CANON.get(method, method or ""),
            _BASIS_CANON.get(basis, basis or ""))


def _steps_of(spec: dict) -> List[dict]:
    steps = list(spec.get("steps") or [])
    for branch in (spec.get("branches") or {}).values():
        steps.extend(branch.get("steps") if isinstance(branch, dict) else branch or [])
    return steps


# ---------------------------------------------------------------------------
# The physics rules that live across steps or in the request
# ---------------------------------------------------------------------------

def _mixed_level_reason(steps: List[dict]) -> Optional[str]:
    """More than one distinct level of theory inside one spec that compares
    energies is a defect: a total energy compares only within ONE method."""
    levels = {lvl for s in steps if (lvl := _level_of(s))}
    if len(levels) > 1:
        shown = sorted(f"{m}/{b}" for m, b in levels)
        return (f"one energy difference, one level of theory: the spec mixes "
                f"{shown} — all compared steps must share method AND basis")
    return None


def _unverified_ts_reason(steps: List[dict]) -> Optional[str]:
    """A TS optimisation must be followed by a frequency check (one
    imaginary mode) or an IRC — an unverified stationary point is not a
    transition state."""
    if not any(_norm(s.get("step")) == "orca_ts_opt" for s in steps):
        return None
    ok = any(_norm(s.get("step")) in ("orca_ts_freq_check", "orca_ts_irc")
             for s in steps)
    if not ok:
        return ("transition state without verification: orca_ts_opt needs a "
                "following orca_ts_freq_check (exactly one imaginary mode) "
                "or an orca_ts_irc")
    return None


def _anion_reason(spec: dict, steps: List[dict]) -> Optional[str]:
    """An anion needs diffuse functions: the extra electron is bound
    weakly and a plain def2-SVP/TZVP tail cannot describe it."""
    for s in steps:
        kw = s.get("kwargs", {}) or {}
        if int(kw.get("charge") or 0) < 0:
            basis = _norm(kw.get("basis"))
            if basis and not any(m in basis for m in _DIFFUSE_MARKERS):
                return (f"anion without diffuse functions: basis '{basis}' "
                        f"cannot bind the extra electron; use an augmented "
                        f"basis (aug-/ma- or a -d diffuse variant)")
    return None


def _relativistic_reason(request: str, steps: List[dict]) -> Optional[str]:
    """A 5d (or 6p) element needs a relativistic Hamiltonian.

    Element matching is on WHOLE words in both spellings (s25 finding
    A1-A4): the word boundary must wrap the whole alternation —
    ``\\b(?:au|gold)\\b`` — or a prefix hits inside unrelated words
    ('redox' contains 're', 'IR' contains 'ir', 'water' contains 'w').
    The SYMBOL is matched case-sensitively (a chemist writes formulas
    capitalised: 'IR' is infrared, 'Ir' is iridium) and the ELEMENT NAME
    case-insensitively."""
    if (not any(re.search(rf"\b(?:{sym}|{name})\b", request)
                for sym, name in _5D_SYMBOLS)
            and not any(re.search(rf"\b{name}\b", request.lower())
                        for _, name in _5D_SYMBOLS)):
        return None
    for s in steps:
        kw = s.get("kwargs", {}) or {}
        if _norm(kw.get("relativity")):
            return None
    return ("a 5d metal needs a relativistic Hamiltonian: set relativity "
            "to zora/dkh/x2c with a matching recontracted basis")


def _word_in(text: str, word: str) -> bool:
    """Whole-word match (s25 finding B5): a substring test fires inside
    unrelated words ('frequencies of the fixed chair conformation' would
    otherwise demand a conformer search for a rigid molecule)."""
    return re.search(rf"\b{re.escape(word)}\b", text) is not None


def _conformer_reason(request: str, steps: List[dict]) -> Optional[str]:
    """A flexible molecule asked for its conformer landscape needs a cheap
    conformer search BEFORE the expensive level. Whole-word matching: a
    rigid-molecule request that merely mentions 'conformation' (as in
    'the fixed chair conformation') must not trigger it."""
    text = request.lower()
    asks_landscape = (
        _word_in(text, "conformer") or _word_in(text, "conformers")
        or _word_in(text, "flexible") or _word_in(text, "rotamer"))
    if not asks_landscape:
        return None
    if any(_norm(s.get("step")) in ("crest_conformers", "xtb_goat")
           for s in steps):
        return None
    return ("a flexible molecule goes to the conformer search first: add a "
            "crest_conformers or xtb_goat step before the DFT level")


def _coverage_reason(request: str, steps: List[dict]) -> Optional[str]:
    """Every capability the request names must be covered by a step.
    Keywords match as whole words (s25 finding B5)."""
    text = request.lower()
    names = {_norm(s.get("step")) for s in steps}
    missing = []
    for keywords, satisfying in _COVERAGE:
        if (any(_word_in(text, k.strip()) for k in keywords)
                and not (names & set(satisfying))):
            missing.append(keywords[0])
    if missing:
        return (f"the request asks for {missing} but no step covers it "
                f"(steps present: {sorted(names)})")
    return None


def _lint_verdict(findings: List[dict], answer: str) -> (List[dict], List[str]):
    """Science-lint warnings are tolerated ONLY when the answer justifies
    them in the same breath; errors never pass (science findings are
    advisory by design, _science.py:14-15 — the tolerance here is the
    grader's part of that contract)."""
    reasons: List[str] = []
    out: List[dict] = []
    for f in findings:
        entry = dict(f)
        if f.get("level") == "error":
            entry["justified"] = False
            reasons.append(f"lint error in step '{f.get('step')}': "
                           f"{f.get('message')}")
            out.append(entry)
            continue
        needle = (f.get("rule") or "").split("-")[0].lower()
        # The justification must mention the substance, not the rule id:
        # a reason counts where the deviation stands (wave-20 lesson A4).
        words = [w for w in re.split(r"\W+", f.get("message", ""))
                 if len(w) > 3]
        hit = any(w.lower() in answer.lower() for w in words)
        entry["justified"] = bool(hit)
        if not hit:
            reasons.append(
                f"lint warning '{f.get('rule')}' on step '{f.get('step')}' "
                f"is not justified in the answer: {f.get('message')}")
        out.append(entry)
    return out, reasons


def _cost_estimate(resolved: dict) -> float:
    total = 0.0
    for s in _steps_of(resolved):
        total += _CORE_HOUR_ESTIMATE.get(_norm(s.get("step")), 0.1)
    return total


# ---------------------------------------------------------------------------
# grade — the one entry point
# ---------------------------------------------------------------------------

def grade(spec: Dict[str, Any], request: str, *,
          answer: str = "", budget_core_hours: float = 4.0) -> dict:
    """Grade a built workflow spec against a chemist's request.

    Returns ``{"status": "pass"|"fail", "reasons": [..], "findings": {..}}``.
    ``findings`` carries the per-check detail: ``resolve``, ``validate``
    (ok + summary), ``lint`` (per finding, with ``justified``), ``coverage``,
    ``dry_run``, ``cost``. Raises ``ValueError`` for a spec with no name —
    the same error :func:`delfin.tools.platform.validate_spec` raises, on
    purpose: broken input is refused, not half-graded.
    """
    from delfin.tools import platform as _platform

    if not isinstance(spec, dict) or not spec.get("name"):
        raise ValueError(
            "spec has no 'name' — validate_spec refuses it and so does the "
            "grader (broken input is never half-graded)")

    reasons: List[str] = []
    steps = _steps_of(spec)

    # 1. resolve
    try:
        resolved = _platform.resolve_spec(spec)
        resolve_note = "resolved"
    except Exception as exc:  # noqa: BLE001 — a resolver failure IS the verdict
        resolved = spec
        resolve_note = f"resolve failed: {exc}"
        reasons.append(resolve_note)

    # 2. structural validation through the real call site.
    # geometry=True: a run hands the pipeline an input geometry (the
    # fixtures do), so a step consuming `geometry` is satisfied — the
    # validator's documented way to say so (pipeline.py:1637-1639).
    findings: Dict[str, Any] = {"resolve": resolve_note}
    try:
        rep = _platform.validate_spec(spec, geometry=True)
        findings["validate"] = {"ok": rep.ok, "summary": rep.summary()}
        if not rep.ok:
            reasons.append(f"validate_spec failed: {rep.summary()}")
    except Exception as exc:  # noqa: BLE001
        findings["validate"] = {"ok": False, "summary": str(exc)}
        reasons.append(f"validate_spec raised: {exc}")

    # 3. science lint, justified in the answer
    lint_findings = _platform.scientific_lint(spec)
    lint_entries, lint_reasons = _lint_verdict(lint_findings, answer or "")
    findings["lint"] = lint_entries
    reasons.extend(lint_reasons)

    # 4. the cross-step physics rules
    physics = (
        _mixed_level_reason(steps),
        _unverified_ts_reason(steps),
        _anion_reason(spec, steps),
        _relativistic_reason(request, steps),
        _conformer_reason(request, steps),
    )
    findings["physics"] = [r for r in physics if r]
    reasons.extend(r for r in physics if r)

    # 5. request coverage
    cov = _coverage_reason(request, steps)
    findings["coverage"] = cov or "covered"
    if cov:
        reasons.append(cov)

    # 6. dry run of the resolved pipeline (builds, never executes)
    try:
        from delfin.tools._serialize import from_dict
        from delfin.tools.pipeline import Pipeline
        built = from_dict(resolved)
        if not isinstance(built, Pipeline):
            built = built.build()
        built.dry_run()
        findings["dry_run"] = "passed"
    except Exception as exc:  # noqa: BLE001
        findings["dry_run"] = f"failed: {exc}"
        reasons.append(f"dry run failed: {exc}")

    # 7. cost estimate (documented as an estimate, never a measurement)
    est = _cost_estimate(resolved)
    over = est > budget_core_hours
    findings["cost"] = {
        "core_hours_estimate": est, "budget_core_hours": budget_core_hours,
        "within_budget": not over,
        "note": "order-of-magnitude estimate from step weights; real cost "
                "is EF1's cluster_efficiency.account after the run"}
    if over:
        reasons.append(
            f"cost estimate {est:.2f} core-hours exceeds the budget "
            f"{budget_core_hours:.2f}")

    return {"status": "fail" if reasons else "pass",
            "reasons": reasons, "findings": findings}


