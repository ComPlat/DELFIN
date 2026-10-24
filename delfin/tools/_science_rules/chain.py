"""P2 chain-consistency realism rules for the scientific lint.

Called from delfin.tools._science.scientific_lint — the single entry
point. These rules are SPEC-LEVEL: they need the sequence of steps
(level consistency, freq after opt, one solvent through the chain), so
they run over the walked step list, not on one step's kwargs.

Sources:
- Same level of theory per energy difference: DELFIN wave-20 rule
  "compare methods, not programs" (a total energy compares only within
  one functional/basis/dispersion/solvent); ZPE/thermal corrections
  valid on a strongly converged geometry at the optimisation level —
  J. Chem. Theory Comput. 2021, 17, DOI 10.1021/acs.jctc.1c00474.
- Thermochemistry at explicit T/P: ORCA 6.1.1 manual §4.7 (ideal-gas
  statistical mechanics; defaults 298.15 K, 1 atm).
- TS needs a frequency check and an IRC to both minima: DELFIN's own TS
  adapters delfin/tools/adapters/orca_ts_adapters.py — orca_ts_opt,
  orca_ts_freq_check, orca_ts_irc; the IRC fails closed with
  irc_endpoints=None when no endpoints are reached (:294-298).
"""

from __future__ import annotations

from typing import Any, Dict, List

from delfin.tools._science import ScienceFinding, _is_placeholder

# Steps that establish the level of theory of an energy difference.
_LEVEL_STEP_NAMES = frozenset({"orca_opt", "orca_freq", "orca_sp",
                               "orca_opt_freq", "orca_tddft", "orca_ts_opt",
                               "orca_ts_freq_check", "orca_ts_irc"})
# Steps that produce a geometry (an optimisation of any kind).
_OPT_STEP_NAMES = frozenset({"orca_opt", "orca_opt_freq"})
_FREQ_STEP_NAMES = frozenset({"orca_freq", "orca_opt_freq"})
_TS_OPT = "orca_ts_opt"
_TS_CHECK = "orca_ts_freq_check"
_TS_IRC = "orca_ts_irc"


def _norm(v: Any) -> str:
    return str(v).strip().lower() if v is not None and not _is_placeholder(v) else ""


def _level(kwargs: Dict[str, Any]) -> tuple:
    """The (method, basis, dispersion, solvent, solvent_model) tuple that
    identifies one level of theory; empty entries mean 'inherited at
    runtime', i.e. not comparable."""
    return tuple(_norm(kwargs.get(k)) for k in
                 ("method", "basis", "dispersion", "solvent", "solvent_model"))


def lint_chain(steps: List[Dict[str, Any]]) -> List[ScienceFinding]:
    """Rules 1, 2, 8, 9 over one walked step sequence (trunk or branch)."""
    out: List[ScienceFinding] = []
    seen_level: tuple = ()
    seen_level_label = ""
    seen_solvent = ""
    seen_model = ""
    opt_pending: List[str] = []       # labels of orca_opt without a freq yet
    ts_pending: List[str] = []        # labels of orca_ts_opt
    await_irc: List[str] = []         # TS labels frequency-checked, IRC owed
    ts_check_seen = False

    for st in steps:
        if not isinstance(st, dict):
            continue
        name = _norm(st.get("step") or st.get("name") or "")
        label = str(st.get("label") or name or "?")
        kw = st.get("kwargs", {}) or {}

        # ── Rule 1: one level of theory per energy difference ──
        if name in _LEVEL_STEP_NAMES:
            lv = _level(kw)
            if any(lv) and seen_level and lv != seen_level:
                # compare only the entries BOTH steps state concretely
                differing = [k for k, a, b in zip(
                    ("method", "basis", "dispersion", "solvent", "solvent_model"),
                    seen_level, lv) if a and b and a != b]
                if differing:
                    out.append(ScienceFinding(
                        label, "warning", "level-consistency",
                        f"step changes the level of theory within one chain "
                        f"({' / '.join(differing)} {seen_level_label} vs here); "
                        "one energy difference must use one level of theory "
                        "— ZPE/thermal corrections belong to the level the "
                        "geometry was optimised at"))
            if any(lv):
                seen_level, seen_level_label = lv, label

        # ── Rule 8: one solvent model + solvent through the chain ──
        solvent = _norm(kw.get("solvent"))
        model = _norm(kw.get("solvent_model"))
        if solvent and seen_solvent and solvent != seen_solvent:
            out.append(ScienceFinding(
                label, "warning", "solvent-inconsistent",
                f"solvent changes within one chain ({seen_solvent} -> "
                f"{solvent}); keep solvent and solvent model consistent "
                "through the whole workflow"))
        if model and seen_model and model != seen_model:
            out.append(ScienceFinding(
                label, "warning", "solvent-inconsistent",
                f"solvent model changes within one chain ({seen_model} -> "
                f"{model}); CPCM and SMD free energies are not comparable"))
        if solvent:
            seen_solvent = solvent
        if model:
            seen_model = model

        # ── Rule 2: frequencies after every optimisation ──
        if name in _FREQ_STEP_NAMES:
            opt_pending.clear()
        elif name in _OPT_STEP_NAMES:
            # orca_opt_freq IS its own frequency calculation — nothing pending
            opt_pending.append(label)
        if name == "orca_freq":
            # ── Rule 9: explicit conditions for thermochemistry ──
            if not _norm(kw.get("freq_temp")):
                out.append(ScienceFinding(
                    label, "warning", "thermochemistry-conditions-missing",
                    "thermochemistry (ZPE, thermal corrections) is "
                    "computed at ORCA's default 298.15 K / 1 atm unless "
                    "stated — set freq_temp explicitly so the conditions "
                    "of the reported Gibbs energy are unambiguous"))
        if name in _TS_OPT:
            ts_pending.append(label)
            ts_check_seen = False
        if name in _TS_CHECK:
            # the frequency check CLOSES its preceding TS optimisation
            # (the finding only fires if the chain ends without a check);
            # the chain still owes it an IRC to both minima
            await_irc.extend(ts_pending)
            ts_pending.clear()
            ts_check_seen = True
        if name in _TS_IRC:
            for plabel in ts_pending:
                out.append(ScienceFinding(
                    plabel, "warning", "ts-freq-missing",
                    f"TS optimisation '{plabel}' has no frequency check "
                    "before its IRC — a TS must show exactly 1 imaginary "
                    "mode before following the reaction path"))
            ts_pending.clear()
            await_irc.clear()
            if not ts_check_seen:
                out.append(ScienceFinding(
                    label, "warning", "ts-irc-missing",
                    "IRC without a prior TS frequency check — confirm the "
                    "stationary point is a TS (exactly 1 imaginary mode) "
                    "before following the reaction path to both minima"))
            ts_check_seen = False

    # End of chain: pending optimisations without a frequency calculation.
    for plabel in opt_pending:
        out.append(ScienceFinding(
            plabel, "warning", "freq-after-opt-missing",
            f"geometry optimisation '{plabel}' has no frequency calculation "
            "after it — a minimum has 0 imaginary modes; run orca_freq "
            "on the optimised geometry"))
    for plabel in ts_pending:
        out.append(ScienceFinding(
            plabel, "warning", "ts-freq-missing",
            f"TS optimisation '{plabel}' has no frequency check and no IRC — "
            "a TS must show exactly 1 imaginary mode plus an IRC to both "
            "minima"))
        out.append(ScienceFinding(
            plabel, "warning", "ts-irc-missing",
            f"TS optimisation '{plabel}' has no IRC — a TS must be followed "
            "to both minima (orca_ts_irc) to confirm the reaction path"))
    for plabel in await_irc:
        out.append(ScienceFinding(
            plabel, "warning", "ts-irc-missing",
            f"TS optimisation '{plabel}' was frequency-checked but has no "
            "IRC — a TS must be followed to both minima (orca_ts_irc) to "
            "confirm the reaction path"))
    return out


__all__ = ["lint_chain"]
