"""Rank DELFIN's existing workflows against a calculation request.

Rule 2 of the on-demand contract: for a calculation request the order is
(a) an existing DELFIN workflow that fits, (b) extend it if it lacks a
step, (c) only then compose one from the ORCA building blocks.  This
module is the deterministic, model-free first half of (a): keyword rules
over the registered workflow descriptors, every candidate carrying the
reason it matched and the CONTROL keys to set (from
``delfin/tools/workflows.py``, which derives them from the
control_validator).  No model call, no ORCA run.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

__all__ = ["workflow_candidates"]


# (workflow name, [triggers], reason).  Triggers are lowercase substrings
# matched against the task text; the order of the list is the tie-break
# order (first list wins a tie).  Only workflows that exist in
# ``delfin.tools.workflows`` are ranked — the classic CONTROL.txt
# workflows are listed, never changed.
_RULES: list[tuple[str, tuple[str, ...], str]] = [
    ("redox",
     ("redox", "oxidation potential", "reduction potential",
      "oxidized", "reduced species", "e_oxidation", "electron transfer"),
     "redox workflow: neutral opt/freq plus oxidised/reduced species; "
     "returns the Gibbs energies the potential is derived from."),
    ("OCCUPIER",
     ("spin state", "spin states", "multiplicit", "ground state",
      "which spin", "high spin", "low spin"),
     "OCCUPIER: optimises every plausible spin state and selects the "
     "preferred one."),
    ("ESD",
     ("excited state", "emission", "phosphorescence", "fluorescence",
      "isc ", "intersystem", "tadf", "tddft", "absorption spectrum"),
     "ESD: TDDFT per state plus emission/ISC/IC rates."),
    ("co2",
     ("co2", "carbon dioxide", "adduct", "co2 binding"),
     "CO2 coordinator: chains ORCA runs across the CO2 adduct species."),
    ("IMAG",
     ("imaginary freq", "imaginary mode", "saddle point",
      "drive the saddle", "saddle to a minimum"),
     "IMAG: reoptimises along the imaginary mode until the structure is "
     "a true minimum."),
    ("classic",
     ("rotational barrier", "barrier", "transition state", " ts ",
      "ts for", "ts of", "optimi", "geometry optimisation",
      "geometry optimization", "frequency", "frequencies"),
     "classic CONTROL workflow: optimise, frequency-check and screen the "
     "structure."),
    ("conformer_dG",
     ("conformer", "rotamer", "ensemble"),
     "conformer_dG: CREST ensemble with CENSO free-energy ranking."),
]

# method_family keywords that pull the TS building blocks as the
# compose-from-blocks hint (rule (c)) when no workflow covers the job.
_TS_TRIGGERS = ("ts ", "ts for", "ts of", "transition state",
                "neb", "irc", "scan")


def _descriptors(available: Optional[List[str]]):
    from delfin.tools.workflows import list_workflows

    flows = list_workflows()
    if available is None:
        return flows
    keep = {str(a).strip().lower() for a in available}
    return [w for w in flows if w.name.lower() in keep]


# Generic classic triggers that a more specific workflow subsumes: a text
# mentioning "imaginary frequency" also matches "frequency" and "optimi",
# but the classic opt/freq screen is not the answer there — IMAG is.
_CLASSIC_SUBSUMED = ("frequency", "frequencies", "optimi",
                     "geometry optimisation", "geometry optimization",
                     "barrier")
_IMAG_MARKERS = ("imaginary", "saddle")


def workflow_candidates(task_text: str,
                        available: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """Ranked list of workflows that fit ``task_text``.

    Each candidate is ``{workflow, kind, reason, control_keys, steps}``;
    ``control_keys`` entries carry ``name``/``type``/``default`` and are
    the validated CONTROL keys of the descriptor (never invented here).
    When nothing matches, the returned list is empty — which is the
    explicit signal to extend an existing workflow (b) or compose from
    the ORCA building blocks (c).
    """
    if not isinstance(task_text, str):
        return []
    text = f" {task_text.lower()} "
    flows = _descriptors(available)
    # A more specific workflow subsumes classic's generic opt/freq wording.
    classic_shrunk = any(m in text for m in _IMAG_MARKERS)
    scored: List[tuple[int, int, Dict[str, Any]]] = []
    for order, (name, triggers, reason) in enumerate(_RULES):
        if name == "classic" and classic_shrunk:
            triggers = tuple(t for t in triggers if t not in _CLASSIC_SUBSUMED)
        hits = [t for t in triggers if t in text]
        if not hits:
            continue
        flow = next((w for w in flows
                     if w.name.lower() == name.lower()), None)
        if flow is None:
            continue
        scored.append((len(hits), -order, _candidate(flow, reason, hits)))
    scored.sort(key=lambda s: (-s[0], s[1]))
    return [c for _, _, c in scored]


def _candidate(flow, reason: str, hits) -> Dict[str, Any]:
    return {
        "workflow": flow.name,
        "kind": getattr(flow, "kind", ""),
        "reason": f"{reason} (matched: {', '.join(sorted(hits))})",
        "control_keys": [
            {"name": c.name, "type": c.type,
             **({"default": c.default} if c.default is not None else {})}
            for c in flow.control_keys
        ],
        "steps": [str(s.get("label", s.get("step", "")))
                  for s in (flow.steps or []) if isinstance(s, dict)],
    }
