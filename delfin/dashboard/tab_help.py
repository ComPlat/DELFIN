"""Help tab: a feature overview and a generated CONTROL-key reference.

The tab has two halves:

* **Features** — a curated overview of what DELFIN can do (the workflows and
  modules), each with a short "what it does" and "which CONTROL keys switch it
  on".  Curated text lives in ``_FEATURES`` below; it is documentation, not
  behaviour, so a wrong line here can be fixed by editing the list.
* **CONTROL keys** — generated from ``CONTROL_FIELD_SPECS`` (the same table
  ``validate_control_config`` validates against), so the reference is always in
  sync with the code: every key, its default, whether it is required, and its
  ``help`` text (the same text ``KEY=?`` prints).  Keys without help text are
  listed with an explicit "no explanation yet" instead of silently vanishing.

The tab is read-only: it mutates nothing and reads no calculation state.
"""

from __future__ import annotations

from typing import Any

import ipywidgets as W

# ---------------------------------------------------------------------------
# Curated feature overview
# ---------------------------------------------------------------------------

#: (feature, what it does, CONTROL keys that switch it on).  Order = display
#: order.  Keep each entry short; the CONTROL-key reference below has details.
_FEATURES: tuple[tuple[str, str, str], ...] = (
    (
        "Standard workflow (SPE / OPT / FREQ)",
        "Reads a geometry or SMILES and runs single-point, optimization and "
        "frequency jobs in order, then reports energies and thermochemistry "
        "into DELFIN.txt and DELFIN_Data.json.",
        "method, FUNCTIONAL, BASIS, charge, multiplicity_global_opt, PAL, "
        "FREQUENCY_CALCULATION (or frequency_calculation), "
        "CONSTRAINTS, SOLVENT",
    ),
    (
        "GUPPY (continuum solvation variants)",
        "Legacy module that runs a grid of solvation / dispersion variants of "
        "the same job and compares them.  Kept so archived CONTROL files still "
        "validate; prefer the standard workflow with SOLVENT today.",
        "GUPPY, GUPPY_run, GUPPY_solvent, GUPPY_functional",
    ),
    (
        "ESD (excited states / spectra)",
        "Excited-state calculations and spectra (UV/Vis, fluorescence, "
        "phosphorescence) with rate and energy reporting.",
        "ESD, ESD_run, ESD_functional, ESD_basis",
    ),
    (
        "IMAG (imaginary-frequency follow-up)",
        "Detects imaginary frequencies after a Hessian and displaces along "
        "them automatically to find the intended minimum or transition state.",
        "IMAG, IMAG_run",
    ),
    (
        "OCCUPIER (adduct / coordination occupancy)",
        "Screening of charge-adduct and coordination states: generates "
        "candidate structures, ranks them, and reports which species is "
        "favoured.  Comparison on FSPE or Gibbs (adds a frequency run).",
        "OCCUPIER_method (auto|manually), OCCUPIER_compare "
        "(FSPE|G), frequency_calculation_OCCUPIER (legacy), "
        "OCCUPIER_generate, OCCUPIER_adducts",
    ),
    (
        "CO2 coordination",
        "Coordination of CO2 to a substrate: candidate structures and ranking, "
        "often combined with the OCCUPIER machinery.",
        "CO2_coordination, CO2_substrate, method=OCCUPIER",
    ),
    (
        "MANTA (structure generation funnel)",
        "SMILES-to-3D funnel: screening, optimization and refinement stages "
        "that turn input SMILES into a ranked set of 3D structures.",
        "MANTA_SCREEN, MANTA_OPT, MANTA_REFINE, MANTA_REFINE_TOPK, "
        "MANTA_SCREEN_KEEP, MANTA_QUALITY, MANTA_MULTIPLICITIES, "
        "smiles_converter (QUICK|NORMAL|MANTA|ARCHITECTOR|MOLSIMPLIFY|MACE)",
    ),
    (
        "Surface scans and relaxed scans",
        "1D/nD relaxed surface scans along chosen coordinates, with energy "
        "profiles in the output and follow-up TS searches.",
        "RSS, SCAN",
    ),
    (
        "MLP (machine-learned potentials)",
        "Preoptimization and sampling with machine-learned potentials (e.g. "
        "MACE) before or instead of DFT levels.",
        "MLP, mlp_method",
    ),
    (
        "SMILES converter",
        "Builds 3D coordinates from SMILES with the configured converter "
        "backend; placeholders in the value fall back to NORMAL.",
        "smiles_converter, thdy_smiles_converter, sc_smiles_converter",
    ),
    (
        "QMMM region splitting",
        "Splits the structure into QM and MM regions with a ``$`` separator "
        "line in the input file; the split is preserved for follow-up jobs.",
        "(input file convention, no CONTROL key)",
    ),
    (
        "Preoptimization",
        "Cheap preoptimization (e.g. xTB) before the DFT optimization.",
        "preopt, thdy_preopt",
    ),
    (
        "Runs / unattended operation",
        "Multi-run CONTROL files and unattended mode let one file drive many "
        "systems; see docs/unattended_mode.md.",
        "RUNS",
    ),
    (
        "Agent tab",
        "The DELFIN Agent tab: ask questions about calculations, errors, "
        "methods and CONTROL keys in chat; read-only unless you approve.",
        "(no CONTROL key)",
    ),
)


# ---------------------------------------------------------------------------
# Widget construction
# ---------------------------------------------------------------------------


def _feature_block() -> W.DOMWidget:
    parts: list[str] = [
        "<h3 style='margin-top:0'>Features</h3>",
        "<p style='color:#666'>What DELFIN can do, and which CONTROL.txt keys "
        "switch it on.  Detailed key semantics: write <code>KEY=?</code> in "
        "CONTROL.txt, or see the key reference below.</p>",
    ]
    for name, what, keys in _FEATURES:
        parts.append(
            "<div style='border-left:3px solid #4a86c8; padding:6px 12px; "
            "margin:10px 0'>"
            f"<b>{name}</b><br>{what}<br>"
            f"<span style='color:#666'><b>Keys:</b> {keys}</span></div>"
        )
    return W.HTML("".join(parts))


def _key_rows() -> list[tuple[str, str, str, str]]:
    """One row per CONTROL key: (name, default, required, help)."""
    from delfin.common.control_validator import CONTROL_FIELD_SPECS

    rows: list[tuple[str, str, str, str]] = []
    for spec in CONTROL_FIELD_SPECS:
        help_text = (spec.help or "").strip()
        if not help_text:
            help_text = "<i style='color:#b00'>no explanation yet — ask the agent or write KEY=?</i>"
        default = "" if spec.default is None else str(spec.default)
        rows.append(
            (spec.name, default, "yes" if spec.required else "", help_text)
        )
    return sorted(rows, key=lambda r: r[0].lower())


def _key_reference_block() -> tuple[W.DOMWidget, W.DOMWidget, Any]:
    """Searchable table of all CONTROL keys.  Returns (block, html, state)."""
    state: dict[str, Any] = {"query": "", "rows": _key_rows()}

    out = W.HTML()

    def _render() -> None:
        query = state["query"].strip().lower()
        esc = lambda s: s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        head = (
            "<table style='border-collapse:collapse; width:100%'>"
            "<tr style='background:#eef'><th style='text-align:left; padding:4px 8px'>Key</th>"
            "<th style='text-align:left; padding:4px 8px'>Default</th>"
            "<th style='text-align:left; padding:4px 8px'>Required</th>"
            "<th style='text-align:left; padding:4px 8px'>What it does</th></tr>"
        )
        shown = 0
        rows_html: list[str] = []
        for name, default, required, help_text in state["rows"]:
            if query and query not in name.lower() and query not in help_text.lower():
                continue
            shown += 1
            style = "border-top:1px solid #ddd; padding:4px 8px; vertical-align:top"
            rows_html.append(
                "<tr>"
                f"<td style='{style}'><code>{esc(name)}</code></td>"
                f"<td style='{style}'>{esc(default)}</td>"
                f"<td style='{style}'>{esc(required)}</td>"
                f"<td style='{style}'>{help_text}</td>"
                "</tr>"
            )
        total = len(state["rows"])
        if not shown:
            rows_html.append(
                f"<tr><td colspan='4' style='padding:8px'>no key matches "
                f"'{esc(query)}' (of {total} keys)</td></tr>"
            )
        out.value = head + "".join(rows_html) + "</table>"
        state["shown"] = shown
        state["total"] = total

    _render()

    search = W.Text(
        placeholder="filter keys or help text…",
        layout=W.Layout(width="40%"),
    )

    counter = W.HTML()

    def _on_search(change: Any) -> None:
        state["query"] = change["new"] or ""
        _render()
        counter.value = (
            f"<span style='color:#666'>{state['shown']} of "
            f"{state['total']} keys shown</span>"
        )

    search.observe(_on_search, names="value")
    _on_search({"new": ""})

    box = W.VBox(
        [
            W.HTML(
                "<h3 style='margin-top:0'>CONTROL.txt key reference</h3>"
                "<p style='color:#666'>Generated from "
                "<code>CONTROL_FIELD_SPECS</code> — the same table the "
                "validator checks your CONTROL.txt against, so it cannot drift "
                "from the code.  The help column is the same text that "
                "<code>KEY=?</code> prints.</p>"
            ),
            W.HBox([search, counter]),
            out,
        ]
    )
    return box, out, state


def create_tab(ctx: Any):
    """Create the Help tab. Returns ``(widget, refs)``."""
    feature = _feature_block()
    keys_box, _keys_html, keys_state = _key_reference_block()
    header = W.HTML(
        "<h2 style='margin-top:0'>DELFIN Help</h2>"
        "<p>An overview of what DELFIN can do and which CONTROL.txt keys "
        "configure it.  Read-only: nothing here starts a calculation.</p>"
        "<hr>"
    )
    widget = W.VBox([header, feature, W.HTML("<hr>"), keys_box])
    return widget, {"help_key_state": keys_state}


# Register additively with the dashboard tab registry.
try:
    from delfin.dashboard.tab_registry import register_tab
    register_tab("help", "Help", create_tab, order=8600)
except Exception:  # pragma: no cover
    pass


__all__ = ["_FEATURES", "create_tab"]
