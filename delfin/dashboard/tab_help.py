"""Help tab: a feature overview and a generated CONTROL-key reference.

The tab is a two-level discovery surface:

* **Features** (subtab 1) — generalised DELFIN functions as an accordion:
  one collapsed section per major capability (structure generation,
  optimisation, redox workflows, excited states, thermodynamics, ...),
  opening to the sub-functions, what they do, and the CONTROL keys that
  switch them on.  A user who does not know the features yet sees the
  big picture first, details on demand.
* **CONTROL keys** (subtab 2) — the full, searchable reference generated
  from ``CONTROL_FIELD_SPECS`` (the same table ``validate_control_config``
  validates against), so it cannot drift from the code: every key, its
  default, whether it is required, and its ``help`` text (the same text
  ``KEY=?`` prints).  Grouped by key family; keys without help text would
  be listed with an explicit marker instead of silently vanishing.

The tab is read-only: it mutates nothing and reads no calculation state.
"""

from __future__ import annotations

from typing import Any

import ipywidgets as W


# ---------------------------------------------------------------------------
# Curated feature overview: generalised functions -> sub-functions
# ---------------------------------------------------------------------------
# (major function, one-line summary, [(sub-function, what it does, keys), ...])
_FEATURE_GROUPS: tuple[tuple[str, str, tuple[tuple[str, str, str], ...]], ...] = (
    (
        "Structure input & generation",
        "From a formula, a geometry file or a SMILES to 3D coordinates.",
        (
            ("Input file / SMILES",
             "input.txt with a geometry, or SMILES= in CONTROL.txt; DELFIN never "
             "overwrites your input.",
             "SMILES, input_file"),
            ("SMILES converter",
             "Builds 3D coordinates from a SMILES: QUICK (one embedding, fast), "
             "NORMAL (multi-seed embedding + force-field refinement), MANTA "
             "(metal-complex ensemble), or the external builders ARCHITECTOR, "
             "MOLSIMPLIFY, MACE.",
             "smiles_converter, thdy_smiles_converter"),
            ("MANTA structure funnel",
             "Isomer-conformer ensemble of a metal complex: build, deduplicate, "
             "screen (single points), optimise (xTB), refine the best (GOAT/CREST); "
             "the winner becomes start.txt.",
             "MANTA_SCREEN, MANTA_OPT, MANTA_REFINE, MANTA_QUALITY, MANTA_SEEDS, "
             "MANTA_MAX_ISOMERS, MANTA_MULTIPLICITIES (legacy: GUPPY_*)"),
            ("Explicit solvent shells",
             "xTB ALPB solvation adds explicit solvent molecules to the geometry.",
             "XTB_SOLVATOR, number_explicit_solv_molecules"),
        ),
    ),
    (
        "Preparation & optimisation",
        "Cheap pre-treatment before the DFT level, then the real optimisation.",
        (
            ("Preoptimisation (xTB)",
             "A quick xTB geometry optimisation before the DFT optimisation.",
             "XTB_preOPT (XTB_OPT), thdy_preopt, hyperpol_xTB_preopt"),
            ("Global conformer search",
             "GOAT or CREST finds the global minimum / conformer ensemble.",
             "global_optimizer, MANTA_REFINE, MANTA_REFINE_TOPK"),
            ("Geometry optimisation & frequencies",
             "DFT optimisation, optional frequency calculation per optimisation "
             "(thermochemistry, Gibbs, imaginary-mode check).",
             "geom_opt, freq_type, frequency_calculation, geom_opt_OCCUPIER, "
             "multiplicity_global_opt"),
            ("Imaginary-mode handling (IMAG)",
             "Displaces along imaginary frequencies to reach the intended minimum; "
             "noise floor and round limits are configurable.",
             "IMAG, IMAG_scope, IMAG_option, IMAG_max_rounds, "
             "allow_imaginary_freq, IMAG_sp_energy_window"),
        ),
    ),
    (
        "Method & environment",
        "The quantum-chemical level everything runs at.",
        (
            ("Functional, basis set, relativistics",
             "DFT functional, main/metal basis sets (relativistic variants "
             "separately), RI approximations, SCF initial guess.",
             "functional, main_basisset, metal_basisset, main_basisset_rel, "
             "metal_basisset_rel, aux_jk, aux_jk_rel, ri_jkx, initial_guess, "
             "relativity"),
            ("Solvation",
             "Implicit solvation model and solvent for every DFT job.",
             "implicit_solvation_model, solvent"),
            ("Parallelism & resources",
             "Core budget, job parallelism and how ORCA is started on the node.",
             "PAL, pal_jobs, orca_parallel_strategy"),
            ("xTB settings",
             "Which xTB Hamiltonian the pre-optimisation and xTB workflows use.",
             "xTB_method"),
        ),
    ),
    (
        "Redox workflows",
        "Redox potentials across charges: the classic DELFIN use case.",
        (
            ("classic",
             "Plain SPE/OPT/FREQ chain over the requested redox steps.",
             "method=classic, oxidation_steps, reduction_steps"),
            ("manually",
             "You define the chain via the step lists; DELFIN schedules it.",
             "method=manually"),
            ("OCCUPIER",
             "Automatic occupancy search: candidate configurations per redox "
             "state are generated, ranked (by electronic energy or Gibbs) and "
             "the winner seeds the next step.",
             "method=OCCUPIER, OCCUPIER_method, OCCUPIER_compare, "
             "OWN_progressive_from, OWN_TREE_PURE_WINDOW, "
             "approximate_spin_projection_APMethod"),
            ("Redox potentials",
             "Potentials from the computed states; several derivation methods.",
             "calc_potential_method, properties_of_interest, "
             "reorganisation_energy"),
            ("CO2 coordination",
             "How CO2 coordinates to the substrate: candidates built and ranked "
             "across the redox states, on top of the OCCUPIER machinery.",
             "co2_coordination, co2_species_delta"),
        ),
    ),
    (
        "Excited states (ESD)",
        "Photophysics and photochemistry from TD-DFT or delta-SCF.",
        (
            ("States & rates",
             "Electronic states, then internal conversion, intersystem crossing "
             "and emission rates between them.",
             "ESD_modul, ESD_modus, ESD_T1_opt, states, ICs, ISCs, "
             "emission_rates"),
            ("TD-DFT settings",
             "NRoots, Davidson space, TDA, spin-orbit coupling, extra keywords "
             "for every TD-DFT job DELFIN writes.",
             "TDDFT_nroots, TDDFT_maxdim, TDDFT_maxiter, TDDFT_TDA, TDDFT_SOC, "
             "TDDFT_additions"),
            ("xTB hyperpolarizability",
             "First hyperpolarizability beta with the sTD2/sTDA xTB response "
             "engine, static or at chosen wavelengths.",
             "hyperpol_xTB, hyperpol_xTB_engine, hyperpol_xTB_wavelengths, "
             "hyperpol_xTB_bfw, hyperpol_xTB_energy_window"),
            ("xTB TADF screening",
             "Singlet/triplet screening for thermally activated delayed "
             "fluorescence.",
             "tadf_xTB, tadf_xTB_excited_method, tadf_xTB_energy_window, "
             "tadf_xTB_run_t1_opt"),
        ),
    ),
    (
        "Thermodynamics",
        "Stability constants from a Born-Haber cycle.",
        (
            ("log K computation",
             "Gibbs energies of complex, ligands and displaced solvent give the "
             "stability constant; an experimental value can be compared.",
             "stability_constant, stability_constant_mode, stability_reaction, "
             "n_explicit_solvent, logK_exp (aliases: thermodynamics, "
             "thermodynamics_mode, thermodynamics_reaction)"),
        ),
    ),
    (
        "Automation & operation",
        "Driving many systems and reading the results.",
        (
            ("Multi-run CONTROL (RUNS)",
             "One CONTROL file drives many systems; unattended mode runs without "
             "interaction. See docs/unattended_mode.md.",
             "RUNS"),
            ("QMMM region splitting",
             "A ``$`` separator line splits QM and MM regions; preserved for "
             "follow-up jobs. No CONTROL key needed.",
             "(input-file convention)"),
            ("Outputs",
             "DELFIN.txt (potentials), DELFIN_Data.json (structured data), "
             "OCCUPIER.txt, delfin_run.log.",
             "—"),
            ("Agent tab",
             "Ask questions about calculations, errors, methods and CONTROL "
             "keys in chat; read-only unless you approve.",
             "—"),
        ),
    ),
)


# ---------------------------------------------------------------------------
# Widget construction
# ---------------------------------------------------------------------------


def _feature_block(on_open_key) -> tuple[W.DOMWidget, W.Accordion]:
    """The accordion of major functions with sub-function sections.

    ``on_open_key(name)`` is called when a key name is clicked: the tab
    switches to the CONTROL-keys subtab and filters to that key.

    Returns ``(widget, accordion)`` -- the accordion handle so the key
    table's "Used in" buttons can open a specific section.
    """
    from html import escape as esc

    header = W.HTML(
        "<p style='color:#666'>What DELFIN can do, grouped by major function. "
        "Open a section for the sub-functions and the CONTROL.txt keys that "
        "switch them on.  Key names are clickable: they open the "
        "<b>CONTROL keys</b> subtab filtered to that key.  In CONTROL.txt, "
        "write <code>KEY=?</code> for the same explanation.</p>"
    )
    no_match = W.HTML("", layout=W.Layout(display="none"))

    def _feature_filter(query: str) -> None:
        """Cross-tab search: hide sections that match neither their group
        title, a sub-function name nor any key; open the first match."""
        q = (query or "").lower()
        hits = []
        for child, (title, summary, subs) in zip(
            accordion.children, _FEATURE_GROUPS
        ):
            match = (
                not q or q in title.lower() or q in summary.lower()
                or any(
                    q in sub.lower() or q in what.lower() or q in keys.lower()
                    for sub, what, keys in subs
                )
            )
            child.layout.display = "" if match else "none"
            if match:
                hits.append(child)
        if q and hits:
            accordion.selected_index = (
                accordion.children.index(hits[0])
                if hasattr(accordion.children, "index") else 0
            )
        no_match.layout.display = "none" if hits else ""

    children = []
    # Help text per key, for hover tooltips on the key buttons.
    help_of = {row[0]: row[3] for row in _key_rows()}
    for title, summary, subs in _FEATURE_GROUPS:
        group_parts = [W.HTML(f"<p style='color:#666'>{summary}</p>")]
        for sub, what, keys in subs:
            key_names = [k.strip() for k in keys.split(",") if k.strip()]
            buttons = []
            for k in key_names:
                btn = W.Button(
                    description=k, button_style="",
                    style=W.ButtonStyle(
                        text_color="#4a86c8", font_weight="normal",
                        padding="0", background="none",
                    ),
                    tooltip=help_of.get(k, k),
                    layout=W.Layout(margin="0 4px 0 0", padding="0"),
                )
                btn.on_click(lambda _b, name=k: on_open_key(name))
                buttons.append(btn)
            group_parts.append(W.VBox([
                W.HTML(
                    "<div style='border-left:3px solid #4a86c8; "
                    "padding:6px 12px; margin:8px 0'>"
                    f"<b>{esc(sub)}</b><br>{esc(what)}</div>"
                ),
                W.HBox(
                    [W.HTML("<span style='color:#666'><b>Keys:</b></span>")]
                    + buttons
                ),
            ]))
        children.append(W.VBox(group_parts))
    # Tag each section with its group title so "Used in" buttons can open it.
    for (title, _summary, _subs), child in zip(_FEATURE_GROUPS, children):
        child._feature_group = title
    # Accordion titles are PLAIN TEXT (ipywidgets renders them verbatim,
    # not as HTML) -- the screenshot showed literal <b>/<span> tags.
    titles = [f"{title} — {summary}" for title, summary, _ in _FEATURE_GROUPS]
    accordion = W.Accordion(children=children, titles=tuple(titles))
    # Start collapsed except the first group, so the page opens compact.
    accordion.selected_index = 0

    return W.VBox([header, accordion]), accordion, _feature_filter


def _key_to_group() -> dict[str, str]:
    """Map every CONTROL key to the feature group that mentions it.

    Keys no feature section names return "" -- they get no Used-in
    button (the key reference stays complete regardless).
    """
    used: dict[str, str] = {}
    for title, _summary, subs in _FEATURE_GROUPS:
        for _sub, _what, keys in subs:
            for k in keys.split(","):
                k = k.strip()
                if k and k not in used:
                    used[k] = title
    return used


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


def _group_of(name: str) -> str:
    """Key family for grouping: the first MANTA_/GUPPY_/TDDFT_ prefix, else General."""
    for prefix, label in (
        ("GUPPY_", "GUPPY (legacy)"),
        ("MANTA_", "MANTA"),
        ("TDDFT_", "TD-DFT"),
        ("hyperpol_xTB", "xTB hyperpolarizability"),
        ("tadf_xTB", "xTB TADF"),
        ("OCCUPIER_", "OCCUPIER"),
        ("OWN_", "OCCUPIER"),
        ("IMAG", "IMAG"),
        ("stability_", "Thermodynamics"),
        ("thdy_", "Thermodynamics"),
        ("co2_", "CO2"),
    ):
        if name.startswith(prefix):
            return label
    return "General"


def _key_reference_block(
    on_open_feature=None,
) -> tuple[W.DOMWidget, W.DOMWidget, Any]:
    """Searchable key table, grouped by key family with per-group toggles.

    ``on_open_feature(group_title)`` is called when a row's "Used in"
    button is clicked: the tab switches to the Features subtab and opens
    that accordion section.
    """
    from html import escape as esc

    state: dict[str, Any] = {"rows": _key_rows(), "query": ""}
    used_in = _key_to_group()

    head = W.HTML(
        "<div style='display:flex; gap:12px; padding:4px 0; "
        "border-bottom:2px solid #333; font-weight:bold'>"
        "<div style='flex:0 0 220px'>Key</div>"
        "<div style='flex:0 0 90px'>Default</div>"
        "<div style='flex:0 0 70px'>Required</div>"
        "<div style='flex:1'>What it does</div>"
        "<div style='flex:0 0 170px'>Used in</div>"
        "</div>"
    )

    out = W.VBox()

    groups: dict[str, list[tuple[str, str, str, str]]] = {}
    for row in state["rows"]:
        groups.setdefault(_group_of(row[0]), []).append(row)
    group_names = sorted(groups)
    # General first, the rest alphabetical.
    group_names.remove("General") if "General" in group_names else None
    order = ["General"] + group_names

    group_toggles = W.HBox()
    state["active_groups"] = set(order)

    def _make_toggle(label: str, count: int) -> W.ToggleButton:
        btn = W.ToggleButton(
            value=True, description=f"{label} ({count})",
            button_style="", layout=W.Layout(margin="0 4px 4px 0"),
        )

        def _on_toggle(change: Any) -> None:
            if change["new"]:
                state["active_groups"].add(label)
            else:
                state["active_groups"].discard(label)
            _render()

        btn.observe(_on_toggle, names="value")
        return btn

    group_toggles.children = tuple(
        _make_toggle(g, len(groups[g])) for g in order
    )

    def _render() -> None:
        query = state["query"].lower()
        total = len(state["rows"])
        shown = 0
        children: list[W.DOMWidget] = []
        for group in order:
            if group not in state["active_groups"]:
                continue
            group_header = W.HTML(
                f"<div style='padding:10px 0 2px; color:#4a86c8'>"
                f"<b>{esc(group)}</b></div>"
            )
            children.append(group_header)
            for name, default, required, help_text in groups[group]:
                if query and query not in name.lower() and query not in help_text.lower():
                    continue
                shown += 1
                used_group = used_in.get(name, "")
                if used_group and on_open_feature:
                    btn = W.Button(
                        description=f"Used in: {used_group}",
                        button_style="",
                        style=W.ButtonStyle(
                            text_color="#4a86c8", font_size="12px",
                            padding="0", background="none",
                        ),
                        layout=W.Layout(margin="0", padding="0"),
                    )
                    btn.on_click(
                        lambda _b, g=used_group: on_open_feature(g)
                    )
                    used_slot = btn
                else:
                    used_slot = W.HTML(
                        f"<div style='flex:0 0 170px; color:#999'>"
                        f"{esc(used_group or '—')}</div>"
                    )
                cells = W.HBox(
                    [
                        W.HTML(
                            f"<div style='flex:0 0 220px'><code "
                            f"title='{esc(help_text, quote=True)}'>"
                            f"{esc(name)}</code></div>"
                        ),
                        W.HTML(
                            f"<div style='flex:0 0 90px'>{esc(default)}</div>"
                        ),
                        W.HTML(
                            f"<div style='flex:0 0 70px'>{esc(required)}</div>"
                        ),
                        W.HTML(f"<div style='flex:1'>{help_text}</div>"),
                        used_slot,
                    ],
                    layout=W.Layout(
                        display="flex", flex_flow="row", width="100%",
                        border_top="1px solid #ddd", padding="4px 0",
                    ),
                )
                cells.add_class("help-key-row")
                children.append(cells)
        if not shown:
            children.append(
                W.HTML(
                    f"<div style='padding:8px'>no key matches "
                    f"'{esc(query)}' (of {total} keys)</div>"
                )
            )
        out.children = tuple(children)
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

    def rerender() -> None:
        """Re-filter the table (used when the Features subtab opens a key)."""
        search.value = state["query"]

    state["rerender"] = rerender

    box = W.VBox(
        [
            W.HTML(
                "<h3 style='margin-top:0'>CONTROL.txt key reference</h3>"
                "<p style='color:#666'>Generated from "
                "<code>CONTROL_FIELD_SPECS</code> — the same table the "
                "validator checks your CONTROL.txt against, so it cannot drift "
                "from the code.  The help column is the same text that "
                "<code>KEY=?</code> prints.  Toggle a family below to hide it; "
                "the search filters everything.</p>"
            ),
            W.HBox([search, counter]),
            group_toggles,
            out,
        ]
    )
    return box, out, state


def create_tab(ctx: Any):
    """Create the Help tab. Returns ``(widget, refs)``."""
    header = W.HTML(
        "<h2 style='margin-top:0'>DELFIN Help</h2>"
        "<p>An overview of what DELFIN can do and which CONTROL.txt keys "
        "configure it.  Read-only: nothing here starts a calculation.</p>"
    )
    tabs = W.Tab()

    def _open_feature(group_title: str) -> None:
        """A "Used in" button was clicked: show Features, open the section."""
        for idx, child in enumerate(accordion.children):
            if getattr(child, "_feature_group", "") == group_title:
                accordion.selected_index = idx
                break
        tabs.selected_index = 0

    keys_box, _keys_html, keys_state = _key_reference_block(
        on_open_feature=_open_feature
    )

    def _open_key(name: str) -> None:
        """A key button was clicked: show the key subtab, filtered to it."""
        keys_state["query"] = name
        keys_state["rerender"]()
        tabs.selected_index = 1

    feature, accordion, feature_filter = _feature_block(_open_key)

    global_search = W.Text(
        placeholder="search features and keys…",
        layout=W.Layout(width="40%", margin="4px 0"),
    )
    _search_hint = W.HTML(
        "<span style='color:#999; font-size:12px'>searches both subtabs"
        "</span>"
    )

    def _on_global_search(change: Any) -> None:
        q = change["new"] or ""
        feature_filter(q)
        # Only touch the key subtab's filter when the user is looking at it
        # (or the query is empty), so a Features search does not clobber a
        # key the user had pinned there via a "Keys:" button.
        if tabs.selected_index == 1 or not q:
            keys_state["query"] = q
            keys_state["rerender"]()

    global_search.observe(_on_global_search, names="value")

    tabs.children = [feature, keys_box]
    tabs.set_title(0, "Features")
    tabs.set_title(1, "CONTROL keys")
    widget = W.VBox([header, W.HBox([global_search, _search_hint]), tabs])
    return widget, {"help_key_state": keys_state}


# Register additively with the dashboard tab registry.
try:
    from delfin.dashboard.tab_registry import register_tab
    register_tab("help", "Help", create_tab, order=8600)
except Exception:  # pragma: no cover
    pass


__all__ = ["_FEATURE_GROUPS", "create_tab"]
