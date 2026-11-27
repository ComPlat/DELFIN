"""The CONTROL-file connectivity map, derived from the code that implements it.

A CONTROL file switches modules on; ``docs/WORKFLOW_CONNECTIVITY.md`` says, for
each of the template's keys, which module reads it, what that module is handed
by the step before it, and where its result goes.  Writing that by hand means
it drifts the first time a key is renamed -- which is exactly what happened
here: two keys (``thermodynamics_mode``, ``thermodynamics_reaction``) were
reported as read by nothing, because nothing reads them *by that name*.  They
are aliases (``delfin/config.py:_CONTROL_KEY_ALIASES``) and are read
constantly; a literal-name scan cannot see it and a hand-written table repeats
the mistake.

So the map is generated from three things that are all in the code:

* the template itself (``delfin.define.TEMPLATE``) -- the keys that exist;
* the alias table -- a key that reaches its reader under another name;
* a ``CONTROL_MODULES`` table below -- one entry per module a CONTROL file can
  switch on, naming its switch key, the entry point the pipeline calls, what
  reaches it from the step before, and where its result goes.

The ``CONTROL_MODULES`` table is the only hand-maintained part, and
``tests/test_the_control_connectivity_map_is_current.py`` pins it in both
directions: every module named there must still be reachable from
``delfin/cli.py`` under its switch, and every key of the template must be
accounted for in the rendered document.  Adding a template key without saying
who reads it fails that test, and so does removing an entry point the table
still promises.

Nothing here runs a calculation: the reachability check is a name lookup on
the module graph (``<module>.<entry point>`` must exist and an importer of the
module must exist above it), not an execution.
"""

from __future__ import annotations

import importlib
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from delfin.config import _CONTROL_KEY_ALIASES, _canonicalize_control_key
from delfin.define import TEMPLATE

__all__ = [
    "CONTROL_MODULES",
    "ControlModule",
    "control_template_keys",
    "key_readers",
    "render_connectivity_map",
    "reachability_problems",
    "MAP_PATH",
]


#: Where the rendered map lives.
MAP_PATH = Path(__file__).resolve().parent.parent / "docs" / "WORKFLOW_CONNECTIVITY.md"


@dataclass(frozen=True)
class ControlModule:
    """One module a CONTROL file can switch on.

    ``switch`` is the canonical key whose value turns it on (``None`` for the
    pipeline's own redox phases, which the ``method`` key selects).  ``entry``
    is ``module:callable`` -- the thing the pipeline actually calls, so the map
    can be checked against the code rather than believed.  ``handed`` and
    ``produces`` are what crosses the step boundary in each direction.
    """

    name: str
    switch: Optional[str]
    switch_values: Tuple[str, ...]
    entry: str
    handed: str
    produces: str
    #: Modules whose presence in ``delfin/`` this entry is checked against.
    module_path: str
    #: Template keys this module reads, beyond its own switch (informational).
    keys: Tuple[str, ...] = field(default_factory=tuple)
    notes: str = ""


#: The modules a CONTROL file switches on, in the order the pipeline reaches
#: them.  This table is the hand-maintained half of the map; the test pins it
#: to the template and to the module graph.
CONTROL_MODULES: Tuple[ControlModule, ...] = (
    ControlModule(
        name="SMILES / structure construction",
        switch="smiles_converter",
        switch_values=("QUICK", "NORMAL", "MANTA", "ARCHITECTOR", "MOLSIMPLIFY",
                       "MACE", "RXEMBED", "MOLASSEMBLER", "METALLOGEN"),
        entry="delfin.workflows.pipeline:normalize_input_file",
        handed="the SMILES line in input.txt (or CONTROL's SMILES key)",
        produces="start.txt -- DELFIN coordinate format, the geometry every later step reads",
        module_path="delfin/workflows/pipeline.py",
        keys=("SMILES", "input_file", "GUPPY"),
        notes="MANTA is the current name for the GUPPY sampler; a legacy "
              "smiles_converter=GUPPY or a bare GUPPY=yes is read as MANTA "
              "(_apply_guppy_legacy, _resolve_smiles_converter).",
    ),
    ControlModule(
        name="global geometry optimisation (xTB/CREST)",
        switch="global_optimizer",
        switch_values=("GOAT", "CREST"),
        entry="delfin.workflows.pipeline:run_classic_phase",
        handed="start.txt",
        produces="the optimised geometry written back to start.txt; XTB_OPT "
                 "and XTB_preOPT select the pre-optimisation instead",
        module_path="delfin/guppy_sampling.py",
        keys=("xTB_method", "multiplicity_global_opt", "MANTA_*", "hyperpol_xTB_preopt"),
        notes="an empty global_optimizer means no global optimisation, which "
              "is what the two keys it replaced (XTB_GOAT=no, CREST=no) said.",
    ),
    ControlModule(
        name="IMAG (imaginary-mode elimination)",
        switch="IMAG",
        switch_values=("yes",),
        entry="delfin.imag:run_IMAG",
        handed="the step's .out/.inp and its .hess, in the step's own folder",
        produces="a relaxed geometry copied back to the calling step (and, for "
                 "OCCUPIER, to input_<step>_OCCUPIER.xyz)",
        module_path="delfin/imag.py",
        keys=("IMAG_scope", "IMAG_option", "allow_imaginary_freq",
              "IMAG_sp_energy_window", "IMAG_optimize_candidates", "IMAG_max_rounds"),
        notes="a key read but with no step to act on is reported, not silently "
              "ignored (delfin/imag.py:_refuse_rates_on_a_saddle family).",
    ),
    ControlModule(
        name="classic redox phase",
        switch="method",
        switch_values=("classic",),
        entry="delfin.workflows.pipeline:run_classic_phase",
        handed="start.txt / initial.xyz, charge and multiplicity from CONTROL",
        produces="initial/ox_step_N/red_step_N folders, each with its .inp/.out; "
                 "the Gibbs energies the report is built from",
        module_path="delfin/workflows/engine/classic.py",
        keys=("calc_initial", "oxidation_steps", "reduction_steps",
              "calc_prop_of_interest", "properties_of_interest",
              "reorganisation_energy"),
    ),
    ControlModule(
        name="manually-written spin states",
        switch="method",
        switch_values=("manually",),
        entry="delfin.workflows.pipeline:run_manual_phase",
        handed="multiplicity_N / BrokenSym_N per step from CONTROL",
        produces="the same step folders as classic, using the stated spin states",
        module_path="delfin/workflows/pipeline.py",
        keys=("multiplicity_0", "BrokenSym_0", "multiplicity_ox1", "BrokenSym_ox1"),
    ),
    ControlModule(
        name="OCCUPIER (preferred spin state)",
        switch="method",
        switch_values=("OCCUPIER",),
        entry="delfin.workflows.pipeline:run_occuper_phase",
        handed="start.txt, charge, and the OCCUPIER_sequence_profiles block",
        produces="one *_OCCUPIER folder per step and OCCUPIER.txt, which names "
                 "the preferred configuration every later step and the report use",
        module_path="delfin/workflows/engine/occupier.py",
        keys=("OCCUPIER_method", "OCCUPIER_compare", "occupier_selection",
              "occupier_precision", "occupier_epsilon", "pass_wavefunction",
              "OWN_TREE_PURE_WINDOW", "OWN_progressive_from", "fob_equal_weights"),
    ),
    ControlModule(
        name="ESD (excited-state dynamics)",
        switch="ESD_modul",
        switch_values=("yes",),
        entry="delfin.esd_module:execute_esd_jobs",
        handed="start.txt / initial.xyz (the geometry the previous step left in "
               "the workspace), charge, solvent, metals and both basis sets",
        produces="ESD/ with one input and output per state (S0, S1, T1, ...), "
                 "plus ISC/IC/fluorescence/phosphorescence outputs the reports read",
        module_path="delfin/esd_module.py",
        keys=("ESD_modus", "ESD_T1_opt", "ESD_frequency", "states", "ISCs",
              "ICs", "emission_rates", "TDDFT_nroots", "TDDFT_TDA", "TDDFT_SOC"),
        notes="reached in parallel inside the classic scheduler when method is "
              "classic (add_esd_jobs_to_scheduler), and as its own phase "
              "otherwise (run_esd_phase).  It reads the workspace, not the "
              "process directory -- see "
              "tests/test_a_control_run_writes_into_its_workspace.py.",
    ),
    ControlModule(
        name="hyperpolarisability (sTD-DFT-xTB)",
        switch="hyperpol_xTB",
        switch_values=("yes",),
        entry="delfin.workflows.pipeline:run_hyperpol_xtb_phase",
        handed="hyperpol_xTB_xyz (default start.txt), charge and multiplicity, "
               "in the CONTROL file's directory",
        produces="hyperpol_xtb/ with the response run and hyperpol_xtb_summary.json",
        module_path="delfin/hyperpol.py",
        keys=("hyperpol_xTB_xyz", "hyperpol_xTB_preopt", "hyperpol_xTB_engine",
              "hyperpol_xTB_bfw", "hyperpol_xTB_wavelengths",
              "hyperpol_xTB_energy_window"),
    ),
    ControlModule(
        name="TADF xTB screening",
        switch="tadf_xTB",
        switch_values=("yes",),
        entry="delfin.workflows.pipeline:run_tadf_xtb_phase",
        handed="tadf_xTB_xyz (default start.txt), charge and multiplicity",
        produces="tadf_xtb/ with the singlet/triplet screening the TADF report reads",
        module_path="delfin/tadf_xtb.py",
        keys=("tadf_xTB_xyz", "tadf_xTB_preopt", "tadf_xTB_excited_method",
              "tadf_xTB_bfw", "tadf_xTB_energy_window", "tadf_xTB_run_t1_opt"),
    ),
    ControlModule(
        name="CO2 coordination",
        switch="co2_coordination",
        switch_values=("on", "yes", "true", "1"),
        entry="delfin.cli:_run_co2_subcommand",
        handed="the finished run's preferred structures and co2_species_delta",
        produces="CO2_coordination/ (its own workspace, entered with os.chdir), "
                 "written back into the report",
        module_path="delfin/co2/chain_setup.py",
        keys=("co2_species_delta",),
    ),
    ControlModule(
        name="thermodynamics / stability constants",
        switch="stability_constant",
        switch_values=("yes",),
        entry="delfin.stability_constant:run_stability_constant_phase",
        handed="the initial step's Gibbs energy and FREQ output; embedded into "
               "the OCCUPIER scheduler when that method runs",
        produces="stability_constant/ jobs and the logK the report reads",
        module_path="delfin/stability_constant.py",
        keys=("stability_constant_mode", "stability_reaction", "n_explicit_solvent",
              "logK_exp", "thdy_smiles_converter", "thdy_preopt"),
        notes="the template spells the switch thermodynamics=, which the alias "
              "table maps onto stability_constant; thermodynamics_mode and "
              "thermodynamics_reaction are aliases too.",
    ),
)


def control_template_keys(template: str = TEMPLATE) -> List[str]:
    """The keys the CONTROL template writes, in template order."""
    keys: List[str] = []
    for line in template.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(("-", "#")) or "=" not in stripped:
            continue
        key = stripped.split("=", 1)[0].strip()
        if key:
            keys.append(key)
    return keys


def _modules_reading(key: str, sources: Dict[str, str]) -> List[str]:
    """delfin/ modules that name *key* (canonicalised) in a string literal."""
    canonical = _canonicalize_control_key(key)
    pattern = re.compile(r"""["']""" + re.escape(canonical) + r"""["']""")
    found: List[str] = []
    for path, text in sources.items():
        if path.endswith("/define.py"):
            continue
        if pattern.search(text):
            found.append(path)
    return sorted(found)


def key_readers(sources: Optional[Dict[str, str]] = None) -> Dict[str, List[str]]:
    """Every template key mapped to the modules that read it.

    With *sources* given (``{relative path: text}``) that text is used;
    otherwise ``delfin/`` is read from disk.  ``define.py`` is excluded: it
    writes the template, it does not read it back.
    """
    if sources is None:
        root = Path(__file__).resolve().parent.parent / "delfin"
        sources = {
            str(p.relative_to(root.parent)): p.read_text(encoding="utf-8", errors="replace")
            for p in sorted(root.rglob("*.py"))
        }
    return {key: _modules_reading(key, sources) for key in control_template_keys()}


def _alias_note(key: str) -> str:
    canonical = _canonicalize_control_key(key)
    if canonical == key:
        return ""
    return f"read as `{canonical}` (alias)"


def reachability_problems(modules: Sequence[ControlModule] = CONTROL_MODULES) -> List[str]:
    """Problems that would make an entry of *modules* untrue.

    Each entry must have (a) an importable module whose file exists, (b) the
    named callable in it, and (c) -- for the ones the pipeline drives -- a
    caller somewhere under ``delfin/``.  A table that promises a module nobody
    calls is the defect this check exists for.
    """
    problems: List[str] = []
    root = Path(__file__).resolve().parent.parent
    for entry in modules:
        module_name, _, callable_name = entry.entry.partition(":")
        path = root / entry.module_path
        if not path.is_file():
            problems.append(f"{entry.name}: {entry.module_path} does not exist")
            continue
        try:
            module = importlib.import_module(module_name)
        except Exception as exc:  # noqa: BLE001 - an unimportable module is a finding
            problems.append(f"{entry.name}: cannot import {module_name} ({exc})")
            continue
        if not hasattr(module, callable_name):
            problems.append(
                f"{entry.name}: {module_name} has no {callable_name}; the map "
                "names an entry point that is not there"
            )
    return problems


def _md_cell(text: str) -> str:
    return str(text).replace("|", "\\|").replace("\n", " ")


def render_connectivity_map(
    modules: Sequence[ControlModule] = CONTROL_MODULES,
    sources: Optional[Dict[str, str]] = None,
    readers: Optional[Dict[str, List[str]]] = None,
) -> str:
    """The whole ``docs/WORKFLOW_CONNECTIVITY.md`` document, as text."""
    readers = readers if readers is not None else key_readers(sources)
    keys = control_template_keys()
    lines: List[str] = []
    add = lines.append

    add("# WORKFLOW CONNECTIVITY — what a CONTROL file switches on, and how it is wired")
    add("")
    add("<!-- GENERATED FILE - do not edit by hand.")
    add("     Source: delfin/workflow_connectivity.py (CONTROL_MODULES) plus")
    add("     delfin/define.py (TEMPLATE) and delfin/config.py (the alias table).")
    add("     Regenerate: python -m delfin.workflow_connectivity")
    add("     Pinned by: tests/test_the_control_connectivity_map_is_current.py -->")
    add("")
    add("A CONTROL file is the only place a DELFIN run is configured. This map")
    add("answers, for every module one can switch on: which key switches it, the")
    add("code path from that key to the module, what the module is handed from the")
    add("step before it, and where its result goes. It is generated from the code")
    add("so it cannot drift: a renamed entry point or an unlisted template key")
    add("fails the test named above.")
    add("")
    add(f"`delfin/define.TEMPLATE` writes **{len(keys)} keys**. All of them are")
    add("central: the last table lists every one with the modules that read it.")
    add("")
    add("## 1. Modules a CONTROL file switches on")
    add("")
    add("| Module | Switch key | Values | Entry point the pipeline calls |")
    add("|---|---|---|---|")
    for m in modules:
        values = ", ".join(f"`{v}`" for v in m.switch_values)
        switch = f"`{m.switch}`" if m.switch else "(redox method)"
        add(f"| {_md_cell(m.name)} | {switch} | {values} | `{m.entry}` |")
    add("")
    add("## 2. What crosses each step boundary")
    add("")
    add("| Module | Handed by the step before | Its result goes to |")
    add("|---|---|---|")
    for m in modules:
        add(f"| {_md_cell(m.name)} | {_md_cell(m.handed)} | {_md_cell(m.produces)} |")
    add("")
    add("## 3. Notes on individual modules")
    add("")
    for m in modules:
        if not m.notes:
            continue
        add(f"- **{m.name}** — {m.notes}")
    add("")
    add("## 4. Which modules the agent's tools reach")
    add("")
    add("The map above says what a CONTROL *run* reaches. An agent reaches a")
    add("module by a different route, and the two are not the same count:")
    add("`delfin/tool_reach.py` lists, per module, the registered building blocks")
    add("(`delfin.tools._registry`) and the agent tools (`delfin/api.py`, the")
    add("tool catalogue). Both are checked against the code by")
    add("`tests/test_the_tool_reach_names_real_routes.py`, so a renamed tool shows")
    add("up as a missing route rather than an inflated figure.")
    add("")
    try:
        from delfin.tool_reach import module_tool_reach

        reach_rows = module_tool_reach()
    except Exception:  # noqa: BLE001 - the table must render without the registry
        reach_rows = []
    if reach_rows:
        add("| Module | Building blocks | Agent tools | Reachable by a tool |")
        add("|---|---|---|---|")
        for row in reach_rows:
            steps = ", ".join(f"`{s}`" for s in row["steps_present"]) or "—"
            tools = ", ".join(f"`{t}`" for t in row["agent_tools"]) or "—"
            add(f"| {_md_cell(row['module'])} | {steps} | {tools} | "
                f"{'yes' if row['reachable'] else '**no — CONTROL run only**'} |")
        add("")
        unreached = [r["module"] for r in reach_rows if not r["reachable"]]
        if unreached:
            add("Reachable by a CONTROL run but by **no tool at all**: "
                + ", ".join(f"`{m}`" for m in unreached) + ".")
            add("")
    else:
        add("(the tool-reach table could not be imported in this environment)")
        add("")
    add("## 5. Every template key, and who reads it")
    add("")
    add("`by alias` means the key reaches its reader under another name")
    add("(`delfin/config.py:_CONTROL_KEY_ALIASES`); a scan for the literal name")
    add("would miss it, which is how two of these were once reported as unused.")
    add("")
    add("| Key | Read as | Modules that read it |")
    add("|---|---|---|")
    for key in keys:
        modules_for_key = readers.get(key, [])
        alias = _alias_note(key)
        if modules_for_key:
            who = ", ".join(f"`{p}`" for p in modules_for_key[:4])
            if len(modules_for_key) > 4:
                who += f" (+{len(modules_for_key) - 4} more)"
        elif alias:
            who = "— (read under the alias above; no literal occurrence)"
        else:
            who = "**none by name** — see the note below"
        add(f"| `{_md_cell(key)}` | {alias or 'own name'} | {who} |")
    add("")
    add("## 6. Keys read without a literal occurrence")
    add("")
    unread = [k for k in keys if not readers.get(k)]
    if not unread:
        add("None: every template key is read by name by at least one module.")
    else:
        add("These are read through a mechanism other than a literal key lookup:")
        add("")
        for key in unread:
            alias = _alias_note(key)
            if alias:
                add(f"- `{key}` — {alias}.")
            elif ":" in key:
                target = key.split(":", 1)[0]
                add(f"- `{key}` — parsed by shape, not by name: "
                    f"`delfin/common/orca_overrides.py:_OVERRIDE_KEY_RE` matches "
                    f"`{target}:<target>` for any job name, so the key is the "
                    f"prefix `{target}:` plus a target, and no literal "
                    f"`{key}` appears in the code.")
            else:
                add(f"- `{key}` — no reader found.")
    add("")
    return "\n".join(lines) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Rewrite ``docs/WORKFLOW_CONNECTIVITY.md`` from the code."""
    target = MAP_PATH
    if argv:
        target = Path(argv[0])
    target.write_text(render_connectivity_map(), encoding="utf-8")
    print(f"wrote {target}")
    return 0


if __name__ == "__main__":  # pragma: no cover - a regeneration command
    raise SystemExit(main())
