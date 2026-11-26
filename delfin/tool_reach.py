"""Which CONTROL-switched modules the agent's tools can reach.

Package o5 item (4): the measurement package (wave 53) counts reachability, and
that count is only meaningful if it says which *route* was counted.  A CONTROL
key, a registered building block, and a tool the agent's dispatch actually
offers are three different levels, and a module can be reachable at one while
being unreachable at the next -- the "Unknown tool" defect is exactly a module
present in the registry but absent from the dispatch.

The table below is derived from the registry (``delfin.tools._registry``) and
from ``delfin/api.py``'s catalogue, so it cannot claim a route the product does
not have.  ``tests/test_the_control_connectivity_map_is_current.py`` checks the
routes named here against both.
"""

from __future__ import annotations

import json
from typing import Dict, List, Sequence, Tuple

__all__ = [
    "TOOL_ROUTES",
    "ToolRoute",
    "module_tool_reach",
    "registered_step_names",
]


class ToolRoute:
    """One CONTROL module and the building blocks / tools that reach it."""

    __slots__ = ("module", "control_switch", "steps", "agent_tools", "note")

    def __init__(self, module: str, control_switch: str,
                 steps: Sequence[str], agent_tools: Sequence[str], note: str = "") -> None:
        self.module = module
        self.control_switch = control_switch
        self.steps = tuple(steps)
        self.agent_tools = tuple(agent_tools)
        self.note = note

    def as_dict(self) -> Dict[str, object]:
        return {
            "module": self.module,
            "control_switch": self.control_switch,
            "steps": list(self.steps),
            "agent_tools": list(self.agent_tools),
            "note": self.note,
        }


#: CONTROL module -> the registered building blocks and agent tools that reach
#: it.  ``steps`` are entries of ``delfin.tools._registry.list_steps()``;
#: ``agent_tools`` are the names the agent's dispatch offers (delfin/api.py
#: catalogue and the ops server), which is what the measurement counts.
TOOL_ROUTES: Tuple[ToolRoute, ...] = (
    ToolRoute(
        module="ESD (excited-state dynamics)", control_switch="ESD_modul",
        steps=("esd", "esd_opi", "orca_tddft"),
        agent_tools=("pipeline_run",),
        note="reachable as a building block (delfin/tools/adapters/esd.py wraps "
             "execute_esd_jobs in a subprocess whose cwd is work_dir), and "
             "through a CONTROL run started with pipeline_run.",
    ),
    ToolRoute(
        module="IMAG (imaginary-mode elimination)", control_switch="IMAG",
        steps=("imag_fix", "imag_opi"),
        agent_tools=(),
        note="as a building block only; there is no separate ops-server tool "
             "name for IMAG, so an agent reaches it through a CONTROL run or "
             "the imag_fix block.",
    ),
    ToolRoute(
        module="OCCUPIER (preferred spin state)", control_switch="method=OCCUPIER",
        steps=("occupier", "occupier_opi"),
        agent_tools=(),
        note="building block only; the agent drives it by running a CONTROL "
             "file with method=OCCUPIER.",
    ),
    ToolRoute(
        module="SMILES / structure construction", control_switch="smiles_converter",
        steps=("smiles_to_xyz", "obabel_from_smiles", "guppy_sampling"),
        agent_tools=("geometry_build_smiles",),
    ),
    ToolRoute(
        module="global geometry optimisation (xTB/CREST)", control_switch="global_optimizer",
        steps=("crest_conformers", "xtb_goat", "orca_goat_opi"),
        agent_tools=(),
    ),
    ToolRoute(
        module="hyperpolarisability (sTD-DFT-xTB)", control_switch="hyperpol_xTB",
        steps=(),
        agent_tools=("hyperpol",),
        note="no registered building block; the ops server's hyperpol tool runs "
             "the workflow, which is the same entry the pipeline calls.",
    ),
    ToolRoute(
        module="TADF xTB screening", control_switch="tadf_xTB",
        steps=(),
        agent_tools=("tadf_xtb",),
        note="no registered building block; the ops server's tadf_xtb tool runs it.",
    ),
    ToolRoute(
        module="CO2 coordination", control_switch="co2_coordination",
        steps=("co2_coordinator_opi",),
        agent_tools=("co2",),
    ),
    ToolRoute(
        module="thermodynamics / stability constants", control_switch="stability_constant",
        steps=(),
        agent_tools=(),
        note="reached only by a CONTROL run: no building block and no agent tool "
             "name it. This is the gap the measurement should count, not hide.",
    ),
    ToolRoute(
        module="classic redox phase", control_switch="method=classic",
        steps=("orca_opt_freq", "orca_sp", "orca_freq"),
        agent_tools=("pipeline_run", "run_orca_input", "submit_calculation"),
    ),
)


def registered_step_names() -> List[str]:
    """The building blocks the registry currently holds."""
    from delfin.tools._registry import list_steps

    return sorted(list_steps())


def module_tool_reach() -> List[Dict[str, object]]:
    """Per CONTROL module: its routes, plus which of them currently resolve.

    ``steps_present`` is the subset of ``steps`` that the registry really has
    now, so a renamed building block shows up as a missing route instead of a
    stale claim.
    """
    present = set(registered_step_names())
    rows: List[Dict[str, object]] = []
    for route in TOOL_ROUTES:
        row = route.as_dict()
        row["steps_present"] = [s for s in route.steps if s in present]
        row["steps_missing"] = [s for s in route.steps if s not in present]
        row["reachable"] = bool(row["steps_present"] or route.agent_tools)
        rows.append(row)
    return rows


if __name__ == "__main__":  # pragma: no cover - a reporting command
    print(json.dumps(module_tool_reach(), indent=2))
