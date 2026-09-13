"""Everything DELFIN can install, and the one place that says how.

What DELFIN installs used to be listed wherever something got installed: the
shell installer, the on-demand installs in :mod:`delfin.qm_health`, the tool
registry the Tools tab reads, the CENSO auto-install of the browser workflows.
The lists drifted, and each drift was a bug somebody met -- g-xTB installable
from the structure editor and absent from Settings, a request for cclib that
installed CENSO and Packmol besides, a switch spelt ``INSTALL_TORCHANI`` on one
side and ``INSTALL_ANI2X`` on the other, the Tools tab installing into a
directory the resolver never reads.

This module is the list, and :func:`install` is the way through it. The family
installers under ``qm_tools/``, ``analysis_tools/``, ``mlp_tools/``,
``csp_tools/`` and ``ai_tools/`` still do the work; what is asked of them is
decided here.

Nothing outside the standard library is imported at module level, so the shell
installer can read the list with any Python before DELFIN's venv exists::

    python -m delfin.installer --list
    python -m delfin.installer --profile standard
    python -m delfin.installer --install crest gxtb ketcher
    python -m delfin.installer --update          # every installed tool
    python -m delfin.installer --repair          # check, and fix what is broken
"""

from __future__ import annotations

import os
import subprocess
import sys
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Tuple

#: In the order they are installed and listed.
GROUPS: Tuple[str, ...] = ("qm", "analysis", "mlp", "csp", "ai", "ketcher")

PROFILES: Tuple[str, ...] = ("core", "standard", "all")


@dataclass(frozen=True)
class Tool:
    """One installable thing.

    ``switch`` is the ``INSTALL_*`` variable of a family installer that does
    everything whose switch is not turned off; the QM installer takes names as
    arguments instead and has none. ``modules`` are the imports that prove a
    Python package is there. ``licensed`` tools are never part of a profile or
    an automatic install: they are installed when somebody names them.
    """

    name: str
    group: str
    label: str
    switch: str = ""
    modules: Tuple[str, ...] = ()
    aliases: Tuple[str, ...] = ()
    licensed: bool = False
    size: str = ""


_TORCH = "it brings PyTorch with it, which is a large download"

TOOLS: Tuple[Tool, ...] = (
    # QM programs. install_qm_tools.sh takes these names, aliases included.
    Tool("xtb", "qm", "xtb"),
    Tool("gxtb", "qm", "g-xTB", aliases=("g-xtb",)),
    Tool("mopac", "qm", "MOPAC"),
    Tool("crest", "qm", "CREST"),
    Tool("dftb+", "qm", "DFTB+", aliases=("dftbplus",)),
    Tool("xtb4stda", "qm", "xtb4stda/sTDA", aliases=("stda",)),
    Tool("std2", "qm", "std2"),
    # Analysis: install_analysis_tools.sh.
    Tool("censo", "analysis", "CENSO", switch="INSTALL_CENSO", aliases=("c2anmr", "nmrplot")),
    Tool("anmr", "analysis", "anmr", switch="INSTALL_ANMR"),
    Tool("cclib", "analysis", "cclib", switch="INSTALL_CCLIB", modules=("cclib",)),
    Tool("morfeus", "analysis", "morfeus", switch="INSTALL_MORFEUS", modules=("morfeus",),
         aliases=("morfeus-ml",)),
    Tool("nglview", "analysis", "nglview", switch="INSTALL_NGLVIEW", modules=("nglview",)),
    Tool("packmol", "analysis", "Packmol", switch="INSTALL_PACKMOL"),
    Tool("multiwfn", "analysis", "Multiwfn", switch="INSTALL_MULTIWFN", licensed=True),
    # Machine-learning potentials: install_mlp_tools.sh.
    Tool("ani2x", "mlp", "TorchANI", switch="INSTALL_ANI2X", modules=("torchani",),
         aliases=("torchani",), size=_TORCH),
    Tool("aimnet2", "mlp", "AIMNet2", switch="INSTALL_AIMNET2", modules=("aimnet2calc",),
         aliases=("aimnet2calc",), size=_TORCH),
    Tool("mace", "mlp", "MACE", switch="INSTALL_MACE", modules=("mace",), size=_TORCH),
    Tool("chgnet", "mlp", "CHGNet", switch="INSTALL_CHGNET", modules=("chgnet",), size=_TORCH),
    Tool("m3gnet", "mlp", "M3GNet", switch="INSTALL_M3GNET", modules=("matgl",),
         aliases=("matgl",), size=_TORCH),
    Tool("schnetpack", "mlp", "SchNetPack", switch="INSTALL_SCHNETPACK", modules=("schnetpack",),
         size=_TORCH),
    Tool("nequip", "mlp", "NequIP", switch="INSTALL_NEQUIP", modules=("nequip",), size=_TORCH),
    Tool("alignn", "mlp", "ALIGNN", switch="INSTALL_ALIGNN", modules=("alignn",), size=_TORCH),
    # Crystal structure prediction: install_csp_tools.sh.
    Tool("genarris", "csp", "Genarris", aliases=("gnrs",)),
    # AI tools: install_ai_tools.sh.
    Tool("molformer", "ai", "MoLFormer", switch="INSTALL_MOLFORMER", modules=("transformers",)),
    Tool("chemberta", "ai", "ChemBERTa", switch="INSTALL_CHEMBERTA", modules=("transformers",)),
    Tool("unimol", "ai", "Uni-Mol", switch="INSTALL_UNIMOL", modules=("unimol_tools",)),
    Tool("reinvent", "ai", "REINVENT4", switch="INSTALL_REINVENT", modules=("reinvent",)),
    Tool("synthemol", "ai", "SyntheMol", switch="INSTALL_SYNTHEMOL", modules=("synthemol",)),
    Tool("geomol", "ai", "GeoMol", switch="INSTALL_GEOMOL", modules=("geomol",)),
    Tool("torsional_diffusion", "ai", "torsional-diffusion", switch="INSTALL_TORSIONAL_DIFFUSION",
         modules=("torsional_diffusion",), aliases=("torsional-diffusion",)),
    Tool("mattergen", "ai", "MatterGen", switch="INSTALL_MATTERGEN", modules=("mattergen",)),
    Tool("cdvae", "ai", "CDVAE", switch="INSTALL_CDVAE", modules=("cdvae",)),
    Tool("aizynthfinder", "ai", "AiZynthFinder", switch="INSTALL_AIZYNTHFINDER",
         modules=("aizynthfinder",)),
    Tool("localretro", "ai", "LocalRetro", switch="INSTALL_LOCALRETRO", modules=("localretro",)),
    Tool("rxnmapper", "ai", "RXNMapper", switch="INSTALL_RXNMAPPER", modules=("rxnmapper",)),
    Tool("deepchem", "ai", "DeepChem", switch="INSTALL_DEEPCHEM", modules=("deepchem",)),
    Tool("admetlab", "ai", "ADMETlab", switch="INSTALL_ADMETLAB", modules=("admetlab3",)),
    Tool("molsimplify", "ai", "molSimplify", switch="INSTALL_MOLSIMPLIFY", modules=("molSimplify",)),
    Tool("architector", "ai", "architector", switch="INSTALL_ARCHITECTOR", modules=("architector",)),
    Tool("plotly", "ai", "plotly", switch="INSTALL_PLOTLY", modules=("plotly",)),
    # The structure editor, fetched into DELFIN's published directory.
    Tool("ketcher", "ketcher", "Ketcher"),
)

_BY_NAME: Dict[str, Tool] = {}
for _tool in TOOLS:
    for _key in (_tool.name, *_tool.aliases):
        _BY_NAME[_key.lower()] = _tool


class UnknownTool(KeyError):
    """A name that is not in :data:`TOOLS`."""


def find(name: str) -> Optional[Tool]:
    """The tool a name or alias stands for, or None."""
    return _BY_NAME.get(str(name or "").strip().lower())


def names(group: Optional[str] = None) -> List[str]:
    return [tool.name for tool in TOOLS if group is None or tool.group == group]


def profile(which: str) -> List[str]:
    """What a profile installs. Licensed tools are in none of them."""
    if which == "core":
        return []
    if which == "standard":
        return [tool.name for tool in TOOLS
                if not tool.licensed and tool.group in ("qm", "analysis", "ketcher")]
    if which == "all":
        return [tool.name for tool in TOOLS if not tool.licensed]
    raise KeyError(f"unknown profile: {which} (one of {', '.join(PROFILES)})")


def plan(requested: Iterable[str]) -> Dict[str, List[Tool]]:
    """The tools asked for, by group and in install order, each once."""
    wanted: List[Tool] = []
    unknown: List[str] = []
    for name in requested:
        tool = find(name)
        if tool is None:
            unknown.append(str(name))
        elif tool not in wanted:
            wanted.append(tool)
    if unknown:
        raise UnknownTool(", ".join(unknown))
    grouped: Dict[str, List[Tool]] = {}
    for group in GROUPS:
        members = [tool for tool in TOOLS if tool in wanted and tool.group == group]
        if members:
            grouped[group] = members
    return grouped


def switches(group: str) -> Tuple[str, ...]:
    """Every ``INSTALL_*`` switch the installer of *group* has."""
    return tuple(tool.switch for tool in TOOLS if tool.group == group and tool.switch)


def switch_env(group: str, wanted: Iterable[object]) -> Dict[str, str]:
    """Each switch of *group* off, and the wanted ones on.

    All of them, not only the wanted ones: these installers do everything whose
    switch is not turned off, so asking for one means naming all the others.
    """
    chosen = set()
    for item in wanted:
        tool = item if isinstance(item, Tool) else find(str(item))
        if tool is not None and tool.switch:
            chosen.add(tool.switch)
    return {switch: ("1" if switch in chosen else "0") for switch in switches(group)}


def qm_installable() -> Tuple[str, ...]:
    """Every name the QM installer accepts."""
    out: List[str] = []
    for tool in TOOLS:
        if tool.group == "qm":
            out.extend((tool.name, *tool.aliases))
    return tuple(out)


def packages() -> Dict[str, Dict[str, str]]:
    """Python modules that can be installed on demand, and what installs each."""
    out: Dict[str, Dict[str, str]] = {}
    for tool in TOOLS:
        if tool.group not in ("analysis", "mlp", "ai") or tool.licensed:
            continue
        for module in tool.modules:
            if module in out:
                continue
            spec = {"family": tool.group, "flag": tool.switch, "label": tool.label, "tool": tool.name}
            if tool.size:
                spec["size"] = tool.size
            out[module] = spec
    return out


# ---- Running the installers --------------------------------------------------

_SCRIPTS: Dict[str, Tuple[str, str, Tuple[str, ...]]] = {
    "qm": ("stage_packaged_qm_tools", "install_qm_tools.sh", ("DELFIN_QM_ROOT", "DELFIN_QM_TOOLS_ROOT")),
    "analysis": ("stage_packaged_analysis_tools", "install_analysis_tools.sh", ("DELFIN_ANALYSIS_TOOLS_ROOT",)),
    "mlp": ("stage_packaged_mlp_tools", "install_mlp_tools.sh", ("DELFIN_MLP_TOOLS_ROOT",)),
    "csp": ("stage_packaged_csp_tools", "install_csp_tools.sh", ("DELFIN_CSP_TOOLS_ROOT",)),
    "ai": ("stage_packaged_ai_tools", "install_ai_tools.sh", ("DELFIN_AI_TOOLS_ROOT",)),
}


def _stage(group: str) -> Path:
    """The user's copy of a family's installer -- where the Settings buttons install."""
    from delfin import runtime_setup

    return Path(getattr(runtime_setup, _SCRIPTS[group][0])())


def _run(command: List[str], *, cwd: Path, env: Dict[str, str],
         on_line: Optional[Callable[[str], None]], timeout: Optional[float]) -> Tuple[bool, List[str]]:
    """Run an installer, handing on each line as it is printed."""
    lines: List[str] = []
    try:
        process = subprocess.Popen(
            command, cwd=str(cwd), env=env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
        )
    except OSError as exc:
        return False, [f"the installer could not be started: {exc}"]
    timer = threading.Timer(timeout, process.kill) if timeout else None
    if timer is not None:
        timer.start()
    try:
        assert process.stdout is not None
        for raw in process.stdout:
            text = raw.rstrip("\n")
            lines.append(text)
            if on_line is not None:
                try:
                    on_line(text)
                except Exception:
                    pass
        process.wait()
    finally:
        if timer is not None:
            timer.cancel()
    if timer is not None and process.returncode is not None and process.returncode < 0:
        lines.append(f"stopped after {timeout:.0f} s")
    return process.returncode == 0, lines


def _core_constraints(root: Path) -> Optional[Path]:
    """DELFIN's own requirements as a pip constraints file, or None.

    Every pip the family installers start is held to them. Installed without,
    MACE pulled numpy 2.4 and Uni-Mol numpy 2.2 into the environment DELFIN
    runs in, which requires numpy<2 -- rdkit, pymol and mendeleev with it. A
    tool that cannot live with DELFIN's pins now fails to install and is
    reported missing, instead of breaking what everything else runs on.
    """
    try:
        from importlib.metadata import requires

        wanted = [line for line in (requires("delfin-complat") or []) if "extra ==" not in line]
    except Exception:
        return None
    if not wanted:
        return None
    target = root / "delfin_core_constraints.txt"
    try:
        target.write_text("\n".join(wanted) + "\n", encoding="utf-8")
    except OSError:
        return None
    return target


def _install_family(group: str, tools: List[Tool], *, on_line, env, timeout) -> Tuple[bool, List[str]]:
    import sysconfig

    root = _stage(group)
    script = root / _SCRIPTS[group][1]
    run_env = os.environ.copy()
    # Into the interpreter that will import what is installed.
    run_env.setdefault("DELFIN_PYTHON", sys.executable)
    # And its commands on the PATH, as they are for anybody using DELFIN. Run
    # from a shell whose venv is not activated, an installer that had just put
    # c2anmr into the venv looked for it on the PATH and said it was missing.
    scripts = sysconfig.get_path("scripts")
    if scripts:
        run_env["PATH"] = scripts + os.pathsep + run_env.get("PATH", "")
    for variable in _SCRIPTS[group][2]:
        run_env[variable] = str(root)
    constraints = _core_constraints(root)
    if constraints is not None:
        run_env["PIP_CONSTRAINT"] = str(constraints)
    run_env.update(switch_env(group, tools))
    if group == "qm" and any(tool.name == "std2" for tool in tools):
        # std2 has no binary release; the Settings button builds it too.
        run_env.setdefault("INSTALL_STD2_FROM_SOURCE", "1")
    if env:
        run_env.update({str(key): str(value) for key, value in env.items()})
    command = ["bash", str(script)]
    if group == "qm":
        command.extend(tool.name for tool in tools)
    return _run(command, cwd=root, env=run_env, on_line=on_line, timeout=timeout)


def _install_ketcher(on_line, force: bool, timeout: Optional[float]) -> Tuple[bool, List[str]]:
    from delfin.dashboard import ketcher

    lines: List[str] = []

    def say(text: str) -> None:
        lines.append(text)
        if on_line is not None:
            try:
                on_line(text)
            except Exception:
                pass

    kept = ketcher.stored_version()
    if kept and not force:
        say(f"Ketcher {kept} is already installed at {ketcher.stored_directory()}")
        return True, lines
    result = ketcher.install(on_line=say, timeout=timeout or 900.0, folder=ketcher.stored_directory())
    say(str(result.get("status") or ""))
    return bool(result.get("ok")), lines


#: What makes each family's installer update rather than keep what is there --
#: the same switches the update buttons in Settings have always used.
_UPDATE_ENV: Dict[str, Dict[str, str]] = {
    "qm": {"FORCE_REDOWNLOAD": "1", "FORCE_CONDA_UPDATE": "1"},
    "analysis": {"FORCE_REINSTALL": "1", "CENSO_PREFER_LATEST": "1"},
    "mlp": {"FORCE_REINSTALL": "1"},
    "csp": {"FORCE_REINSTALL": "1"},
    "ai": {"FORCE_REINSTALL": "1"},
}


def install(requested: Iterable[str], *, on_line: Optional[Callable[[str], None]] = None,
            env: Optional[Dict[str, str]] = None, timeout: Optional[float] = None,
            force: bool = False, update: bool = False) -> Dict[str, object]:
    """Install the named tools, each family's installer run once.

    ``update`` fetches and reinstalls what is already there instead of keeping
    it. Returns ``{'ok': bool, 'results': [{'group', 'tools', 'ok', 'lines'}]}``.
    One family failing does not stop the next: a missing compiler for Genarris
    is no reason to leave CREST uninstalled.
    """
    results: List[Dict[str, object]] = []
    for group, tools in plan(requested).items():
        if group == "ketcher":
            ok, lines = _install_ketcher(on_line, force or update, timeout)
        else:
            group_env = dict(_UPDATE_ENV.get(group, {})) if update else {}
            group_env.update(env or {})
            ok, lines = _install_family(group, tools, on_line=on_line, env=group_env, timeout=timeout)
        # Asked afterwards, not taken from the exit code: these installers say
        # "Packmol installation requires conda" and still return 0.
        missing = [tool.name for tool in tools if not present(tool)]
        if missing:
            ok = False
            lines = list(lines) + [f"not installed after its installer ran: {' '.join(missing)}"]
            _say(on_line, lines[-1])
        results.append({"group": group, "tools": [tool.name for tool in tools], "ok": ok,
                        "missing": missing, "lines": lines})
    return {"ok": all(result["ok"] for result in results), "results": results}


# ---- What is there, updating it, putting it right ----------------------------

def _say(on_line: Optional[Callable[[str], None]], text: str) -> None:
    if on_line is not None:
        try:
            on_line(text)
        except Exception:
            pass


def _module_present(module: str) -> bool:
    import importlib.util

    try:
        return importlib.util.find_spec(module) is not None
    except (ImportError, ValueError):
        return False


def _probe_name(tool: Tool) -> str:
    """The name qm_health checks a program under."""
    return "gnrs" if tool.name == "genarris" else tool.name


def present(tool: Tool) -> bool:
    """Whether *tool* is installed where DELFIN looks. Found, nothing run."""
    if tool.group == "ketcher":
        from delfin.dashboard import ketcher

        return bool(ketcher.stored_version())
    if tool.modules:
        return any(_module_present(module) for module in tool.modules)
    import sysconfig

    scripts = sysconfig.get_path("scripts")
    if scripts and os.access(os.path.join(scripts, _probe_name(tool)), os.X_OK):
        return True
    from delfin import qm_health

    return bool(qm_health.check_tool(_probe_name(tool), depth="present").present)


def installed() -> List[str]:
    return [tool.name for tool in TOOLS if present(tool)]


def update(requested: Optional[Iterable[str]] = None, *,
           on_line: Optional[Callable[[str], None]] = None,
           timeout: Optional[float] = None) -> Dict[str, object]:
    """Update the named tools, or every tool that is installed.

    Only what is there when nothing is named: an update that installed forty
    tools nobody had asked for would be an install.
    """
    wanted = list(requested) if requested else installed()
    if not wanted:
        return {"ok": True, "results": []}
    return install(wanted, on_line=on_line, timeout=timeout, update=True)


def _imports(module: str, timeout: float = 300.0) -> Tuple[bool, str]:
    """Whether *module* imports in a fresh interpreter -- found, not guessed."""
    try:
        done = subprocess.run([sys.executable, "-c", f"import {module}"],
                              capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return False, f"importing it did not finish in {timeout:.0f} s"
    said = (done.stderr or "").strip().splitlines()
    return done.returncode == 0, (said[-1] if said else "")


def _checked(tool: Tool, ok: bool, action: str, status: str,
             lines: Optional[List[str]] = None) -> Dict[str, object]:
    return {"group": tool.group, "tools": [tool.name], "ok": ok, "action": action,
            "status": status, "lines": list(lines or [])}


def _repair_one(tool: Tool, named: bool, *, on_line, timeout) -> Optional[Dict[str, object]]:
    def reinstall(why: str) -> Dict[str, object]:
        _say(on_line, f"{tool.label}: {why}; installing it again")
        outcome = install([tool.name], on_line=on_line, timeout=timeout, update=True)
        lines = [line for result in outcome["results"] for line in result["lines"]]
        ok = bool(outcome["ok"]) and present(tool)
        return _checked(tool, ok, "reinstalled", why, lines)

    if tool.group == "ketcher":
        if present(tool):
            return _checked(tool, True, "checked", f"{tool.label} is installed")
        return reinstall("not installed") if named else None

    if tool.modules:
        found = [module for module in tool.modules if _module_present(module)]
        if not found:
            return reinstall("not installed") if named else None
        ok, why = _imports(found[0])
        if ok:
            return _checked(tool, True, "checked", f"{tool.label} imports")
        return reinstall(f"{found[0]} does not import" + (f" ({why})" if why else ""))

    from delfin import qm_health

    probe = _probe_name(tool)
    health = qm_health.check_tool(probe, depth="answer")
    if not health.present and health.level != "fail":
        return reinstall("not installed") if named else None
    if health.repair:
        _say(on_line, f"{tool.label}: {health.why}; "
                      f"{qm_health.REPAIRS.get(health.repair, health.repair)}")
        answer = qm_health.repair_tool(probe, health.repair, on_line=on_line,
                                       timeout=timeout or 1800.0)
        return _checked(tool, bool(answer.get("ok")), "repaired",
                        str(answer.get("status") or ""), answer.get("lines"))
    if health.healthy or health.level in ("ok", "warn"):
        return _checked(tool, True, "checked", health.why or f"{tool.label} works")
    if tool.licensed and not named:
        return _checked(tool, False, "", f"{tool.label}: {health.why}. It is licensed, so it "
                                         f"is only reinstalled when named ({tool.name}).")
    return reinstall(health.why or "it does not work")


def repair(requested: Optional[Iterable[str]] = None, *,
           on_line: Optional[Callable[[str], None]] = None,
           timeout: Optional[float] = None) -> Dict[str, object]:
    """Check each tool and put right what does not work.

    The checks and repairs of the programs are qm_health's -- a link into a
    deleted environment repointed, an xtb that cannot optimise replaced -- and
    each is proven afterwards there. A Python package that is present and does
    not import is reinstalled. A tool that is simply not installed is left
    alone unless it was named: repairing is not installing everything.
    """
    grouped = plan(requested) if requested else None
    scope = [tool for tools in grouped.values() for tool in tools] if grouped else list(TOOLS)
    results: List[Dict[str, object]] = []
    for tool in scope:
        outcome = _repair_one(tool, grouped is not None, on_line=on_line, timeout=timeout)
        if outcome is not None:
            results.append(outcome)
    return {"ok": all(result["ok"] for result in results), "results": results}


# ---- Command line ------------------------------------------------------------

def _print_plan(grouped: Dict[str, List[Tool]]) -> None:
    for group, tools in grouped.items():
        print(f"{group}: {' '.join(tool.name for tool in tools)}")


def main(argv: Optional[List[str]] = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(prog="python -m delfin.installer",
                                     description="List and install what DELFIN uses.")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--list", action="store_true", help="every installable tool, by group")
    action.add_argument("--plan", metavar="TOOLS", help="the groups the named tools fall into")
    action.add_argument("--profile", choices=PROFILES, help="the tools a profile installs")
    action.add_argument("--install", nargs="+", metavar="TOOL", help="install these tools")
    action.add_argument("--update", nargs="*", metavar="TOOL",
                        help="update these tools, or every installed one")
    action.add_argument("--repair", nargs="*", metavar="TOOL",
                        help="check these tools, or every installed one, and fix what is broken")
    action.add_argument("--status", action="store_true", help="which tools are installed")
    parser.add_argument("--force", action="store_true", help="fetch Ketcher again even if it is there")
    args = parser.parse_args(argv)

    if args.list:
        for group in GROUPS:
            print(f"{group}: {' '.join(names(group))}")
        licensed = [tool.name for tool in TOOLS if tool.licensed]
        if licensed:
            print(f"licensed, installed only when named: {' '.join(licensed)}")
        return 0
    if args.status:
        for tool in TOOLS:
            state = "installed" if present(tool) else "missing  "
            print(f"{state}  {tool.group:<9} {tool.name}", flush=True)
        return 0
    try:
        if args.profile:
            _print_plan(plan(profile(args.profile)))
            return 0
        if args.plan is not None:
            _print_plan(plan(args.plan.replace(",", " ").split()))
            return 0
        if args.update is not None:
            outcome = update(args.update, on_line=print)
        elif args.repair is not None:
            outcome = repair(args.repair, on_line=print)
        else:
            outcome = install(args.install, on_line=print, force=args.force)
    except UnknownTool as exc:
        print(f"unknown tool: {exc.args[0]} (see --list)", file=sys.stderr)
        return 2
    if not outcome["results"]:
        print("[delfin-install] nothing to do", flush=True)
    for result in outcome["results"]:
        mark = "ok    " if result["ok"] else "FAILED"
        said = f" ({result['action']}: {result['status']})" if result.get("action") else ""
        if result.get("missing"):
            said += f" (missing: {' '.join(result['missing'])})"
        print(f"[delfin-install] {mark} {result['group']}: {' '.join(result['tools'])}{said}", flush=True)
    return 0 if outcome["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
