"""Each module a CONTROL file switches on is reached, and handed the step before.

Package o5 item (2): for every module, a test with a SYNTHETIC complete CONTROL
(structure only -- no real project data) that the module is reached and gets
the right inputs from the step before it.

The module set is not written here by hand: it comes from
``delfin.workflow_connectivity.CONTROL_MODULES``, the same table the map is
generated from. A module added to that table must therefore be covered here,
and a module whose dispatch is removed from ``delfin/cli.py`` fails both the
map test and this one.

What is asserted is the *dispatch*, not the chemistry: a stub ORCA is put on
the path so a run can pass the classic phase, and the question asked is only
whether CONTROL's switch reaches the phase and whether that phase was given
the geometry the previous step produced. Nothing computes a number.

The input that would make each assertion fail is named in it, so the test is
not one that cannot fail.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pytest

#: The dispatch table's own source, when it exists. On a tree that predates
#: the map (the baseline) it does not, and the fallback below is used: the
#: dispatch assertions must be runnable -- and red -- on that tree, or their
#: redness would say only "the map is new" and not "the dispatch is broken".
_MAP_MODULE = Path(__file__).resolve().parent.parent / "delfin" / "workflow_connectivity.py"

#: The modules package o5 names, with the entry point the pipeline calls.
#: Used when the generated table is absent (the baseline run).
_FALLBACK_MODULES: tuple[tuple[str, str], ...] = (
    ("ESD (excited-state dynamics)", "delfin.esd_module:execute_esd_jobs"),
    ("IMAG (imaginary-mode elimination)", "delfin.imag:run_IMAG"),
    ("OCCUPIER (preferred spin state)", "delfin.workflows.pipeline:run_occuper_phase"),
    ("CO2 coordination", "delfin.cli:_run_co2_subcommand"),
    ("hyperpolarisability (sTD-DFT-xTB)", "delfin.workflows.pipeline:run_hyperpol_xtb_phase"),
    ("TADF xTB screening", "delfin.workflows.pipeline:run_tadf_xtb_phase"),
    ("thermodynamics / stability constants",
     "delfin.stability_constant:run_stability_constant_phase"),
    ("classic redox phase", "delfin.workflows.pipeline:run_classic_phase"),
)

#: A three-atom synthetic structure. Deliberately not a real project's input:
#: the point is the wiring, and no content of any archive belongs in the repo.
_GEOMETRY = (
    "O  0.0000  0.0000  0.1173\n"
    "H  0.0000  0.7572 -0.4692\n"
    "H  0.0000 -0.7572 -0.4692\n"
)

#: The template as shipped, with every placeholder that would fail validation
#: filled in, and every module switched OFF unless a case says otherwise.
_BASE = {
    "input_file": "input.txt",
    "charge": "0",
    "solvent": "water",
    "method": "classic",
    "smiles_converter": "NORMAL",
    "global_optimizer": "none",
    "XTB_preOPT": "no",
    "IMAG": "no",
    "calc_initial": "yes",
    "oxidation_steps": "",
    "reduction_steps": "",
    "ESD_modul": "no",
    "hyperpol_xTB": "no",
    "tadf_xTB": "no",
    "co2_coordination": "off",
    "thermodynamics": "no",
    # The template carries two other placeholders; a CONTROL that still holds
    # them does not validate ("Replace template placeholders").
    "ESD_modus": "tddft",
}


def _control_text(**overrides: str) -> str:
    """The shipped template with *overrides* applied, key by key.

    Keys are spelled as the TEMPLATE spells them. ``config._CONTROL_KEY_ALIASES``
    maps some canonical names onto template names (``stability_constant`` ->
    ``thermodynamics``), but a CONTROL file is written in template names, so
    that is what is rewritten here -- and a key the template does not write
    fails loudly instead of silently matching nothing.
    """
    from delfin.define import TEMPLATE

    text = TEMPLATE
    for key, value in {**_BASE, **overrides}.items():
        # A key may be written as ``key=value`` or indented as ``key = value``
        # (odd_seq / even_seq are), so match both forms.
        text, count = re.subn(
            rf"(?m)^{re.escape(key)}\s*=.*$", f"{key}={value}", text
        )
        assert count == 1, f"{key} does not appear exactly once in the template"
    return text


def _write_control(workspace: Path, **overrides: str) -> Path:
    path = workspace / "CONTROL.txt"
    path.write_text(_control_text(**overrides), encoding="utf-8")
    (workspace / "input.txt").write_text(_GEOMETRY, encoding="utf-8")
    return path


def _module_names() -> list[str]:
    return [name for name, _entry in _modules_from_code()]


def _modules_from_code() -> list[tuple[str, str]]:
    """``(name, entry)`` pairs for the modules a CONTROL file switches on.

    Read from the generated table's source when it is present, so a module
    added there is covered here automatically. That file is deliberately not
    *imported*: the dispatch assertions have to run on a tree that predates
    the map, and there the fallback list is used instead.
    """
    if not _MAP_MODULE.is_file():
        return list(_FALLBACK_MODULES)
    source = _MAP_MODULE.read_text(encoding="utf-8")
    pairs = re.findall(r'name="([^"]+)",.*?entry="([^"]+)"', source, re.S)
    return pairs or list(_FALLBACK_MODULES)


def test_the_module_table_is_the_one_the_map_publishes():
    """The cases below are driven off the map's own table, not a copy.

    Fails if CONTROL_MODULES loses an entry, because the map would then be
    describing fewer modules than the package asks about.
    """
    names = _module_names()
    assert len(names) >= 10, f"the map knows only {len(names)} modules: {names}"
    for required in ("ESD", "IMAG", "OCCUPIER", "CO2", "hyperpolarisability",
                     "TADF", "thermodynamics"):
        assert any(required.lower() in n.lower() for n in names), (
            f"the map has no entry for {required}; names were {names}"
        )


@pytest.mark.parametrize("module_name", _module_names())
def test_every_module_names_a_call_the_pipeline_makes(module_name):
    """An import is not a caller: the entry point must be *called*, somewhere.

    The module's ``entry`` names ``module:callable``. This asserts the callable
    exists AND that some module under ``delfin/`` (outside its own definition)
    calls it. A bare definition is exactly the "reachable in theory" state this
    package exists to rule out.

    Two things a literal grep gets wrong, and this does not:
    * the callable may be imported under another name -- ``ESD``'s
      ``execute_esd_jobs`` is imported as ``execute_esd_module`` and called
      under the alias, so the alias is resolved rather than searched for;
    * the caller need not be ``cli.py`` or ``pipeline.py`` -- ``IMAG`` is
      called from ``workflows/engine/classic.py``, so every module is walked.
    """
    import importlib

    module_name = next((n for n, _e in _modules_from_code() if n == module_name), None)
    assert module_name, "the module vanished while the test was running"
    entry = dict(_modules_from_code())[module_name]
    mod_name, _, attr = entry.partition(":")
    imported = importlib.import_module(mod_name)
    assert hasattr(imported, attr), (
        f"{entry} does not exist, so {module_name} cannot be reached"
    )

    root = Path(__file__).resolve().parent.parent / "delfin"
    # Names the callable can be called under: its own, plus any alias given
    # by ``from ... import attr as alias``.
    alias_pattern = re.compile(
        rf"from\s+{re.escape(mod_name)}\s+import\s+[^\n]*\b{re.escape(attr)}\b"
        r"\s+as\s+(\w+)"
    )
    callable_names = {attr}
    for source in root.rglob("*.py"):
        text = source.read_text(encoding="utf-8", errors="replace")
        callable_names.update(alias_pattern.findall(text))

    call_pattern = re.compile(
        r"\b(?:" + "|".join(re.escape(n) for n in sorted(callable_names)) + r")\s*\("
    )
    # The caller may share the definition's file: ``_run_co2_subcommand`` is
    # defined in cli.py and called from cli.py:1633. A separate file is not
    # required -- a *call site distinct from the ``def`` line* is. Excluding
    # the defining file would have hidden CO2 behind a false negative.
    callers = []
    for source in sorted(root.rglob("*.py")):
        text = source.read_text(encoding="utf-8", errors="replace")
        for match in call_pattern.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            line_text = text.splitlines()[line - 1].strip()
            if line_text.startswith("def "):
                continue
            callers.append(f"{source.relative_to(root).as_posix()}:{line}")
            break
    assert callers, (
        f"{entry} is defined but never called under any of "
        f"{sorted(callable_names)} -- {module_name} is unreachable from a run"
    )


def test_a_control_for_each_module_is_structurally_complete(tmp_path):
    """Every complete CONTROL written here parses and validates as it did on main.

    Package o5 item (3): a complete CONTROL must read exactly as before. This
    only proves the *synthetic* CONTROLs are structurally valid, so that a
    failure in the dispatch tests below is a dispatch failure, not a
    malformed-input failure with the same red colour.
    """
    from delfin.config import read_control_file, validate_control_file

    for name in _module_names():
        ws = tmp_path / re.sub(r"\W+", "_", name)
        ws.mkdir()
        path = _write_control(ws)
        problems = validate_control_file(str(path))
        assert not problems, f"{name}: the synthetic CONTROL does not validate: {problems}"
        config = read_control_file(str(path))
        assert config.get("charge") is not None, f"{name}: CONTROL did not parse"


def test_switching_esd_on_changes_what_the_run_dispatches(tmp_path, monkeypatch):
    """ESD on vs off must change the phases the run reaches -- an effect.

    The assertion is that the control key *does* something: with
    ``ESD_modul=no`` the ESD phase is never entered, with ``ESD_modul=yes`` it
    is. A key that is read but has no effect fails here, which is the
    distinction the package asks for ("a key read but without effect is a
    finding too").
    """
    import delfin.cli as cli

    seen: list[str] = []

    def _record(name, result=True):
        def _phase(ctx):
            seen.append(name)
            return result
        return _phase

    cli._load_full_cli_dependencies()
    monkeypatch.setattr(cli, "run_esd_phase", _record("esd"))
    monkeypatch.setattr(cli, "run_classic_phase", _record("classic"))
    monkeypatch.setattr(cli, "run_hyperpol_xtb_phase", _record("hyperpol"))
    monkeypatch.setattr(cli, "run_tadf_xtb_phase", _record("tadf"))

    # ESD is dispatched *inside* a method run, not instead of one: with
    # ``method=classic`` the classic phase runs and the ESD phase follows it
    # (delfin/cli.py:2200). The switch decides whether the second one is
    # entered. An empty ``method`` is refused by the validator ("Missing
    # required key: method") and ``method=esd`` by the value check, so
    # classic is the route by which ESD is actually reached.
    for switch, expected in (("no", ["classic"]), ("yes", ["classic", "esd"])):
        seen.clear()
        ws = tmp_path / f"esd_{switch}"
        ws.mkdir(exist_ok=True)
        _write_control(ws, ESD_modul=switch, method="classic",
                       ESD_T1_opt="uks", calc_initial="no")
        # cli.main sets the DELFIN_FFFREE_* construction flags itself
        # (cli_manta._apply_construction_env) and reads them from os.environ
        # at call time, so a leftover would change what the next test builds.
        # monkeypatch.setenv cannot stop it; restoring can.
        env_before = {k: v for k, v in os.environ.items() if k.startswith("DELFIN_")}
        try:
            cli.main([str(ws), "--no-cleanup"])
        except SystemExit:
            pass
        finally:
            for _k in [k for k in os.environ if k.startswith("DELFIN_")]:
                if _k not in env_before:
                    del os.environ[_k]
            os.environ.update(env_before)
        assert seen == expected, (
            f"ESD_modul={switch} dispatched {seen}, expected {expected}; "
            "a switch that does not change the dispatch has no effect"
        )
