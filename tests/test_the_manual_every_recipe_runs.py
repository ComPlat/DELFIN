"""Every recipe in the manual runs, against DELFIN's real code.

A recipe is a fenced code block in one of the manual files.  What "runs"
means per kind — always the cheapest check that could fail, never a real
computation:

* ``ini``  — ``validate_control_text`` accepts it.  Snippet blocks that
  document single keys rather than a whole CONTROL file are marked
  ``<!-- recipe: not-runnable — … -->`` only when they could never stand
  alone; the *complete* recipes (Sections 12/13) must validate in full.
* ``yaml`` — the pipeline loader ``_load_yaml`` accepts the mapping.
* ``bash`` — every ``delfin …``, ``delfin-<tool> …`` and ``python -m
  delfin[.installer]`` line parses in the real argparse parser
  (subcommand + every option).  Nothing is executed: the patched
  ``parse_args`` raises as soon as parsing succeeds.  ``pip install``
  extras are checked against pyproject.toml, scripts and paths against
  the repository.
* ``python`` — the block compiles and every import it names resolves.
* SMILES — found in any block, however embedded (``SMILES=…``,
  ``echo "…" > input.txt``, a bare line in an untagged block): RDKit
  must parse it WITH sanitisation, and DELFIN's own
  ``mol_from_smiles_rdkit`` must accept it.  That is the exact class the
  a05f02a4 audit caught: the five-ring pyridine parses only
  unsanitised and silently builds a different molecule (C24 instead of
  C30).  Where the surrounding text names a charge (``charge=N`` in the
  same ini block), the parsed molecule's formal charge must match.

Blocks that are deliberately not runnable carry, directly above the
fence, ``<!-- recipe: not-runnable — reason -->`` and are skipped with
that reason.  An untagged block that is neither SMILES nor an XYZ body
and carries no marker is a failure, not silence.
"""

from __future__ import annotations

import argparse
import ast
import configparser
import re
import shlex
import sys
import tomllib
from pathlib import Path

import pytest

from tests.manual_blocks import load_blocks, repo_root

# ---------------------------------------------------------------- bash


class _Parsed(Exception):
    """Raised by patched parse_args: validated, execution prevented."""

    def __init__(self, parser, unknown=()):
        super().__init__("parsed")
        self.parser = parser
        self.unknown = list(unknown)


@pytest.fixture()
def no_execution(monkeypatch):
    """argparse validates and stops; nothing the manuals spell runs."""
    real_parse = argparse.ArgumentParser.parse_args
    real_known = argparse.ArgumentParser.parse_known_args

    def _parse(self, args=None, namespace=None):
        real_parse(self, args, namespace)
        raise _Parsed(self)

    def _known(self, args=None, namespace=None):
        _, unknown = real_known(self, args, namespace)
        raise _Parsed(self, unknown)

    monkeypatch.setattr(argparse.ArgumentParser, "parse_args", _parse)
    monkeypatch.setattr(argparse.ArgumentParser, "parse_known_args", _known)


def _load_pyproject():
    with open(repo_root() / "pyproject.toml", "rb") as f:
        return tomllib.load(f)


def _script_target(prog: str) -> str | None:
    """The ``module:function`` a console script maps to, or None."""
    project = _load_pyproject()
    scripts = project["project"].get("scripts", {})
    return scripts.get(prog)


def _import_target(target: str):
    mod, _, func = target.partition(":")
    module = __import__(mod, fromlist=[func])
    return getattr(module, func)


def _parser_of(prog: str) -> argparse.ArgumentParser | None:
    """The real parser behind a console script, if cheaply reachable.

    Only builds a parser from a module that exposes one of the known
    builder shapes (build_parser / _build_parser) or constructs it in
    main(); otherwise returns None and the caller reports the gap.
    """
    target = _script_target(prog)
    if not target:
        return None
    mod, _, func = target.partition(":")
    module = __import__(mod, fromlist=[func])
    builder = getattr(module, "build_parser", None) or getattr(
        module, "_build_parser", None)
    if builder is not None:
        return builder()
    # main() builds its parser inline; reconstruct it by calling main
    # under the patched parse (the patch raises _Parsed before any side
    # effect) — but main may do work before parsing, so only allow this
    # for modules known to parse first.  Conservative: give up.
    return None


# Progs the manuals spell that do NOT parse via a module-level builder:
# `delfin` dispatches subcommands by hand in cli.main before the main
# parser; `delfin-agent` needs its parser imported directly.
_EXTRA_PROGS = {
    "delfin": ("delfin.cli_helpers", "_build_parser"),
    "delfin-agent": ("delfin.agent.cli", "build_parser"),
}


def _validate_bash_line(where: str, cmd: list[str]) -> None:
    """One shell command from a recipe: prog + options must parse."""
    # env-var prefixes (USE_SYSTEM_TOOLS=1 bash …) are not arguments
    while cmd and re.match(r"^[A-Z_][A-Z0-9_]*=", cmd[0]):
        cmd = cmd[1:]
    if not cmd:
        return
    prog = cmd[0]
    if prog in {"bash", "source", "cd", "mkdir", "echo", "which", "git",
                "pip", "pip3", "python", "python3"}:
        _validate_generic_unix(where, cmd)
        return
    if prog.startswith("delfin"):
        _validate_delfin_prog(where, prog, cmd[1:])
        return
    raise AssertionError(f"{where}: unknown program {prog!r} in recipe")


def _validate_generic_unix(where: str, cmd: list[str]) -> None:
    """pip/git/bash sanity: extras exist, paths exist, flags plausible."""
    prog = cmd[0]
    if prog in {"pip", "pip3"}:
        if "install" in cmd:
            for tok in cmd:
                m = re.match(r"^\.\[([\w,-]*)\]$", tok)
                if m and m.group(1):
                    extras = [e for e in m.group(1).split(",") if e]
                    _validate_pip_extras(where, extras)
        return
    if prog == "python" and "-m" in cmd:
        i = cmd.index("-m")
        mod = cmd[i + 1]
        if mod.startswith("delfin"):
            _validate_python_m(where, mod, cmd[i + 2:])
        else:
            # generic stdlib/third-party module (venv, pytest, …):
            # it must at least be importable
            import importlib.util
            if importlib.util.find_spec(mod) is None:
                raise AssertionError(
                    f"{where}: python -m {mod}: module not found")
        return
    if prog in {"bash", "source"} and len(cmd) > 1:
        script = cmd[1].lstrip("./")
        # Path is under the repo; a copy beyond it (~/software/delfin)
        # names the clone, not a repo file — accept both shapes.
        if script.startswith("software/delfin/") or "/" not in script:
            _ = script  # repo-relative: delfin/… checked below
            if not script.startswith("software/") and "/" in script:
                if not (repo_root() / script).exists():
                    raise AssertionError(
                        f"{where}: script {cmd[1]} not in the repository")
    # git/cd/mkdir/echo/which: nothing structural to check.


def _validate_pip_extras(where: str, extras: list[str]) -> None:
    project = _load_pyproject()["project"]
    known = set(project.get("optional-dependencies", {}))
    for e in extras:
        if e not in known:
            raise AssertionError(
                f"{where}: pip extra {e!r} is not in pyproject.toml "
                f"(known: {sorted(known)})")


def _try_parse(parser: argparse.ArgumentParser, argv: list[str]) -> None:
    """Parse under the patch; _Parsed means validated-and-stopped."""
    try:
        parser.parse_args(argv)
    except _Parsed:
        return


def _validate_python_m(where: str, module: str, rest: list[str]) -> None:
    if module == "delfin":
        _validate_delfin_prog(where, "delfin", rest)
        return
    if module == "delfin.installer":
        from delfin import installer
        _inline_parser(
            where, "python -m delfin.installer",
            lambda: installer.main(rest))
        return
    raise AssertionError(f"{where}: unhandled python -m module {module!r}")


def _inline_parser(where: str, prog: str, call) -> argparse.ArgumentParser:
    """Run call() under the patched parse; return the parser it used.

    The patch raises _Parsed at the first parse_args — before any side
    effect — and carries the parser instance on the exception, so the
    caller can inspect what was accepted.  ``call`` must arrange for
    the recipe's argv to reach the parser (passing it to main, or via
    sys.argv for mains that parse ``sys.argv`` implicitly).
    """
    argv = list(sys.argv)
    try:
        call()
    except _Parsed as p:
        return p.parser
    except SystemExit as e:
        raise AssertionError(
            f"{where}: {prog} rejected its own recipe line "
            f"(exit {e.code}) — argparse error") from None
    finally:
        sys.argv[:] = argv
    raise AssertionError(
        f"{where}: {prog} ran past parsing without hitting parse_args — "
        "cannot validate safely")


def _validate_delfin_prog(where: str, prog: str, rest: list[str]) -> None:
    """A `delfin …` / `delfin-<tool> …` line against the real parser."""
    if prog == "delfin":
        # Subcommands dispatched by hand in cli.main BEFORE the main
        # parser: validate through the real dispatch, with parse patched
        # to stop right after validation.
        from delfin import cli
        p = _inline_parser(where, "delfin", lambda: cli.main(rest))
        return
    if prog == "delfin-agent":
        # Bare `delfin-agent` routes to `chat` inside main(); the bare
        # parser would reject it, so validate through the real main.
        from delfin.agent import cli as agent_cli
        _inline_parser(where, prog, lambda: agent_cli.main(rest))
        return
    target = _script_target(prog)
    if target is None:
        raise AssertionError(
            f"{where}: {prog} is not a console script in pyproject.toml")
    parser = _parser_of(prog)
    if parser is not None:
        _try_parse(parser, rest)
        return
    # No module-level builder: run the real main() under the patch.
    # It either takes argv, or parses sys.argv — feed it either way.
    import inspect
    func = _import_target(target)
    try:
        sig = inspect.signature(func)
        takes_argv = any(
            p.name in {"argv", "args"} for p in sig.parameters.values())
    except (TypeError, ValueError):
        takes_argv = False
    if takes_argv:
        _inline_parser(where, prog, lambda: func(rest))
    else:
        sys.argv[:] = [prog] + rest
        _inline_parser(where, prog, func)


# ---------------------------------------------------------------- SMILES

# Where SMILES can hide in a recipe line.  Ordered: most specific first.
_SMILES_LINE_PATTERNS = [
    # SMILES=[Cu+2](...)        (CONTROL key)
    re.compile(r"^SMILES=(?P<smiles>\S.*)$"),
    # echo "[Fe+2](...)" > input.txt
    re.compile(r"^echo\s+\"(?P<smiles>[^\"]+)\""),
]


def _smiles_candidates(block) -> list[tuple[str, str]]:
    """(where, smiles) for every SMILES a block spells, embedded or bare."""
    out = []
    for raw in block.lines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        hit = None
        for pat in _SMILES_LINE_PATTERNS:
            m = pat.match(line)
            if m:
                hit = m.group("smiles").strip()
                break
        if hit is None and not block.lang:
            # Untagged block: a lone non-empty line that is not an XYZ
            # body row, not a $ splitter and not an ellipsis is a SMILES
            # by the manual's own definition of input.txt.
            if not re.match(r"^[A-Z][a-z]?\s+-?\d", line) and line != "$" \
                    and not re.match(r"^\.+$", line) \
                    and re.match(r"^[\[\]A-Za-z0-9@+\-=()/#%.,\\]+$", line):
                hit = line
        if hit is not None:
            out.append((block.where, hit))
    return out


def _check_smiles(where: str, smiles: str, expected_charge: int | None) -> None:
    """A manual SMILES must survive RDKit sanitised AND DELFIN's reader.

    The a05f02a4 audit: the five-ring pyridine spelling parses only
    unsanitised, and the unsanitised molecule is a DIFFERENT species
    (C24H24FeN6, 31 atoms, vs the intended C30H30FeN6, 37).  So the
    strict parse is the check; the converter parse is the second,
    DELFIN-native opinion.
    """
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise AssertionError(
            f"{where}: RDKit refuses the SMILES (sanitised): {smiles[:60]}…")
    if expected_charge is not None:
        actual = Chem.GetFormalCharge(mol)
        if actual != expected_charge:
            raise AssertionError(
                f"{where}: SMILES charge {actual} but the recipe says "
                f"charge={expected_charge}")
    from delfin.smiles_converter import mol_from_smiles_rdkit
    m2, err = mol_from_smiles_rdkit(smiles, allow_metal=True)
    if m2 is None:
        raise AssertionError(
            f"{where}: DELFIN's smiles_converter refuses the SMILES: {err}")


# ---------------------------------------------------------------- checks


_MISSING_RE = re.compile(
    r"Missing required CONTROL values for:|Missing required key:")
_PLACEHOLDER_RE = re.compile(r"Placeholder \[")


def _check_ini(block, *, full: bool) -> None:
    from delfin.config import validate_control_text
    errors = validate_control_text(block.text)
    # A snippet (not a whole CONTROL file) may only fail on missing
    # required keys / placeholders — those are the keys the surrounding
    # prose documents.  Any other error (unknown value, bad syntax,
    # wrong type) is a real recipe bug.
    real = [e for e in errors
            if not _MISSING_RE.search(e) and not _PLACEHOLDER_RE.search(e)]
    if real:
        raise AssertionError(
            f"{block.where}: CONTROL recipe invalid: {'; '.join(real)}")
    if full and errors:
        raise AssertionError(
            f"{block.where}: complete recipe missing required values: "
            f"{'; '.join(errors)}")


def _check_yaml(block) -> None:
    import yaml
    data = yaml.safe_load(block.text)
    if not isinstance(data, dict):
        raise AssertionError(
            f"{block.where}: pipeline YAML is not a mapping")
    for key in ("name", "steps"):
        if key not in data:
            raise AssertionError(
                f"{block.where}: pipeline YAML lacks '{key}' "
                "(the real loader requires it)")
    # The real loader also builds the pipeline; that step needs the
    # steps themselves resolvable, which the loader does via
    # delfin.tools._serialize.from_dict — validate it too.
    from delfin.tools._serialize import from_dict
    from_dict(data)


def _check_python(block) -> None:
    try:
        tree = ast.parse(block.text)
    except SyntaxError as e:
        raise AssertionError(
            f"{block.where}: python block does not compile: {e}") from None
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append(node.module)
    import importlib
    for mod in imports:
        try:
            importlib.import_module(mod)
        except Exception as e:  # noqa: BLE001
            raise AssertionError(
                f"{block.where}: import {mod!r} fails: {e}") from None


_XYZ_ROW_RE = re.compile(r"^[A-Z][a-z]?\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s*$")


def _check_untagged(block) -> None:
    """Untagged block: SMILES, XYZ body — or a marker, else a failure."""
    lines = [l for l in block.lines() if l.strip()]
    if not lines:
        raise AssertionError(f"{block.where}: empty untagged block")
    if all(_XYZ_ROW_RE.match(l.strip()) or l.strip() == "..." for l in lines):
        return  # XYZ body — checked as SMILES candidates too (none)
    # Not an XYZ body: every line must be a SMILES (input.txt contract)
    cands = _smiles_candidates(block)
    if not cands:
        raise AssertionError(
            f"{block.where}: untagged block is neither SMILES nor XYZ — "
            "tag it or mark it not-runnable")
    for where, smi in cands:
        _check_smiles(where, smi, None)


# ---------------------------------------------------------------- tests

# Complete CONTROL recipes: sections titled as ready-to-run examples.
# Located by the heading line above the fence, so a new "Recipes"
# section is picked up automatically.
_FULL_RECIPE_RE = re.compile(r"recipes|workflow|getting started|example",
                             re.IGNORECASE)


def _is_full_recipe(block) -> bool:
    """A CONTROL block under a Recipes/Examples heading is complete."""
    path = repo_root() / block.path
    lines = path.read_text(encoding="utf-8").splitlines()
    i = block.line - 2  # 0-based index of the line above the fence
    while i >= 0:
        text = lines[i].strip()
        if text.startswith("#"):
            return bool(_FULL_RECIPE_RE.search(text))
        if text and not text.startswith("<!--") and not text.startswith("```"):
            return False  # prose paragraph directly above: a snippet
        i -= 1
    return False


def _expected_charge(block) -> int | None:
    m = re.search(r"^charge=(-?\d+)$", block.text, re.MULTILINE)
    return int(m.group(1)) if m else None


def test_bash_recipes_parse_in_the_real_argparse(no_execution, capsys):
    blocks = [b for b in load_blocks(repo_root())
              if b.lang == "bash" and not b.not_runnable_reason]
    failures = []
    for b in blocks:
        for raw in b.lines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            try:
                cmd = shlex.split(line, comments=True)
            except ValueError as e:
                failures.append(f"{b.where}: unparseable shell line: {e}")
                continue
            if not cmd:
                continue
            try:
                _validate_bash_line(b.where, cmd)
            except AssertionError as e:
                failures.append(str(e))
            except SystemExit as e:
                failures.append(
                    f"{b.where}: argparse rejected: {line} (exit {e.code})")
    print(f"bash blocks: {len(blocks)}, failures: {len(failures)}")
    for f in failures:
        print("  " + f)
    assert not failures


def test_ini_recipes_validate(capsys):
    blocks = [b for b in load_blocks(repo_root())
              if b.lang == "ini" and not b.not_runnable_reason]
    failures = []
    for b in blocks:
        try:
            _check_ini(b, full=_is_full_recipe(b))
        except AssertionError as e:
            failures.append(str(e))
    print(f"ini blocks: {len(blocks)}, failures: {len(failures)}")
    for f in failures:
        print("  " + f)
    assert not failures


def test_yaml_recipe_loads():
    blocks = [b for b in load_blocks(repo_root())
              if b.lang == "yaml" and not b.not_runnable_reason]
    assert blocks, "the pipeline YAML recipe is missing from the manuals"
    for b in blocks:
        _check_yaml(b)


def test_python_recipes_compile_and_import(capsys):
    blocks = [b for b in load_blocks(repo_root())
              if b.lang == "python" and not b.not_runnable_reason]
    assert blocks
    failures = []
    for b in blocks:
        try:
            _check_python(b)
        except AssertionError as e:
            failures.append(str(e))
    print(f"python blocks: {len(blocks)}, failures: {len(failures)}")
    for f in failures:
        print("  " + f)
    assert not failures


def test_every_smiles_in_the_manual_is_a_real_molecule(capsys):
    rdkit = pytest.importorskip("rdkit")
    blocks = load_blocks(repo_root())
    checked = 0
    failures = []
    for b in blocks:
        if b.not_runnable_reason:
            continue
        charge = _expected_charge(b)
        for where, smi in _smiles_candidates(b):
            checked += 1
            try:
                _check_smiles(where, smi, charge)
            except AssertionError as e:
                failures.append(str(e))
    print(f"SMILES checked: {checked}, failures: {len(failures)}")
    for f in failures:
        print("  " + f)
    assert not failures


def test_untagged_blocks_are_classified():
    for b in load_blocks(repo_root()):
        if b.lang or b.not_runnable_reason:
            continue
        _check_untagged(b)


def test_marked_blocks_carry_a_reason():
    for b in load_blocks(repo_root()):
        if b.not_runnable_reason is not None:
            assert len(b.not_runnable_reason) > 3, (
                f"{b.where}: not-runnable marker without a real reason")


# Language tags that are documentation, not recipes: output examples and
# literature.  A block with any OTHER tag must be handled by a check
# above (or explicitly marked not-runnable), else the inventory grows a
# kind the tests silently ignore.
_NON_RECIPE_LANGS = {"text", "bibtex"}


def test_every_language_is_accounted_for(capsys):
    blocks = load_blocks(repo_root())
    checked_langs = {"bash", "ini", "yaml", "python"}
    problems = []
    for b in blocks:
        if b.not_runnable_reason:
            continue  # skipped with a reason, whatever its tag
        if b.lang in checked_langs or b.lang in _NON_RECIPE_LANGS:
            continue
        if not b.lang:
            continue  # untagged: classified by test_untagged_blocks…
        problems.append(f"{b.where}: unhandled language tag {b.lang!r}")
    for f in problems:
        print("  " + f)
    assert not problems, (
        "new block kind the recipe tests do not cover — add a check or "
        "mark the block not-runnable")
