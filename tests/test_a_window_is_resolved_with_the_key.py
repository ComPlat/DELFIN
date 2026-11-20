"""Every capability resolve in the product sends the provider key.

``model_capabilities.resolve`` sends ``api_key`` as a Bearer token on the
``/v1/models`` probe. Without it the probe is refused and the resolver falls
through to a heuristic window, so the call returns a plausible object with a
wrong window and a ``supports_tools`` that was guessed rather than read.

Measured on 2026-11-20, same key, same endpoints, with and without the key:

    kit.glm-5.3-flash     no key    131,072 static     with key    524,288 live
    alias-code            no key     32,768 heuristic  with key    262,144 live
    alias-huge            no key     32,768 heuristic  with key  1,048,576 live

A window read a factor of 32 too small is not a visible failure: the agent
keeps answering and only the compaction count shows it. The two call sites
this file was written for decided more than a number -- one fed a model's
fitness warning, the other decided whether a cheaper model may carry a
sub-agent at all, and a heuristic ``supports_tools`` answers "no tools" for
a model that has them.

Input: every ``.py`` file under ``delfin/``. Output: the resolve calls that
pass no ``api_key``. Semantics: a tree scan rather than a test of the two
known sites, because the defect is reintroduced by the next call site and an
importer grep does not see one that has not been written yet. Aliases are
resolved per file from the import statements, so renaming the import does
not hide a call.
"""
from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRODUCT = ROOT / "delfin"

#: The attribute of ``model_capabilities`` under test.
_FUNC = "resolve"
#: Module names that carry it.
_MODULE_TAILS = ("model_capabilities",)


def _aliases(tree: ast.Module) -> tuple[set[str], set[str]]:
    """Names bound to the resolver in one module.

    Returns (direct, module) -- direct names called as ``name(...)`` and
    module names called as ``name.resolve(...)``.
    """
    direct: set[str] = set()
    module: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            # Two spellings reach the same function and both must be seen.
            # ``from .model_capabilities import resolve`` names the module in
            # node.module, while ``from . import model_capabilities as mc``
            # leaves node.module empty and names it among node.names -- a
            # scan that reads only node.module misses the second, which is
            # how the second of the two live call sites first went unseen.
            mod = node.module or ""
            from_module = any(mod.endswith(t) for t in _MODULE_TAILS)
            for a in node.names:
                if a.name in _MODULE_TAILS:
                    module.add(a.asname or a.name)
                elif from_module and a.name == _FUNC:
                    direct.add(a.asname or a.name)
        elif isinstance(node, ast.Import):
            for a in node.names:
                if any(a.name.endswith(t) for t in _MODULE_TAILS):
                    module.add(a.asname or a.name.split(".")[-1])
    return direct, module


def _is_resolver_call(call: ast.Call, direct: set[str],
                      module: set[str]) -> bool:
    f = call.func
    if isinstance(f, ast.Name):
        return f.id in direct
    if isinstance(f, ast.Attribute) and f.attr == _FUNC:
        return isinstance(f.value, ast.Name) and f.value.id in module
    return False


def _passes_key(call: ast.Call) -> bool:
    """True when the call needs no key, or sends one.

    A key is needed only for the ``/v1/models`` probe, so a call that
    switches the probe off with ``allow_live=False`` is exempt: it reads the
    static table and reaches no endpoint. Two dashboard call sites do exactly
    that to read ``is_reasoning``, which is not a live property. The
    exemption is the literal ``False``; ``allow_live=True`` without a key is
    the defect this file is about.

    Otherwise the call names ``api_key``, forwards ``**kwargs``, or passes a
    fourth positional argument, the signature being
    ``resolve(provider, model, base_url, api_key="")``.
    """
    for k in call.keywords:
        if (k.arg == "allow_live" and isinstance(k.value, ast.Constant)
                and k.value.value is False):
            return True
    if any(k.arg == "api_key" for k in call.keywords):
        return True
    if any(k.arg is None for k in call.keywords):      # **kwargs
        return True
    return len(call.args) >= 4


def _offenders() -> list[str]:
    out: list[str] = []
    for path in sorted(PRODUCT.rglob("*.py")):
        try:
            tree = ast.parse(path.read_text(errors="replace"))
        except SyntaxError:
            continue
        direct, module = _aliases(tree)
        if not direct and not module:
            continue
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call)
                    and _is_resolver_call(node, direct, module)
                    and not _passes_key(node)):
                rel = path.relative_to(ROOT)
                out.append(f"{rel}:{node.lineno}")
    return out


def test_the_scan_finds_the_resolver_at_all():
    """A scan that matches nothing would make the assertion below pass
    while checking no call site."""
    found = 0
    for path in PRODUCT.rglob("*.py"):
        try:
            tree = ast.parse(path.read_text(errors="replace"))
        except SyntaxError:
            continue
        direct, module = _aliases(tree)
        if not direct and not module:
            continue
        found += sum(
            1 for n in ast.walk(tree)
            if isinstance(n, ast.Call) and _is_resolver_call(n, direct, module))
    assert found >= 3, (
        f"the scan sees {found} resolver call(s); it is not looking at the "
        "product any more")


def test_every_resolve_sends_the_provider_key():
    """Without the key the probe is refused and the window is guessed."""
    bad = _offenders()
    assert bad == [], (
        "capability resolve without api_key at: " + ", ".join(bad)
        + " -- the /v1/models probe is then refused and the window comes "
          "from the heuristic, which was 32,768 against a measured 262,144 "
          "to 1,048,576")
