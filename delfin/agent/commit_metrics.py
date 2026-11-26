"""What one commit is worth, per agent, as one JSON row.

The question this answers is "did the agents get better", and it needs three
kinds of number that existed separately and were never written down together
for the same commit:

* **Reachability** -- how many capabilities in the product nothing reaches,
  against the ceiling the ratchet allows. A commit that orphans a mechanism
  costs something even when every test passes.
* **What each agent is OFFERED** -- the tool surface and the system-prompt
  size of the generalist and of EVERY specialist separately. One count for
  "the agent" answers a question about no agent that runs.
* **Quality with repeats** -- read from the benchmark rounds
  ``bench_history`` already keeps (``delfin-agent bench nightly`` writes them),
  overall and per capability tag, with the min/max of the per-round pass
  rates as the spread. Nothing here stores a benchmark result a second time.

The fourth part is the comparison: a quality difference is ``better`` or
``worse`` only when ``bench_ab.newcombe_interval`` on the pooled pass counts
excludes zero, and rounds that differ on a ``bench_history.comparable`` axis
are not compared at all.

Why the reachability scan is here rather than imported
------------------------------------------------------
``tests/test_every_capability_has_a_caller.py`` is the instrument for the
first number, and this module does not import it:

1. ``tests/`` is not an importable package, so there is no import to write;
2. it is a pytest module with a module-scoped fixture, and importing it
   would drag pytest collection into every command that records a row;
3. worse, it would make the PRODUCT depend on the SUITE -- a metrics command
   that fails on an installation without pytest is not a measurement.

So the scan is a faithful copy here, and
``tests/test_the_metrics_scan_agrees_with_the_ratchet.py`` pins it against
the original: the definitions, the unmounted SET (not only its size) and each
name's tier must agree with what the ratchet computes, so the two cannot
drift apart while both suites stay green. Divergence is a defect of this
module, not of the ratchet. The same pattern and the same three reasons are
already recorded in ``scripts/ab_task_pools.py``.

What the OFFERED count includes, and why that is not
``role_tool_surface_report``
----------------------------------------------------
``api_client.role_tool_surface_report`` reads the built-in catalogue only.
Production also appends the tools the MCP servers expose -- ``delfin-tools``
(26 tools) and ``delfin-ops`` (102) -- and those are where an editor preset's
tools actually live. Measured on 41e15fe9f: ``role_tool_surface_report``
reports ``xyz-editor`` as 5 tools while the session offers it 23, and
``method-reviewer`` as 17 against 24. A row built on the catalogue alone
would have said the structure specialists are offered almost nothing.

The MCP half is derived statically, from the modules that register the
tools, so record it needs no server started and no key present.
"""

from __future__ import annotations

import ast
import io
import json
import re
import tokenize
from pathlib import Path
from typing import Any, Iterable, Optional

#: The ratchet's ceiling, and the number the fixture was measured at.
#: A row records the MEASURED count; the ceiling is here so a row can say
#: whether the commit is inside it without reading the test file.
_REACHABILITY_CEILING = 86

#: A root that yields fewer definitions than this was not measured. Set two
#: orders of magnitude below the 5368 this tree carried on 2026-11-25, so only
#: a broken or foreign root trips it, never a tree that lost code.
_MIN_DEFINITIONS = 300

#: The trees the ratchet censuses (``_SCAN_DIRS`` in
#: ``tests/test_every_capability_has_a_caller.py``), first-wins in this order.
_SCAN_DIRS = ("delfin/agent", "delfin/tools", "delfin/agent/pack")

_BINARY_SUFFIXES = frozenset({
    ".png", ".jpg", ".jpeg", ".gif", ".pdf", ".ico", ".woff", ".woff2",
    ".zip", ".gz", ".so", ".pyc", ".xlsx", ".docx", ".pptx", ".odt", ".ods",
})

_WORD = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_DOTTED = re.compile(r"[A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)+")

# ---------------------------------------------------------------------------
# Reachability
# ---------------------------------------------------------------------------

def definitions(tree_root: Path, definition_root: Optional[Path] = None
                ) -> dict[str, tuple[str, int]]:
    """``{name: (file, line)}`` for every module-level callable.

    Module level only, as the ratchet defines it: a method is reached through
    its class, and asking the same question of it would need the class's own
    reachability first.

    *definition_root* defaults to the ratchet's own ``_SCAN_DIRS``, the scope
    its ceiling is measured on. It is a
    parameter so the same question can be asked of another root without
    copying the function again.
    """
    roots = ([Path(definition_root)] if definition_root
             else [Path(tree_root) / d for d in _SCAN_DIRS])
    out: dict[str, tuple[str, int]] = {}
    for root in roots:
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*.py")):
            try:
                tree = ast.parse(path.read_text(errors="replace"))
            except SyntaxError:
                continue
            for node in tree.body:
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) \
                        and not node.name.startswith("__"):
                    out.setdefault(node.name,
                                   (str(path.relative_to(tree_root)), node.lineno))
    return out


def mentions(tree_root: Path, defs: dict, *, skip_own_definition: bool,
             skip_paths: Iterable[Path] = (),
             repo_root: Optional[Path] = None) -> dict[str, int]:
    """How often each name is mentioned under *tree_root*.

    Python files are tokenised, so a name in a comment does not count while a
    name in a string literal does -- that is name-based dispatch, which is a
    real reference. Everything else readable as text is scanned for words,
    which is how a module named from YAML is found.

    *repo_root* is what ``defs``' paths are relative to, and it defaults to
    *tree_root*. It has to be passed when the scan covers a SUBTREE -- the
    suite half walks ``<repo>/tests`` while the definitions are recorded as
    ``delfin/agent/<file>.py``. Deriving the relative path from *tree_root*
    in that case raised for every file, the tokeniser's ``except`` swallowed
    it, and the suite half came out empty: measured as ``unmounted: 0`` on a
    tree where the ratchet counts 86.

    *skip_paths* must contain the unmounted baseline: it lists every name it
    holds, so scanning it made all of them read as "reached by the suite" and
    the instrument measured itself.
    """
    skip = {Path(p) for p in skip_paths}
    hits = {name: 0 for name in defs}
    root = Path(tree_root)
    anchor = Path(repo_root) if repo_root else root
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix.lower() in _BINARY_SUFFIXES:
            continue
        if path in skip:
            continue
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        # Relative to the REPOSITORY (see *repo_root* above), which is what
        # ``defs`` holds.
        try:
            rel = path.relative_to(anchor).as_posix()
        except ValueError:
            rel = path.as_posix()
        if path.suffix == ".py":
            try:
                for tok in tokenize.generate_tokens(io.StringIO(text).readline):
                    if tok.type == tokenize.NAME and tok.string in hits:
                        name = tok.string
                        if (skip_own_definition and rel == defs[name][0]
                                and tok.start[0] == defs[name][1]):
                            continue
                        hits[name] += 1
                    elif tok.type == tokenize.STRING:
                        for word in _WORD.findall(tok.string):
                            if word in hits:
                                hits[word] += 1
                continue
            except (tokenize.TokenError, IndentationError, SyntaxError):
                pass        # unparseable: fall through to the word scan
        for word in _WORD.findall(text):
            if word in hits:
                hits[word] += 1
    count_dotted_dispatch(tree_root, hits, skip_paths=skip)
    return hits


def count_dotted_dispatch(tree_root: Path, hits: dict, *,
                          skip_paths: Iterable[Path] = ()) -> None:
    """Also count a name spelled with dots instead of underscores.

    ``supervisor`` dispatches with ``getattr(swarm_rpc, method.replace(".", "_"))``,
    so a protocol declares ``"swarm.observe.feed"`` and the underscored name
    appears nowhere. Two of the names the scan called unreachable are reached
    that way.
    """
    skip = {Path(p) for p in skip_paths}
    dotted = {n.replace("_", "."): n for n in hits if "_" in n}
    if not dotted:
        return
    for path in sorted(Path(tree_root).rglob("*")):
        if not path.is_file() or path.suffix.lower() in _BINARY_SUFFIXES:
            continue
        if path in skip:
            continue
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        for run in set(_DOTTED.findall(text)):
            name = dotted.get(run)
            if name:
                hits[name] += 1


def read_baseline(path: Path) -> dict[str, str]:
    """``{name: tier}`` from the unmounted-capability fixture."""
    rows: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            rows[parts[1]] = parts[0]
    return rows


def reachability(tree_root: Path) -> dict[str, Any]:
    """The unmounted capabilities, their tiers, and the ratchet's ceiling.

    Reported as counts over a STATED population: module-level callables under
    the ratchet's ``_SCAN_DIRS``. A count from another root answers another question and
    is not quoted against the 86 ceiling.
    """
    tree_root = Path(tree_root)
    defs = definitions(tree_root)
    # A scan that read nothing must not report "nothing unreachable". Measured
    # 2026-11-25: an unreadable root gave `unmounted: 0, within_ceiling: true`
    # -- a row that reads as a clean bill of health for a tree it never
    # opened, which is worse than no row. The ratchet guards its own scan with
    # the same idea (`assert len(defs) > 3000`); a row cannot assert, so it
    # says what happened instead and leaves the numbers absent rather than
    # zero. The floor is deliberately far below any tree that has agents: on
    # 2026-11-25 this tree defined 5368.
    if len(defs) < _MIN_DEFINITIONS:
        return {
            "population": "module-level callables under " + ", ".join(_SCAN_DIRS[:2]),
            "callables_defined": len(defs),
            "unmounted": None,
            "unmounted_names": [],
            "tiers": None,
            "ceiling": _REACHABILITY_CEILING,
            "within_ceiling": None,
            "newly_unreachable": [],
            "stale_in_baseline": [],
            "error": (f"the scan found {len(defs)} definitions under "
                      f"{tree_root}/delfin/agent (< {_MIN_DEFINITIONS}), so "
                      "this root was not measured -- not a count of zero"),
        }
    baseline_path = tree_root / "tests" / "fixtures" / "unmounted_capabilities.txt"
    baseline = read_baseline(baseline_path) if baseline_path.is_file() else {}
    # The product half scans ``delfin/`` ONLY, exactly as the ratchet does.
    # Scanning the whole tree here counted the suite's own mentions as
    # product reachability, so every test-only capability looked connected
    # and the count came out 0 against the ratchet's 86.
    in_product = mentions(tree_root / "delfin", defs, skip_own_definition=True,
                          skip_paths=[baseline_path], repo_root=tree_root)
    unmounted = sorted(n for n in defs if in_product[n] == 0)
    in_suite = mentions(tree_root / "tests", defs, skip_own_definition=False,
                        skip_paths=[baseline_path], repo_root=tree_root)
    tiers: dict[str, int] = {"nowhere": 0, "test-only": 0}
    for name in unmounted:
        tiers["test-only" if in_suite.get(name) else "nowhere"] += 1
    newly = sorted(n for n in unmounted if n not in baseline)
    stale = sorted(n for n in baseline if n not in unmounted)
    return {
        "population": "module-level callables under " + ", ".join(_SCAN_DIRS[:2]),
        "callables_defined": len(defs),
        "unmounted": len(unmounted),
        "unmounted_names": unmounted,
        "tiers": tiers,
        "ceiling": _REACHABILITY_CEILING,
        "within_ceiling": len(unmounted) <= _REACHABILITY_CEILING,
        # A name in neither direction is the failure the ratchet exists for,
        # so a row carries both lists rather than only the count.
        "newly_unreachable": newly,
        "stale_in_baseline": stale,
    }


# ---------------------------------------------------------------------------
# Per agent: offered tools and system-prompt size
# ---------------------------------------------------------------------------

def _tools_server_names(tree_root: Path) -> frozenset[str]:
    """Tool names ``delfin/tools/mcp_server.py`` registers, by AST.

    Read from the decorators rather than by importing the server: the module
    needs the MCP SDK at run time, and a metrics row must be recordable
    without it. ``@mcp.tool()`` and ``@mcp.tool(name="x")`` are both seen;
    the bare form names the tool after the function it decorates.
    """
    path = tree_root / "delfin" / "tools" / "mcp_server.py"
    if not path.is_file():
        return frozenset()
    try:
        tree = ast.parse(path.read_text(errors="replace"))
    except SyntaxError:
        return frozenset()
    found: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for dec in node.decorator_list:
            text = ast.unparse(dec)
            if ".tool(" not in text and not text.endswith(".tool"):
                continue
            named = re.search(r"""name\s*=\s*['"]([a-z_][a-z0-9_]*)['"]""", text)
            found.add(named.group(1) if named else node.name)
    return frozenset(found)


def _ops_server_names(tree_root: Optional[Path] = None) -> frozenset[str]:
    """Tool names the ops server exposes, read from an AST of its own file.

    Derived from the registrations rather than hand-listed, so a tool the
    server registers later is counted automatically.

    Read as SOURCE from the measured tree when one is given, not imported.
    Measured 2026-11-25: the import form read the *installed* package, so the
    PR table -- which claims that one instrument reads both trees -- counted
    a PR's own added or removed ops tool on both sides and cancelled it out,
    reporting an unchanged surface for a change it exists to show. The import
    also drags the MCP SDK in, which a metrics row must not need.

    Both forms are accepted because both are used in the server: the
    decorator (``@mcp.tool(name="x")``) and the functional registration
    (``mcp.tool(name="x")(fn)``). Only ``name=`` is read, since a tool named
    after a wrapped function is not how this server registers anything --
    if that ever changes, the count drops and
    ``tests/test_a_metrics_row_measures_what_the_product_writes.py``
    (``test_the_shipped_ops_surface_is_still_read_in_full``) fails rather
    than the surface quietly shrinking.
    """
    if tree_root is not None:
        path = Path(tree_root) / "delfin" / "ops_server" / "server.py"
        if not path.is_file():
            # A tree without the ops server is not a tree whose ops server
            # exposes nothing. Measured 2026-11-25 through build_row: the
            # empty set published `delfin-ops: 0`, which reads as a
            # regression against any tree that has the server.
            raise OpsSurfaceUnreadable(f"no ops server at {path}")
        try:
            tree = ast.parse(path.read_text(errors="replace"))
        except (SyntaxError, OSError) as exc:
            raise OpsSurfaceUnreadable(
                f"cannot parse {path}: {type(exc).__name__}: {exc}") from exc
        return frozenset(
            m.group(1) for node in ast.walk(tree)
            for m in [re.search(r"""\.tool\(\s*name\s*=\s*['"]([A-Za-z0-9_]+)['"]""",
                                ast.unparse(node))]
            if m is not None)
    try:
        from delfin.ops_server import server as ops_server
    except Exception as exc:                                # noqa: BLE001
        # Reached only when no tree was given. Returning an empty set here
        # would be the "absent reads as zero" defect: the caller cannot tell
        # "the server exposes no tools" from "the server could not be
        # imported". Raised instead, so `mcp_tool_names` can record which.
        raise OpsSurfaceUnreadable(
            f"cannot import the ops server: {type(exc).__name__}: {exc}"
        ) from exc
    return frozenset(n[len("tool_"):] for n in dir(ops_server)
                     if n.startswith("tool_"))


class OpsSurfaceUnreadable(RuntimeError):
    """The ops server's tool surface could not be read at all."""


def mcp_tool_names(tree_root: Path,
                   errors: Optional[dict[str, str]] = None
                   ) -> dict[str, frozenset[str]]:
    """The MCP servers production starts, and the tool names each exposes.

    Only the two DELFIN ships. A user's own server is discovered inside a
    turn and has nothing to measure at record time, so it is not guessed at.

    Both halves are read from ``tree_root``, so a row measures ONE tree.

    An unreadable server is reported into *errors* rather than raised: a
    checkout without the ops server must not fail the whole row, and must
    not report that server as exposing zero tools either. Measured
    2026-11-25 -- a bare ``except Exception`` made "could not read" and
    "registers nothing" the same answer.
    """
    sink: dict[str, str] = {} if errors is None else errors
    out: dict[str, frozenset[str]] = {}
    for server, read in (("delfin-tools", _tools_server_names),
                         ("delfin-ops", _ops_server_names)):
        try:
            out[server] = read(Path(tree_root))
        except OpsSurfaceUnreadable as exc:
            out[server] = frozenset()
            sink[server] = str(exc)
    return out


def mcp_surface_row(tree_root: Path) -> dict[str, Any]:
    """Each MCP server's tool count, and why it is missing when it is.

    The row records an unreadable MCP server as an error with a reason,
    not as a server that exposes no tools. Measured 2026-11-25 -- a bare
    ``except Exception`` made the two indistinguishable.
    """
    errors: dict[str, str] = {}
    names = mcp_tool_names(Path(tree_root), errors)
    servers = {server: len(found) for server, found in names.items()}
    return {"servers": servers, "errors": errors, "complete": not errors}


def _role_sources() -> list[dict[str, Any]]:
    """The three places a runnable role can be declared, lazily imported.

    Each entry resolves its names at call time, so an import that fails is
    caught by the caller and reported by ``role_sources_row`` rather than
    silently reducing the row's agent count.
    """

    def _presets() -> list[str]:
        from .subagents import subagent_type_names
        return sorted(subagent_type_names())

    def _whitelisted() -> list[str]:
        from .engine import _ROLE_TOOL_WHITELIST
        return sorted(_ROLE_TOOL_WHITELIST)

    def _allowlisted() -> list[str]:
        from .api_client import _ROLE_EXEC_ALLOWLIST
        return sorted(_ROLE_EXEC_ALLOWLIST)

    return [
        {"source": "subagent presets", "names": _presets},
        {"source": "engine role whitelist", "names": _whitelisted},
        {"source": "api_client exec allowlist", "names": _allowlisted},
    ]


def role_sources_row() -> dict[str, Any]:
    """Which role sources this row could read, and which it could not.

    The row's agent count is a claim about "every agent", and a source that
    failed to import makes that claim false while leaving the count looking
    complete. Measured 2026-11-25: ten of twenty-one roles were absent and
    nothing in the row said so. This is the field that says so.
    """
    read: dict[str, int] = {}
    unreadable: dict[str, str] = {}
    for source in _role_sources():
        try:
            read[source["source"]] = len(list(source["names"]()))
        except Exception as exc:                            # noqa: BLE001
            unreadable[source["source"]] = f"{type(exc).__name__}: {exc}"
    return {"read": read, "unreadable": unreadable,
            "complete": not unreadable}


def role_names() -> list[str]:
    """Every agent the product can run, generalist first.

    Three sources, because a role can run without being a delegation preset
    and this row's question is "what is EACH agent offered", not "what can
    the generalist delegate to":

    * the empty role and ``solo_agent`` -- the generalist, spelled two ways;
    * the subagent presets (``subagents.subagent_type_names``), so a preset
      added to the pack appears without an edit here;
    * every role the engine restricts with its own tool whitelist and every
      role ``api_client`` restricts with an execution allow-list. These are
      the SESSION roles -- ``dashboard_agent``, ``office_agent``,
      ``chief_agent`` and the rest. Measured 2026-11-25: the preset registry
      alone left ten of them out of every row, including the two with the
      tightest allow-lists in the product, so a row that reported "every
      agent" measured eleven of twenty-one. The question this answers for
      ``dashboard_agent`` -- what is it offered? -- had no answer at all.

    Anything a source cannot supply is skipped rather than guessed: a role
    this row cannot describe is missing from the row, never reported with
    the generalist's surface.
    """
    names = ["", "solo_agent"]
    for source in _role_sources():
        try:
            names.extend(sorted(source["names"]()))
        except Exception:                                   # noqa: BLE001
            # A source that cannot be imported must not shrink the row in
            # silence: "21 agents" and "21 agents, one source unreadable"
            # are different claims, and only the second is honest. The
            # reason is recorded by `role_sources_row` into the row itself.
            continue
    seen: set[str] = set()
    out: list[str] = []
    for name in names:
        if name not in seen:
            seen.add(name)
            out.append(name)
    return out


def offered_tools(role: str, tree_root: Path) -> dict[str, Any]:
    """The tool names one agent is offered, catalogue AND MCP servers.

    The two halves are reported apart because they are derived differently
    and a reader debugging a surface needs to know which half moved. The
    total is what the model's request carries.
    """
    from .api_client import (
        ToolSurfaceContext, _role_allowlist, role_tool_surface_report,
        tool_unavailable_reason,
    )

    allow = _role_allowlist(role)
    ctx = ToolSurfaceContext(role=role,
                             session_allowed=frozenset(allow or ()))
    # The catalogue half comes from the product's own per-role report rather
    # than a second filter written here: it already applies the role's
    # allowlist through ``tool_unavailable_reason``, and two readings of the
    # same rules is how a count and the surface it claims to describe drift
    # apart. What it does NOT cover is the MCP servers below, which is why
    # this module exists at all.
    surface = role_tool_surface_report([role])[role]
    catalogue_names = list(surface["names"])
    mcp: list[str] = []
    mcp_errors: dict[str, str] = {}
    for server, names in mcp_tool_names(tree_root, mcp_errors).items():
        for name in sorted(names):
            namespaced = f"mcp__{server}__{name}"
            if tool_unavailable_reason(namespaced, ctx) is None:
                mcp.append(namespaced)
    result = {
        "declared_tools": len(allow) if allow else 0,
        "no_preset_restriction": allow is None,
        "catalogue_count": int(surface["count"]),
        "catalogue_tokens": int(surface["total_tokens"]),
        "mcp_count": len(mcp),
        "count": int(surface["count"]) + len(mcp),
        "catalogue_names": sorted(catalogue_names),
        "mcp_names": sorted(mcp),
    }
    if mcp_errors:
        result["mcp_errors"] = mcp_errors
    return result


def prompt_size(role: str) -> dict[str, Any]:
    """Chars and estimated tokens of one agent's composed system prompt.

    The house estimator (``chars // 4``) is the loader's own, so this number
    and the enforced prompt budget are the same currency. A role the loader
    cannot compose says so instead of reporting a zero that reads like a
    measurement.
    """
    try:
        from .prompt_loader import PromptLoader
    except Exception as exc:                                # noqa: BLE001
        return {"error": f"prompt loader unavailable: {exc}"}
    loader = PromptLoader()
    role_id = role or "solo_agent"
    if not loader.role_exists(role_id):
        return {"error": f"unknown role {role_id!r}"}
    try:
        report = loader.prompt_size_report(
            role_id=role_id, mode_id="solo", task_text="prompt size report",
            permission_mode="default")
    except Exception as exc:                                # noqa: BLE001
        return {"error": f"{type(exc).__name__}: {exc}"}
    return {
        "prompt_chars": report["total_chars"],
        "prompt_tokens": report["total_tokens"],
        "sections": len(report["sections"]),
    }


def agents_row(tree_root: Path) -> dict[str, Any]:
    """One entry per agent: what it is offered and what its prompt costs."""
    out: dict[str, Any] = {}
    for role in role_names():
        entry = {"role": role or "(generalist)"}
        entry.update(offered_tools(role, tree_root))
        entry.update(prompt_size(role))
        out[role or "(generalist)"] = entry
    return out


# ---------------------------------------------------------------------------
# Quality: read from the kept benchmark rounds, never stored a second time
# ---------------------------------------------------------------------------
#
# The benchmark rounds are already kept by ``bench_history`` (one line per
# round, one per instance, each instance carrying its capability tag), and
# the interval on a difference of pass rates is ``bench_ab.newcombe_interval``.
# This section only READS those: a metrics row does not copy a quality number
# into its own file, so the two can never disagree.

def _rounds_for(rounds: list, commit: str) -> list:
    """The kept rounds of one commit, matched by prefix in either direction."""
    c = str(commit or "").split("-dirty")[0]
    if not c:
        return []
    return [r for r in rounds
            if r.commit and (r.commit.startswith(c) or c.startswith(r.commit))]


def _tally(rounds: list) -> dict[str, dict[str, Any]]:
    """Pass counts overall and per capability tag over *rounds*.

    The capability tag is ``bench_history.inst_capability`` -- the field the
    benchmark writer stamps on every instance. An instance without a tag is
    counted under ``(untagged)`` rather than dropped.
    """
    from .bench_history import inst_capability

    out: dict[str, dict[str, Any]] = {}
    for rnd in rounds:
        per_round: dict[str, list[int]] = {}
        for inst in rnd.graded():
            ok = 1 if inst.get("success") else 0
            for key in ("all", f"capability:{inst_capability(inst) or '(untagged)'}"):
                per_round.setdefault(key, [0, 0])
                per_round[key][0] += ok
                per_round[key][1] += 1
        for key, (p, n) in per_round.items():
            entry = out.setdefault(key, {"passed": 0, "n": 0, "rounds": 0,
                                         "round_rates": []})
            entry["passed"] += p
            entry["n"] += n
            entry["rounds"] += 1
            entry["round_rates"].append(p / n)
    for entry in out.values():
        rates = entry["round_rates"]
        entry["rate"] = entry["passed"] / entry["n"] if entry["n"] else None
        entry["min"] = min(rates) if rates else None
        entry["max"] = max(rates) if rates else None
    return out


def quality_for(commit: str, history_root: Optional[Path] = None
                ) -> dict[str, Any]:
    """Quality of one commit from the kept benchmark rounds, with spread.

    Repeats are the kept rounds of that commit; the spread is the min/max of
    the per-round pass rates. Per AGENT quality is not available: a kept
    round records the model, not which agent produced each answer, and the
    row says so instead of relabelling the model as an agent.
    """
    from . import bench_history

    try:
        rounds = _rounds_for(bench_history.load(history_root), commit)
    except Exception as exc:                                # noqa: BLE001
        return {"rounds": 0, "error": f"{type(exc).__name__}: {exc}"}
    measured = [r for r in rounds if not bench_history.not_measurable_reason(r)]
    return {
        "source": "bench history (kept rounds)",
        "rounds": len(measured),
        "unmeasured_rounds": len(rounds) - len(measured),
        "models": sorted({r.model for r in measured}),
        "numbers": _tally(measured),
        "per_agent": None,
        "note": ("per-agent quality is not recorded by a benchmark round; "
                 "the numbers are per capability tag and overall"),
    }


def compare_quality(old_commit: str, new_commit: str,
                    history_root: Optional[Path] = None) -> list[dict[str, Any]]:
    """better / worse / within noise for overall and each capability.

    The two commits' rounds must agree on every ``bench_history.comparable``
    axis (seed mode, task-set version, model); otherwise the pair is
    ``not comparable`` and no direction is claimed. The verdict is read off
    ``bench_ab.newcombe_interval`` on the pooled pass counts: an interval
    that contains zero is ``within noise``.
    """
    from . import bench_history
    from .bench_ab import newcombe_interval

    rounds = bench_history.load(history_root)
    old_r = [r for r in _rounds_for(rounds, old_commit)
             if not bench_history.not_measurable_reason(r)]
    new_r = [r for r in _rounds_for(rounds, new_commit)
             if not bench_history.not_measurable_reason(r)]
    if not old_r or not new_r:
        return [{"number": "quality.all", "verdict": "not measured",
                 "reason": (f"kept rounds: {len(old_r)} for the old commit, "
                            f"{len(new_r)} for the new one")}]
    for r in old_r[1:] + new_r:
        ok, why = bench_history.comparable(r, old_r[0])
        if not ok:
            return [{"number": "quality.all", "verdict": "not comparable",
                     "reason": f"round {r.round_id}: {why}"}]
    o_t, n_t = _tally(old_r), _tally(new_r)
    out: list[dict[str, Any]] = []
    for key in sorted(set(o_t) | set(n_t), key=lambda k: (k != "all", k)):
        o, n = o_t.get(key), n_t.get(key)
        entry: dict[str, Any] = {
            "number": f"quality.{key}",
            "old": None if not o else o["rate"],
            "new": None if not n else n["rate"],
            "old_rounds": 0 if not o else o["rounds"],
            "new_rounds": 0 if not n else n["rounds"],
        }
        if not o or not n:
            entry.update(verdict="not measured", reason="one side has no rounds")
            out.append(entry)
            continue
        lo, hi = newcombe_interval(o["passed"], o["n"], n["passed"], n["n"])
        entry["delta"] = n["rate"] - o["rate"]
        entry["interval"] = [lo, hi]
        entry["verdict"] = ("better" if lo > 0 else "worse" if hi < 0
                            else "within noise")
        entry["reason"] = (f"{o['passed']}/{o['n']} -> {n['passed']}/{n['n']}, "
                           f"95% interval on the difference [{lo:+.2f}, {hi:+.2f}]")
        out.append(entry)
    return out


# ---------------------------------------------------------------------------
# Comparison: better / worse / within noise
# ---------------------------------------------------------------------------

def compare_rows(old: dict[str, Any], new: dict[str, Any], *,
                 bench_history_root: Optional[Path] = None) -> dict[str, Any]:
    """Every number the two rows agree on, classified one by one.

    Reachability is a count that should not rise; quality comes from the
    kept benchmark rounds of the two commits; the per-agent offered-tool
    count and prompt size are reported (a larger surface is not by itself
    better or worse).
    """
    comparisons: list[dict[str, Any]] = []

    o_reach = (old.get("reachability") or {})
    n_reach = (new.get("reachability") or {})
    if o_reach.get("unmounted") is not None and n_reach.get("unmounted") is not None:
        d = n_reach["unmounted"] - o_reach["unmounted"]
        comparisons.append({
            "number": "reachability.unmounted",
            "old": o_reach["unmounted"], "new": n_reach["unmounted"],
            "delta": d,
            "verdict": ("better" if d < 0 else "worse" if d > 0
                        else "within noise"),
            "reason": ("counted, not sampled: every capability is enumerated, "
                       "so there is no spread to exceed"),
            "never_worse": d <= 0,
        })

    comparisons.extend(compare_quality(old.get("commit", ""),
                                       new.get("commit", ""),
                                       bench_history_root))

    o_surf = (old.get("agents") or {})
    n_surf = (new.get("agents") or {})
    for agent in sorted(set(o_surf) & set(n_surf)):
        for key in ("count", "prompt_tokens"):
            o_val = o_surf[agent].get(key)
            n_val = n_surf[agent].get(key)
            if o_val is None or n_val is None:
                continue
            comparisons.append({
                "number": f"agents.{agent}.{key}",
                "old": o_val, "new": n_val, "delta": n_val - o_val,
                "verdict": "reported" if o_val != n_val else "within noise",
                "reason": ("a surface size is not better or worse by itself; "
                           "it is reported so a change in it is visible"),
            })

    verdicts = [c["verdict"] for c in comparisons]
    return {
        "old_commit": old.get("commit", ""),
        "new_commit": new.get("commit", ""),
        "comparisons": comparisons,
        "counts": {v: verdicts.count(v) for v in sorted(set(verdicts))},
        "never_worse_violations": [
            c["number"] for c in comparisons
            if c.get("never_worse") is False],
        "any_worse": "worse" in verdicts,
    }


def format_comparison(result: dict[str, Any]) -> str:
    """The before/after table, for a person or a PR comment."""
    lines = [
        f"metrics {result.get('old_commit', '?')[:9]} -> "
        f"{result.get('new_commit', '?')[:9]}",
        f"{'number':<38} {'old':>10} {'new':>10} {'delta':>10}  verdict",
    ]
    for c in result.get("comparisons", []):
        def _fmt(v):
            if v is None:
                return "-"
            if isinstance(v, float):
                return f"{v:.4g}"
            return str(v)
        delta = c.get("delta")
        lines.append(
            f"{c['number']:<38} {_fmt(c.get('old')):>10} "
            f"{_fmt(c.get('new')):>10} {_fmt(delta):>10}  {c['verdict']}")
    counts = result.get("counts") or {}
    lines.append("")
    lines.append("summary: " + ", ".join(f"{k}={v}" for k, v in counts.items()))
    for name in result.get("never_worse_violations") or []:
        lines.append(f"NEVER-WORSE VIOLATION: {name}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# The row, and the kept history
# ---------------------------------------------------------------------------

def git_commit(tree_root: Path) -> str:
    """The commit this row describes, read from the tree itself.

    ``-dirty`` is appended when the working tree has changes: a row that
    describes a commit while measuring uncommitted code is a row about
    neither.
    """
    from . import safe_git
    try:
        sha = safe_git.run(["rev-parse", "--short=9", "HEAD"], str(tree_root),
                           timeout=15).stdout.strip()
        dirty = safe_git.run(["status", "--porcelain"], str(tree_root),
                             timeout=30).stdout.strip()
    except Exception:                                       # noqa: BLE001
        return "unknown"
    return f"{sha or 'unknown'}{'-dirty' if dirty else ''}"


def build_row(tree_root: Path, *, commit: str = "") -> dict[str, Any]:
    """The kept row for one commit: reachability and what each agent is offered.

    Needs no model call and no key. Quality is not copied into the row: it is
    read from the kept benchmark rounds (``bench_history``) when a row is
    shown or compared, so there is one store of benchmark results.
    """
    tree_root = Path(tree_root)
    row: dict[str, Any] = {
        "schema": 1,
        "commit": commit or git_commit(tree_root),
        "reachability": reachability(tree_root),
        "agents": agents_row(tree_root),
        "agent_sources": role_sources_row(),
        "mcp_surface": mcp_surface_row(tree_root),
        "packs": _packs_row(tree_root),
    }
    return row


def _packs_row(tree_root: Path) -> dict[str, Any]:
    """How much shipped material the agents are given, as counts.

    Reported because a prompt can grow without a tool changing: skills,
    agent prompts and knowledge notes all land in the system prompt. Counting
    them separately is what makes such a growth visible in a diff of rows.

    ``knowledge/operating/`` is counted too. It was 0 before, because the
    glob was ``knowledge/*.md`` while ten notes live one level down -- and
    those ten DO reach an agent's turn, so a row that counted only the top
    level reported the note layer as ten files smaller than it is. The two
    levels are kept apart rather than summed: the taxonomy is flat by
    design (``knowledge.global_dir`` reads ``*.md``) and the operating notes
    are a separate layer, so a reader can see which one moved.

    A missing directory reports ``None``, not 0: "this tree has no skills"
    and "this tree has no ``pack/`` I could read" are different answers, and
    the second one is a defect in the measurement.
    """
    pack = tree_root / "delfin" / "agent" / "pack"
    out: dict[str, Any] = {}
    for name, sub in (("agents", "agents"), ("skills", "skills"),
                      ("knowledge", "knowledge"), ("shared", "shared")):
        directory = pack / sub
        out[name] = (len(list(directory.glob("*.md")))
                     if directory.is_dir() else None)
    operating = pack / "knowledge" / "operating"
    out["knowledge_operating"] = (len(list(operating.glob("*.md")))
                                  if operating.is_dir() else None)
    return out


def history_path(root: Optional[Path] = None) -> Path:
    """Where the kept rows live: one append-only JSONL file per root.

    ``~/.delfin/metrics_history`` by default -- outside any checkout, because
    a measurement committed with the code it grades is a measurement someone
    can edit to agree with a claim.
    """
    if root is not None:
        return Path(root) / "rows.jsonl" if Path(root).suffix != ".jsonl" \
            else Path(root)
    return Path.home() / ".delfin" / "metrics_history" / "rows.jsonl"


def record_row(row: dict[str, Any], root: Optional[Path] = None) -> Path:
    """Append one row to the kept history and return the file written.

    Append-only and one line per row: a kept history that can be rewritten
    in place cannot answer "what did this commit measure" after the fact.
    """
    path = history_path(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(row, sort_keys=True) + "\n")
    return path


def read_rows(path: Path) -> list[dict[str, Any]]:
    """Every kept row in one file, oldest first."""
    rows: list[dict[str, Any]] = []
    if not Path(path).is_file():
        return rows
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return rows


def latest_for(rows: list[dict[str, Any]], commit: str) -> Optional[dict]:
    """The newest row for a commit, matched by prefix so a short SHA works."""
    for row in reversed(rows):
        if str(row.get("commit", "")).startswith(commit):
            return row
    return None


def to_json(row: dict[str, Any]) -> str:
    """The row as indented JSON, the form a PR comment or a diff carries."""
    return json.dumps(row, indent=2, sort_keys=True)
