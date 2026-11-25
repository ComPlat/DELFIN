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
* **Quality with repeats** -- the benchmark result with its spread. A single
  run is not a result: measured on the committed reference, one task scored
  45, 93, 93 across three runs of unchanged code, so a comparison built on
  one run per side would call that a regression or an improvement at random.

The fourth part is the comparison, and its rule is the point of the whole
module: a difference is ``better`` or ``worse`` only when it is larger than
the spread of the repeats says it could be by chance. Otherwise it is
``within noise``, including when the means moved a lot.

Why the reachability scan is here rather than imported
------------------------------------------------------
``tests/test_every_capability_has_a_caller.py`` is the instrument for the
first number, and this module does not import it:

1. ``tests/`` is not an importable package, so there is no import to write;
2. it is a pytest module with a module-scoped fixture, and importing it
   would drag pytest collection into every command that records a row;
3. worse, it would make the PRODUCT depend on the SUITE -- a metrics command
   that fails on an installation without pytest is not a measurement.

So the scan is a faithful copy here, and ``tests/
test_a_commit_row_measures_every_agent.py`` pins it against the original:
where both exist, ``definitions()`` and ``mentions()`` must answer
identically, and the unmounted count must equal the ratchet's. Divergence is
a defect of this module, not of the ratchet. The same pattern and the same
three reasons are already recorded in ``scripts/ab_task_pools.py``.

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
import math
import re
import statistics
import tokenize
from dataclasses import asdict, dataclass, field
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

#: Below this many repeats per agent per side, no comparison is made. Two
#: runs cannot tell a spread from a coincidence, and one run has no spread at
#: all -- the difference between "measured" and "looked at once".
_MIN_REPEATS = 3

#: How many standard errors a difference must clear before it is called
#: better or worse. Two, because the row is read by a person deciding whether
#: to look further, and the cost of a false "worse" is a wasted bisect while
#: the cost of a false "better" is a regression merged.
_NOISE_SIGMAS = 2.0


def _zero_spread_floor(n: int) -> float:
    """The smallest difference that cannot be a repeat coincidence, for *n*.

    Used only when every repeat on both sides agreed and the measured spread
    is therefore exactly 0. The score is a rate over a task set, so one task
    flipping moves the rate by ``1/k`` for a set of ``k`` tasks; the set size
    is not in the row, so the conservative reading is one repeat's worth of
    resolution, ``1/n``. Conservative in the direction that matters: it makes
    a zero-spread comparison HARDER to call, never easier.

    A project whose repeats are genuinely deterministic (an exact classifier,
    say) pays for this with repeats -- ``n >= 20`` gets the floor below 0.05
    -- which is the honest price of a claim that the change moved something.
    """
    return 1.0 / max(int(n), 1)


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

    *definition_root* defaults to ``<tree_root>/delfin/agent``, which is the
    scope the ratchet and its 86-name ceiling were measured on. It is a
    parameter so the same question can be asked of another root without
    copying the function again.
    """
    root = Path(definition_root) if definition_root else Path(tree_root) / "delfin" / "agent"
    out: dict[str, tuple[str, int]] = {}
    if not root.is_dir():
        return out
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
    ``delfin/agent/``. A count from another root answers another question and
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
            "population": "module-level callables under delfin/agent/",
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
        "population": "module-level callables under delfin/agent/",
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
    if that ever changes, the count drops and `test_every_ops_tool_is_counted`
    fails rather than the surface quietly shrinking.
    """
    if tree_root is not None:
        path = Path(tree_root) / "delfin" / "ops_server" / "server.py"
        if not path.is_file():
            return frozenset()
        try:
            tree = ast.parse(path.read_text(errors="replace"))
        except SyntaxError:
            return frozenset()
        return frozenset(
            m.group(1) for node in ast.walk(tree)
            for m in [re.search(r"""\.tool\(\s*name\s*=\s*['"]([A-Za-z0-9_]+)['"]""",
                                ast.unparse(node))]
            if m is not None)
    try:
        from delfin.ops_server import server as ops_server
    except Exception:                                       # noqa: BLE001
        return frozenset()
    return frozenset(n[len("tool_"):] for n in dir(ops_server)
                     if n.startswith("tool_"))


def mcp_tool_names(tree_root: Path) -> dict[str, frozenset[str]]:
    """The MCP servers production starts, and the tool names each exposes.

    Only the two DELFIN ships. A user's own server is discovered inside a
    turn and has nothing to measure at record time, so it is not guessed at.

    Both halves are read from ``tree_root``, so a row measures ONE tree.
    """
    return {
        "delfin-tools": _tools_server_names(Path(tree_root)),
        "delfin-ops": _ops_server_names(Path(tree_root)),
    }


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
    try:
        from .subagents import subagent_type_names
        names.extend(sorted(subagent_type_names()))
    except Exception:                                       # noqa: BLE001
        pass
    try:
        from .engine import _ROLE_TOOL_WHITELIST
        names.extend(sorted(_ROLE_TOOL_WHITELIST))
    except Exception:                                       # noqa: BLE001
        pass
    try:
        from .api_client import _ROLE_EXEC_ALLOWLIST
        names.extend(sorted(_ROLE_EXEC_ALLOWLIST))
    except Exception:                                       # noqa: BLE001
        pass
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
    for server, names in mcp_tool_names(tree_root).items():
        for name in sorted(names):
            namespaced = f"mcp__{server}__{name}"
            if tool_unavailable_reason(namespaced, ctx) is None:
                mcp.append(namespaced)
    return {
        "declared_tools": len(allow) if allow else 0,
        "no_preset_restriction": allow is None,
        "catalogue_count": int(surface["count"]),
        "catalogue_tokens": int(surface["total_tokens"]),
        "mcp_count": len(mcp),
        "count": int(surface["count"]) + len(mcp),
        "catalogue_names": sorted(catalogue_names),
        "mcp_names": sorted(mcp),
    }


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
# Quality: repeats, spread, and the noise rule
# ---------------------------------------------------------------------------

def spread(values: list[float]) -> dict[str, Any]:
    """mean, min, max and the standard error of a set of repeats.

    The standard error is computed from the sample's own spread, so a row
    measured with more repeats tightens every comparison made from it by
    itself, rather than a percentage fixed in advance that is too loose at
    three runs and absurd at thirty.
    """
    n = len(values)
    if n == 0:
        return {"n": 0, "mean": None, "min": None, "max": None, "stderr": None}
    mean = statistics.fmean(values)
    if n < 2:
        return {"n": n, "mean": mean, "min": min(values), "max": max(values),
                "stderr": None}
    stdev = statistics.stdev(values)
    return {
        "n": n,
        "mean": mean,
        "min": min(values),
        "max": max(values),
        "stdev": stdev,
        "stderr": stdev / math.sqrt(n),
    }


def quality_row(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """Per-agent quality from repeated benchmark runs.

    *runs* is one dict per repeat. Two shapes are accepted, and the reason is
    a defect measured on 2026-11-25: the first shape below was the only one
    this module understood, and NOTHING in the product writes it, so feeding
    it the file ``delfin-agent bench run`` actually produces recorded
    ``agents: {}`` -- a quality half that had never read a real run.

    * what ``benchmark.write_run`` writes, one result row per (task, repeat):
      ``{"task_id", "quality_0_100", "success", "profile_name", ...}``. Grouped
      here by ``profile_name``, which is the MODEL profile the result was
      stamped with. That is not an agent role, and this module does not
      pretend it is: per-agent quality needs a run in which the agent is
      known, which is why the summary says which column it grouped by.
    * the agent-labelled summary this command was designed for:
      ``{"agent": {name: score}, "delegation": {...}, "delegation_expected":
      {...}}``. A caller that knows which agent produced each score -- a
      harness that runs the generalist and a specialist and labels the rows
      -- writes this, and then the row IS per agent.

    An agent measured fewer than ``_MIN_REPEATS`` times is reported with its
    mean and marked ``is_result: false``: the number is what was seen, and
    saying it is a result would be the lie the repeats exist to prevent.
    """
    labelled = [r for r in runs if isinstance(r.get("agent"), dict)]
    if labelled:
        per_agent: dict[str, list[float]] = {}
        for run in labelled:
            for agent, value in (run.get("agent") or {}).items():
                if value is None:
                    continue
                per_agent.setdefault(agent, []).append(float(value))
        grouped_by = "agent"
    else:
        per_agent, grouped_by = _per_agent_from_bench_rows(runs)

    agents: dict[str, Any] = {}
    for agent, values in sorted(per_agent.items()):
        stats = spread(values)
        stats["is_result"] = stats["n"] >= _MIN_REPEATS
        stats["values"] = values
        agents[agent] = stats

    delegation = delegation_accuracy(runs)
    out: dict[str, Any] = {
        "repeats": len(runs),
        "min_repeats": _MIN_REPEATS,
        "grouped_by": grouped_by,
        "agents": agents,
        "delegation": delegation,
    }
    if not agents:
        out["note"] = (
            f"{len(runs)} row(s) were read but none carried a score to group: "
            "expected benchmark result rows (quality_0_100 / success) or an "
            "`agent` mapping. This is not a quality of zero.")
    elif grouped_by == "profile_name":
        out["note"] = (
            "grouped by profile_name, which is the MODEL profile stamped on "
            "each result -- not an agent role. Per-agent quality needs a run "
            "that records which agent produced each score.")
    return out


def _per_agent_from_bench_rows(
        runs: list[dict[str, Any]]) -> tuple[dict[str, list[float]], str]:
    """Group real benchmark result rows into per-profile quality samples.

    A "repeat" here is one *file*: one complete pass of the suite. The score
    for that repeat is the mean ``quality_0_100`` over the tasks it covered,
    so a repeat that covered more tasks cannot outweigh one that covered
    fewer, and a suite half-run is reported as the rate it is.

    Unmeasured rows are skipped, exactly as ``aggregate_replicates`` skips
    them: a run that never reached the model is not a failing model. When
    every row is unmeasured the profile gets no sample at all rather than a
    zero -- ``is_result`` stays false and the caller sees ``n: 0``.
    """
    per_profile: dict[str, list[float]] = {}
    for row in runs:
        if not isinstance(row, dict) or row.get("quality_0_100") is None:
            continue
        if row.get("unmeasured"):
            continue
        who = str(row.get("profile_name") or row.get("model") or "default")
        per_profile.setdefault(who, []).append(float(row["quality_0_100"]))
    # One sample per profile: this file is one repeat of the suite.
    return ({who: [sum(values) / len(values)]
             for who, values in per_profile.items() if values},
            "profile_name")


def delegation_accuracy(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """Did the generalist hand each task to the specialist that owned it.

    Counted over the runs supplied, per (task, expected specialist, actual).
    A task with no expected specialist in the run file is not counted at
    all -- an unlabelled task is not a wrong delegation, and counting it as
    one would move the rate with the size of the task set instead of with the
    routing.
    """
    hits = 0
    total = 0
    wrong: dict[str, dict[str, int]] = {}
    per_run: list[float] = []
    for run in runs:
        labels = run.get("delegation_expected") or {}
        if not labels:
            continue
        handed = run.get("delegation") or {}
        run_hits = 0
        run_total = 0
        for task, expected in labels.items():
            actual = handed.get(task)
            run_total += 1
            if actual == expected:
                hits += 1
                run_hits += 1
            else:
                key = str(actual or "(none)")
                wrong.setdefault(str(expected), {})
                wrong[str(expected)][key] = wrong[str(expected)].get(key, 0) + 1
        total += run_total
        per_run.append(run_hits / run_total if run_total else 0.0)
    if total == 0:
        return {"tasks": 0, "accuracy": None, "note":
                "no run carried delegation_expected labels, so routing was "
                "not measured -- this is not a zero"}
    stats = spread([r for r in per_run])
    return {
        "tasks": total,
        "hits": hits,
        "accuracy": hits / total,
        "per_run": stats,
        "misroutes": wrong,
    }


# ---------------------------------------------------------------------------
# Comparison: better / worse / within noise
# ---------------------------------------------------------------------------

def compare_number(name: str, old: Optional[dict], new: Optional[dict],
                   *, higher_is_better: bool = True) -> dict[str, Any]:
    """Classify one number between two commits, using its own spread.

    The rule, and the reason the row needs repeats at all:

    * either side with fewer than ``_MIN_REPEATS`` samples -> ``not a result``
      and no direction is claimed, whatever the means say;
    * ``|delta| > _NOISE_SIGMAS * sqrt(stderr_old**2 + stderr_new**2)`` ->
      ``better`` or ``worse`` by the sign and the direction of merit;
    * otherwise -> ``within noise``.

    A one-run difference on each side can therefore never be called an
    improvement, which is the whole point: on unchanged code this project's
    own benchmark moved one task from 45 to 93.
    """
    out: dict[str, Any] = {
        "number": name,
        "old": None if not old else old.get("mean"),
        "new": None if not new else new.get("mean"),
        "old_n": 0 if not old else old.get("n", 0),
        "new_n": 0 if not new else new.get("n", 0),
    }
    if not old or not new or old.get("mean") is None or new.get("mean") is None:
        out.update(verdict="not measured", reason="one side has no value")
        return out
    delta = new["mean"] - old["mean"]
    out["delta"] = delta
    out["old_spread"] = [old.get("min"), old.get("max")]
    out["new_spread"] = [new.get("min"), new.get("max")]
    if out["old_n"] < _MIN_REPEATS or out["new_n"] < _MIN_REPEATS:
        out.update(
            verdict="not a result",
            reason=(f"{out['old_n']} and {out['new_n']} repeats; "
                    f"{_MIN_REPEATS} are needed before a difference can be "
                    "told from the spread"))
        return out
    # A MISSING spread is not a zero spread. A row read back from the kept
    # history, or written by a caller that supplied a mean and an n but no
    # stderr, has no spread to compare against; reading the absent value as
    # 0.0 would let it through the zero-spread branch below as if its
    # repeats had agreed, which is the opposite of what an absent value
    # means. Measured 2026-11-25: `stderr: None` and `stderr: 0.0` took the
    # same path.
    missing = [side for side, entry in (("old", old), ("new", new))
               if entry.get("stderr") is None]
    if missing and not all(
            entry.get("n", 0) >= _MIN_REPEATS and (entry.get("max") == entry.get("min"))
            for entry in (old, new)):
        out.update(
            verdict="not a result",
            reason=(f"{' and '.join(missing)} side carries no spread, so no "
                    "difference can be told from noise; record the repeats "
                    "rather than a mean"))
        return out
    se = math.sqrt((old.get("stderr") or 0.0) ** 2
                   + (new.get("stderr") or 0.0) ** 2)
    threshold = _NOISE_SIGMAS * se
    out["threshold"] = threshold
    if se == 0.0:
        # Every repeat agreed on BOTH sides, so the spread measured here is
        # zero. That is NOT the same as "no noise": a per-run score is a
        # rate over a task set, and a binomial process gives identical
        # repeats by chance -- measured on this project's own reference, one
        # task scored 45, 93, 93 and an afternoon of re-running the same
        # code would eventually produce three equal runs and a false
        # direction. So a zero spread falls back to the range the repeats
        # must have covered, which cannot be zero for a rate in [0, 1]: at
        # least 1/n of it can move by chance.
        floor = _zero_spread_floor(max(out["old_n"], out["new_n"]))
        out["threshold"] = floor
        out["reason_no_spread"] = (
            "every repeat agreed on both sides, so the measured spread is 0; "
            f"the 1/n floor for {max(out['old_n'], out['new_n'])} repeats "
            "bounds it instead of treating three coincidences as certainty")
        if abs(delta) <= floor:
            out.update(verdict="within noise", reason=out["reason_no_spread"])
            return out
        out["verdict"] = ("better" if (delta > 0) == higher_is_better
                          else "worse")
        out["reason"] = out["reason_no_spread"]
        return out
    if abs(delta) <= threshold:
        out.update(verdict="within noise",
                   reason=(f"|{delta:+.4g}| <= {threshold:.4g} "
                           f"({_NOISE_SIGMAS:g} standard errors)"))
        return out
    improved = (delta > 0) == higher_is_better
    out.update(verdict="better" if improved else "worse",
               reason=(f"|{delta:+.4g}| > {threshold:.4g} "
                       f"({_NOISE_SIGMAS:g} standard errors)"))
    return out


def compare_rows(old: dict[str, Any], new: dict[str, Any]) -> dict[str, Any]:
    """Every number the two rows agree on, classified one by one.

    Reachability is compared as a count that should not rise; quality per
    agent; and the per-agent offered-tool count and prompt size, which are
    only reported (a larger surface is not by itself better or worse).
    """
    comparisons: list[dict[str, Any]] = []

    o_reach = (old.get("reachability") or {})
    n_reach = (new.get("reachability") or {})
    if "unmounted" in o_reach and "unmounted" in n_reach:
        d = n_reach["unmounted"] - o_reach["unmounted"]
        comparisons.append({
            "number": "reachability.unmounted",
            "old": o_reach["unmounted"], "new": n_reach["unmounted"],
            "delta": d,
            # Fewer unreached capabilities is the better direction, and this
            # is a count of defects, not a sample: no spread applies.
            "verdict": ("better" if d < 0 else "worse" if d > 0
                        else "within noise"),
            "reason": ("counted, not sampled: every capability is enumerated, "
                       "so there is no spread to exceed"),
            "never_worse": d <= 0,
        })

    o_agents = ((old.get("quality") or {}).get("agents") or {})
    n_agents = ((new.get("quality") or {}).get("agents") or {})
    for agent in sorted(set(o_agents) & set(n_agents)):
        comparisons.append(compare_number(
            f"quality.{agent}", o_agents[agent], n_agents[agent]))
    # An agent measured on ONE side only is compared too, and comes out
    # `not measured` rather than vanishing. Measured 2026-11-25: the
    # intersection-only loop dropped such an agent silently, so a specialist
    # that disappeared between two commits produced no row at all -- the
    # absence that this whole command exists to make visible.
    for agent in sorted(set(o_agents) ^ set(n_agents)):
        comparisons.append(compare_number(
            f"quality.{agent}", o_agents.get(agent), n_agents.get(agent)))

    o_del = ((old.get("quality") or {}).get("delegation") or {})
    n_del = ((new.get("quality") or {}).get("delegation") or {})
    if o_del.get("per_run") and n_del.get("per_run"):
        comparisons.append(compare_number(
            "quality.delegation_accuracy", o_del["per_run"], n_del["per_run"]))

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
    import subprocess
    try:
        sha = subprocess.run(["git", "rev-parse", "--short=9", "HEAD"],
                             cwd=str(tree_root), capture_output=True,
                             text=True, timeout=15).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"],
                               cwd=str(tree_root), capture_output=True,
                               text=True, timeout=30).stdout.strip()
    except Exception:                                       # noqa: BLE001
        return "unknown"
    return f"{sha or 'unknown'}{'-dirty' if dirty else ''}"


def build_row(tree_root: Path, *, runs: Optional[list[dict[str, Any]]] = None,
              commit: str = "") -> dict[str, Any]:
    """The whole row for one commit: reachability, agents, quality.

    *runs* is optional because the reachability half needs no model call and
    no key: a CI job with none still records the half it can measure, and the
    quality section then says why it is empty instead of reporting a zero.
    """
    tree_root = Path(tree_root)
    row: dict[str, Any] = {
        "schema": 1,
        "commit": commit or git_commit(tree_root),
        "reachability": reachability(tree_root),
        "agents": agents_row(tree_root),
        "packs": _packs_row(tree_root),
    }
    if runs:
        row["quality"] = quality_row(runs)
    else:
        row["quality"] = {
            "repeats": 0, "min_repeats": _MIN_REPEATS, "agents": {},
            "delegation": {"tasks": 0, "accuracy": None,
                           "note": "no benchmark runs supplied; quality was "
                                   "not measured for this commit"},
            "note": "no benchmark runs supplied (--bench-runs)",
        }
    row["costs"] = {"model_calls": 0 if not runs else len(runs)}
    return row


def _packs_row(tree_root: Path) -> dict[str, Any]:
    """How much shipped material the agents are given, as counts.

    Reported because a prompt can grow without a tool changing: skills,
    agent prompts and knowledge notes all land in the system prompt. Counting
    them separately is what makes such a growth visible in a diff of rows.
    """
    pack = tree_root / "delfin" / "agent" / "pack"
    out: dict[str, int] = {}
    for name, sub in (("agents", "agents"), ("skills", "skills"),
                      ("knowledge", "knowledge"), ("shared", "shared")):
        directory = pack / sub
        out[name] = (len(list(directory.glob("*.md")))
                     if directory.is_dir() else 0)
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


def load_bench_runs(path: Path) -> list[dict[str, Any]]:
    """Read a ``--bench-runs`` input: one path, or several separated by commas.

    Each path is either a JSON list or JSONL. This is the seam that lets the
    quality half be tested without a model: a file of measured per-run results
    is a first-class input, and a fake one proves the noise rule exactly as
    the real one would be measured by it.

    SEVERAL paths, because one file is one repeat. ``delfin-agent bench run``
    writes one file per invocation, and the noise rule needs three per side;
    a single path cannot express that, which would have made ``is_result``
    unreachable for every real benchmark. Each file's ``run_header`` line is
    skipped -- it is provenance, not a result -- and unparseable lines are
    named rather than dropped: a truncated run file that silently contributes
    fewer rows is a measurement of the wrong population.
    """
    runs: list[dict[str, Any]] = []
    for one in str(path).split(","):
        one = one.strip()
        if not one:
            continue
        file_path = Path(one)
        text = file_path.read_text(encoding="utf-8")
        stripped = text.strip()
        if not stripped:
            continue
        if stripped.startswith("["):
            data = json.loads(stripped)
            runs.extend(d for d in data if isinstance(d, dict))
            continue
        for number, line in enumerate(stripped.splitlines(), 1):
            line = line.strip()
            if not line:
                continue
            try:
                parsed = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"{file_path}:{number}: not JSON ({exc.msg}) -- a truncated "
                    "run file is a different population, not a shorter one"
                ) from exc
            if not isinstance(parsed, dict):
                continue
            if "run_header" in parsed and len(parsed) == 1:
                continue
            runs.append(parsed)
    return runs


def to_json(row: dict[str, Any]) -> str:
    """The row as indented JSON, the form a PR comment or a diff carries."""
    return json.dumps(row, indent=2, sort_keys=True)
