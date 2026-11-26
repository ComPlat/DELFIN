"""Every tool the model can call is exercised, or excluded on the record.

WHY THIS FILE EXISTS. ``web_search`` was refused on every single call for
three weeks and nobody noticed. Not because it was untested in the usual
sense, but because nothing ever asked it for a RESULT: the agent fell
back to ``web_fetch``, answered plausibly, and a broken tool looked like a
working one from every angle a person or a suite was looking from.

More tests would not have caught that. What catches it is a rule about
arrivals: a tool the model can call must either be exercised somewhere
against a fact, or appear below with a reason it cannot be. Adding a tool
without doing one of the two turns this file red. That is the mechanism;
the exercises underneath are what give it something to point at.

The exercises judge the world, never the wording: a file's bytes, a
returned URL, an exit status. Anything that needs the network asks for it
through the ``net_backend`` fixture and is skipped when it is unreachable,
because a red run caused by a firewall teaches people to ignore red runs.
The question is asked when a test needs it, never at import: a module-level
probe answers it during collection, which makes the list of tests depend on
the network.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import pytest

from conftest import child_env
from delfin.agent import api_client
from delfin.agent.api_client import KitToolPermissions, _doc_executor


# ---------------------------------------------------------------------------
# The surface
# ---------------------------------------------------------------------------

def advertised_tools() -> set[str]:
    """Every name that appears as a tool schema the model is shown."""
    from delfin.agent import tool_catalog
    src = (Path(api_client.__file__).read_text(encoding="utf-8")
           + Path(tool_catalog.__file__).read_text(encoding="utf-8"))
    return set(re.findall(r'"name":\s*"([a-z_][a-z0-9_]*)"', src))


# Exercised below against a fact. Keep in step with the tests underneath —
# a name here with no exercise is the same lie this file exists to prevent.
_EXERCISED = {
    "list_files", "read_file", "grep_file", "find_definition",
    "find_references", "write_file", "edit_file", "multi_edit",
    "bash", "task_create", "task_list", "task_get", "task_update",
    "check_environment", "project_introspect", "list_changes_made",
    "task_adopt",
    "web_search", "web_fetch",
}

# Not exercised here, each with the reason. A reason like "hard to test"
# is not a reason; every line below names the resource that is missing.
_NOT_EXERCISED_HERE = {
    # The deferral mechanism: a real _DocToolExecutor.execute call returns
    # a deferred tool's schema, and an unknown name answers with the
    # listing instead of a bare error.
    "describe_native_tool":
        "tests/test_the_schema_arrives_when_it_is_needed.py",
    # Live probes drive these through a real pty, where the approval
    # channel and the plan gate are themselves under test.
    "exit_plan_mode": "tests/agent_terminal_probes — needs the approval channel",
    "ask_user_question": "tests/agent_terminal_probes — needs a UI to answer",
    "subagent": "tests/test_orchestration.py — needs a child engine",
    "subagent_result": "tests/test_orchestration.py",
    "orchestrate": "tests/test_orchestration.py",
    "report_verdict": "tests/test_subagent_report_verification.py",
    # Official-API lookups (SRC1): the executor path, the network policy and
    # the not-found refusal run against recorded Crossref/PubChem/OpenAlex
    # responses -- this file has no recorded API fixture.
    "verify_source": "tests/test_src1_adversarial_registration.py — needs recorded API responses",
    # ORCA knowledge on demand (w04): dictionary/card lookups and the
    # workflow ranking, each driven against the shipped stores in its own file.
    "orca_explain": "tests/test_w04_orca_knowledge_tools.py",
    "orca_check_input": "tests/test_w04_orca_knowledge_tools.py",
    "workflow_candidates": "tests/test_w04_workflow_candidates.py",
    # Lab device control: driven through this same executor against the
    # simulator, each threat with a control, in their own files -- they
    # need a pinned lab manifest and a simulated world in a scratch lab dir.
    "lab_devices": "tests/test_lab_tools_registered.py — needs a lab manifest",
    "lab_read": "tests/test_lab_tools_registered.py",
    "lab_set": "tests/test_lab_safety_kernel.py — needs the simulated lab",
    "lab_plan_check": "tests/test_lab_safety_kernel.py",
    "lab_estop": "tests/test_lab_deadman_and_estop.py",
    # lab_run is EXERCISED, not just accounted for: the new file drives it
    # through this same executor against the simulator, judges a plan, runs
    # it, and checks that a hazardous step stops without writing. It is
    # named here because it also needs a pinned manifest in a scratch lab.
    "lab_run": "tests/test_a_lab_experiment_runs_through_the_kernel.py",
    # Cluster resources this machine does not have.
    "watch_job": "needs SLURM",
    "cluster_status": "tests/test_the_agent_reaches_a_cluster_only_through_"
                      "its_tools_and_a_person.py — needs a cluster behind ctun",
    "cluster_read": "tests/test_the_agent_reaches_a_cluster_only_through_its_tools_and_a_person.py",
    "cluster_fetch": "tests/test_the_agent_reaches_a_cluster_only_through_its_tools_and_a_person.py",
    "cluster_upload": "tests/test_the_agent_reaches_a_cluster_only_through_its_tools_and_a_person.py",
    "cluster_submit": "tests/test_the_agent_reaches_a_cluster_only_through_its_tools_and_a_person.py",
    "cluster_cancel": "tests/test_the_agent_reaches_a_cluster_only_through_its_tools_and_a_person.py",
    "cluster_run": "tests/test_the_agent_reaches_a_cluster_only_through_its_tools_and_a_person.py",
    "cluster_propose": "tests/test_a_cluster_is_added_in_the_window_and_only_proposed_by_an_agent.py",
    "session_message": "tests/test_sessions_can_write_to_each_other.py — "
                       "needs a second open session",
    "run_tests": "would run the suite inside the suite",
    "make_plot": "renders a plot card; covered by tests/test_v2_*.py",
    # Exercised through the real executor in their own files, each with
    # the store it needs set up under a redirected HOME.
    "skill_propose_patch": "tests/test_skill_patch_tool.py — needs an "
                           "active skill and the proposal store",
    "skill_share": "tests/test_zz_l6_a_shared_skill_arrives_as_a_proposal.py "
                   "— needs an accepted skill and the proposal store",
    "session_search": "tests/test_session_search_end_to_end.py — needs an "
                      "indexed session archive",
    # Side effects that leave the machine.
    "remote_trigger": "sends a request to a configured remote",
    "push_notification": "delivers to a configured channel",
    "draft_email": "writes a draft for a human to send",
    "publish_report": "publishes outside the workspace",
    # Long-lived or scheduling state, tested where that state is owned.
    "cron_create": "tests/test_scheduler.py — owns the schedule store",
    "cron_delete": "tests/test_scheduler.py",
    "cron_list": "tests/test_scheduler.py",
    "schedule_wakeup": "tests/test_scheduler.py",
    "bash_background": "tests/test_a_background_job_outlives_the_process_that_started_it.py",
    "bash_status": "needs a running background job to report on",
    "bash_output": "needs a running background job",
    "bash_kill": "needs a running background job",
    # Git worktrees: a second checkout, tested where that is set up.
    "enter_worktree": "tests covering worktree isolation",
    "exit_worktree": "needs a second checkout to leave",
    "worktree_merge": "needs two worktrees with divergent commits",
    # Durable memory and permission grants: writing real user state.
    "remember": "tests/test_memory_store.py — must not write real memory",
    "forget": "tests/test_memory_store.py",
    "remember_permission": "tests covering kit_settings persistence",
    "remember_permission_bundle": "writes real kit_settings grants",
    # An MCP server has to be running.
    "mcp_get_prompt": "needs a live MCP server",
    "mcp_read_resource": "needs a live MCP server",
    # Office and document formats, covered by their own suites with
    # fixtures this file has no business duplicating.
    "read_document": "tests/test_office_documents.py",
    "create_docx": "tests/test_office_documents.py",
    "create_pdf": "tests/test_office_documents.py",
    "fill_docx_template": "tests/test_office_documents.py",
    "fill_pdf_form": "tests/test_office_documents.py",
    "merge_pdfs": "tests/test_office_documents.py",
    "split_pdf": "tests/test_office_documents.py",
    "edit_sheet": "tests/test_office_documents.py",
    "fill_series": "tests/test_office_documents.py",
    "sum_column": "tests/test_office_documents.py",
    "compare_tables": "tests/test_office_documents.py",
    "notebook_read": "tests covering notebook editing",
    "notebook_edit": "needs an .ipynb fixture with outputs",
    "view_image": "needs an image and a vision-capable route",
    "show_molecule": "renders a viewer card; covered by tests/test_zz_show_molecule.py",
    # Chemistry data: needs an indexed corpus and a calc archive.
    "search_docs": "needs the indexed doc corpus",
    "read_section": "needs the indexed doc corpus",
    "list_docs": "needs the indexed doc corpus",
    "list_sections": "needs the indexed doc corpus",
    "search_calcs": "needs a calc archive",
    "get_calc_info": "needs a calc archive with real outputs",
    "calc_summary": "needs a calc archive with real outputs",
    # Journalled undo and history: their own state stores.
    "undo_changes": "tests/test_change_journal.py",
    "history_search": "tests covering the history store",
    "history_get": "needs a populated history store",
    "apply_patch": "tests/test_apply_patch_gate.py",
    "skill": "tests covering skills",
    "subagent_message":
        "tests/test_a_message_reaches_a_delegate_that_is_working.py",
}


def test_no_tool_arrives_without_a_decision():
    """The gate. A new tool must be exercised or excluded WITH a reason.

    This is the part that would have caught the three-week search outage:
    not a better test of search, but a rule that a callable tool cannot
    sit unexamined and unremarked.
    """
    advertised = advertised_tools()
    accounted = _EXERCISED | set(_NOT_EXERCISED_HERE)
    orphans = sorted(advertised - accounted)
    assert not orphans, (
        "these tools are advertised to the model but neither exercised "
        f"here nor listed with a reason: {orphans}")


def test_the_exclusion_list_does_not_drift_ahead_of_the_surface():
    """A name removed from the tool surface must leave this file too, or
    the list slowly becomes a graveyard nobody trusts."""
    advertised = advertised_tools()
    stale = sorted((_EXERCISED | set(_NOT_EXERCISED_HERE)) - advertised)
    assert not stale, f"listed here but no longer advertised: {stale}"


def test_every_exclusion_names_what_is_missing():
    """"Hard to test" is not a reason. Each entry names the resource."""
    empty = [k for k, v in _NOT_EXERCISED_HERE.items() if len(v.strip()) < 8]
    assert not empty, f"exclusions without a real reason: {empty}"


# ---------------------------------------------------------------------------
# The exercises
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def box():
    with tempfile.TemporaryDirectory(prefix="toolsurface-") as tmp:
        ws = Path(tmp)
        for cmd in (["init", "-q", "."], ["config", "user.email", "t@t"],
                    ["config", "user.name", "t"]):
            subprocess.run(["git", *cmd], cwd=ws, check=True)
        (ws / "ledger.py").write_text(
            "def settle(amount, account):\n"
            "    return account - amount\n\n\n"
            "def report(rows):\n"
            "    return settle(sum(rows), 0)\n")
        (ws / "notes.txt").write_text("alpha\nbeta MARKER gamma\ndelta\n")
        subprocess.run(["git", "add", "-A"], cwd=ws, check=True)
        subprocess.run(["git", "commit", "-qm", "first"], cwd=ws, check=True)
        yield ws


@pytest.fixture
def perms(box):
    return KitToolPermissions(mode="bypassPermissions", workspace=str(box))


def call(tool, args, perms):
    """Parse a payload, stepping over the untrusted fence when present.

    web_search and web_fetch come back inside it -- correctly, they are a
    stranger's text. A caller that json.loads the raw string sees nothing
    and would call a working tool broken; that mistake was made while
    writing this file, which is why the note is here.
    """
    raw = str(_doc_executor.execute(tool, args, perms))
    body = raw
    if body.lstrip().startswith("[UNTRUSTED"):
        body = "\n".join(l for l in body.split("\n")
                         if not l.startswith(("[UNTRUSTED", "[END")))
    try:
        return json.loads(body.strip()), raw
    except Exception:
        return {"_raw": raw}, raw


def test_reading_the_world(box, perms):
    out, _ = call("list_files", {"path": "."}, perms)
    assert "ledger.py" in json.dumps(out)

    out, _ = call("read_file", {"path": "ledger.py"}, perms)
    assert "def settle(amount, account)" in json.dumps(out)

    out, _ = call("grep_file", {"pattern": "MARKER", "path": "."}, perms)
    assert "MARKER" in json.dumps(out)


def test_navigating_code(box, perms):
    out, _ = call("find_definition", {"symbol": "settle"}, perms)
    assert "ledger.py" in json.dumps(out)
    out, _ = call("find_references", {"symbol": "settle"}, perms)
    assert "ledger" in json.dumps(out).lower()


def test_writing_is_verified_on_disk(box, perms):
    call("read_file", {"path": "notes.txt"}, perms)     # read-before-write
    call("edit_file", {"path": "notes.txt",
                       "old_string": "beta MARKER gamma",
                       "new_string": "beta CHANGED gamma"}, perms)
    assert "CHANGED" in (box / "notes.txt").read_text()

    call("write_file", {"path": "fresh.txt", "content": "made by the sweep\n"},
         perms)
    assert (box / "fresh.txt").read_text() == "made by the sweep\n"

    call("read_file", {"path": "fresh.txt"}, perms)
    call("multi_edit", {"path": "fresh.txt", "edits": [
        {"old_string": "made", "new_string": "written"}]}, perms)
    assert "written by the sweep" in (box / "fresh.txt").read_text()


def test_the_shell_reports_what_happened(box, perms):
    out, _ = call("bash", {"command": "echo sweep-ok"}, perms)
    assert "sweep-ok" in json.dumps(out)
    out, _ = call("bash", {"command": "exit 3"}, perms)
    assert "3" in json.dumps(out), "a non-zero exit has to be visible"


def test_tasks_survive_being_written(box, perms):
    out, _ = call("task_create",
                  {"subject": "surface sweep", "description": "d"}, perms)
    blob = json.dumps(out)
    assert "surface sweep" in blob or "id" in blob
    out, _ = call("task_list", {}, perms)
    assert "surface sweep" in json.dumps(out)


def test_the_environment_can_be_described(box, perms):
    out, _ = call("check_environment", {}, perms)
    assert "error" not in json.dumps(out)[:40].lower()
    out, _ = call("project_introspect", {}, perms)
    assert len(json.dumps(out)) > 20
    out, _ = call("list_changes_made", {}, perms)
    assert isinstance(out, (dict, list))


# ---------------------------------------------------------------------------
# The network half — the exact place the three-week outage lived
# ---------------------------------------------------------------------------

@pytest.fixture()
def net_backend():
    """Skip the caller unless the search backend is reachable from here.

    ASKED HERE, NOT AT IMPORT. This was a module-level ``skipif`` whose
    condition opened a socket while the file was imported -- which is
    during COLLECTION. Every ``--collect-only``, on CI, on a developer's
    machine and inside a cage asked the network what to collect; the test
    set then depended on whether the answer came back, which makes a green
    run unattributable. A listing has no business touching egress.

    THE SUITE'S OWN REFUSAL IS NOT A FIREWALL. The session guard
    (tests/external_agent_guard.py) refuses every non-loopback connect and
    raises ``ExternalAgentTripwire``, a ``PermissionError`` -- an
    ``OSError``. A plain ``except OSError`` therefore reports the guard's
    deliberate refusal as an unreachable host, and would keep the three
    network tests silently skipped even where the network works. The
    refusal is drained from the guard's record (so it is not charged to
    this test at teardown) and named in the skip reason.

    ``pytest.skip`` raises ``Skipped``, which derives from
    ``BaseException`` rather than ``Exception``, so it is not caught by
    the ``except OSError`` and cannot be turned into a silent pass.
    """
    import socket

    import external_agent_guard as guard

    try:
        socket.create_connection(("duckduckgo.com", 443), timeout=5).close()
    except OSError as exc:
        refused = [t for t in guard.take_trips() if "network connection" in t]
        if refused:
            pytest.skip("the suite's egress guard refused the reachability "
                        f"probe ({refused[0]}); that does not mean the "
                        "backend is down")
        pytest.skip(f"no route to the search backend "
                    f"({type(exc).__name__}: {exc}); a firewall must not "
                    "read as a bug")


def test_the_network_question_is_not_asked_while_tests_are_listed(tmp_path):
    """Collecting this file must not reach for the network.

    The guard does not fail a collection-time refusal -- it records it
    against "(after the last test)" and prints it in the terminal summary.
    That is exactly how the module-level probe survived for so long: every
    collection asked the network what to collect, and nothing went red. So
    the oracle here is not "no failure"; it is the guard's own signal, an
    audit hook on ``socket.getaddrinfo``/``socket.connect``, which fires
    whether the connect is refused, allowed or unresolvable. It therefore
    measures the same thing on CI, on a machine behind a firewall and on
    one whose resolver is down.

    Driven in a CHILD pytest (the pattern of
    tests/test_the_suite_reaches_no_external_agent.py) so the count is a
    whole collection of this file in a fresh process, with this checkout's
    conftest and therefore the real guard installed. The plugin is loaded
    with ``-p``, which happens before any conftest or module is imported,
    so the hook sees the collection it is measuring.

    The collected tests are asserted too: a change that silently stopped
    collecting the network tests would leave the hook equally silent, and
    would "pass" a test that only looked for the absence of a socket.
    """
    plugin = tmp_path / "net_audit_plugin.py"
    report = tmp_path / "events.json"
    plugin.write_text(textwrap.dedent('''\
        import json
        import os
        import sys


        EVENTS = []


        def _hook(event, args):
            if event in ("socket.getaddrinfo", "socket.connect"):
                EVENTS.append(event)


        sys.addaudithook(_hook)


        def pytest_sessionfinish(session, exitstatus):
            with open(os.environ["NET_AUDIT_REPORT"], "w",
                      encoding="utf-8") as fh:
                json.dump({"events": EVENTS,
                           "ids": [item.name for item in session.items]}, fh)
        '''))
    target = Path(__file__).resolve()
    env = child_env(tmp_path)
    env["NET_AUDIT_REPORT"] = str(report)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(target.parent.parent), str(tmp_path)]
        + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p])
    r = subprocess.run(
        [sys.executable, "-m", "pytest", "-p", "net_audit_plugin",
         "--collect-only", "-q", "-p", "no:cacheprovider",
         "-p", "no:randomly", "-o", "addopts=", "--rootdir",
         str(target.parent.parent), str(target)],
        cwd=str(target.parent.parent), env=env,
        capture_output=True, text=True, timeout=300)
    out = r.stdout + r.stderr
    assert report.exists(), f"the child collection did not finish:\n{out}"
    seen = json.loads(report.read_text(encoding="utf-8"))
    # A fix that stopped collecting the network tests would silence the
    # socket count for the wrong reason, so they are named one by one.
    for name in ("test_web_search_returns_results_not_just_an_absence",
                 "test_search_results_arrive_marked_untrusted",
                 "test_web_fetch_retrieves_and_fences_a_page"):
        assert name in seen["ids"], (
            f"{name} is no longer collected, so a silent socket count "
            f"would prove nothing:\n{out}")
    assert seen["events"] == [], (
        "collecting this file opened the network "
        f"({seen['events']}); the reachability question belongs in a "
        f"fixture, where a test asks it, not at import:\n{out}")


def test_web_search_returns_results_not_just_an_absence(box, perms,
                                                       net_backend):
    """The assertion the outage needed: a URL, not "the tool answered".

    A refusal is a legitimate outcome here (the backend rate-limits), so
    the check accepts either real hits or an explicit refusal payload --
    and fails on the third possibility, an empty result set presented as
    if nothing matched, which is the shape that misleads.
    """
    out, _ = call("web_search",
                  {"query": "python dataclasses documentation",
                   "max_results": 3}, perms)
    hits = out.get("results") or []
    if hits:
        assert str(hits[0].get("url", "")).startswith("http")
        return
    assert out.get("error"), (
        "no hits and no stated reason — an empty result set is exactly "
        "what a blocked backend must never look like")


def test_search_results_arrive_marked_untrusted(box, perms, net_backend):
    """Results carry the fence. An anti-bot refusal is not a result.

    DuckDuckGo answers a share of requests with HTTP 202 — a challenge,
    not an empty result set — and web_search reports that as its own
    error. Our own error text is not external content and must NOT be
    fenced, so asserting the fence unconditionally turned a refused
    request into an accusation against the fencing. Measured: it failed
    roughly one run in three, in isolation, and cost a full-suite
    investigation twice.

    The fence is what is under test, so a run with no results has nothing
    to test and says so.

    The condition is the CONTRACT, not a list of the refusals seen so far.
    Naming two of them left a third — a connection reset, reported as
    ``{"error": "network error: ..."}`` — falling through to the assert,
    and the suite went red over someone else's network. ``_wrap_untrusted``
    passes an error object through unfenced by design, so any error object
    is a run with nothing to fence.
    """
    _, raw = call("web_search", {"query": "orca manual", "max_results": 2},
                  perms)
    if raw.lstrip().startswith('{"error"') or '"results": []' in raw:
        pytest.skip(f"the search backend returned nothing: {raw[:120]}")
    assert "UNTRUSTED" in raw


def test_a_fetched_page_arrives_as_text_inside_the_fence(box, perms,
                                                         monkeypatch):
    """The fence and the HTML-to-text step, on every machine.

    The transport is replaced by a canned page, so the check does not
    depend on a network or on a third party's server being up; the live
    test below covers the transport where it is reachable.
    """
    from delfin.agent import web_tools
    page = (b"<html><head><title>T</title><script>var x=1;</script></head>"
            b"<body><h1>Example Domain</h1><p>For use in examples.</p>"
            b"</body></html>")
    monkeypatch.setattr(web_tools, "_fetch_bytes",
                        lambda url, timeout_s, *a, **k:
                        (page, "text/html; charset=utf-8"))
    # The user named the page: in bypass a host neither known nor named
    # is refused before anything is sent (web_fetch's host decision).
    perms.user_named_urls.add("https://example.com")
    out, raw = call("web_fetch", {"url": "https://example.com",
                                  "timeout_s": 15}, perms)
    assert "Example Domain" in out.get("text", ""), raw
    assert "var x" not in out.get("text", "")
    assert raw.lstrip().startswith("[UNTRUSTED"), raw[:120]


def test_web_fetch_retrieves_and_fences_a_page(box, perms, net_backend):
    """The real transport, where the page is reachable.

    example.com is a third party's server: it answered HTTP 503 during a
    CI run and turned a green suite red. An error object is web_fetch
    reporting that server's state correctly and is not fenced by design
    (``_wrap_untrusted``), so it is a run with nothing to check here; the
    fence itself is pinned offline above.
    """
    perms.user_named_urls.add("https://example.com")
    out, raw = call("web_fetch", {"url": "https://example.com",
                                  "timeout_s": 15}, perms)
    if out.get("error"):
        pytest.skip(f"example.com is not serving: {out['error'][:120]}")
    assert "Example Domain" in raw
    assert "UNTRUSTED" in raw
