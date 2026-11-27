"""The two places the CO2 coordinator is meant to grow, marked so that an
extension cannot land half-wired.

The group extends the CO2 coordinator next, so the workflow has to fail a
future change that plugs a piece in at one end and leaves the other end
unaware of it. There are exactly two such joints, and each gets a test that
a real extension must keep satisfied:

JOINT 1 — a new CONTROL option (or a new ``job``) the coordinator reads.
    ``CO2_Coordinator6.main()`` reads its own CONTROL.txt, and
    ``setup_co2_from_delfin`` writes it from the main CONTROL.txt. A key in
    the read-set but not in the transfer set is an option a user can write
    and see nothing happen — measured before this file existed:
    ``perform_orientation_scan`` was in exactly that hole, and because the
    coordinator treats an ABSENT key as "run the scan", turning the scan off
    from the main CONTROL silently ran it anyway.

JOINT 2 — a new step adapter for the coordinator.
    ``delfin-step co2_coordinator_opi --job=sp|scan`` is the product entry
    point to the OPI input writers. An extension that adds a third job must
    also make it discoverable in the contract and keep the existing two
    reachable; a new job that is accepted but writes nothing, or one that
    silently takes over the default, is the defect this test refuses.

Neither test names a fixed set of keys or jobs — both derive the set from the
code — so a legitimate extension passes without editing them, and an
incomplete one fails with the name of what it forgot.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
COORDINATOR = REPO / "delfin" / "co2" / "CO2_Coordinator6.py"


def _coordinator_read_keys() -> set[str]:
    """Every CONTROL key ``CO2_Coordinator6.main()`` reads, from its source.

    Derived by AST rather than from a hand-kept list, so the test follows the
    code instead of a copy of it that drifts.
    """
    tree = ast.parse(COORDINATOR.read_text(encoding="utf-8"))
    main = next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "main")
    keys: set[str] = set()
    for node in ast.walk(main):
        # args.get("key", ...) and args["key"]
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "get"
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "args"
                and node.args
                and isinstance(node.args[0], ast.Constant)
                and isinstance(node.args[0].value, str)):
            keys.add(node.args[0].value)
        if (isinstance(node, ast.Subscript)
                and isinstance(node.value, ast.Name)
                and node.value.id == "args"
                and isinstance(node.slice, ast.Constant)
                and isinstance(node.slice.value, str)):
            keys.add(node.slice.value)
    return keys


# ---------------------------------------------------------------------------
# JOINT 1 — a CONTROL option the coordinator reads must be reachable
# ---------------------------------------------------------------------------


def test_every_control_key_the_coordinator_reads_can_be_set_from_the_main_control():
    """No option may be readable by the coordinator and untransferable.

    Red on the base commit: ``perform_orientation_scan``, ``orca_keywords``
    and ``rot_orca_keywords`` were read by ``main()`` and absent from
    ``_TRANSFER_KEYS``, so a user could write them in the main CONTROL.txt and
    the hand-over dropped them without a word.

    A future extension that adds an option to the coordinator inherits this
    gate automatically: it fails here until the option is transferable.
    """
    from delfin.co2 import chain_setup

    read = _coordinator_read_keys()
    transferable = set(chain_setup._TRANSFER_KEYS) | {
        # Set by the hand-over itself from the OCCUPIER result, never
        # copied from the main CONTROL (see _TRANSFER_KEYS' comment).
        "charge", "multiplicity", "broken_sym",
    }
    unreachable = sorted(read - transferable)

    assert not unreachable, (
        "these keys are read by CO2_Coordinator6.main() but cannot be passed "
        f"on from the main CONTROL.txt: {unreachable}. Add them to "
        "delfin/co2/chain_setup.py's _TRANSFER_KEYS (and to a section in "
        "_build_co2_control if they should be written) or the option is "
        "silently ignored.")


def test_a_transferable_option_is_written_into_a_section():
    """The transfer set alone is not enough — the writer must emit the key.

    ``_build_co2_control`` writes only the keys its section list names, so a
    key added to ``_TRANSFER_KEYS`` and forgotten here is still dropped. That
    is exactly how the first attempt at this fix failed.
    """
    from delfin.co2 import chain_setup
    import inspect

    source = inspect.getsource(chain_setup._build_co2_control)
    missing = sorted(
        key for key in chain_setup._TRANSFER_KEYS
        if f'"{key}"' not in source
    )
    assert not missing, (
        f"these transferable keys are never written by _build_co2_control, so "
        f"the hand-over drops them: {missing}")


def test_an_option_with_no_shipped_default_is_not_invented_for_the_user():
    """An option with no shipped default must not be given one.

    The rule that keeps every existing complete CONTROL reading as it did on
    origin/main: a key with no shipped default is written ONLY when the user
    set it, so the coordinator's own default applies otherwise. Inventing a
    value (or writing an empty ``key=`` line) would change what an untouched
    file means.

    A control case: it holds on the base commit and must keep holding after
    any extension. The base commit satisfies it by having no such keys at all,
    which is why the two tests above — not this one — are the red ones.
    """
    from delfin.co2 import chain_setup

    optional = getattr(chain_setup, "_OPTIONAL_KEYS", frozenset())
    defaults = set(chain_setup._CO2_DEFAULTS)

    # A key cannot be both: a shipped default is always emitted, an optional
    # key only when the user wrote it.
    assert optional.isdisjoint(defaults), (
        "a key cannot be both a shipped default (always emitted) and optional "
        "(emitted only when set): "
        f"{sorted(optional & defaults)}")
    for key in sorted(optional):
        assert key in chain_setup._TRANSFER_KEYS, (
            f"{key} is declared optional but not transferable, so it could "
            "never be read from the main CONTROL at all")


# ---------------------------------------------------------------------------
# JOINT 2 — a new job for the OPI step adapter
# ---------------------------------------------------------------------------


def _adapter_class():
    from delfin.tools.adapters.co2_coordinator_opi import Co2CoordinatorOpiAdapter

    return Co2CoordinatorOpiAdapter


def _declared_jobs() -> tuple[str, ...]:
    """The jobs the adapter's contract advertises."""
    from delfin.tools._spec import ParamSpec  # noqa: F401  (documents the type)

    for spec in _adapter_class().params:
        if spec.name == "job":
            return tuple(spec.enum or ())
    return ()


def test_the_step_adapter_is_reachable_in_the_registry():
    """The adapter is registered and discoverable, not merely importable.

    A module that registers on import is only reachable if the registry walk
    actually loaded it; this is the "an import is not a caller" check.
    """
    from delfin.tools._registry import list_steps

    steps = list_steps()
    assert "co2_coordinator_opi" in steps, (
        "the CO2 input-builder step is not in the registry, so no pipeline "
        f"can run it. Registered: {sorted(steps)[:10]}...")
    assert steps["co2_coordinator_opi"].params, (
        "the registered step declares no params, so a caller cannot discover "
        "the job choice")


def test_every_declared_job_is_actually_offered(capsys):
    """Each job the contract advertises must be named by ``--describe``.

    This is the extension gate: adding a job to the enum without making it
    real is caught here, and adding it without advertising it is caught by
    the next test.
    """
    from delfin.cli_step import main as step_main

    assert step_main(["co2_coordinator_opi", "--describe"]) == 0
    out = capsys.readouterr().out
    jobs = _declared_jobs()
    assert jobs, "the step advertises no job choices at all"
    for job in jobs:
        assert f"'{job}'" in out, (
            f"job {job!r} is declared in the adapter's enum but the contract "
            f"does not show it, so a caller cannot find it:\n{out}")


def _accepted_jobs() -> frozenset[str]:
    """The jobs the adapter's own gate in ``execute`` lets through.

    Read from the AST of the refusal — ``if job not in (...):`` — rather than
    from a copy of the list, and rather than by string-matching the source.
    A string match was the first attempt at this and silently found only
    ``sp`` of ``("sp", "scan")``, which made the comparison below pass no
    matter what the gate said.
    """
    import inspect

    from delfin.tools.adapters import co2_coordinator_opi as module

    tree = ast.parse(
        inspect.getsource(module.Co2CoordinatorOpiAdapter.execute).lstrip())
    for node in ast.walk(tree):
        if (isinstance(node, ast.Compare)
                and isinstance(node.left, ast.Name) and node.left.id == "job"
                and len(node.ops) == 1 and isinstance(node.ops[0], ast.NotIn)):
            comparator = node.comparators[0]
            values = getattr(comparator, "elts", None) or getattr(
                comparator, "values", None)
            if values is None:
                raise AssertionError(
                    "the job gate no longer compares against a literal tuple "
                    "or set, so this test cannot derive the accepted jobs")
            return frozenset(
                v.value for v in values
                if isinstance(v, ast.Constant) and isinstance(v.value, str))
    raise AssertionError(
        "no ``if job not in (...)`` gate found in "
        "Co2CoordinatorOpiAdapter.execute, so an unknown job is not refused")


def test_the_jobs_the_adapter_accepts_and_advertises_are_the_same_set():
    """The accepted set and the advertised set must be the SAME set.

    The adapter refuses an unknown job with its own hard gate inside
    ``execute`` and advertises the choice in ``params``. A job accepted but
    not advertised is undiscoverable; one advertised but refused is a lie in
    ``--describe``. Both are extension bugs, and this test must fail on
    EITHER — a one-directional subset check passes the first of them, which
    is how the first version of this test let a half-wired ``orient`` job
    through.
    """
    accepted = _accepted_jobs()
    advertised = set(_declared_jobs())

    assert advertised, "the adapter advertises no jobs"
    assert accepted, "the adapter's gate accepts no jobs at all"
    # Equality, not a subset with the expected answers spelled out: an
    # extension that adds a job to BOTH passes without editing this test, and
    # one that adds it to only one fails naming which side is missing it.
    assert accepted == advertised, (
        "the jobs the adapter ACCEPTS and the jobs it ADVERTISES must be the "
        f"same set. Accepted only: {sorted(accepted - advertised)} (callers "
        f"cannot discover these). Advertised only: {sorted(advertised - accepted)} "
        "(``--describe`` offers a job the adapter then refuses).")
