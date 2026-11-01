"""SW1 adversarial tests, phase 1 (reviewer nacht-s12).

Attacks against delfin/agent/swarm/spec.py glob-overlap and protection
logic, formed BEFORE reading the implementation (G1-G8 in
.gate/HANDOFF.md) plus package-name and chemistry-core cases.

Green cases pin correct behaviour; red cases are findings reported to
the builder with file:line.
"""

from __future__ import annotations

import pytest

from delfin.agent.swarm.spec import SpecError, validate_spec


def _spec(*packages: dict, max_sessions: int = 8) -> dict:
    return {"goal": "x", "packages": list(packages),
            "max_sessions": max_sessions}


def _pkg(name: str, *owned: str) -> dict:
    return {"name": name, "owned": list(owned), "task": "t"}


def test_star_crossing_slash_globs_must_conflict():
    # fnmatch '*' crosses '/', so 'delfin/dashboard/*.py' covers
    # 'delfin/dashboard/sw/spec.py' at write time: two packages with
    # these scopes race. Neither glob touches a protected path, so the
    # failure must be the OVERLAP check, not the protected check.
    with pytest.raises(SpecError, match="overlap"):
        validate_spec(_spec(
            _pkg("a", "delfin/dashboard/*.py"),
            _pkg("b", "delfin/dashboard/sw/**"),
        ))


def test_broad_glob_covering_a_protected_file_is_refused():
    # The protected check must work on the EXPANDED match, not on a
    # string comparison against the glob text.
    for broad in ("delfin/agent/*", "delfin/**", "**"):
        with pytest.raises(SpecError, match="protected"):
            validate_spec(_spec(_pkg("a", broad)))


def test_duplicate_package_names_are_an_error():
    with pytest.raises(SpecError):
        validate_spec(_spec(
            _pkg("core", "delfin/agent/swarm/**"),
            _pkg("core", "delfin/agent/dashboard/**"),
        ))


def test_unicode_lookalike_path_is_not_a_false_overlap():
    # Cyrillic 'а' in 'agent' is a different path: neither package
    # conflicts with the other and the ASCII one stays valid.
    spec = _spec(
        _pkg("ascii", "delfin/agent/swarm/**"),
        _pkg("cyr", "delfin/\u0430gent/swarm/**"),
    )
    validate_spec(spec)


def test_case_difference_is_not_a_false_overlap():
    # Case-sensitive repo: 'Delfin/agent/**' is a different tree from
    # 'delfin/agent/**' (documented choice). Both are valid here: the
    # upper one is far enough from every protected file not to be
    # refused for that reason instead.
    validate_spec(_spec(
        _pkg("lower", "delfin/agent/swarm/**"),
        _pkg("upper", "Delfin/agent/**"),
    ))


def test_overlap_within_one_package_is_legal():
    # Only BETWEEN packages is an overlap an error; one package naming
    # two shapes of the same tree is one writer, no race.
    validate_spec(_spec(
        _pkg("a", "delfin/agent/", "delfin/agent/swarm/**"),
    ))


def test_empty_or_none_owned_list_is_an_error():
    with pytest.raises(SpecError):
        validate_spec(_spec({"name": "a", "owned": [], "task": "t"}))
    with pytest.raises(SpecError):
        validate_spec(_spec({"name": "a", "task": "t"}))


def test_max_sessions_counts_across_packages_with_roles():
    # Two packages with 2 declared roles each = 4 members: a
    # max_sessions of 3 must refuse, not count per package.
    spec = _spec(
        {**_pkg("a", "delfin/a/**"), "roles": {"builder": True,
                                              "reviewer": True}},
        {**_pkg("b", "delfin/b/**"), "roles": {"builder": True,
                                              "reviewer": True}},
        max_sessions=3,
    )
    with pytest.raises(SpecError, match="max_sessions"):
        validate_spec(spec)


def _roles(builder: bool = True, reviewer: bool = True) -> dict:
    return {"builder": builder, "reviewer": reviewer}


def test_max_sessions_across_packages_with_roles():
    spec = _spec(
        {**_pkg("a", "delfin/a/**"), "roles": _roles()},
        {**_pkg("b", "delfin/b/**"), "roles": _roles()},
        max_sessions=4,
    )
    validate_spec(spec)
    spec["max_sessions"] = 3
    with pytest.raises(SpecError, match="max_sessions"):
        validate_spec(spec)


def test_chemistry_core_globs_are_refused():
    # The chemistry core is read-only for every swarm member: a spec
    # that claims it as a write scope is refused like the safety frame.
    for owned in ("delfin/api.py", "delfin/smiles_converter.py",
                  "delfin/smiles_converter/**", "delfin/manta/**"):
        with pytest.raises(SpecError, match="protected"):
            validate_spec(_spec(_pkg("a", owned)))
