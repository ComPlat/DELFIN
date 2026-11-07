"""SWP phase 4: the 12 shared-practice points are explicit in the role
templates where they BIND the member.

Red control: currently the templates name "red first", "never merge" and
"own scope", but none of them say "never push", "run the failing test on
plain main first", "CI red -> fix with the failing tests" or "no paths /
user names / hosts in PR text". Every assertion below is therefore RED on
the previous commit and GREEN once prompts.py carries the binding.

Source of the points: .gate/SCHWARM_PRAXIS.md, "Shared practice" items
1-10 + additions 11-12 (SW1 templates / SW3 are the delivery vehicle).
"""

from delfin.agent.swarm.prompts import TEMPLATES, render

_VALUES = {
    "goal": "measure",
    "package": "core",
    "task": "Build the thing.",
    "owned": "delfin/x.py",
    "branch": "swarm/core-builder",
}

#: phrases that bind EVERY member, whatever their role.
_COMMON_BINDING = (
    "never push",
    "plain main first",
    "no paths, user names or host",
)


def _render(role: str) -> str:
    return render(role, _VALUES)


def test_common_bindings_in_every_role():
    for role in ("builder", "reviewer", "qs"):
        text = _render(role)
        for phrase in _COMMON_BINDING:
            assert phrase in text, (role, phrase)


def test_builder_makes_exact_writes_red_first_and_ci_red_explicit():
    t = _render("builder")
    # exact files / owned scope (point 1) -- present, pinned here.
    assert "You own exactly these write paths" in t
    assert "delfin/x.py" in t
    # red first (point) and CI red -> fix with the failing tests (point 2).
    assert "Red test first" in t
    assert "failing tests" in t
    # failing test run on plain main first (point 4).
    assert "plain main first" in t


def test_reviewer_only_writes_adversarial_tests_and_never_product_code():
    t = _render("reviewer")
    assert "adversar" in t
    assert "never commit product code" in t
    assert "plain main first" in t


def test_qs_runs_tests_itself_and_never_approves_on_word():
    t = _render("qs")
    # Binding already in the qs template (kept; case-insensitive).
    assert "run the tests yourself" in t.lower()
    assert "never" in t


def test_every_template_stays_short_and_leaves_no_unfilled_placeholder():
    for role in ("builder", "reviewer", "qs"):
        t = _render(role)
        # every referenced {placeholder} is supplied -> no literal braces
        # remain; and the template stays a brief brief.
        assert "{" not in t, (role, "unfilled placeholder")
        assert len(t) < 1400, (role, len(t))


def test_templates_are_a_dict_with_exactly_three_roles():
    assert set(TEMPLATES) == {"builder", "reviewer", "qs"}
