"""The prompt teaches the form in inline code, and the parser rejected it.

``dashboard_agent.md`` lists the five accepted ACTION forms and prints
every one of them in backticks, because that is how a markdown document
shows a command. A model writes the shape back the way it was shown —
`ACTION: /tab calc` — and the parser, anchored at the start of the line,
found a backtick there and returned nothing. The action was dropped in
silence, in the dashboard as well as in the benchmark.

Measured 2026-09-08 on kit.deepseek-v4-flash: workflow_verify_after_modify
emitted both required actions, inside backticks, and scored 28 with the
signals recorded as missing.
"""

import pytest

from delfin.agent.benchmark_runner import extract_actions


@pytest.mark.parametrize("line, expected", [
    ("ACTION: /tab calc", ["/tab calc"]),
    ("`ACTION: /tab calc`", ["/tab calc"]),
    ("**ACTION: /tab calc**", ["/tab calc"]),
    ("*ACTION: /tab calc*", ["/tab calc"]),
    ("`ACTION:/tab calc`", ["/tab calc"]),
    ("`Action /tab calc`", ["/tab calc"]),
    ("```\nACTION: /tab calc\n```", ["/tab calc"]),
])
def test_the_wrappers_a_model_writes_are_unwrapped(line, expected):
    assert extract_actions(line) == expected


def test_a_bare_slash_in_code_is_still_a_command():
    assert extract_actions("`/tab calc`") == ["/tab calc"]


def test_a_path_in_prose_is_still_not_a_command():
    """The whitelist guard is what keeps /home/user/x.py out; unwrapping
    must not reach past it."""
    assert extract_actions("Siehe `/home/user/x.py` im Ordner.") == []
    assert extract_actions("Die Datei `/etc/hosts` ist gemeint.") == []


def test_two_wrapped_actions_are_both_found():
    text = ("Ich setze das Funktional:\n\n"
            "`ACTION: /orca set method B3LYP`\n\n"
            "Und zeige danach den Zustand:\n\n"
            "`ACTION: /orca show`\n")
    assert extract_actions(text) == ["/orca set method B3LYP", "/orca show"]


def test_the_dashboard_parser_unwraps_the_same_way():
    """The benchmark mirrors the dashboard; a fix in one that is not in
    the other means the benchmark stops measuring the product."""
    import inspect

    from delfin.dashboard import tab_agent

    src = inspect.getsource(tab_agent.create_tab)
    assert "_ACTION_WRAPPERS" in src
    assert 'strip(_ACTION_WRAPPERS)' in src


@pytest.mark.parametrize("line", [
    "ACTION: /tab calc", "ACTION:/tab calc", "ACTION /tab calc",
    "Action: /tab calc", "Action /tab calc", "action: /tab calc",
    "ACTION. /tab calc",
])
def test_the_mirror_accepts_every_form_the_dashboard_does(line):
    """The dashboard matches one IGNORECASE pattern; the benchmark spelled
    out three case-sensitive ones and missed "Action /tab calc". A mirror
    stricter than the product under-reports the product."""
    assert extract_actions(line) == ["/tab calc"], line
