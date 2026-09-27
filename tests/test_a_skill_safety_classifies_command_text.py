"""A skill proposal's commands go through DELFIN's own classifiers.

Paket 3 of the skill-learning wave: ``skill_safety.check(text)`` reads a
draft SKILL.md and reports every finding. The core idea (per the wave
design): whatever the LIVE bash gate would ask about or refuse must be a
finding in a skill text — an unsafe success must not become a rule
("Skill Misevolution", arXiv 2608.12851).

This file covers phase 1: extracting commands from the text and judging
them with DELFIN's own functions (imported from api_client, never
copied, never changed).
"""
from __future__ import annotations



from delfin.agent import skill_safety


def test_the_module_exists_and_exposes_check():
    assert callable(skill_safety.check)


def test_a_clean_command_in_a_bash_block_is_clean(tmp_path):
    text = (
        "# ORCA input check\n\n"
        "Run the lint before committing:\n\n"
        "```bash\n"
        "pytest -q tests/test_x.py\n"
        "```\n"
    )
    assert skill_safety.check(text) == []


def test_a_deny_listed_command_in_a_bash_block_is_a_finding():
    text = (
        "```bash\n"
        "rm -rf /tmp/scratch\n"
        "```\n"
    )
    findings = skill_safety.check(text)
    assert any("rm -rf" in f or "deny" in f.lower() for f in findings), findings


def test_a_command_the_gate_would_ask_about_is_a_finding():
    # curl is neither on the deny list nor auto-allowed: live it asks.
    text = (
        "```bash\n"
        "curl -o data.json https://example.invalid/data\n"
        "```\n"
    )
    findings = skill_safety.check(text)
    assert findings, "a command that would ask in live operation must be a finding"


def test_dollar_lines_are_commands(tmp_path):
    text = "Install the tool:\n\n$ pip install --user numpy\n"
    findings = skill_safety.check(text)
    assert findings


def test_inline_code_spans_are_commands():
    text = "Clean up with `chmod 000 /etc/hosts` when it locks."
    findings = skill_safety.check(text)
    assert findings


def test_a_command_wrapped_across_lines_is_one_command():
    # A line-wrapped delete: the classifier must see it whole.
    text = (
        "```bash\n"
        "rm -rf \\\n"
        "  build/\n"
        "```\n"
    )
    findings = skill_safety.check(text)
    assert findings


def test_bash_c_payload_is_checked_too():
    text = "```bash\nbash -c 'rm -rf build'\n```"
    findings = skill_safety.check(text)
    assert findings


def test_plain_prose_without_commands_is_clean():
    text = (
        "# Reading ORCA output\n\n"
        "Use the extract_energy_table tool for energies. Compare only\n"
        "within one method; a total energy across functionals is not\n"
        "comparable. Report imaginary frequencies even when the headline\n"
        "number looks fine.\n"
    )
    assert skill_safety.check(text) == []


def test_findings_are_named_strings():
    text = "```bash\nsudo apt-get install jq\n```"
    findings = skill_safety.check(text)
    assert findings and all(isinstance(f, str) and f for f in findings)
