"""Adversarial tests for package U5 (reviewer s21): DSML leak detection.

The builder adds ``text_sanitize.leaked_tool_call(text)`` returning ``(name,
arguments)`` when the text is a leaked tool call written as text, else None.
These tests pin the contract from the reviewer's side:

- real DSML shapes (from .gate/dsml_samples.txt) MUST be classified as calls,
  with the tool name and arguments recovered;
- prose, code or XML that merely MENTIONS the tags must NOT be classified as a
  call (no false positives), and must survive sanitise.

``leaked_tool_call`` does not exist yet: that import is wrapped so collection
succeeds and every test here reports red on the unfixed code (the control).
Each test goes green only when the builder's implementation satisfies the
contract. A red one that stays red after the fix is a finding.
"""
from __future__ import annotations

from delfin.agent.text_sanitize import parse_leaked_tool_calls, sanitize_agent_text

try:
    from delfin.agent.text_sanitize import leaked_tool_call
    _HAS_CALL = True
except ImportError:  # builder has not landed it yet -> every test below is red
    _HAS_CALL = False


def _call(text):
    """leaked_tool_call helper; explicit failure tells us the red reason."""
    assert _HAS_CALL, "leaked_tool_call not implemented yet (phase-2 control)"
    return leaked_tool_call(text)


# --- real DSML shapes from .gate/dsml_samples.txt ---------------------------

MARRIED_BASH = (
    "<invoke name=\"bash\">\n"
    "<parameter name=\"command\" string=\"true\">git status</parameter>\n"
    "<parameter name=\"description\" string=\"true\">Check repo state</parameter>\n"
    "</invoke>"
)

MARRIED_BASH_PLAIN = (
    "<invoke name=\"bash\">\n"
    "<parameter name=\"command\">git status</parameter>\n"
    "</invoke>"
)

MARRIED_WRITE_FILE = (
    "<invoke name=\"write_file\">\n"
    "<parameter name=\"path\">notes.md</parameter>\n"
    "<parameter name=\"content\">hello world</parameter>\n"
    "</invoke>"
)


def test_detects_married_bash_call():
    got = _call(MARRIED_BASH)
    assert got is not None
    name, args = got
    assert name == "bash"
    assert args.get("command") == "git status"


def test_detects_plain_without_string_attr():
    got = _call(MARRIED_BASH_PLAIN)
    assert got is not None
    name, args = got
    assert name == "bash"
    assert args.get("command") == "git status"


def test_recovers_all_parameters():
    got = _call(MARRIED_WRITE_FILE)
    assert got is not None
    name, args = got
    assert name == "write_file"
    assert args.get("path") == "notes.md"
    assert args.get("content") == "hello world"


def test_leaked_tool_call_recovers_married_dsml():
    # The detection seam is leaked_tool_call (phase 2); the phase-3 wiring
    # (engine.py, protected) calls it on the whole answer and re-requests.
    got = leaked_tool_call(MARRIED_BASH)
    assert got is not None
    assert got[0] == "bash"
    assert got[1]["command"] == "git status"
    assert got[1]["description"] == "Check repo state"


# --- false positives: prose / code that merely mentions the tags -----------

PROSE_MENTION = (
    "The DSML file uses <parameter> and <invoke> elements; see the schema. "
    "Do not confuse </invoke> with a real call."
)

CODE_SNIPPET = (
    "def render(tags):\n"
    "    return '<parameter name=\"x\">' + str(tags) + '</invoke>'\n"
)

XML_DOC = (
    "<schema>\n"
    "  <element name=\"invoke\"/>\n"
    "  <element name=\"parameter\"/>\n"
    "</schema>"
)

XSS_SHAPE = "<parameter name='command'>alert(1)</parameter>"


def test_prose_mention_is_not_a_call():
    assert _call(PROSE_MENTION) is None
    res = sanitize_agent_text(PROSE_MENTION)
    assert res.text == PROSE_MENTION


def test_code_snippet_is_not_a_call():
    assert _call(CODE_SNIPPET) is None
    res = sanitize_agent_text(CODE_SNIPPET)
    # the tag content survives sanitise (only indentation is normalised —
    # a pre-existing behaviour, not a DSML concern); no new tags introduced.
    assert "<parameter name=\"x\">" in res.text
    assert "</invoke>" in res.text


def test_xml_doc_is_not_a_call():
    assert _call(XML_DOC) is None
    res = sanitize_agent_text(XML_DOC)
    # the element tags survive (only indentation is normalised); the doc is
    # not misread as a call.
    assert "<schema>" in res.text
    assert "invoke" in res.text


def test_unclosed_invoke_is_not_a_call():
    assert _call("<invoke name=\"bash\">") is None


def test_escaped_entity_is_not_a_call():
    assert _call("&lt;invoke name=&quot;bash&quot;&gt;") is None


# --- edge cases ---------------------------------------------------------------

def test_empty_string_is_none():
    assert _call("") is None


def test_none_is_none():
    assert leaked_tool_call(None) is None


def test_unicode_params_preserved():
    got = _call("<invoke name=\"bash\">"
                "<parameter name=\"command\">ls „daten“</parameter>"
                "</invoke>")
    assert got is not None
    assert got[1]["command"] == "ls „daten“"


def test_nested_brackets_in_value():
    got = _call("<invoke name=\"bash\">"
                "<parameter name=\"command\">echo {a:{b:1}}</parameter>"
                "</invoke>")
    assert got is not None
    assert got[1]["command"] == "echo {a:{b:1}}"


def test_multiline_json_value():
    got = _call("<invoke name=\"write_file\">"
                "<parameter name=\"content\">{'a': 1,\n 'b': 2}</parameter>"
                "</invoke>")
    assert got is not None
    assert got[1]["content"] == "{'a': 1,\n 'b': 2}"


def test_leading_prose_still_detects_declared_call():
    # A model may narrate before leaking a call; that must not hide the call.
    got = _call("Let me check the repo state.\n" + MARRIED_BASH)
    assert got is not None
    assert got[0] == "bash"


def test_fragmented_DSMLparameter_form_is_not_a_call():
    # .gate/dsml_samples.txt includes a streamed/fragmented attribute form
    # (<DSMLparameter name="path … </parameter> truncated mid-attribute). This
    # is an artifact of a cut-off stream, not a complete leak; it must be None.
    frag = '<DSMLparameter name="path\n'
    assert _call(frag) is None


def test_trailing_prose_still_detects_declared_call():
    got = _call(MARRIED_BASH + "\nI'll run that now.")
    assert got is not None
    assert got[0] == "bash"


def test_single_parameter_with_string_attr_detected():
    got = _call('<invoke name="read_file">'
                '<parameter name="path" string="true">/tmp/a.txt</parameter>'
                "</invoke>")
    assert got is not None
    assert got[0] == "read_file"
    assert got[1]["path"] == "/tmp/a.txt"
