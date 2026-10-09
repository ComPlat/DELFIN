"""Phase 2 of U5: detect a DSML tool call written as text.

deepseek-v4-flash sometimes emits its tool calls as literal DSL markup
(``<invoke name=\"bash\">…</parameter>…</invoke>``) in the text channel instead
of calling the tool.  The call never runs and the turn stands still.  Unlike
the Harmony path (which re-dispatches the parseable JSON), the remedy here is
to notice that the WHOLE answer is a leaked call and ask the model once to
actually call the tool.

The detector must be strict: only a COMPLETE invocation block — an opening
``<invoke name=\"TOOL\">`` AND its matching ``</invoke>`` — is a leak.  Prose
or code that merely *mentions* the tags (an XML-documentation snippet, a bare
``<parameter>`` without an enclosing invoke, an unmatched open) is not, and
must yield None so no real turn is ever falsely cut short.
"""

from __future__ import annotations

from delfin.agent.text_sanitize import leaked_tool_call


# The reconstructed whole-call shapes from the wave-13 logs (the fragments
# in the stream, joined, are exactly these).  ``string=\"true\"`` is the DSML
# type serialization; the s21 log shows the same call without it.
_BASH_CALL = (
    '<invoke name="bash">'
    '<parameter name="command" string="true">git log --oneline -3</parameter>'
    '<parameter name="description" string="true">Read the last commits</parameter>'
    '</invoke>'
)

_WRITE_CALL = (
    '<invoke name="write_file">'
    '<parameter name="path">.gate/commitmsg.txt</parameter>'
    '<parameter name="content">Fix the sort order</parameter>'
    '</invoke>'
)

_PLAIN_S21_CALL = (
    '<invoke name="bash">'
    '<parameter name="command">git status</parameter>'
    '</invoke>'
)


def test_returns_tool_name_and_arguments():
    got = leaked_tool_call(_BASH_CALL)
    assert got is not None
    name, args = got
    assert name == "bash"
    assert args["command"] == "git log --oneline -3"
    assert args["description"] == "Read the last commits"


def test_plain_variant_without_type_attribute():
    got = leaked_tool_call(_PLAIN_S21_CALL)
    assert got is not None
    name, args = got
    assert name == "bash"
    assert args == {"command": "git status"}


def test_multiple_parameters_and_value_with_punctuation():
    got = leaked_tool_call(_WRITE_CALL)
    assert got is not None
    name, args = got
    assert name == "write_file"
    assert args == {
        "path": ".gate/commitmsg.txt",
        "content": "Fix the sort order",
    }


def test_arguments_survive_angle_brackets_and_unusual_chars():
    # A command that legitimately contains < > and & inside its value.
    text = (
        '<invoke name="bash">'
        '<parameter name="command">grep -n "<invoke" delfin/agent/x.py</parameter>'
        '</invoke>'
    )
    got = leaked_tool_call(text)
    assert got is not None
    assert got[1]["command"] == 'grep -n "<invoke" delfin/agent/x.py'


def test_tool_name_with_underscore():
    text = (
        '<invoke name="read_file">'
        '<parameter name="path">delfin/agent/text_sanitize.py</parameter>'
        '</invoke>'
    )
    got = leaked_tool_call(text)
    assert got is not None
    assert got[0] == "read_file"


# --- No false positives: prose / code that merely mentions the tags ---

def test_prose_mentioning_invoke_without_close_is_not_a_leak():
    text = (
        "Check whether <invoke name=\"bash\"> appears; if it does, "
        "the model wrote its call as text."
    )
    assert leaked_tool_call(text) is None


def test_unmatched_open_invoke_is_not_a_leak():
    # Streaming cut off mid-call: an opening tag but no closing </invoke>.
    text = '<invoke name="read_file">\n<parameter name="path">delfin/agent'
    assert leaked_tool_call(text) is None


def test_bare_parameter_without_invoke_wrapper_is_not_a_leak():
    # Documentation of the low-level tag alone — no invoke, no tool name.
    text = (
        "Each argument is serialised as <parameter name=\"path\">VALUE</parameter> "
        "inside an <invoke> block."
    )
    assert leaked_tool_call(text) is None


def test_closing_tag_alone_is_not_a_leak():
    assert leaked_tool_call("The fragment ends with </invoke> but came alone.") is None


def test_clean_prose_with_word_invoke_is_not_a_leak():
    text = (
        "Noch bevor ich das Tool invoke, prüfe ich, ob eine Rechnung läuft. "
        "Die Antwort ist 42."
    )
    assert leaked_tool_call(text) is None


def test_xml_documentation_snippet_of_a_full_call_is_a_leak():
    # A complete <invoke name="x">…</invoke> pair IS a leaked call even when
    # it looks like documentation — real leaked calls and the documentation
    # of them are indistinguishable, and the whole-answer check below decides.
    doc = (
        "For example: <invoke name=\"search_docs\">"
        "<parameter name=\"query\">CASSCF</parameter></invoke> which runs a search."
    )
    got = leaked_tool_call(doc)
    assert got is not None
    assert got[0] == "search_docs"


# --- Property-style checks, both ways ---

_TAG_PIECES = [
    "<invoke", "</invoke>", "<parameter name=\"x\">y</parameter>",
    "<DSMLparameter name=\"path", "<invoke name=\"bash\">",
]

def test_property_no_tag_fragment_alone_triggers():
    # None of these fragments is a complete closed <invoke name=…>…</invoke>
    # block, so none may trigger — however they appear in a sentence.
    surrounds = [
        "Erkläre {}",
        "Schau dir {} an, bitte.",
        "The tag {} is never allowed here.",
    ]
    for piece in _TAG_PIECES:
        for s in surrounds:
            assert leaked_tool_call(s.format(piece)) is None, piece


def test_property_complete_blocks_always_trigger():
    tools = ["bash", "read_file", "search_docs", "grep_file"]
    for tool in tools:
        text = (
            '<invoke name="{t}"><parameter name="arg">1</parameter></invoke>'
        ).format(t=tool)
        got = leaked_tool_call(text)
        assert got is not None
        assert got[0] == tool
