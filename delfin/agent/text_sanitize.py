"""Sanitise corrupted model output (harmony tool-channel leaks + glitch tokens).

Observed in production with azure.gpt-5.x served through the KIT
OpenAI-compatible endpoint: the model signals tool calls in its native
"harmony" channel syntax (``to=<tool> <json>``) which the endpoint passes
through as **text** instead of structured ``tool_calls``, and the channel's
special tokens decode into low-frequency multilingual "glitch" tokens
(e.g. ``手机天天中彩票``, ``출장샵``, ``ացին``).  The user sees Chinese / Korean /
Armenian garbage interleaved with a leaked ``to=search_docs {…}`` fragment.

This module repairs that text: it strips the leaked tool-channel fragments
(reporting which tools the model *intended* to call) and removes runs of
scripts that never legitimately appear in DELFIN's German/English chemistry
output.  It is a pure, dependency-free function — safe to run on every
response (a no-op on clean text).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field


# Leaked harmony tool-channel fragment: ``to=<tool> <junk> {json}``.
# The junk between the tool name and the JSON is the corrupted special
# tokens (``json_schema``, glitch chars).  DOTALL so the JSON can wrap.
_LEAKED_TOOL = re.compile(
    r"to=\s*([A-Za-z_][\w\-]*)[^\{\n]*?(\{.*?\})",
    re.DOTALL,
)

# Bare leftover ``to=<tool>`` with no JSON (defensive — strip the marker).
_BARE_TOOL_MARKER = re.compile(r"to=\s*[A-Za-z_][\w\-]*")

# Qwen/Gemma-family native tool-call markup. When the serving chat template
# (KIT vLLM, Ollama) lacks a tool parser, these models emit their calls as
# literal ``<tool_call>{json}</tool_call>`` blocks in the TEXT channel.
_QWEN_TOOL_CALL = re.compile(
    r"<tool_call>\s*(\{.*?\})\s*</tool_call>", re.DOTALL)

# GLM-family markup that wraps TEXT rather than a call object. Given a tool
# surface, GLM-5.3 on the KIT deployment wraps its ordinary answer in its own
# call tags and the serving parser leaves them in the text channel:
#
#   <tool_call>ACTION: /tab calc</arg_value></tool_call>
#
# The pattern above cannot see it — there is no JSON inside — so the answer
# survived with markup around it and the ACTION line no longer began a line.
# Measured 2026-09-07: the dashboard reported "[empty turn]" for a request
# the model had answered correctly, 0/3 across four benchmark arms, while the
# same prompt sent by hand returned "ACTION: /tab calc" every time. Only the
# tags go; the text between them is the answer.
_GLM_CALL_TAGS = re.compile(
    r"</?(?:tool_call|arg_key|arg_value)\s*>")

# A fenced ```json block whose entire payload is ONE call object. Only
# treated as a call when the object's keys are exactly a call shape
# (see _call_shape below) — ordinary JSON output must never be executed.
_FENCED_JSON_CALL = re.compile(r"```(?:json)?\s*(\{.*?\})\s*```", re.DOTALL)

# DSML invocation block — deepseek-v4-flash writes its tool calls as literal
# ``<invoke name="TOOL">…</invoke>`` XML in the text channel instead of
# calling the tool (wave 13, .gate/dsml_samples.txt: ``<invoke name="bash">``
# wrapping ``<parameter name="command" string="true">…</parameter>`` pairs).
# The discriminator is the PAIR: an opening ``<invoke name="…">`` with a tool
# name AND its matching ``</invoke>``.  Prose or code that merely mentions
# the tags (an XML-documentation snippet, a bare ``<parameter>``, an
# unmatched open from a stream cut short) is not a call — nothing ran, and
# ending a real turn on one would be worse than the leak it guards against.
_DSML_INVOKE = re.compile(
    r"<invoke\s+name\s*=\s*[\"']([A-Za-z_][\w\-]*)[\"']\s*>(.*?)</invoke>",
    re.DOTALL,
)
# A single argument, ``<parameter name="K" …>VALUE</parameter>``.  The ``[^>]*``
# after the name skips the DSML type serialization (``string="true"``) and
# any other attributes.
_DSML_PARAM = re.compile(
    r"<parameter\s+name\s*=\s*[\"']([\w.\-]+)[\"']\s*[^>]*>(.*?)</parameter>",
    re.DOTALL,
)

# Harmony special-token leftovers that decode as literal words.
#
# A leftover sits against what it was announcing: `json_schema{"doc_id"`
# is the recorded shape, and the special-token bar `constrain|>json{` is
# the other. The follow-set used to be "any character that is not a
# space or a letter", which cannot tell that from an ordinary word --
# and a file extension is an ordinary word: `settings.json` written in
# backticks, before a period or before a comma matched, so every JSON
# path an answer named arrived as `settings. `. Seen on 2026-09-08 in an
# answer that named the file eleven times and mangled it every time.
#
# Narrowed rather than widened on purpose. Missing a leftover leaves one
# stray word in prose; catching an extension hands the user a path that
# does not exist, with nothing in the answer to say why.
_HARMONY_TOKENS = re.compile(
    r"\b(?:json_schema|json|constrain)\b(?=\s*[{\[<|])")

# Reasoning-tag models (deepseek-r1, qwq, qwen3-thinking, …) emit their chain
# of thought as <think>…</think> in the visible text channel. Strip the whole
# block; also strip a dangling unterminated <think> (streaming can cut off
# before </think>) so no reasoning leaks into the user-visible answer.
_THINK_BLOCK = re.compile(r"<think>.*?</think>", re.DOTALL | re.IGNORECASE)
_THINK_DANGLING = re.compile(r"<think>.*$", re.DOTALL | re.IGNORECASE)

# Scripts that do not occur in DELFIN's de/en chemistry output — a run of
# these is glitch-token corruption, not content.  (CJK, kana, Hangul,
# Armenian, Cyrillic, Malayalam, Georgian, Thai, Devanagari, fullwidth.)
_GLITCH = re.compile(
    "["
    "　-〿぀-ヿ㐀-䶿一-鿿豈-﫿"
    "가-힯ᄀ-ᇿ"
    "Ѐ-ӿ"
    "԰-֏"
    "ഀ-ൿ"
    "Ⴀ-ჿ"
    "฀-๿"
    "ऀ-ॿ"
    "＀-￯"
    "]+"
)


@dataclass
class SanitizeResult:
    text: str
    leaked_tools: list[str] = field(default_factory=list)
    glitch_chars: int = 0
    think_stripped: bool = False
    # How much text went in. Only interesting next to an empty ``text``:
    # an answer that was entirely think-blocks or tool-call markup leaves
    # nothing behind, and a caller that sees only the empty string cannot
    # tell that apart from a backend that said nothing at all. The two
    # have different causes and different remedies, so the caller has to
    # be able to tell them apart.
    source_chars: int = 0

    @property
    def changed(self) -> bool:
        return (bool(self.leaked_tools) or self.glitch_chars > 0
                or self.think_stripped)

    @property
    def emptied(self) -> bool:
        """The model produced text and none of it survived cleaning."""
        return self.source_chars > 0 and not self.text


def strip_glitch(text: str) -> str:
    """``text`` with glitch-token runs removed and nothing else touched.

    For tool arguments a person reads, where sanitize_agent_text's tool-call
    and think-block repairs have no business. Report 20260915-132613 put
    "Nein, erst отчетen" on a dialog button.
    """
    if not isinstance(text, str) or not _GLITCH.search(text):
        return text
    return re.sub(r"[ \t]{2,}", " ", _GLITCH.sub("", text)).strip()


def sanitize_agent_text(text: str) -> SanitizeResult:
    """Return a cleaned copy of ``text`` plus what was repaired.

    - Strips leaked ``to=<tool> {json}`` fragments (the model's intended
      tool calls that leaked into the text channel) and reports the tool
      names in ``leaked_tools``.
    - Removes runs of non-DELFIN scripts (glitch tokens), counting how many
      characters were dropped in ``glitch_chars``.

    Clean text passes through unchanged (``changed`` is False).
    """
    if not text:
        return SanitizeResult(text=text or "", source_chars=0)

    leaked: list[str] = []
    for m in _LEAKED_TOOL.finditer(text):
        name = m.group(1)
        if name not in leaked:
            leaked.append(name)
    for m in _QWEN_TOOL_CALL.finditer(text):
        try:
            shape = _call_shape(json.loads(m.group(1)))
        except Exception:
            shape = None
        if shape and shape[0] not in leaked:
            leaked.append(shape[0])

    think_stripped = bool(_THINK_BLOCK.search(text) or _THINK_DANGLING.search(text))
    cleaned = _THINK_BLOCK.sub(" ", text)
    cleaned = _THINK_DANGLING.sub(" ", cleaned)
    cleaned = _LEAKED_TOOL.sub(" ", cleaned)
    cleaned = _QWEN_TOOL_CALL.sub(" ", cleaned)
    # After the JSON form above, so a real leaked call is still recognised
    # as a call and removed whole rather than unwrapped into prose. The tag
    # becomes a NEWLINE, not nothing: two wrapped ACTION lines are written
    # back to back, and the dashboard's parser reads one action per line —
    # joining them would trade an empty turn for a mangled one.
    cleaned = _GLM_CALL_TAGS.sub("\n", cleaned)
    cleaned = _BARE_TOOL_MARKER.sub(" ", cleaned)
    cleaned = _HARMONY_TOKENS.sub(" ", cleaned)

    glitch_chars = sum(len(m.group(0)) for m in _GLITCH.finditer(cleaned))
    cleaned = _GLITCH.sub("", cleaned)

    # Collapse the whitespace/newlines left behind by the removals.
    cleaned = re.sub(r"[ \t]{2,}", " ", cleaned)
    cleaned = re.sub(r"\n{3,}", "\n\n", cleaned)
    cleaned = cleaned.strip()

    return SanitizeResult(
        text=cleaned,
        leaked_tools=leaked,
        glitch_chars=glitch_chars,
        think_stripped=think_stripped,
        source_chars=len(text),
    )


def _call_shape(obj) -> tuple[str, dict] | None:
    """Return (name, arguments) when ``obj`` is unambiguously a tool-call
    object: a dict whose keys are exactly a name plus an arguments dict
    (``arguments`` or ``parameters``). Anything looser must be rejected so
    ordinary JSON the model *prints* is never mistaken for a call."""
    if not isinstance(obj, dict):
        return None
    keys = set(obj.keys())
    for args_key in ("arguments", "parameters"):
        if keys == {"name", args_key}:
            name = obj.get("name")
            args = obj.get(args_key)
            if isinstance(name, str) and name and isinstance(args, dict):
                return name, args
    return None


def parse_leaked_tool_calls(text: str) -> list[dict]:
    """Best-effort recovery of the JSON args from leaked tool fragments.

    Recognised grammars, in order:
    - harmony ``to=<tool> {json}`` leaks (gpt-oss / gpt-5 family)
    - ``<tool_call>{json}</tool_call>`` markup (qwen/gemma family on
      serving stacks whose chat template lacks a tool parser)
    - fenced ```json blocks that ARE one strict call object

    Returns ``[{"name": str, "arguments": dict}]`` for each parseable
    fragment. The caller decides whether to re-dispatch; only cleanly
    parsed dict arguments are ever returned for the strict grammars.
    """
    out: list[dict] = []
    for m in _LEAKED_TOOL.finditer(text):
        name = m.group(1)
        try:
            args = json.loads(m.group(2))
        except Exception:
            args = {}
        out.append({"name": name, "arguments": args})
    for m in _QWEN_TOOL_CALL.finditer(text):
        try:
            shape = _call_shape(json.loads(m.group(1)))
        except Exception:
            continue
        if shape:
            out.append({"name": shape[0], "arguments": shape[1]})
    if not out:
        for m in _FENCED_JSON_CALL.finditer(text):
            try:
                shape = _call_shape(json.loads(m.group(1)))
            except Exception:
                continue
            if shape:
                out.append({"name": shape[0], "arguments": shape[1]})
    return out


def leaked_tool_call(text: str) -> tuple[str, dict] | None:
    """Return ``(tool_name, arguments)`` when ``text`` is a leaked DSML tool
    call, else ``None``.

    deepseek-v4-flash sometimes emits its tool calls as literal
    ``<invoke name="…">…</invoke>`` markup in the text channel instead of
    calling the tool, so the call never runs and the turn stands still.
    Only a COMPLETE block — an opening ``<invoke name="TOOL">`` AND its
    matching ``</invoke>`` — is treated as a leak, so prose or code that
    merely mentions the tags (an XML snippet, a bare ``<parameter>``, an
    unmatched open) is never mistaken for one.  Arguments are the
    ``<parameter name="K">VALUE</parameter>`` pairs inside the block.

    Returns a single ``(name, args)`` for the first complete block: the
    caller's remedy is to ask the model once to call the tool, which only
    makes sense for a turn dominated by one leaked call.
    """
    if not isinstance(text, str):
        return None
    m = _DSML_INVOKE.search(text)
    if not m:
        return None
    name = m.group(1)
    args: dict = {}
    for p in _DSML_PARAM.finditer(m.group(2)):
        args[p.group(1)] = p.group(2)
    return name, args


# Leftover (text outside the matched <invoke>…) threshold for DOMINATION.
# The real leaked answers (wave 13, .gate/dsml_samples.txt) are the bare
# <invoke>…</invoke> block and nothing else — empty leftover.  Prose that
# merely CITIES a complete block ("For example: <invoke …>…</invoke> which
# runs a search.") leaves a real sentence outside, well over this bound, so
def leaked_tool_call_dominates(text: str) -> tuple[str, dict] | None:
    """``(name, args)`` when *text* IS one leaked DSML call, else ``None``.

    The re-request remedy only makes sense for a turn that is NOTHING but the
    leaked call: the model wrote the tool call as text and nothing else, so
    nothing ran and the turn should ask again.  Prose that quotes a full block
    is a real answer and must never be cut short or re-requested — and a
    one-word aside ("ok <block>", "<block> done", "Maybe: <block>") is still
    prose.  No length bound can tell a bare call from a short wrapper, so the
    rule is strict: after removing every complete block, any non-whitespace
    text outside the block makes the answer prose (``None``), never the call.
    """
    if not isinstance(text, str):
        return None
    m = _DSML_INVOKE.search(text)
    if not m:
        return None
    leftover = _DSML_INVOKE.sub("", text).strip()
    if leftover:
        return None
    args: dict = {}
    for p in _DSML_PARAM.finditer(m.group(2)):
        args[p.group(1)] = p.group(2)
    return m.group(1), args
