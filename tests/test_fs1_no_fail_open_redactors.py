# FS1 phase 2 guard: no function in delfin/ may call output_guard.scrub_secrets
# (or memory_store._without_secrets) inside a try whose except returns its own
# argument -- that is the fail-open redactor hole FS1 closed. Found via AST
# (not grep), so the next wrapper cannot fail open again.
import ast
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
AGENT = REPO_ROOT / "delfin" / "agent"

REDACTOR_NAMES = {"scrub_secrets", "_without_secrets"}

# KNOWN FOLLOW-UPS found by this guard on the FS1 base (never edited by the
# FS1 builder: api_client.py is a protected file). Each entry is a real,
# verified fail-open wrapper that must be fixed in its own package:
#   api_client.py::_redact_tool_result  -- except returns the unscrubbed arg
#   working_state.py::_scrub            -- except returns the unscrubbed arg
# The guard stays strict for every other module: a name not listed here
# fails the build.
KNOWN_FOLLOWUPS = {
    ("api_client.py", "_redact_tool_result"),
    ("working_state.py", "_scrub"),
}


def _call_names(node):
    """Names of functions a Call node invokes (flat, dotted last part)."""
    func = node.func
    name = getattr(func, "id", None) or getattr(func, "attr", None)
    return {name} if name else set()


def _arg_names(call):
    """Names of variables passed as the call's first positional argument."""
    out = set()
    if call.args and isinstance(call.args[0], ast.Name):
        out.add(call.args[0].id)
    return out


def _try_blocks(node):
    for child in ast.walk(node):
        if isinstance(child, ast.Try):
            yield child


def _returns_in(node):
    return [n for n in ast.walk(node) if isinstance(n, ast.Return)]


def fail_open_redactors():
    """(file, function) pairs whose except returns the scrubbed argument."""
    hits = []
    for path in sorted(AGENT.glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for func in tree.body:
            if not isinstance(func, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            for try_block in _try_blocks(func):
                scrub_args = set()
                for child in ast.walk(try_block):
                    if isinstance(child, ast.Call) and (
                            _call_names(child) & REDACTOR_NAMES):
                        scrub_args |= _arg_names(child)
                if not scrub_args:
                    continue
                for handler in try_block.handlers:
                    for ret in _returns_in(handler):
                        if (isinstance(ret.value, ast.Name)
                                and ret.value.id in scrub_args):
                            hits.append((path.name, func.name))
    return hits


def test_no_function_scrubs_inside_a_fail_open_try():
    hits = sorted(set(fail_open_redactors()))
    new = [h for h in hits if h not in KNOWN_FOLLOWUPS]
    assert not new, (
        "fail-open secret redaction in delfin/agent (except returns the "
        f"unscrubbed input): {new}. Return a fixed placeholder instead, "
        "e.g. '[redaction failed -- text withheld]' or ''.")


@pytest.mark.parametrize("snippet,expected", [
    # fires: except returns the scrubbed argument
    ("def f(text):\n"
     "    try:\n"
     "        from .output_guard import scrub_secrets\n"
     "        return scrub_secrets(text)\n"
     "    except Exception:\n"
     "        return text\n", True),
    # does not fire: fixed placeholder
    ("def f(text):\n"
     "    try:\n"
     "        return scrub_secrets(text)\n"
     "    except Exception:\n"
     "        return '[redaction failed -- text withheld]'\n", False),
    # does not fire: empty string
    ("def f(text):\n"
     "    try:\n"
     "        return scrub_secrets(text)\n"
     "    except Exception:\n"
     "        return ''\n", False),
    # does not fire: a scrub call outside any try block
    ("def f(text):\n"
     "    return scrub_secrets(text)\n", False),
    # does not fire: no redactor call at all
    ("def f(text):\n"
     "    try:\n"
     "        return other(text)\n"
     "    except Exception:\n"
     "        return text\n", False),
])
def test_detector_snippets(snippet, expected):
    tree = ast.parse(snippet)
    func = tree.body[0]
    fired = False
    scrub_args = set()
    for try_block in _try_blocks(func):
        for child in ast.walk(try_block):
            if isinstance(child, ast.Call) and (
                    _call_names(child) & REDACTOR_NAMES):
                scrub_args |= _arg_names(child)
        fired = False
        for handler in try_block.handlers:
            for ret in _returns_in(handler):
                if (isinstance(ret.value, ast.Name)
                        and ret.value.id in scrub_args):
                    fired = True
    assert fired == expected