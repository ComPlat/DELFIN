# FS1 phase 2 guard: no function in delfin/ may let the unscrubbed text
# through when output_guard.scrub_secrets (or a sibling redactor it calls)
# raises -- that is the fail-open redactor hole FS1 closed. Found via AST
# (not grep), so the next wrapper cannot fail open again.
#
# Guard scope (review findings 3-5, nacht-s6): the scan is RECURSIVE over
# delfin/agent/ (service/ included), it walks class METHODS as well as
# module-level functions, and it catches the fail-open idiom in every
# ordinary spelling:
#   * a handler that RETURNS the guarded input, directly, wrapped
#     (str(text), text.strip(), f"{text}") or via an alias (out = text);
#   * a handler that SWALLOWS the exception (no return, no raise) when the
#     guarded name is still returned by the function on the success path --
#     the unscrubbed value escapes later. A swallow that only misses a
#     DETECTION (the guarded value is never returned, e.g. swarm_memory's
#     visible "secret scan unavailable" reason) is not a leak and stays
#     allowed.
import ast
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
AGENT = REPO_ROOT / "delfin" / "agent"

REDACTOR_NAMES = {"scrub_secrets", "_without_secrets", "_redact_secrets",
                  "_shapes"}

# KNOWN FOLLOW-UPS: every entry a real, verified fail-open redaction site
# outside FS1's owned files (output_guard.py is read-only for the builder;
# the others are not owned by any FS1 file). Each was confirmed by reading
# the code, not by the scan alone:
#   api_client.py::_redact_tool_result -- except returns the unscrubbed arg
#   working_state.py::_scrub           -- except returns the unscrubbed arg
#   output_guard.py::scrub_secrets     -- :231-232 replace-loop swallow with
#       `text` returned; :236-238 except returns the raw input
#   output_guard.py::run_output_guards -- :309-310 resets `guarded = text`
#       (the raw, unscrubbed input) and returns it
#   session_store.py::_walk            -- :626-627 except returns the value
#       unredacted into the crash-recovery checkpoint; the scan also
#       attributes the nested closure's hit to the enclosing
#       _redact_checkpoint, so both names are listed
#   session_store.py::_redact_checkpoint -- the nested _walk closure hit
#   turn_metrics.py::_clean_error      -- :85-86 except returns the
#       unscrubbed error string, which is written to disk
#   redact.py::text                    -- service redactor: :54-55 except
#       pass around the shape pass, :43-44 an unreadable credential store
#       empties the held-value list; fix delivered as the FS1 operator
#       patch .gate/fs1_service_redact_fail_closed.patch -- remove this
#       entry once the operator has committed it.
KNOWN_FOLLOWUPS = {
    ("api_client.py", "_redact_tool_result"),
    ("working_state.py", "_scrub"),
    ("output_guard.py", "scrub_secrets"),
    ("output_guard.py", "run_output_guards"),
    ("session_store.py", "_walk"),
    ("session_store.py", "_redact_checkpoint"),
    ("turn_metrics.py", "_clean_error"),
    ("redact.py", "text"),
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


def _has_raise(node):
    return any(isinstance(n, ast.Raise) for n in ast.walk(node))


def _names_in(node):
    if node is None:
        return set()
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}


def _guarded_names(func, try_block):
    """Names holding the redactor's input inside this try, plus aliases
    assigned from them anywhere in the enclosing function (an alias may
    be bound before the try: out = text ... except: return out)."""
    scrub = set()
    for child in ast.walk(func):
        if isinstance(child, ast.Call) and (
                _call_names(child) & REDACTOR_NAMES):
            scrub |= _arg_names(child)
    if not scrub:
        return set()
    # Restrict to what flows through THIS try: the call must sit inside it
    # (or the try must sit inside the function's own guarded path).
    try_end = max((n.lineno for n in ast.walk(try_block)
                   if hasattr(n, "lineno")), default=try_block.lineno)
    call_inside = any(
        isinstance(child, ast.Call) and (_call_names(child) & REDACTOR_NAMES)
        for child in ast.walk(try_block))
    call_after = any(
        isinstance(child, ast.Call) and (_call_names(child) & REDACTOR_NAMES)
        and child.lineno > try_end
        for child in ast.walk(func))
    if not call_inside:
        # Only the explicit-return branch is sound here. The swallow branch
        # stays off: without a scrub inside the try we cannot tell a safe
        # re-assignment from an escaping raw value without dataflow.
        if not call_after:
            return set()
        returns_guarded = any(
            isinstance(ret, ast.Return)
            and (_names_in(ret.value) & scrub)
            for handler in try_block.handlers
            for ret in ast.walk(handler))
        if not returns_guarded:
            return set()
    aliases = set()
    for child in ast.walk(func):
        if isinstance(child, ast.Assign) and isinstance(
                child.value, ast.Name) and child.value.id in scrub:
            for target in child.targets:
                if isinstance(target, ast.Name):
                    aliases.add(target.id)
    return scrub | aliases


def _returned_later(func, try_block, guarded):
    """True when any return statement of ``func`` outside this try's
    handlers hands back one of the guarded names."""
    handler_nodes = {id(h) for tb in _try_blocks(func)
                     for h in tb.handlers}
    for ret in ast.walk(func):
        if not isinstance(ret, ast.Return) or id(ret) in handler_nodes:
            continue
        if _names_in(ret.value) & guarded:
            return True
    return False


def _handler_fails_open(handler, guarded, func, try_block):
    rets = _returns_in(handler)
    for ret in rets:
        if _names_in(ret.value) & guarded:
            return True
    if not rets and not _has_raise(handler):
        # Swallowed: a leak only when the scrubbed value escapes later.
        if _returned_later(func, try_block, guarded):
            return True
    return False


def fail_open_redactors():
    """(file, function) pairs whose except lets the unscrubbed input through."""
    hits = []
    for path in sorted(AGENT.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for func in ast.walk(tree):
            if not isinstance(func, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            for try_block in _try_blocks(func):
                guarded = _guarded_names(func, try_block)
                if not guarded:
                    continue
                for handler in try_block.handlers:
                    if _handler_fails_open(handler, guarded, func,
                                           try_block):
                        hits.append((path.name, func.name))
    return hits


def test_no_function_scrubs_inside_a_fail_open_try():
    hits = sorted(set(fail_open_redactors()))
    new = [h for h in hits if h not in KNOWN_FOLLOWUPS]
    assert not new, (
        "fail-open secret redaction in delfin/agent (except lets the "
        f"unscrubbed input through): {new}. Return a fixed placeholder "
        "instead, e.g. '[redaction failed -- text withheld]' or ''.")


def test_the_guard_scans_the_service_subpackage_files():
    """Finding 3: delfin/agent/service/ holds 20+ files including a
    fail-open redactor; a glob('*.py') scan never saw them."""
    scanned = {path.name for path in AGENT.rglob("*.py")}
    assert "redact.py" in scanned, (
        "the guard must scan delfin/agent/service/, not delfin/agent/*.py")


@pytest.mark.parametrize("snippet,expected", [
    # fires: except returns the scrubbed argument
    ("def f(text):\n"
     "    try:\n"
     "        from .output_guard import scrub_secrets\n"
     "        return scrub_secrets(text)\n"
     "    except Exception:\n"
     "        return text\n", True),
    # fires: except returns the argument wrapped (str(), strip(), f-string)
    ("def f(text):\n"
     "    try:\n"
     "        return scrub_secrets(text)\n"
     "    except Exception:\n"
     "        return str(text)\n", True),
    ("def f(text):\n"
     "    try:\n"
     "        return scrub_secrets(text)\n"
     "    except Exception:\n"
     "        return text.strip()\n", True),
    ("def f(text):\n"
     "    try:\n"
     "        return scrub_secrets(text)\n"
     "    except Exception:\n"
     "        return f'{text}'\n", True),
    # fires: except returns via an alias
    ("def f(text):\n"
     "    out = text\n"
     "    try:\n"
     "        return scrub_secrets(text)\n"
     "    except Exception:\n"
     "        return out\n", True),
    # fires: swallowed exception, the scrubbed value escapes on the
    # success path afterwards (the service redactor's shape pass)
    ("def f(text):\n"
     "    try:\n"
     "        text = scrub_secrets(text)\n"
     "    except Exception:\n"
     "        pass\n"
     "    return text\n", True),
    # fires inside a class METHOD (the Redactor.text shape)
    ("class R:\n"
     "    def text(self, text):\n"
     "        try:\n"
     "            text = self._shapes(text, [])\n"
     "        except Exception:\n"
     "            pass\n"
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
    # does not fire: re-raised
    ("def f(text):\n"
     "    try:\n"
     "        return scrub_secrets(text)\n"
     "    except Exception:\n"
     "        raise\n", False),
    # does not fire: a scrub call outside any try block
    ("def f(text):\n"
     "    return scrub_secrets(text)\n", False),
    # does not fire: no redactor call at all
    ("def f(text):\n"
     "    try:\n"
     "        return other(text)\n"
     "    except Exception:\n"
     "        return text\n", False),
    # does not fire: the redactor runs BEFORE the try, so the handler
    # returns an already-scrubbed value (probe12 false positive, s6)
    ("def f(text):\n"
     "    out = scrub_secrets(text)\n"
     "    try:\n"
     "        record(out)\n"
     "    except Exception:\n"
     "        return out\n", False),
    # does not fire: swallow only misses a DETECTION -- the guarded value
    # is compared, never returned (swarm_memory / skill_safety shape)
    ("def f(text):\n"
     "    try:\n"
     "        changed = scrub_secrets(text) != text\n"
     "    except Exception:\n"
     "        pass\n"
     "    return [] if not changed else ['secret']\n", False),
    # does not fire: swallow whose failure is made visible to the caller
    # (the guarded value never reaches a return)
    ("def f(text):\n"
     "    try:\n"
     "        if scrub_secrets(text) != text:\n"
     "            raise ValueError('secret')\n"
     "    except Exception:\n"
     "        raise\n"
     "    return text\n", False),
])
def test_detector_snippets(snippet, expected):
    tree = ast.parse(snippet)
    fired = bool(fail_open_redactors_tree(tree))
    assert fired == expected, snippet


def fail_open_redactors_tree(tree):
    """The scanner over an in-memory tree (snippet tests)."""
    hits = []
    for func in ast.walk(tree):
        if not isinstance(func, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for try_block in _try_blocks(func):
            guarded = _guarded_names(func, try_block)
            if not guarded:
                continue
            for handler in try_block.handlers:
                if _handler_fails_open(handler, guarded, func, try_block):
                    hits.append(func.name)
    return hits
