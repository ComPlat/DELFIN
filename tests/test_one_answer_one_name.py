"""The quote-aware redirect pass needs one name, not two.

`refusal_memory._redirect_targets_safe` was introduced (160a5a92) so the
cd/pushd/popd guard could be applied on top of `_redirect_targets`, and was
kept afterwards (86a992624) on the stated ground that the operator-built
hook in `api_client.py` calls it by name. That ground is false and was
false when written: `api_client.py:24681` imports `_redirect_targets`
literally. The guard itself had already moved INSIDE `_redirect_targets`,
so the alias was a second name for one behaviour.

What is claimed here, and how it is measured:

  * the behaviour the alias existed for still happens through the REAL
    route a refused shell write takes -- `_DocToolExecutor._run_permission_gate`
    recording the ledger and `_gate_bash_write_targets` refusing a later
    spelling of the same write (this is what the alias was kept for);
  * nothing in the tracked tree reaches the alias, so its removal changes
    no refusal. That is stated as a measurement over the repository rather
    than as an import that happens to work.

These tests FAIL on a tree that still offers the alias only because the
alias-presence test is written against its absence, and they fail if the
behaviour ever stops being reachable without it. The behaviour tests
themselves pass before and after the deletion -- deliberately: the claim is
"the alias was not carrying this", and a behaviour test that went green
only on the change would prove the opposite.
"""

from __future__ import annotations

import pathlib
import subprocess

import pytest

from delfin.agent import api_client as A
from delfin.agent import refusal_memory as R

_REPO = pathlib.Path(__file__).resolve().parent.parent


def _perms(tmp_path, asked):
    """A permission set whose dialog refuses and records what it was asked."""
    def _deny(name, args, preview):
        asked.append((name, dict(args or {})))
        return False
    return A.KitToolPermissions(workspace=tmp_path, mode="default",
                                confirm_callback=_deny)


def test_the_quoted_path_is_recorded_and_later_refused_without_the_alias(
        tmp_path):
    """The record path and the check path, through the real executor.

    A quoted redirect/destination word with a space is what the quote-aware
    pass exists for: the gate's own tokenizer splits it at the space. After
    refusal, a DIFFERENT spelling of the same write must be refused from the
    ledger rather than reaching the dialog a second time.
    """
    assert not hasattr(R, "_redirect_targets_safe"), (
        "the alias is back; this test measures the single-name route")
    ex = A._DocToolExecutor()
    asked: list = []
    perms = _perms(tmp_path, asked)

    f = tmp_path / "engine with space.py"
    cmd1 = f'cp src.txt "{f}" && git checkout abc123 -- engine.py'
    verdict1 = ex._run_permission_gate("bash", {"command": cmd1}, perms)
    assert verdict1 is not None, "setup: command should be denied"
    assert len(asked) == 1
    assert perms.action_denied_earlier("write", str(f.resolve())), (
        f"the whole quoted path was not recorded; ledger keys: "
        f"{sorted(perms.denied_actions)}")

    cmd2 = f'printf x > "{f}"'
    verdict2 = ex._run_permission_gate("bash", {"command": cmd2}, perms)
    assert verdict2 is not None and "already refused a write" in verdict2, (
        "the same write in another spelling reached the dialog again; "
        f"ledger keys: {sorted(perms.denied_actions)}")
    assert len(asked) == 1


def test_the_directory_change_guard_still_refuses_to_resolve_redirect_words():
    """The guard the alias was introduced for, on the surviving function.

    A command that changes directory resolves a relative redirect word
    against a base this module cannot know, so the quote-aware pass takes
    nothing from it. Taking nothing only ever removes candidates from the
    ADDING direction, so it cannot cause a refusal that was not there.
    """
    assert R._redirect_targets('cat x > "my file.txt"') == ["my file.txt"]
    assert R._redirect_targets("cd sub && ./tool > f.txt") == []
    assert R._redirect_targets("pushd sub && ./tool > f.txt") == []
    assert R._command_changes_directory("cd sub && ./tool > f.txt") is True


def _production_mentions(symbol: str) -> list[str]:
    """Whole-word mentions of *symbol* in the shipped package, ``delfin/``.

    The population is production code on purpose. A scan over the whole tree
    cannot answer "does anything call this" for a name a test names: this
    file talks about both deleted names, so once it is committed it reports
    itself. That is not a corner case -- it is the difference between the
    measurement and the thing measured, and an earlier version of this file
    passed only because it was still untracked when it was run.

    ``-w`` is required, not cosmetic: `expected_digest` is a prefix of the
    surviving `expected_digests` parameter and of `EXPECTED_DIGEST`, and a
    substring search reports those as callers of the deleted accessor. The
    positive control is in `test_the_scanner_finds_a_name_that_is_there`.
    """
    out = subprocess.run(
        ["git", "-c", "grep.threads=1", "grep", "-nw", symbol,
         "--", "delfin/"],
        cwd=_REPO, capture_output=True, text=True, check=False,
    )
    if out.returncode not in (0, 1):  # 1 = no match, which is the answer
        pytest.fail(f"git grep failed ({out.returncode}): {out.stderr[:400]}")
    lines = [ln for ln in out.stdout.splitlines() if ln.strip()]
    return [ln for ln in lines if _is_a_reference(ln.split(":", 1)[0], symbol)]


def _is_a_reference(rel: str, symbol: str) -> bool:
    """Whether *symbol* appears in *rel* as CODE, not only as prose.

    The grep above answers "do these letters occur", and a comment that
    says what a name used to be occurs. Measured: after #834 deleted
    ``_redirect_targets_safe``, the only mention left in ``delfin/`` was
    ``refusal_memory.py``'s own note -- "an alias, _redirect_targets_safe,
    used to sit above" -- and this scan read that sentence as a caller. The
    same false alarm appeared in #834's own CI, from a docstring naming the
    alias it had just removed.

    So the rule, and it is deliberately not "ignore every string":

    * an identifier (``Name``/``Attribute``) IS a reference;
    * a string literal IS one too -- ``getattr(mod, "name")`` dispatches
      by name, and the reachability ratchet counts those on purpose;
    * a module, class or function DOCSTRING is not, and neither is a
      comment, because neither can dispatch anything.

    A file that will not parse counts as a reference: the question is
    whether anything still calls the name, and an unparseable file is not
    an answer.
    """
    import ast

    path = _REPO / rel
    try:
        tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"))
    except (OSError, SyntaxError):
        return True
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)):
            body = getattr(node, "body", None) or []
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                docstrings.add(id(body[0].value))
    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and node.id == symbol:
            return True
        if isinstance(node, ast.Attribute) and node.attr == symbol:
            return True
        if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                and id(node) not in docstrings and symbol in node.value):
            return True
    return False


def test_the_scanner_finds_a_name_that_is_there():
    """A negative result from a scanner nobody checked is not evidence.

    The deleted accessors are reported as "no mention"; that only means
    something if the same call finds a name that is present, and if it still
    tells a deleted name from the surviving one it prefixes.
    """
    hits = _production_mentions("expected_digests")
    assert hits, "the scanner found nothing for a name that must be there"
    assert not _production_mentions("expected_digest"), (
        "the scanner does not distinguish expected_digest from "
        "expected_digests, so its 'no caller' result means nothing")


@pytest.mark.parametrize("symbol", ["_redirect_targets_safe",
                                    "expected_digest"])
def test_the_removed_predicate_has_no_caller(symbol):
    """No production code reaches the name -- measured, not assumed.

    A hit here means the shipped package still refers to the deleted name,
    so removing it removed something reachable.
    """
    hits = _production_mentions(symbol)
    assert hits == [], (
        f"{symbol} is still referenced in delfin/, so removing it removes "
        f"something reachable:\n" + "\n".join(hits[:12]))


def test_a_reference_is_told_from_prose_about_one(tmp_path, monkeypatch):
    """The filter above, in both directions, on files written for it.

    Without this, "no caller" could mean the filter swallows callers. The
    four cases are the ones that actually occur: a call, a dynamic
    dispatch by name, a docstring that explains what a name used to be,
    and a comment that does the same.
    """
    import test_one_answer_one_name as mod

    cases = {
        "a_call.py": "from x import gone\ngone()\n",
        "an_attribute.py": "import x\nx.gone()\n",
        "a_dispatch.py": 'import x\ngetattr(x, "gone")()\n',
        "a_docstring.py": '"""An alias, gone, used to sit above."""\n',
        "a_comment.py": "# An alias, gone, used to sit above.\npass\n",
        "unparseable.py": "def (\n",
    }
    for name, text in cases.items():
        (tmp_path / name).write_text(text, encoding="utf-8")
    monkeypatch.setattr(mod, "_REPO", tmp_path)

    assert mod._is_a_reference("a_call.py", "gone")
    assert mod._is_a_reference("an_attribute.py", "gone")
    assert mod._is_a_reference("a_dispatch.py", "gone"), (
        "a name dispatched by string is a caller; the reachability ratchet "
        "counts those on purpose")
    assert not mod._is_a_reference("a_docstring.py", "gone")
    assert not mod._is_a_reference("a_comment.py", "gone")
    assert mod._is_a_reference("unparseable.py", "gone"), (
        "a file that will not parse is not an answer")
