"""``stop --help`` must not promise a removal that no longer happens.

``stop --cleanup`` still said, in its ``--help`` text,
"remove scratch/OCCUPIER artifacts". That was true before this wave and is
false after it: ``*_OCCUPIER`` stage directories are resume state, kept
deliberately so a resumed stage finds its ``OCCUPIER.txt``, its per-stage
``CONTROL.txt`` and its restart ``.gbw``
(``delfin/workflows/engine/occupier.py:646``,
``delfin/copy_helpers.py:153``).

The report dict was guarded by
``tests/test_stop_cleanup_keeps_what_a_resume_reads.py:203``, but the help
text is a separate surface and was missed. An operator who reads ``--help``
concludes the resume state is gone and re-runs what did not need re-running.

Driven through ``cli.main``: the help comes out of the real argument parser,
not out of a string constant read directly.

The help text is WRAPPED by argparse -- "--cleanup" sits on one line and its
description continues on the next -- so the option's text is joined before it
is read. A selector that grabbed only the first line would return "" and every
"not in" assertion on it would pass without measuring anything.
"""
from __future__ import annotations

import contextlib
import io

from delfin.cli import main


def _help_text(argv: list[str]) -> str:
    """Run ``main`` and return everything it printed."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        try:
            main(argv)
        except SystemExit:
            pass
    return buf.getvalue()


def _option_help(help_text: str, option: str) -> str:
    """Return the help text of one option, joining its continuation lines.

    ``--cleanup`` is defined at one indentation and its description wraps at a
    deeper one; the block ends at the next option or a blank line.
    """
    lines = help_text.splitlines()
    start = next(
        (n for n, line in enumerate(lines) if line.strip().startswith(option + " ")
         or line.strip() == option),
        None,
    )
    assert start is not None, f"{option} is not in the help output at all"

    block = [lines[start]]
    for line in lines[start + 1:]:
        if not line.strip():
            break
        # A new option starts with two spaces and a dash at the option column.
        if line.lstrip().startswith("--"):
            break
        block.append(line)
    return " ".join(part.strip() for part in block)


def test_the_stop_cleanup_help_does_not_claim_stage_dirs_are_removed():
    """The help must not say the OCCUPIER stage dirs are removed.

    Fails on a tree whose help still says "remove scratch/OCCUPIER
    artifacts", and fails loudly rather than silently if the option text
    cannot be located at all.
    """
    text = _option_help(_help_text(["stop", "--help"]), "--cleanup")
    lowered = text.lower()

    # Guard: the joined block must be substantive. Without this, a failed
    # selector yields "" and the assertions below measure nothing.
    assert len(lowered) > 40, f"could not read the --cleanup help text; got {text!r}"
    assert "artifacts" in lowered, f"unexpected --cleanup help text: {text!r}"

    assert "scratch/occupier" not in lowered, (
        "--help still names scratch/OCCUPIER as something cleanup removes; "
        f"got: {text!r}"
    )
    assert "resume" in lowered or "keep" in lowered, (
        "the --cleanup help should say the OCCUPIER stage dirs are kept as "
        f"resume state; got: {text!r}"
    )


def test_the_cleanup_option_is_reachable_in_the_help():
    """Control: the option is parsed and displayed on both trees.

    Passes before and after the fix; it proves the option still reaches the
    help output, so the test above cannot pass because ``--cleanup`` vanished.
    """
    out = _help_text(["stop", "--help"])

    assert "--cleanup" in out
    assert "artifacts" in out
