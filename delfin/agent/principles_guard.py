"""Guard for the maintainer's principles addendum.

``pack/shared/principles_addendum.md`` is the first shared contract in
every role prompt (prompt_loader.py, layer 0). This module pins its
digest so tampering becomes visible: a file that is missing, an empty
scaffold, or edited text fails :func:`check` with a clear reason.

The digest is pinned in TWO independent places — the constant below and
a copy the maintainer keeps at a protected location. ``check`` accepts
an optional list of additional expected digests so callers can verify
against both.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

# SHA-256 of the normalised principles text (LF line endings, exactly
# one trailing newline). If the maintainer edits the file, this constant
# moves with it — tests/test_principles_guard.py names it in the failure
# message.
EXPECTED_DIGEST = (
    "d1d488ddbc0317ce31ffe3c67425e4641f87b0181782b4b8939b574e11d8f25f"
)

_PRINCIPLES_REL = Path("shared") / "principles_addendum.md"


@dataclass(frozen=True)
class GuardResult:
    ok: bool
    reason: str


def _normalise(text: str) -> str:
    """Line endings unified to LF, trailing blank lines collapsed to one
    final newline — so a checkout on CRLF or an editor's extra blank
    line does not read as tampering."""
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    return text.rstrip("\n") + "\n"


def _has_body(text: str) -> bool:
    """True when the file holds text outside its heading and comments —
    the same bar the prompt loader applies before injecting it."""
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("<!--") and stripped.endswith("-->"):
            continue
        if stripped.startswith("<!--"):
            continue  # comment body line
        if stripped.endswith("-->"):
            continue  # comment end line
        return True
    return False


def expected_digest() -> str:
    """The pinned digest of the expected principles text."""
    return EXPECTED_DIGEST


def check(
    pack_dir: Path | str,
    expected_digests: list[str] | None = None,
) -> GuardResult:
    """Verify ``pack_dir``'s principles addendum against the pinned digest.

    ``pack_dir`` is the directory containing ``shared/`` (the prompt
    pack). ``expected_digests`` optionally lists additional digests the
    maintainer pinned elsewhere; the file may match any of them.
    """
    path = Path(pack_dir) / _PRINCIPLES_REL
    if not path.is_file():
        return GuardResult(
            ok=False,
            reason=f"principles_addendum.md is missing: {path}",
        )
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        return GuardResult(
            ok=False,
            reason=f"principles_addendum.md is unreadable: {exc}",
        )
    if not _has_body(text):
        return GuardResult(
            ok=False,
            reason=("principles_addendum.md holds only a heading or "
                    "scaffold comment — the body is gone"),
        )
    digest = hashlib.sha256(_normalise(text).encode("utf-8")).hexdigest()
    allowed = {EXPECTED_DIGEST, *(expected_digests or [])}
    if digest not in allowed:
        return GuardResult(
            ok=False,
            reason=(f"principles_addendum.md digest mismatch: got "
                    f"{digest}, expected one of {len(allowed)} pinned "
                    f"digests — the text was changed or tampered with"),
        )
    return GuardResult(ok=True, reason="principles addendum matches the pinned digest")
