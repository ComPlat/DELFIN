"""Guard for the maintainer's principles addendum.

``pack/shared/principles_addendum.md`` is the first shared contract in
every role prompt (prompt_loader.py, layer 0). This module pins its
digest so tampering becomes visible: a file that is missing, an empty
scaffold, or edited text fails :func:`check` with a clear reason.

The digest is pinned in independent places -- the constant below and a
copy in protected code (api_client.PRINCIPLES_DIGEST) -- and the file
must match EVERY copy: removing or rewriting the principles means
finding and changing all of them, and a copy left behind stops the agent.
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
    """The prompt loader's own bar for "this addendum says something"."""
    from .prompt_loader import _has_body as _loader_has_body
    return _loader_has_body(text)


def expected_digest() -> str:
    """The pinned digest of the expected principles text."""
    return EXPECTED_DIGEST


def check(
    pack_dir: Path | str,
    expected_digests: list[str] | None = None,
) -> GuardResult:
    """Verify ``pack_dir``'s principles addendum against the pinned digest.

    ``pack_dir`` is the directory containing ``shared/`` (the prompt
    pack). ``expected_digests`` lists the copies pinned elsewhere; the
    file must match the module constant AND every one of them.
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
    pinned = [EXPECTED_DIGEST, *(expected_digests or [])]
    stale = [d for d in pinned if d != digest]
    if stale:
        return GuardResult(
            ok=False,
            reason=(f"principles_addendum.md digest mismatch: got "
                    f"{digest}, {len(stale)} of {len(pinned)} pinned "
                    f"digests differ — the text was changed or tampered "
                    f"with, or a pinned copy was not moved with it"),
        )
    return GuardResult(ok=True, reason="principles addendum matches the pinned digest")


class PrinciplesTampered(RuntimeError):
    """The principles file does not match every pinned digest."""


def _shipped_pack() -> Path:
    return Path(__file__).resolve().parent / "pack"


def enforce(pack_dir: Path | str | None = None) -> None:
    """Raise :class:`PrinciplesTampered` unless the shipped principles
    match the module constant AND the protected copy in api_client.

    Called where every agent is born (``AgentEngine.__init__``), so no
    entry point -- terminal, dashboard, scheduler, monitor, benchmark --
    runs an agent without them.
    """
    from .api_client import PRINCIPLES_DIGEST
    result = check(pack_dir or _shipped_pack(),
                   expected_digests=[PRINCIPLES_DIGEST])
    if not result.ok:
        raise PrinciplesTampered(result.reason)
