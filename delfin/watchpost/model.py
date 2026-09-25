"""The one shape every check returns.

A finding names a location (path, optional line), what was seen, and
why it matters. Values of secrets are never carried -- `what` and `why`
describe the pattern that matched, not the matched text.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

SEVERITIES = ("info", "warn", "alert")


@dataclass(frozen=True)
class Finding:
    check: str          # which check produced it, e.g. "ssh"
    severity: str       # info | warn | alert
    path: str           # file or pseudo-path the finding refers to
    line: int | None    # 1-based line when the finding is line-bound
    what: str           # what was seen (never a secret value)
    why: str            # why it matters, one sentence

    def as_dict(self) -> dict:
        return asdict(self)

    @property
    def key(self) -> tuple:
        """Identity for baseline/diff: same place, same story."""
        return (self.check, self.path, self.line, self.what)
