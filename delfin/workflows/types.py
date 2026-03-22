"""Base interfaces for DELFIN workflows.

Uses :class:`typing.Protocol` (structural subtyping) so existing workflow
implementations can conform without changing their inheritance hierarchy.
"""

from __future__ import annotations

from typing import Any, Dict, List, Protocol, runtime_checkable


@runtime_checkable
class Workflow(Protocol):
    """Protocol that all DELFIN workflows implement.

    Minimal contract:
      - *name*        -- unique identifier used by the registry.
      - *description* -- human-readable, for CLI ``--help`` / dashboard.
      - *run()*       -- execute the workflow.
      - *run_cli()*   -- CLI entry point.
    """

    name: str
    description: str

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> Any:
        """Execute the workflow.  Returns a workflow-specific result."""
        ...

    def run_cli(self, argv: List[str]) -> int:
        """CLI entry point.  Returns exit code."""
        ...
