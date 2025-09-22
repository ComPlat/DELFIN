"""Helpers for assembling ORCA input sections."""
from __future__ import annotations

from typing import Any, Dict, Iterable, List


def normalize_bang(bang: str) -> str:
    """Ensure the method line ends with a newline for consistent serialization."""
    return bang if bang.endswith("\n") else bang + "\n"


def collect_output_blocks(config: Dict[str, Any]) -> List[str]:
    """Build optional %output blocks based on configuration flags."""
    blocks: List[str] = []
    if str(config.get("print_MOs", "no")).lower() == "yes":
        blocks.append("%output\nprint[p_mos] 1\nprint[p_basis] 2\nend\n")
    if str(config.get("print_Loewdin_population_analysis", "no")).lower() == "yes":
        blocks.append("%output\nprint[P_ReducedOrbPopMO_L] 1\nend\n")
    return blocks


class OrcaInputBuilder:
    """Utility for assembling ordered ORCA input sections without altering behaviour."""

    def __init__(self, bang: str) -> None:
        self._lines: List[str] = [normalize_bang(bang)]

    def add_resources(self, maxcore: int, pal: int) -> None:
        self._lines.append(f"%maxcore {maxcore}\n%pal nprocs {pal} end\n")

    def add_additions(self, additions: str) -> None:
        if additions:
            stripped = additions.strip()
            if stripped:
                self._lines.append(additions if additions.endswith("\n") else f"{additions}\n")

    def add_block(self, block: str) -> None:
        if block:
            self._lines.append(block)

    def add_blocks(self, blocks: Iterable[str]) -> None:
        self._lines.extend(blocks)

    @property
    def lines(self) -> List[str]:
        return self._lines
