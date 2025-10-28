"""Helpers to parse results produced by the ESD module."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Optional

from delfin.common.logging import get_logger
from delfin.energies import find_electronic_energy

logger = get_logger(__name__)

_ISC_RATE_RE = re.compile(
    r"The\s+calculated\s+ISC\s+rate\s+constant\s+is\s+([0-9.+-Ee]+)\s*s(?:-1|\^-1)",
    flags=re.IGNORECASE,
)

_IC_RATE_RE = re.compile(
    r"The\s+calculated\s+internal\s+conversion\s+rate\s+constant\s+is\s+([0-9.+-Ee]+)\s*s(?:-1|\^-1)",
    flags=re.IGNORECASE,
)


@dataclass
class ESDSummary:
    """Structured results parsed from ESD output files."""

    states_fspe: Dict[str, Optional[float]] = field(default_factory=dict)
    isc_rates: Dict[str, Optional[float]] = field(default_factory=dict)
    ic_rates: Dict[str, Optional[float]] = field(default_factory=dict)

    @property
    def has_data(self) -> bool:
        return any(
            (
                any(value is not None for value in self.states_fspe.values()),
                any(value is not None for value in self.isc_rates.values()),
                any(value is not None for value in self.ic_rates.values()),
            )
        )


def _read_text(path: Path) -> Optional[str]:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except FileNotFoundError:
        logger.info("File %s not found; skipping ESD parsing.", path)
        return None


def _parse_rate_constant(text: Optional[str], pattern: re.Pattern[str]) -> Optional[float]:
    if not text:
        return None
    match = pattern.search(text)
    if not match:
        return None
    try:
        value = float(match.group(1))
    except ValueError:
        logger.warning("Failed to convert '%s' to float during ESD parsing.", match.group(1))
        return None
    return value


def collect_esd_results(
    esd_dir: Path,
    states: Iterable[str],
    iscs: Iterable[str],
    ics: Iterable[str],
) -> ESDSummary:
    """Collect FSPE values and ISC/IC rate constants from ESD outputs."""

    summary = ESDSummary()

    if not esd_dir.exists():
        logger.info("ESD directory %s missing; skipping ESD result aggregation.", esd_dir)
        return summary

    # Final single point energies for requested states
    for state in states:
        state_key = state.strip().upper()
        if not state_key:
            continue
        output_path = esd_dir / f"{state_key}.out"
        if output_path.exists():
            summary.states_fspe[state_key] = find_electronic_energy(str(output_path))
        else:
            logger.info("ESD state output %s missing; skipping FSPE extraction.", output_path)
            summary.states_fspe[state_key] = None

    # ISC rate constants
    for isc in iscs:
        isc_key = isc.strip().upper()
        if not isc_key or ">" not in isc_key:
            continue
        init_state, final_state = (part.strip() for part in isc_key.split(">", 1))
        filename = esd_dir / f"{init_state}_{final_state}_ISC.out"
        text = _read_text(filename)
        summary.isc_rates[f"{init_state}>{final_state}"] = _parse_rate_constant(text, _ISC_RATE_RE)

    # IC rate constants
    for ic in ics:
        ic_key = ic.strip().upper()
        if not ic_key or ">" not in ic_key:
            continue
        init_state, final_state = (part.strip() for part in ic_key.split(">", 1))
        filename = esd_dir / f"{init_state}_{final_state}_IC.out"
        text = _read_text(filename)
        summary.ic_rates[f"{init_state}>{final_state}"] = _parse_rate_constant(text, _IC_RATE_RE)

    return summary
