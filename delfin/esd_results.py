"""Helpers to parse results produced by the ESD module."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.energies import find_electronic_energy, find_ZPE, find_gibbs_energy

logger = get_logger(__name__)

# Regular expressions for parsing ESD outputs
_ISC_RATE_RE = re.compile(
    r"The\s+calculated\s+ISC\s+rate\s+constant\s+is\s+([0-9.+-Ee]+)\s*s(?:-1|\^-1)",
    flags=re.IGNORECASE,
)
_IC_RATE_RE = re.compile(
    r"The\s+calculated\s+internal\s+conversion\s+rate\s+constant\s+is\s+([0-9.+-Ee]+)\s*s(?:-1|\^-1)",
    flags=re.IGNORECASE,
)
_TEMP_RE = re.compile(
    r"Temperature\s+used:\s*([0-9.+-Ee]+)\s*K",
    flags=re.IGNORECASE,
)
_DELE_RE = re.compile(
    r"0-0\s+energy\s+difference:\s*([0-9.+-Ee]+)\s*cm-1",
    flags=re.IGNORECASE,
)
_SOC_RE = re.compile(
    r"Reference\s+SOC\s+\(Re\s+and\s+Im\):\s*([0-9.+-Ee]+),\s*([0-9.+-Ee]+)",
    flags=re.IGNORECASE,
)
_FC_HT_RE = re.compile(
    r"with\s+([0-9.+-Ee]+)\s+from\s+FC\s+and\s+([0-9.+-Ee]+)\s+from\s+HT",
    flags=re.IGNORECASE,
)


def _safe_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        logger.warning("Failed to convert '%s' to float during ESD parsing.", value)
        return None


@dataclass
class StateResult:
    """Result information for a single electronic state."""

    fspe: Optional[float]
    zpe: Optional[float]
    gibbs: Optional[float]
    source: Path


@dataclass
class ISCResult:
    """Parsed information for an ISC transition."""

    rate: Optional[float]
    temperature: Optional[float]
    delta_cm1: Optional[float]
    soc: Optional[Tuple[Optional[float], Optional[float]]]
    fc_percent: Optional[float]
    ht_percent: Optional[float]
    source: Path


@dataclass
class ICResult:
    """Parsed information for an IC transition."""

    rate: Optional[float]
    temperature: Optional[float]
    delta_cm1: Optional[float]
    source: Path


@dataclass
class ESDSummary:
    """Structured results parsed from ESD output files."""

    states: Dict[str, StateResult] = field(default_factory=dict)
    isc: Dict[str, ISCResult] = field(default_factory=dict)
    ic: Dict[str, ICResult] = field(default_factory=dict)

    @property
    def states_fspe(self) -> Dict[str, Optional[float]]:
        return {key: result.fspe for key, result in self.states.items()}

    @property
    def isc_rates(self) -> Dict[str, Optional[float]]:
        return {key: result.rate for key, result in self.isc.items()}

    @property
    def ic_rates(self) -> Dict[str, Optional[float]]:
        return {key: result.rate for key, result in self.ic.items()}

    @property
    def has_data(self) -> bool:
        return any(
            (
                any(result.fspe is not None for result in self.states.values()),
                any(result.rate is not None for result in self.isc.values()),
                any(result.rate is not None for result in self.ic.values()),
            )
        )


def _read_text(path: Path) -> Optional[str]:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except FileNotFoundError:
        logger.info("File %s not found; skipping ESD parsing.", path)
        return None


def _parse_isc_output(path: Path) -> ISCResult:
    text = _read_text(path)
    rate = _safe_float(_ISC_RATE_RE.search(text).group(1)) if text and _ISC_RATE_RE.search(text) else None
    temp = _safe_float(_TEMP_RE.search(text).group(1)) if text and _TEMP_RE.search(text) else None
    delta = _safe_float(_DELE_RE.search(text).group(1)) if text and _DELE_RE.search(text) else None
    if text:
        soc_match = _SOC_RE.search(text)
        if soc_match:
            soc = (_safe_float(soc_match.group(1)), _safe_float(soc_match.group(2)))
        else:
            soc = None
        fc_ht_match = _FC_HT_RE.search(text)
        if fc_ht_match:
            fc = _safe_float(fc_ht_match.group(1))
            ht = _safe_float(fc_ht_match.group(2))
        else:
            fc = ht = None
    else:
        soc = None
        fc = ht = None
    return ISCResult(rate, temp, delta, soc, fc, ht, path)


def _parse_ic_output(path: Path) -> ICResult:
    text = _read_text(path)
    rate = _safe_float(_IC_RATE_RE.search(text).group(1)) if text and _IC_RATE_RE.search(text) else None
    temp = _safe_float(_TEMP_RE.search(text).group(1)) if text and _TEMP_RE.search(text) else None
    delta = _safe_float(_DELE_RE.search(text).group(1)) if text and _DELE_RE.search(text) else None
    return ICResult(rate, temp, delta, path)


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

    for state in states:
        state_key = state.strip().upper()
        if not state_key:
            continue
        output_path = esd_dir / f"{state_key}.out"
        if output_path.exists():
            fspe = find_electronic_energy(str(output_path))
            zpe = find_ZPE(str(output_path))
            gibbs = find_gibbs_energy(str(output_path))
        else:
            logger.info("ESD state output %s missing; skipping FSPE/ZPE/Gibbs extraction.", output_path)
            fspe = None
            zpe = None
            gibbs = None
        summary.states[state_key] = StateResult(fspe, zpe, gibbs, output_path)

    for isc in iscs:
        isc_key = isc.strip().upper()
        if not isc_key or ">" not in isc_key:
            continue
        init_state, final_state = (part.strip() for part in isc_key.split(">", 1))
        filename = esd_dir / f"{init_state}_{final_state}_ISC.out"
        summary.isc[f"{init_state}>{final_state}"] = _parse_isc_output(filename)

    for ic in ics:
        ic_key = ic.strip().upper()
        if not ic_key or ">" not in ic_key:
            continue
        init_state, final_state = (part.strip() for part in ic_key.split(">", 1))
        filename = esd_dir / f"{init_state}_{final_state}_IC.out"
        summary.ic[f"{init_state}>{final_state}"] = _parse_ic_output(filename)

    return summary
