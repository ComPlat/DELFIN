"""Helpers to parse results produced by the ESD module."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.energies import find_electronic_energy, find_ZPE, find_gibbs_energy


def _format_ms_suffix(trootssl: int) -> str:
    """Format TROOTSSL value as ms suffix (e.g., -1 -> 'msm1', 0 -> 'ms0', 1 -> 'msp1')."""
    if trootssl < 0:
        return f"msm{abs(trootssl)}"
    elif trootssl > 0:
        return f"msp{trootssl}"
    else:
        return "ms0"


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

_FLUOR_RATE_RE = re.compile(
    r"(?:calculated\s+fluorescence\s+rate\s+constant\s+is|fluorescence\s+rate\s+constant\s+is|k[_\s-]*f)\s*=?\s*([0-9.+-Ee]+)\s*s(?:-1|\^-1)",
    flags=re.IGNORECASE,
)
_PHOSP_RATE_RE = re.compile(
    r"(?:calculated\s+phosphorescence\s+rate\s+constant\s+is|phosphorescence\s+rate\s+constant\s+is|k[_\s-]*p)\s*=?\s*([0-9.+-Ee]+)\s*s(?:-1|\^-1)",
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
class FluorResult:
    """Parsed information for a fluorescence transition."""

    rate: Optional[float]
    temperature: Optional[float]
    delta_cm1: Optional[float]
    source: Path


@dataclass
class PhospResult:
    """Parsed information for a phosphorescence transition (multiple sublevels + mean)."""

    sublevel_rates: Tuple[Optional[float], ...]
    rate_mean: Optional[float]
    temperature: Optional[float]
    delta_cm1: Optional[float]
    source: Path

@dataclass
class ESDSummary:
    """Structured results parsed from ESD output files."""

    states: Dict[str, StateResult] = field(default_factory=dict)
    isc: Dict[str, ISCResult] = field(default_factory=dict)
    ic: Dict[str, ICResult] = field(default_factory=dict)
    fluor: Dict[str, FluorResult] = field(default_factory=dict)
    phosp: Dict[str, PhospResult] = field(default_factory=dict)

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
                any(result.rate is not None for result in self.fluor.values()),
                any(result.rate_mean is not None for result in self.phosp.values()),
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


def _parse_fluor_output(path: Path) -> FluorResult:
    text = _read_text(path)
    rate = _safe_float(_FLUOR_RATE_RE.search(text).group(1)) if text and _FLUOR_RATE_RE.search(text) else None
    temp = _safe_float(_TEMP_RE.search(text).group(1)) if text and _TEMP_RE.search(text) else None
    delta = _safe_float(_DELE_RE.search(text).group(1)) if text and _DELE_RE.search(text) else None
    return FluorResult(rate, temp, delta, path)


def _parse_phosp_output(path: Path) -> PhospResult:
    text = _read_text(path)
    rates: Tuple[Optional[float], ...]
    if text:
        matches = _PHOSP_RATE_RE.findall(text)
        rates = tuple(_safe_float(m) for m in matches)
        vals = [r for r in rates if r is not None]
        rate_mean = sum(vals) / len(vals) if vals else None
        temp = _safe_float(_TEMP_RE.search(text).group(1)) if _TEMP_RE.search(text) else None
        delta = _safe_float(_DELE_RE.search(text).group(1)) if _DELE_RE.search(text) else None
    else:
        rates = tuple()
        rate_mean = None
        temp = None
        delta = None
    return PhospResult(rates, rate_mean, temp, delta, path)


def collect_esd_results(
    esd_dir: Path,
    states: Iterable[str],
    iscs: Iterable[str],
    ics: Iterable[str],
    config: Optional[Dict[str, Any]] = None,
) -> ESDSummary:
    """Collect FSPE values and ISC/IC rate constants from ESD outputs.

    Args:
        esd_dir: ESD working directory
        states: List of electronic states
        iscs: List of ISC transitions
        ics: List of IC transitions
        config: Configuration dictionary (for TROOTSSL and ESD_modus parsing)
    """
    from delfin.esd_input_generator import _resolve_state_filename

    summary = ESDSummary()

    if not esd_dir.exists():
        logger.info("ESD directory %s missing; skipping ESD result aggregation.", esd_dir)
        return summary

    # Get ESD mode to resolve correct file names for hybrid1
    esd_mode = 'tddft'  # default
    if config is not None:
        esd_mode = str(config.get('ESD_modus', 'tddft')).strip().lower()
        if "|" in esd_mode:
            esd_mode = esd_mode.split("|")[0].strip()

    for state in states:
        state_key = state.strip().upper()
        if not state_key:
            continue

        # Resolve output filename for hybrid1 mode
        output_filename = _resolve_state_filename(state_key, 'out', esd_mode)
        output_path = esd_dir / output_filename

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

    # Parse TROOTSSL values from config
    trootssl_values = [0]  # default
    if config is not None:
        trootssl_raw = config.get('TROOTSSL', '0')

        # Handle list/tuple directly
        if isinstance(trootssl_raw, (list, tuple)):
            trootssl_values = [int(x) for x in trootssl_raw]
        else:
            # Handle string format
            trootssl_str = str(trootssl_raw).strip()
            # Remove brackets if present (e.g., "['-1', '0', '1']" or "[-1, 0, 1]")
            trootssl_str = trootssl_str.strip('[]').strip()

            # If empty string after stripping, keep default [0]
            if not trootssl_str:
                pass  # Keep default trootssl_values = [0]
            elif ',' in trootssl_str:
                # Split and clean each value (remove quotes), filter out empty strings
                trootssl_values = [
                    int(x.strip().strip("'\""))
                    for x in trootssl_str.split(',')
                    if x.strip().strip("'\"")  # Skip empty strings
                ]
                # If all values were empty, keep default [0]
                if not trootssl_values:
                    trootssl_values = [0]
            else:
                # Single value - only convert if not empty
                cleaned = trootssl_str.strip("'\"")
                if cleaned:
                    trootssl_values = [int(cleaned)]
                else:
                    trootssl_values = [0]

    for isc in iscs:
        isc_key = isc.strip().upper()
        if not isc_key or ">" not in isc_key:
            continue
        init_state, final_state = (part.strip() for part in isc_key.split(">", 1))

        # Collect results for each TROOTSSL value
        for trootssl in trootssl_values:
            ms_suffix = _format_ms_suffix(trootssl)
            filename = esd_dir / f"{init_state}_{final_state}_ISC_{ms_suffix}.out"
            ms_str = f"{trootssl:+d}" if trootssl != 0 else "0"
            key = f"{init_state}>{final_state}(Ms={ms_str})"
            summary.isc[key] = _parse_isc_output(filename)

    for ic in ics:
        ic_key = ic.strip().upper()
        if not ic_key or ">" not in ic_key:
            continue
        init_state, final_state = (part.strip() for part in ic_key.split(">", 1))
        filename = esd_dir / f"{init_state}_{final_state}_IC.out"
        summary.ic[f"{init_state}>{final_state}"] = _parse_ic_output(filename)

    # Fluorescence / phosphorescence rates (controlled by emission_rates or auto-detected)
    emission_rates = set()
    if config is not None:
        raw = config.get("emission_rates", "")
        emission_rates = {t.strip().lower() for t in str(raw).replace(";", ",").replace(" ", ",").split(",") if t.strip()}

    # Auto-detect fluorescence if FLUOR output exists
    fluor_out = esd_dir / "S1_S0_FLUOR.out"
    if "f" in emission_rates or fluor_out.exists():
        summary.fluor["S1>S0"] = _parse_fluor_output(fluor_out)

    # Auto-detect phosphorescence if PHOSP output exists
    phosp_out = esd_dir / "T1_S0_PHOSP.out"
    if "p" in emission_rates or phosp_out.exists():
        summary.phosp["T1>S0"] = _parse_phosp_output(phosp_out)

    return summary
