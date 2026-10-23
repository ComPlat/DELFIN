"""Intelligent ORCA error recovery system.

This module provides automatic detection and recovery from common ORCA
calculation failures by analyzing output files, identifying error patterns,
and modifying input files with appropriate fixes.

Key features:
- Automatic error classification (SCF convergence, TRAH crashes, etc.)
- Intelligent input file modification with MOREAD for continuation
- Progressive escalation of fixes across retry attempts
- State tracking to prevent infinite loops
- Preservation of basis set specifications and geometry
- Universal preparation of inputs for continuation (xyz + old.gbw)
"""

import json
import re
import shutil
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional

from delfin.common import orca_input
from delfin.common.logging import get_logger
from delfin.common.tddft_settings import ORCA_DEFAULT_MAXDIM

logger = get_logger(__name__)

_OPT_GAVE_UP = "The optimization did not converge"
_JOB_BANNER = re.compile(r"\${4,}\s+JOB NUMBER\s+(\d+)\s+\$")

# How much of the end of an ORCA output is searched for error signatures.
# ORCA reports what went wrong immediately before it aborts, so the tail is
# where the evidence is; reading more only costs memory at the worst moment.
_ERROR_SCAN_TAIL_BYTES = 2 * 1024 * 1024


class OrcaErrorType(Enum):
    """Classification of ORCA error types."""

    SCF_NO_CONVERGENCE = "scf_convergence"
    """SCF failed to converge within maxiter iterations."""

    LEANSCF_NOT_CONVERGED = "leanscf_not_converged"
    """LEANSCF (coupled-perturbed SCF for FREQ) failed to converge."""

    TRAH_SEGFAULT = "trah_crash"
    """Segmentation fault during TRAH-SCF procedure."""

    GEOMETRY_NOT_CONVERGED = "geom_not_converged"
    """Geometry optimization did not converge."""

    FREQUENCY_FAILURE = "freq_failure"
    """Numerical frequency calculation failed."""

    MEMORY_ERROR = "memory_error"
    """Insufficient memory or memory allocation failure."""

    MPI_CRASH = "mpi_crash"
    """MPI communication failure or process crash."""

    CIS_FAILURE = "cis_failure"
    """CIS/TD-DFT module aborted (e.g. RPA diagonalization failure)."""

    TDDFT_ROOT_COLLAPSE = "tddft_root_collapse"
    """TD-DFT excitation energy collapsed to ~0 or went negative.

    Signature of a triplet instability of the closed-shell reference under
    full TD-DFT (RPA). ORCA does not flag this as an error and happily
    carries the garbage root through an entire geometry optimization.
    """

    DIIS_ERROR = "diis_error"
    """DIIS convergence acceleration failed."""

    INPUT_ERROR = "input_error"
    """ORCA refused the input: no retry can change what was asked.

    146 archived OCCUPIER runs stopped on "multiplicity (3) is odd and
    number of electrons (243) is odd -> impossible" -- broken-symmetry pairs
    of an older sequence profile whose unpaired electrons cannot match the
    electron count.
    """

    ESD_RATE_UNPHYSICAL = "esd_rate_unphysical"
    """ORCA's ESD module ended normally with a negative rate constant.

    ORCA says so itself ("negative rates are unphysical ... something went
    wrong with the CorrFunc integration") and suggests more points.
    """

    ESD_WINDOW_TRUNCATED = "esd_window_truncated"
    """The ESD correlation function was cut off before it decayed.

    ORCA does not warn.  A MAXTIME shorter than ORCA's own window leaves the
    rate unconverged: 2.3x too fast for formaldehyde's ISC at 290 fs.
    """

    TRANSIENT_SYSTEM_ERROR = "transient_system"
    """Temporary system errors (disk full, network issues, etc.) - retry without modification."""

    UNKNOWN = "unknown"
    """Unknown or unclassified error."""


class OrcaErrorDetector:
    """Analyzes ORCA output files to detect and classify errors."""

    # Error pattern definitions with priority (lower = higher priority)
    ERROR_PATTERNS = [
        # ORCA refused the input: nothing a retry changes (highest priority)
        {
            "type": OrcaErrorType.INPUT_ERROR,
            "patterns": [
                "and number of electrons",
                "INPUT ERROR",
                "Unknown identifier",
                "Invalid assignment",
                "Unrecognized symbol",
                "UNRECOGNIZED OR DUPLICATED KEYWORD",
                "basis set was either not assigned",
            ],
            "all_required": False,
            "priority": 0,
        },
        # TRAH crashes (highest priority - most specific)
        {
            "type": OrcaErrorType.TRAH_SEGFAULT,
            "patterns": [
                "Signal: Segmentation fault",
                "TRAH",
                "Auto-TRAH",
            ],
            "all_required": True,
            "priority": 1,
        },
        # LEANSCF SCF convergence failures (higher priority than generic MPI crash)
        {
            "type": OrcaErrorType.LEANSCF_NOT_CONVERGED,
            "patterns": [
                "LEANSCF",
                "SCF has not converged",
            ],
            "all_required": True,
            "priority": 2,
        },
        # CIS/TD-DFT module aborts. Must outrank the generic MPI_CRASH
        # block below: an RPA diagonalization failure also prints the mpirun
        # "exited on signal" banner, but it is a deterministic numerical
        # failure, not a parallelization hiccup - retrying it with fewer
        # cores reproduces it exactly.
        {
            "type": OrcaErrorType.CIS_FAILURE,
            "patterns": [
                "error termination in CIS",
            ],
            "all_required": False,
            "priority": 2,
        },
        {
            "type": OrcaErrorType.CIS_FAILURE,
            "patterns": [
                "TRandomPhaseApproximationDiagonalization",
            ],
            "all_required": False,
            "priority": 2,
        },
        # MPI crashes - general process crashes
        {
            "type": OrcaErrorType.MPI_CRASH,
            "patterns": [
                "mpirun noticed that process rank",
                "exited on signal",
            ],
            "all_required": True,
            "priority": 3,
        },
        # MPI crashes - ORCA parallelization bugs (TGeneralVectorSet, etc.)
        {
            "type": OrcaErrorType.MPI_CRASH,
            "patterns": [
                "TGeneralVectorSet",
                "Constructor called with NVecs<=0",
            ],
            "all_required": True,
            "priority": 3,
        },
        # MPI crashes - LEANSCF failures (generic - lower priority than specific LEANSCF SCF convergence)
        {
            "type": OrcaErrorType.MPI_CRASH,
            "patterns": [
                "error termination in LEANSCF",
            ],
            "all_required": False,
            "priority": 4,
        },
        # SCF convergence failures
        {
            "type": OrcaErrorType.SCF_NO_CONVERGENCE,
            "patterns": [
                "SCF NOT CONVERGED",
                "Error : ORCA finished by error termination in SCF",
            ],
            "all_required": False,
            "priority": 3,
        },
        # DIIS errors.
        # NOTE: "DIIS error" must NOT be listed here. Pattern matching is
        # case-insensitive and every ORCA SCF settings block prints
        #     DIIS Error             TolErr          ....  1.000e-06
        #     Last DIIS Error            ...  1.04e-03  Tolerance : 1.0e-06
        # so that substring occurs dozens of times in perfectly healthy
        # output. It made every truncated run (walltime kill, crash) look
        # like a DIIS failure and misrouted its recovery. Only genuine
        # failure phrasings belong here; anything not matched falls through
        # to SCF_NO_CONVERGENCE or UNKNOWN, which is the safer outcome.
        {
            "type": OrcaErrorType.DIIS_ERROR,
            "patterns": [
                "Error in DIIS",
                "DIIS did not converge",
            ],
            "all_required": False,
            "priority": 4,
        },
        # Geometry convergence
        {
            "type": OrcaErrorType.GEOMETRY_NOT_CONVERGED,
            "patterns": [
                "NOT CONVERGED",
                "GEOMETRY OPTIMIZATION",
            ],
            "all_required": True,
            "priority": 5,
        },
        # Frequency failures
        {
            "type": OrcaErrorType.FREQUENCY_FAILURE,
            "patterns": [
                "error termination in SCF RESPONSE",
                "error termination in NUMFREQ",
            ],
            "all_required": False,
            "priority": 6,
        },
        # Memory errors
        # Memory exhaustion as ORCA 6.x words it. This needs to outrank the
        # generic "error termination in LEANSCF" rule below, because ORCA
        # routes an out-of-memory abort through LEANSCF as well: measured
        # against ORCA 6.1.1, a run capped at %maxcore 1 was classified as an
        # MPI crash and answered with MPI fixes, while ORCA had plainly written
        # what it needed. The message also carries the required MaxCore, which
        # _memory_error_fixes reads back instead of guessing.
        {
            "type": OrcaErrorType.MEMORY_ERROR,
            "patterns": [
                "Not enough memory available",
                "Please increase MaxCore to more than",
            ],
            "all_required": False,
            "priority": 2,
        },
        {
            "type": OrcaErrorType.MEMORY_ERROR,
            "patterns": [
                "insufficient memory",
                "cannot allocate memory",
                "memory allocation failed",
                "out-of-memory",
                "out of memory",
                "not a single batch is possible with the present MaxCore",
                "BatchOrganizer",
            ],
            "all_required": False,
            "priority": 7,
        },
        # Transient system errors (lowest priority - retry without modification)
        {
            "type": OrcaErrorType.TRANSIENT_SYSTEM_ERROR,
            "patterns": [
                "cannot create temporary file",
                "no space left on device",
                "disk quota exceeded",
                "stale file handle",
                "connection timed out",
                "broken pipe",
                "input/output error",
                "resource temporarily unavailable",
            ],
            "all_required": False,
            "priority": 99,  # Lowest priority - check this last
        },
    ]

    # Excitation energies below this are unphysical for any real molecule
    # and signal a collapsed TD-DFT root rather than a genuinely low-lying
    # state. Observed values in collapsed runs: 1e-4 .. 7e-4 Eh, and
    # -7.1e-8 Eh (i.e. negative) for an outright triplet instability.
    ROOT_COLLAPSE_THRESHOLD_EH = 1.0e-3

    # A root that drops below this fraction of the median root energy of the
    # same run has collapsed, even if its absolute value looks harmless.
    _RELATIVE_COLLAPSE_FACTOR = 0.2
    _MIN_SAMPLES_FOR_RELATIVE_CHECK = 5
    _MAX_DE_CIS_SAMPLES = 50000

    _DE_CIS_RE = re.compile(
        r"DE\(CIS\)\s*=\s*(-?\d+\.?\d*(?:[eEdD][-+]?\d+)?)"
    )

    # Cheap gate before the full-file scan below. Non-TD-DFT runs - the vast
    # majority - must not pay for it.
    _TDDFT_TAIL_MARKERS = ("DE(CIS)", "TD-DFT", "ORCA-CIS", "orca_cis")

    @classmethod
    def _looks_like_tddft(cls, tail: str) -> bool:
        """True if the output tail shows any sign of a TD-DFT/CIS run."""
        return any(marker in tail for marker in cls._TDDFT_TAIL_MARKERS)

    @classmethod
    def _detect_collapsed_root(
        cls, output_file: Path, terminated_normally: bool
    ) -> Optional[float]:
        """Return the offending DE(CIS) value if a TD-DFT root collapsed.

        A run that ended normally is only rejected when its *final* root is
        collapsed - a transient collapse the optimizer recovered from should
        not invalidate an otherwise converged result. A run that did not end
        normally is rejected on the lowest root it ever produced, because a
        collapse in cycle 12 can leave the optimizer thrashing for another
        140 cycles and a perfectly healthy-looking tail.

        Unlike pattern matching this cannot work on a bounded tail, so the
        file is streamed line by line: constant memory, one pass, and the
        cheap substring guard keeps it fast even on 300k-line outputs.
        """
        minimum = None
        last = None
        samples = []  # bounded reservoir, only needed for the median
        try:
            with open(output_file, 'r', encoding='utf-8', errors='replace') as fh:
                for line in fh:
                    if "DE(CIS)" not in line:
                        continue
                    match = cls._DE_CIS_RE.search(line)
                    if not match:
                        continue
                    token = match.group(1).replace("D", "E").replace("d", "e")
                    try:
                        value = float(token)
                    except ValueError:
                        continue
                    last = value
                    if minimum is None or value < minimum:
                        minimum = value
                    if len(samples) < cls._MAX_DE_CIS_SAMPLES:
                        samples.append(value)
        except OSError:
            return None

        if last is None:
            return None

        candidate = last if terminated_normally else minimum

        # Absolute criterion: an excitation energy at or below ~0.03 eV is
        # unphysical for a closed-shell organic chromophore.
        if candidate <= cls.ROOT_COLLAPSE_THRESHOLD_EH:
            return candidate

        # Relative criterion. The absolute threshold has to stay tight,
        # because a genuinely near-degenerate system (diradical, broken-
        # symmetry metal complex) legitimately has a tiny gap and must not
        # be flagged. A collapse in an otherwise normal run looks different:
        # one root drops far below all the others. Compare against the
        # median, which a handful of collapsed cycles cannot drag down.
        if len(samples) >= cls._MIN_SAMPLES_FOR_RELATIVE_CHECK:
            ordered = sorted(samples)
            median = ordered[len(ordered) // 2]
            if median > 0 and candidate < cls._RELATIVE_COLLAPSE_FACTOR * median:
                return candidate

        return None

    @classmethod
    def unusable_result(cls, output_file: Path, input_file: Optional[Path] = None):
        """A run that ended normally but whose result is not one: (error type, job index, reason) or None.

        ORCA ends with "ORCA TERMINATED NORMALLY" after an optimisation that
        ran out of cycles, a TD-DFT root that collapsed, and an ESD rate that
        is negative or cut off early; judged on that marker, each counted as
        done and the retry for it never ran.
        """
        output_file = Path(output_file)
        if not output_file.exists():
            return None
        collapsed = cls.analyze_output(output_file)
        if collapsed == OrcaErrorType.TDDFT_ROOT_COLLAPSE:
            return OrcaErrorType.TDDFT_ROOT_COLLAPSE, None, "a TD-DFT root collapsed"
        wants_scan = True
        if input_file is not None:
            try:
                head = Path(input_file).read_text(encoding="utf-8", errors="replace").lower()
                wants_scan = "opt" in head or "esd(" in head
            except OSError:
                wants_scan = True
        if not wants_scan:
            return None
        from delfin.common.esd_numerics import rate_problem

        job = 0
        section: List[str] = []
        found = None

        def check(text: str, index: int):
            if _OPT_GAVE_UP in text:
                return OrcaErrorType.GEOMETRY_NOT_CONVERGED, index, "the optimisation ran out of cycles"
            problem = rate_problem(text) if "EXCITED STATE DYNAMICS" in text else None
            if problem:
                kind = (OrcaErrorType.ESD_RATE_UNPHYSICAL if problem.startswith("unphysical")
                        else OrcaErrorType.ESD_WINDOW_TRUNCATED)
                return kind, index, f"the rate is {problem}"
            return None

        try:
            with output_file.open(encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if _JOB_BANNER.search(line):
                        found = check("".join(section), job) if section else None
                        if found:
                            return found
                        job = int(_JOB_BANNER.search(line).group(1)) - 1
                        section = []
                        continue
                    # only the lines the checks read are kept
                    if ("optimization did not converge" in line or "rate constant is" in line
                            or "negative rates" in line or "linewidth is" in line or "Maximum time" in line
                            or "EXCITED STATE DYNAMICS" in line):
                        section.append(line)
        except OSError:
            return None
        return check("".join(section), job) if section else None

    @classmethod
    def analyze_output(cls, output_file: Path) -> Optional[OrcaErrorType]:
        """Analyze ORCA output file and classify the error.

        Args:
            output_file: Path to ORCA .out file

        Returns:
            OrcaErrorType if an error is detected, None if calculation succeeded
        """
        import time

        # Wait for output file with retry logic (race condition fix)
        max_wait_attempts = 10
        wait_interval = 0.5  # seconds

        for attempt in range(max_wait_attempts):
            if output_file.exists():
                # File exists, give it a moment to be fully written
                time.sleep(0.2)
                break
            if attempt < max_wait_attempts - 1:
                time.sleep(wait_interval)
        else:
            # File still doesn't exist after all retries
            logger.warning(f"Output file does not exist after {max_wait_attempts * wait_interval}s: {output_file}")
            return None

        try:
            # Read a bounded tail rather than the whole file. ORCA outputs
            # reach hundreds of MB, and this runs on every failure — often for
            # several parallel jobs at once, which is exactly when memory is
            # already tight. 2 MB is more context than the previous "last
            # 10000 lines" (~1 MB at ORCA's line length) and costs a seek
            # instead of loading the file into a list of lines.
            with open(output_file, 'r', encoding='utf-8', errors='replace') as f:
                size = output_file.stat().st_size
                if size > _ERROR_SCAN_TAIL_BYTES:
                    f.seek(size - _ERROR_SCAN_TAIL_BYTES)
                content = f.read()

            terminated_normally = "ORCA TERMINATED NORMALLY" in content

            # Triplet-instability guard. Must run *before* the
            # terminated-normally shortcut: a collapsed root does not abort
            # ORCA, it silently poisons every downstream energy.
            # The scan below streams the whole file, which the bounded tail
            # read above deliberately avoids. Only pay that cost when the
            # tail actually shows a TD-DFT/CIS run.
            collapsed = (
                cls._detect_collapsed_root(output_file, terminated_normally)
                if cls._looks_like_tddft(content)
                else None
            )
            if collapsed is not None:
                criterion = (
                    f"<= {cls.ROOT_COLLAPSE_THRESHOLD_EH:.1e} Eh absolute"
                    if collapsed <= cls.ROOT_COLLAPSE_THRESHOLD_EH
                    else f"< {cls._RELATIVE_COLLAPSE_FACTOR:g} x median"
                )
                logger.error(
                    "Collapsed TD-DFT root in %s: DE(CIS) = %.3e Eh (%s). "
                    "This is a triplet instability of the reference under "
                    "full TD-DFT; the fix is 'tda true', not fewer cores.",
                    output_file.name,
                    collapsed,
                    criterion,
                )
                return OrcaErrorType.TDDFT_ROOT_COLLAPSE

            # Check for successful termination first
            if terminated_normally:
                return None

            # Check error patterns by priority
            detected_errors = []
            for pattern_def in cls.ERROR_PATTERNS:
                if cls._matches_pattern(content, pattern_def):
                    detected_errors.append(
                        (pattern_def["priority"], pattern_def["type"])
                    )

            if not detected_errors:
                # Generic failure without specific pattern
                if "error termination" in content.lower():
                    return OrcaErrorType.UNKNOWN
                return None

            # Return highest priority error
            detected_errors.sort(key=lambda x: x[0])
            error_type = detected_errors[0][1]

            logger.info(f"Detected ORCA error type: {error_type.value} in {output_file.name}")
            return error_type

        except Exception as e:
            logger.error(f"Error analyzing output file {output_file}: {e}", exc_info=True)
            return None

    @staticmethod
    def _matches_pattern(content: str, pattern_def: dict) -> bool:
        """Check if content matches pattern definition.

        Args:
            content: File content to search
            pattern_def: Pattern definition dict

        Returns:
            True if pattern matches according to definition
        """
        patterns = pattern_def["patterns"]
        all_required = pattern_def.get("all_required", False)

        matches = [pattern.lower() in content.lower() for pattern in patterns]

        if all_required:
            return all(matches)
        else:
            return any(matches)


class RecoveryStrategy:
    """Defines recovery modifications for specific error types and retry attempts."""

    def __init__(self, error_type: OrcaErrorType, attempt: int, config: Dict):
        """Initialize recovery strategy.

        Args:
            error_type: Type of error to recover from
            attempt: Retry attempt number (1-indexed)
            config: DELFIN configuration dict
        """
        self.error_type = error_type
        self.attempt = attempt
        self.config = config
        self.parsed_input = None  # Will be set by apply_recovery before get_modifications
        # Optional: the output that failed. Set by the caller the same way as
        # parsed_input; when present, strategies can read what ORCA itself
        # reported instead of guessing. Absent, every strategy behaves as before.
        self.output_file: Optional[Path] = None
        # Optional: which job of a $new_job file the problem is in, when the
        # caller knows (a run that ended normally ran every job; the last
        # banner is then not the one that failed).
        self.job_index: Optional[int] = None

    def get_modifications(self) -> Dict:
        """Get input file modifications for this error/attempt combination.

        Returns:
            Dictionary of modifications to apply
        """
        if self.error_type == OrcaErrorType.SCF_NO_CONVERGENCE:
            return self._scf_convergence_fixes()

        elif self.error_type == OrcaErrorType.LEANSCF_NOT_CONVERGED:
            return self._leanscf_convergence_fixes()

        elif self.error_type == OrcaErrorType.TRAH_SEGFAULT:
            return self._trah_crash_fixes()

        elif self.error_type == OrcaErrorType.DIIS_ERROR:
            return self._diis_error_fixes()

        elif self.error_type == OrcaErrorType.GEOMETRY_NOT_CONVERGED:
            return self._geometry_convergence_fixes()

        elif self.error_type == OrcaErrorType.MPI_CRASH:
            return self._mpi_crash_fixes()

        elif self.error_type in (
            OrcaErrorType.CIS_FAILURE,
            OrcaErrorType.TDDFT_ROOT_COLLAPSE,
        ):
            return self._tddft_instability_fixes()

        elif self.error_type == OrcaErrorType.FREQUENCY_FAILURE:
            return self._frequency_failure_fixes()

        elif self.error_type == OrcaErrorType.MEMORY_ERROR:
            return self._memory_error_fixes()

        elif self.error_type == OrcaErrorType.TRANSIENT_SYSTEM_ERROR:
            return self._transient_error_fixes()

        elif self.error_type in (OrcaErrorType.ESD_RATE_UNPHYSICAL, OrcaErrorType.ESD_WINDOW_TRUNCATED):
            return self._esd_rate_fixes()

        return {}

    def _printed_esd_grid(self) -> tuple:
        """(number of points, maximum time in a.u.) the failed ESD run printed, or (None, None)."""
        if self.output_file is None:
            return None, None
        try:
            text = Path(self.output_file).read_text(encoding="utf-8", errors="replace")
        except OSError:
            return None, None
        points = re.findall(r"Number of points:\s*(\d+)", text)
        window = re.findall(r"Maximum time:\s*([0-9.]+)\s*fs", text)
        return (int(points[-1]) if points else None,
                float(window[-1]) / 0.02418884326585747 if window else None)

    def _esd_rate_fixes(self) -> Dict:
        """ESD rates that ended normally but are not a number.

        A window cut short (only a MAXTIME set below ORCA's own gives one) is
        handed back to ORCA, which chooses it from the linewidth.  A negative
        rate on ORCA's window gets what ORCA's warning asks for: more points,
        then also a longer window.  The job's own orbitals make the TD-DFT
        part of the rerun cheap; derivatives for Herzberg-Teller terms are
        computed again (ORCA offers no way to read them back).
        """
        points, window_au = self._printed_esd_grid()
        if self.error_type == OrcaErrorType.ESD_WINDOW_TRUNCATED:
            if self.attempt == 1:
                return {"use_moread": True, "esd_block_remove": ["NPOINTS", "MAXTIME"]}
            return {}
        base_points = points or 131072
        if self.attempt == 1:
            return {"use_moread": True, "esd_block_remove": ["MAXTIME"],
                    "esd_block": {"NPOINTS": 4 * base_points}}
        if self.attempt == 2:
            fixes = {"use_moread": True, "esd_block": {"NPOINTS": 16 * base_points}}
            if window_au:
                fixes["esd_block"]["MAXTIME"] = int(2 * window_au)
            return fixes
        return {}

    def _transient_error_fixes(self) -> Dict:
        """Handle transient system errors with exponential backoff + MOREAD.

        For transient errors (disk full, network timeout, etc.), the input
        is fine - we just need to wait and continue from last state.

        Uses MOREAD to continue from last .gbw, plus exponential backoff.
        """
        return {
            "use_moread": True,  # Continue from last state
            "backoff_delay": 2 ** self.attempt,  # Exponential backoff (2s, 4s, 8s, ...)
        }

    def _scf_convergence_fixes(self) -> Dict:
        """Progressive fixes for SCF convergence failures.

        Based on ORCA Manual recommendations.
        DeltaSCF/hybrid1 gets additional aggressive fixes; other jobs keep standard SCF fixes.
        """
        # Check if this is a deltaSCF (or hybrid1 deltaSCF step) calculation
        is_deltascf = self._is_deltascf()

        if is_deltascf:
            if self.attempt == 1:
                # Attempt 1: TightSCF/SlowConv + increased MaxIter
                return {
                    "use_moread": True,
                    "keywords_add": ["TightSCF", "SlowConv"],
                    "scf_block": {
                        "MaxIter": 600,
                    },
                }
            elif self.attempt == 2:
                # Attempt 2: Disable DIIS/TRAH when they are problematic
                return {
                    "use_moread": True,
                    "keywords_add": ["NoDIIS", "NoTRAH"],
                    "scf_block": {
                        "MaxIter": 800,
                    },
                }
            else:
                # Attempt 3+: FreezeAndRelease + GMF + very slow convergence + damping
                logger.info("Adding FreezeAndRelease + GMF for deltaSCF recovery")
                return {
                    "use_moread": True,
                    "keywords_add": ["VerySlowConv", "FreezeAndRelease", "GMF"],
                    "scf_block": {
                        "MaxIter": 1000,
                        "CNVSOSCF": "true",  # %scf switch; a bare SOSCF line opens a sub-block
                        "DampFac": 0.9,  # High damping for pathological cases
                        "DampErr": 0.02,
                        "SOSCFHESSUP": "LBFGS",  # Switch from LSR1 to LBFGS
                    },
                }
        else:
            if self.attempt == 1:
                # Attempt 1: SlowConv keyword + increased MaxIter
                return {
                    "use_moread": True,
                    "keywords_add": ["SlowConv"],
                    "scf_block": {
                        "MaxIter": 300,
                    },
                }
            elif self.attempt == 2:
                # Attempt 2: VerySlowConv + KDIIS + high damping
                return {
                    "use_moread": True,
                    "keywords_add": ["VerySlowConv", "KDIIS"],
                    "scf_block": {
                        "MaxIter": 500,
                        "DampFac": 0.9,  # High damping factor (default: 0.7)
                        "DampErr": 0.02,  # Keep damping longer
                    },
                }
            else:
                # Attempt 3+: SOSCF (second-order) with very high damping
                return {
                    "use_moread": True,
                    "keywords_add": ["VerySlowConv"],
                    "scf_block": {
                        "MaxIter": 800,
                        "CNVSOSCF": "true",  # %scf switch; a bare SOSCF line opens a sub-block
                        "DampFac": 0.95,  # Very high damping for pathological cases
                        "DampErr": 0.001,  # Damp until very converged
                    },
                }

    def _failed_in_hessian(self) -> bool:
        """True when the SCF that failed ran inside the analytic Hessian (the CP-SCF of a frequency step)."""
        if self.output_file is None:
            return False
        try:
            path = Path(self.output_file)
            size = path.stat().st_size
            with path.open(encoding="utf-8", errors="replace") as handle:
                if size > _ERROR_SCAN_TAIL_BYTES:
                    handle.seek(size - _ERROR_SCAN_TAIL_BYTES)
                tail = handle.read()
        except OSError:
            return False
        # the analytic Hessian's own section (ORCA 6.1.1); "CP-SCF" alone also
        # stands in every TD-DFT gradient and does not tell the phase
        response = max(tail.rfind("ORCA SCF RESPONSE CALCULATION"), tail.rfind("SCF HESSIAN"))
        return response >= 0 and response > tail.rfind("SCF ITERATIONS")

    def _leanscf_convergence_fixes(self) -> Dict:
        """Fixes for an SCF that ORCA's LEANSCF program could not converge.

        In ORCA 6 every SCF runs in orca_leanscf, so "error termination in
        LEANSCF" is the ordinary SCF failure: all 472 archived ones were the
        main SCF (426 in the first optimisation cycle), none a CP-SCF of a
        frequency step.  They get the SCF escalation; FREQ is kept -- taking
        it away does not help an SCF converge and leaves IMAG and the
        thermochemistry without a Hessian.  Only an SCF that failed inside the
        analytic Hessian gets the frequency path below.
        """
        if not self._failed_in_hessian():
            return self._scf_convergence_fixes()
        return self._hessian_scf_fixes()

    def _hessian_scf_fixes(self) -> Dict:
        """Progressive fixes for LEANSCF (coupled-perturbed SCF) convergence failures.

        LEANSCF is used for analytical frequencies and can fail to converge for:
        - Highly excited deltaSCF states
        - Difficult electronic structures
        - Unstable wavefunctions

        Strategy:
        - deltaSCF/hybrid1: tighter SCF, NoDIIS/NoTRAH, FreezeAndRelease+GMF+LBFGS
        - other jobs: standard LEANSCF escalation, then skip FREQ
        """
        # Check if this is a deltaSCF (or hybrid1 deltaSCF step) calculation
        is_deltascf = self._is_deltascf()

        if is_deltascf:
            if self.attempt == 1:
                # Attempt 1: Tighter SCF convergence + increased MaxIter
                return {
                    "use_moread": True,
                    "keywords_add": ["TightSCF", "SlowConv"],
                    "scf_block": {
                        "MaxIter": 600,
                        "ConvForced": "true",  # Force convergence for better LEANSCF starting point
                    },
                }
            elif self.attempt == 2:
                # Attempt 2: Disable DIIS/TRAH for difficult deltaSCF
                return {
                    "use_moread": True,
                    "keywords_add": ["NoDIIS", "NoTRAH"],
                    "scf_block": {
                        "MaxIter": 800,
                        "ConvForced": "true",
                        "DampFac": 0.9,  # High damping
                        "DampErr": 0.01,
                    },
                }
            elif self.attempt == 3:
                # Attempt 3: FreezeAndRelease + GMF + strong damping
                logger.info("Adding FreezeAndRelease + GMF for deltaSCF LEANSCF recovery")
                return {
                    "use_moread": True,
                    "keywords_add": ["VerySlowConv", "FreezeAndRelease", "GMF"],
                    "scf_block": {
                        "MaxIter": 1000,
                        "ConvForced": "true",
                        "DampFac": 0.9,
                        "DampErr": 0.02,
                        "SOSCFHESSUP": "LBFGS",
                    },
                }
            else:
                # Attempt 4+: Skip FREQ entirely, just keep the optimized geometry
                logger.info("LEANSCF failed repeatedly - skipping FREQ calculation")
                return {
                    "use_moread": True,
                    "skip_freq": True,
                }
        else:
            if self.attempt == 1:
                # Attempt 1: Tighter SCF convergence + increased MaxIter
                return {
                    "use_moread": True,
                    "keywords_add": ["TightSCF", "SlowConv"],
                    "scf_block": {
                        "MaxIter": 400,
                        "ConvForced": "true",  # Force convergence for better LEANSCF starting point
                    },
                }
            elif self.attempt == 2:
                # Attempt 2: VeryTightSCF + high damping
                return {
                    "use_moread": True,
                    "keywords_add": ["VeryTightSCF", "VerySlowConv"],
                    "scf_block": {
                        "MaxIter": 600,
                        "ConvForced": "true",
                        "DampFac": 0.9,  # High damping
                        "DampErr": 0.01,
                    },
                }
            else:
                # Attempt 3+: Skip FREQ entirely, just keep the optimized geometry
                logger.info("LEANSCF failed repeatedly - skipping FREQ calculation")
                return {
                    "use_moread": True,
                    "skip_freq": True,
                }

    def _is_deltascf(self) -> bool:
        """Return True if input indicates a deltaSCF (including hybrid1 deltaSCF step) job."""
        if not self.parsed_input:
            return False
        keywords = [k.lower() for k in self.parsed_input.get("keywords", [])]
        if "deltascf" in keywords:
            return True
        base_block = self.parsed_input.get("blocks", {}).get("base", "")
        return "deltascf" in base_block.lower()

    def _trah_crash_fixes(self) -> Dict:
        """Fixes for TRAH segmentation faults.

        TRAH crashes are often memory-related or happen with difficult systems.
        Solution: Disable TRAH completely with NoTRAH keyword and use SlowConv for stability.

        ORCA Manual: AutoTRAH is enabled by default since ORCA 5.0.
        For segfaults, use ! NoTRAH to completely disable TRAH-SCF.
        """
        if self.attempt == 1:
            return {
                "use_moread": True,
                "keywords_add": ["SlowConv", "NoTRAH"],  # Disable TRAH completely
                "scf_block": {
                    "MaxIter": 400,
                },
            }
        else:
            # If still failing, add KDIIS and higher damping
            return {
                "use_moread": True,
                "keywords_add": ["VerySlowConv", "KDIIS", "NoTRAH"],
                "scf_block": {
                    "MaxIter": 600,
                    "DampFac": 0.9,
                },
            }

    def _diis_error_fixes(self) -> Dict:
        """Fixes for DIIS errors.

        ORCA Manual recommendations:
        - KDIIS is more robust than standard DIIS
        - SOSCF when DIIS gets stuck at ~0.001
        - TolE relaxation for difficult cases
        """
        if self.attempt == 1:
            # Try KDIIS
            return {
                "use_moread": True,
                "keywords_add": ["KDIIS", "SlowConv"],  # Simple input keywords
                "scf_block": {
                    "MaxIter": 500,
                    "TolE": 5e-6,  # Slightly relaxed from default 1e-6
                },
            }
        elif self.attempt == 2:
            # Use SOSCF with early start
            return {
                "use_moread": True,
                "scf_block": {
                    "CNVSOSCF": "true",  # %scf switch; a bare SOSCF line opens a sub-block
                    "SOSCFStart": 0.00033,  # Start SOSCF earlier
                    "MaxIter": 500,
                    "DampFac": 0.9,
                    "TolE": 1e-5,  # More relaxed
                },
            }
        else:
            # Last resort: Very aggressive damping and relaxed convergence
            return {
                "use_moread": True,
                "scf_block": {
                    "CNVSOSCF": "true",
                    "SOSCFStart": 0.001,
                    "MaxIter": 600,
                    "DampFac": 0.95,  # Strong damping
                    "DampErr": 0.1,  # Damp until error < 0.1
                    "TolE": 5e-5,  # Very relaxed
                },
            }

    def _geometry_convergence_fixes(self) -> Dict:
        """Fixes for geometry optimization convergence.

        ORCA Manual recommendations:
        - Trust < 0: FIXED trust radius (recommended for oscillating energy)
        - Trust > 0: trust radius UPDATE (adaptive)
        - Start with -0.1 and decrease if needed
        - Recalc_Hess helps with difficult PES
        """
        if self.attempt == 1:
            # Fixed small trust radius to prevent overshooting
            return {
                "use_moread": True,
                "geom_block": {
                    "Trust": -0.1,  # NEGATIVE = fixed trust radius!
                    "MaxIter": 250,
                },
            }
        elif self.attempt == 2:
            # Even smaller trust radius + recalc Hessian
            return {
                "use_moread": True,
                "geom_block": {
                    "Trust": -0.05,  # Very small fixed trust radius
                    "Recalc_Hess": 10,  # Recalc Hessian every 10 steps
                    "MaxIter": 300,
                },
            }
        else:
            # Last resort: Loose convergence criteria
            return {
                "use_moread": True,
                "geom_block": {
                    "Trust": -0.05,
                    "Recalc_Hess": 5,  # More frequent Hessian updates
                    "MaxIter": 400,
                },
                "keywords_add": ["LooseOpt"],  # Relaxed convergence
            }

    def _tddft_instability_fixes(self) -> Dict:
        """Fixes for CIS/TD-DFT aborts and collapsed TD-DFT roots.

        Both symptoms share one cause: full TD-DFT (RPA, ``tda false``) on a
        triplet-unstable closed-shell reference. The RPA eigenvalue problem
        then has a zero or negative eigenvalue, which either aborts the
        diagonalizer or yields a garbage root that the optimizer chases
        across a discontinuous surface for hundreds of cycles.

        The Tamm-Dancoff approximation drops the de-excitation block and is
        stable in exactly this situation, so TDA is the fix. Reducing cores
        is not - the failure is deterministic and reproduces exactly.
        """
        if self.attempt == 1:
            return {"tddft_block": {"tda": "true"}}
        if self.attempt == 2:
            return {
                "tddft_block": {"tda": "true"},
                "scf_block": {"Convergence": "Tight"},
            }
        # Last resort: TDA plus a larger Davidson subspace and more
        # iterations, for roots that are merely hard to converge.  MaxDim is
        # a multiplier of NRoots in ORCA, so twice ORCA's own 10 is already
        # a large space; the 50 this used to set meant 750 vectors at 15 roots.
        return {
            "tddft_block": {"tda": "true", "maxdim": 2 * ORCA_DEFAULT_MAXDIM, "maxiter": 1000},
            "scf_block": {"Convergence": "Tight"},
            "reduce_pal": 0.5,
        }

    def _mpi_crash_fixes(self) -> Dict:
        """Fixes for MPI crashes, escalating across attempts.

        MPI crashes often happen early before orbitals are written.
        DO NOT use MOREAD - it will fail with "No orbitals found".

        The escalation matters: this used to return the same modification on
        every attempt, so each retry reproduced the identical crash. Combined
        with the retry budget in orca.py that is three identical runs and
        hours of walltime for nothing.
        """
        env = {
            # Keep vader (SHM transport); only disable its CMA single-copy path,
            # which fails on locked-down kernels (ptrace_scope>=1). Falls back
            # to double-copy SHM — slower than CMA but still faster than ^vader.
            "OMPI_MCA_btl_vader_single_copy_mechanism": "none",
            "OMPI_MCA_btl_base_warn_component_unused": "0",
        }
        if self.attempt == 1:
            return {"reduce_pal": 0.5, "env_vars": env}
        if self.attempt == 2:
            return {"reduce_pal": 0.25, "env_vars": env}
        # Last resort: serial. Slow, but takes every MPI transport and every
        # parallelization bug out of the picture.
        return {"force_pal": 1, "env_vars": env}

    def _frequency_failure_fixes(self) -> Dict:
        """Fixes for frequency calculation failures."""
        if self.attempt == 1:
            return {
                "use_moread": True,
                "freq_method": "NumFreq",  # Try different method
                "scf_block": {
                    "Convergence": "VeryTight",  # Tighter for numerical derivatives
                },
            }
        else:
            return {
                "use_moread": True,
                "skip_freq": True,  # Skip frequency calculation
            }

    def _required_maxcore_mb(self) -> Optional[int]:
        """The MaxCore ORCA asked for, if it said so.

        ORCA does not merely report that memory ran out, it prints the figure
        that would have been enough:

            Error  (ORCA_SCF): Not enough memory available!
                  ====>        Please increase MaxCore to more than:  29.9 MB

        Reading that back beats halving the core count and hoping, which is
        what the fallback below does when the message is absent.
        """
        if self.output_file is None:
            return None
        try:
            size = self.output_file.stat().st_size
            with self.output_file.open(encoding="utf-8", errors="replace") as handle:
                if size > _ERROR_SCAN_TAIL_BYTES:
                    handle.seek(size - _ERROR_SCAN_TAIL_BYTES)
                tail = handle.read()
        except OSError:
            return None

        matches = re.findall(
            r"Please increase MaxCore to more than:\s*([0-9]+(?:\.[0-9]+)?)\s*MB",
            tail,
            re.IGNORECASE,
        )
        if not matches:
            return None
        # Several modules can complain in one run; the largest demand is the
        # one that has to be satisfied. A margin on top keeps the retry from
        # failing again just short of the mark.
        return int(max(float(m) for m in matches) * 1.5) + 1

    def _memory_error_fixes(self) -> Dict:
        """Fixes for memory errors.

        Memory errors often occur with SOSCF (needs Hessian) or high parallelization.
        Strategy:
        1. Raise MaxCore to what ORCA asked for, when it named a figure
        2. Disable SOSCF (most memory-intensive)
        3. Reduce cores significantly (more memory per core)
        """
        fixes: Dict = {
            "keywords_remove": ["SOSCF"],  # Disable SOSCF keyword if present
            "scf_block_remove": ["SOSCF", "CNVSOSCF", "SOSCFStart", "SOSCFConvFactor", "SOSCFMaxStep", "SOSCFHESSUP"],
            # Attempt 1 halves the cores, attempt 2 quarters them: fewer ranks
            # means more memory per rank.
            "reduce_pal": 0.5 if self.attempt == 1 else 0.25,
        }

        required = self._required_maxcore_mb()
        if required is not None:
            fixes["set_maxcore"] = required
            logger.info(
                "ORCA reported it needed more MaxCore; retrying with %d MB", required
            )
        return fixes




class OrcaInputModifier:
    """Writes the retry input for a recovery strategy.

    The input is changed through :mod:`delfin.common.orca_input`, which keeps
    every line it does not change.  Measured on ORCA 6.1.1, the dictionary
    rebuild this replaces wrote retries ORCA refused or misread:

    * one-line blocks (``%scf maxiter 125 end``, in every OCCUPIER input)
      came back nested inside a new ``%scf`` -- an input error;
    * a second ``%scf`` (``BrokenSym``) was written twice, the first lost;
    * a bare ``SOSCF`` line opens a sub-block in ORCA (manual section 2.1),
      so the settings after it were read inside it -- an input error;
    * ``Convergence TightSCF`` is not an ORCA value (``Tight`` is);
    * a metal's ``NewGTO ... end`` on its coordinate line was dropped when the
      coordinates were updated, so the retry ran the metal in another basis;
    * a multi-job file could come back as one job carrying the second job's
      ``%base`` and ``%tddft``.

    Each job of the file is kept apart; the strategy's changes go to the job
    that failed, and what is written must read back as ORCA input.
    """

    def __init__(self, inp_file: Path, config: Dict):
        """Initialize input modifier.

        Args:
            inp_file: Path to original ORCA .inp file
            config: DELFIN configuration dict
        """
        self.inp_file = Path(inp_file)
        self.config = config
        self.work_dir = self.inp_file.parent

    def apply_recovery(self, strategy: RecoveryStrategy) -> Path:
        """Apply recovery strategy and create modified input file.

        Args:
            strategy: Recovery strategy to apply

        Returns:
            Path to the retry input, or the original path when no retry can be written
        """
        try:
            text = self.inp_file.read_text(encoding="utf-8", errors="replace")
            jobs = [orca_input.parse_job(job) for job in orca_input.split_jobs(text)]
        except (OSError, orca_input.OrcaInputError) as exc:
            logger.error("Recovery cannot rewrite %s: %s", self.inp_file.name, exc)
            return self.inp_file

        target = self._failed_job(strategy, len(jobs))
        strategy.parsed_input = self._strategy_view(jobs[target])
        mods = strategy.get_modifications()
        if not mods:
            logger.warning(f"No modifications defined for {strategy.error_type.value}")
            return self.inp_file

        try:
            self._apply(jobs, target, mods)
            new_text = "".join(orca_input.render_job(items) for items in jobs)
            for job in orca_input.split_jobs(new_text):
                orca_input.parse_job(job)
        except orca_input.OrcaInputError as exc:
            logger.error("Recovery for %s would write an input ORCA cannot read (%s); not retrying",
                         self.inp_file.name, exc)
            return self.inp_file

        new_inp = self.inp_file.with_suffix(f'.retry{strategy.attempt}.inp')
        new_inp.write_text(new_text, encoding="utf-8")
        logger.info(
            f"Created recovery input: {new_inp.name} "
            f"(error={strategy.error_type.value}, attempt={strategy.attempt}, job {target + 1} of {len(jobs)})"
        )
        return new_inp

    # ------------------------------------------------------------------ which job

    @staticmethod
    def _failed_job(strategy: RecoveryStrategy, n_jobs: int) -> int:
        """Index of the job ORCA was running when it stopped (its last "JOB NUMBER" banner)."""
        known = getattr(strategy, "job_index", None)
        if known is not None:
            return max(0, min(int(known), n_jobs - 1))
        if n_jobs <= 1 or strategy.output_file is None:
            return 0
        try:
            text = Path(strategy.output_file).read_text(encoding="utf-8", errors="replace")
        except OSError:
            return 0
        numbers = re.findall(r"\${4,}\s+JOB NUMBER\s+(\d+)\s+\$", text)
        if not numbers:
            return 0
        return max(0, min(int(numbers[-1]) - 1, n_jobs - 1))

    @staticmethod
    def _strategy_view(items) -> Dict:
        """What a strategy reads of the failed job: its keywords, and its blocks and settings by name."""
        blocks: Dict[str, str] = {}
        for it in items:
            if it.kind in ("block", "setting"):
                blocks[it.name] = "".join(it.lines).strip()
        return {"keywords": orca_input.job_keywords(items), "blocks": blocks}

    def _original_stem(self) -> str:
        return re.sub(r"(\.retry\d+)+$", "", self.inp_file.stem)

    @staticmethod
    def _base(items) -> Optional[str]:
        value = orca_input.setting_value(items, "base")
        return value.strip().strip('"').strip("'") if value else None

    # ------------------------------------------------------------------ changes

    def _apply(self, jobs, target: int, mods: Dict) -> None:
        job = jobs[target]

        # A job with no %base writes under the input's name; the retry input
        # has another name, so without this its gbw/xyz/hess land beside the
        # failed run's files and every later step reads the failed ones.
        for items in jobs:
            if self._base(items) is None:
                orca_input.set_setting(items, "base", f'"{self._original_stem()}"')

        if mods.get("use_moread"):
            for index in range(target + 1):
                self._continue_from_own_state(jobs[index])

        if "scf_block" in mods:
            orca_input.set_block_values(job, "scf", mods["scf_block"])

        if "geom_block" in mods:
            orca_input.set_block_values(job, "geom", mods["geom_block"])

        esd_holders = [items for items in jobs if orca_input.has_block(items, "esd")]
        if "esd_block_remove" in mods:
            for items in esd_holders:
                if orca_input.remove_block_values(items, "esd", mods["esd_block_remove"]):
                    logger.info("ESD time grid left to ORCA: %s removed", ", ".join(mods["esd_block_remove"]))
        if "esd_block" in mods:
            for items in esd_holders:
                orca_input.set_block_values(items, "esd", mods["esd_block"], create=False)
            logger.info("Set %%esd parameters: %s", ", ".join(f"{k}={v}" for k, v in mods["esd_block"].items()))

        if "tddft_block" in mods:
            # never into a job without one: %tddft in an optimisation makes
            # ORCA optimise an excited state
            holders = [job] if orca_input.has_block(job, "tddft") else \
                [items for items in jobs if orca_input.has_block(items, "tddft")]
            for items in holders:
                orca_input.set_block_values(items, "tddft", mods["tddft_block"], create=False)
            logger.info("Set %%tddft parameters: %s",
                        ", ".join(f"{k}={v}" for k, v in mods["tddft_block"].items()))

        for items in jobs:
            if mods.get("force_pal"):
                self._set_nprocs(items, lambda n: max(1, int(mods["force_pal"])), create=True)
            elif mods.get("reduce_pal"):
                self._set_nprocs(items, lambda n: max(1, int(n * mods["reduce_pal"])), create=False)
            if mods.get("reduce_maxcore"):
                self._set_maxcore(items, lambda mb: max(1000, int(mb * mods["reduce_maxcore"])))
            if mods.get("set_maxcore"):
                self._set_maxcore(items, lambda mb: max(mb, int(mods["set_maxcore"])))

        if mods.get("skip_freq"):
            if orca_input.remove_keywords(job, ["Freq", "NumFreq", "AnFreq"]):
                logger.info("Removed frequency keywords from input")

        if mods.get("freq_method"):
            if orca_input.replace_in_family(job, mods["freq_method"], add_if_absent=False):
                logger.info("Frequencies now computed with %s", mods["freq_method"])

        if "keywords_add" in mods:
            if orca_input.add_keywords(job, mods["keywords_add"]):
                logger.info("Added keyword(s) %s to input", mods["keywords_add"])

        if "keywords_remove" in mods:
            if orca_input.remove_keywords(job, mods["keywords_remove"]):
                logger.info("Removed keyword(s): %s", mods["keywords_remove"])

        if "scf_block_remove" in mods:
            if orca_input.remove_block_values(job, "scf", mods["scf_block_remove"]):
                logger.info("Removed SCF parameter(s): %s", mods["scf_block_remove"])

    @staticmethod
    def _set_nprocs(items, new_value, *, create: bool) -> None:
        current = orca_input.block_value(items, "pal", "nprocs")
        if current is None:
            if create:
                orca_input.set_block_values(items, "pal", {"nprocs": new_value(1)})
            return
        try:
            n = int(current.split()[0])
        except (ValueError, IndexError):
            return
        value = new_value(n)
        if value != n:
            orca_input.set_block_values(items, "pal", {"nprocs": value})
            logger.info(f"PAL nprocs: {n} → {value}")

    @staticmethod
    def _set_maxcore(items, new_value) -> None:
        current = orca_input.setting_value(items, "maxcore")
        if current is None:
            return
        match = re.search(r"\d+", current)
        if not match:
            return
        mb = int(match.group(0))
        value = new_value(mb)
        if value != mb:
            orca_input.set_setting(items, "maxcore", str(value))
            logger.info("MaxCore: %d → %d MB", mb, value)

    def _continue_from_own_state(self, items) -> None:
        """Restart a job from what its failed run left: its own orbitals, and for an optimisation its geometry.

        Only the job's own files count (named by its %base).  A run that left
        none keeps the guess DELFIN gave it, or ORCA's own guess when that
        file is gone -- never the newest gbw of some other job in the folder.
        """
        base = self._base(items) or self._original_stem()
        own_gbw = self.work_dir / f"{base}.gbw"
        if own_gbw.exists() and own_gbw.stat().st_size > 100_000:
            backup = self.work_dir / f"{base}_old.gbw"
            try:
                shutil.copy2(own_gbw, backup)
            except OSError as exc:
                logger.warning("Failed to back up %s: %s", own_gbw.name, exc)
            else:
                orca_input.add_keywords(items, ["MORead"])
                orca_input.set_setting(items, "moinp", f'"{backup.name}"')
                logger.info("Continuing %s from its own orbitals (%s)", base, backup.name)
        else:
            moinp = orca_input.setting_value(items, "moinp")
            source = moinp.strip().strip('"').strip("'") if moinp else ""
            if source and not (self.work_dir / source).exists():
                orca_input.remove_keywords(items, ["MORead"])
                orca_input.remove_setting(items, "moinp")
                logger.info("%s is gone; %s starts from ORCA's own guess", source, base)

        if any(orca_input.keyword_family(k) == "optimisation level" for k in orca_input.job_keywords(items)):
            self._continue_from_own_geometry(items, self.work_dir / f"{base}.xyz")

    @staticmethod
    def _continue_from_own_geometry(items, xyz_file: Path) -> None:
        if not xyz_file.exists():
            return
        try:
            rows = xyz_file.read_text(encoding="utf-8", errors="replace").splitlines()
            n_atoms = int(rows[0].split()[0])
            atoms = []
            for row in rows[2:2 + n_atoms]:
                parts = row.split()
                atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            if len(atoms) != n_atoms:
                raise ValueError("fewer atoms than announced")
            orca_input.with_coordinates(items, atoms)
        except (ValueError, IndexError, orca_input.OrcaInputError) as exc:
            logger.info("Keeping the input's coordinates; %s is not this job's geometry (%s)", xyz_file.name, exc)
            return
        logger.info("Continuing from the geometry in %s (%d atoms)", xyz_file.name, len(atoms))


class RetryStateTracker:
    """Tracks retry attempts to prevent infinite loops and record recovery history."""

    def __init__(self, state_file: Path):
        """Initialize state tracker.

        Args:
            state_file: Path to JSON file for storing state
        """
        self.state_file = state_file
        self.state = self._load_state()

    def _load_state(self) -> Dict:
        """Load state from file."""
        if self.state_file.exists():
            try:
                with open(self.state_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                logger.warning(f"Failed to load retry state: {e}")
        return {}

    def _save_state(self):
        """Save current state to file."""
        try:
            with open(self.state_file, 'w') as f:
                json.dump(self.state, f, indent=2)
        except Exception as e:
            logger.warning(f"Failed to save retry state: {e}")

    def get_attempt(self, job_name: str, error_type: OrcaErrorType) -> int:
        """Get current retry attempt number for job/error combination.

        Args:
            job_name: Name of the job (input file stem)
            error_type: Type of error

        Returns:
            Current attempt number (0 if first attempt)
        """
        key = f"{job_name}_{error_type.value}"
        return self.state.get(key, 0)

    def increment_attempt(self, job_name: str, error_type: OrcaErrorType):
        """Increment retry counter for job/error combination.

        Args:
            job_name: Name of the job
            error_type: Type of error
        """
        key = f"{job_name}_{error_type.value}"
        self.state[key] = self.state.get(key, 0) + 1
        self._save_state()

    def should_retry(
        self,
        job_name: str,
        error_type: OrcaErrorType,
        max_retries: int = 3
    ) -> bool:
        """Check if we should attempt recovery for this job/error.

        Args:
            job_name: Name of the job
            error_type: Type of error
            max_retries: Maximum number of recovery attempts

        Returns:
            True if we should attempt recovery
        """
        return self.get_attempt(job_name, error_type) < max_retries

    def reset_job(self, job_name: str):
        """Reset all retry counts for a specific job.

        Args:
            job_name: Name of the job to reset
        """
        keys_to_remove = [k for k in self.state.keys() if k.startswith(f"{job_name}_")]
        for key in keys_to_remove:
            del self.state[key]
        self._save_state()

    def get_summary(self) -> Dict:
        """Get summary of all retry attempts.

        Returns:
            Dictionary with retry statistics
        """
        summary = {}
        for key, count in self.state.items():
            job_name, error_type = key.rsplit("_", 1)
            if job_name not in summary:
                summary[job_name] = {}
            summary[job_name][error_type] = count
        return summary
