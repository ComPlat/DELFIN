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

from delfin.common.logging import get_logger

logger = get_logger(__name__)


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

    DIIS_ERROR = "diis_error"
    """DIIS convergence acceleration failed."""

    TRANSIENT_SYSTEM_ERROR = "transient_system"
    """Temporary system errors (disk full, network issues, etc.) - retry without modification."""

    UNKNOWN = "unknown"
    """Unknown or unclassified error."""


class OrcaErrorDetector:
    """Analyzes ORCA output files to detect and classify errors."""

    # Error pattern definitions with priority (lower = higher priority)
    ERROR_PATTERNS = [
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
        # DIIS errors
        {
            "type": OrcaErrorType.DIIS_ERROR,
            "patterns": [
                "DIIS error",
                "Error in DIIS",
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
            with open(output_file, 'r', encoding='utf-8', errors='replace') as f:
                # Read last 10000 lines for error detection
                lines = f.readlines()
                content = '\n'.join(lines[-10000:])

            # Check for successful termination first
            if "ORCA TERMINATED NORMALLY" in content:
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

        elif self.error_type == OrcaErrorType.FREQUENCY_FAILURE:
            return self._frequency_failure_fixes()

        elif self.error_type == OrcaErrorType.MEMORY_ERROR:
            return self._memory_error_fixes()

        elif self.error_type == OrcaErrorType.TRANSIENT_SYSTEM_ERROR:
            return self._transient_error_fixes()

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

        Based on ORCA Manual recommendations:
        - SlowConv/VerySlowConv are simple input keywords (! SlowConv)
        - KDIIS as robust DIIS alternative (! KDIIS or %scf KDIIS end)
        - SOSCF when DIIS stucks at ~0.001 (%scf SOSCF end)
        - High damping (0.9+) for pathological cases
        - GMX for deltaSCF calculations (helps with excited state convergence)
        """
        # Check if this is a deltaSCF calculation
        is_deltascf = False
        if self.parsed_input and "keywords" in self.parsed_input:
            keywords_lower = [k.lower() for k in self.parsed_input["keywords"]]
            is_deltascf = "deltascf" in keywords_lower

        if self.attempt == 1:
            # Attempt 1: SlowConv keyword + increased MaxIter
            # For deltaSCF: Add GMX to help with excited state convergence
            keywords_to_add = ["SlowConv"]
            maxiter = 400 if is_deltascf else 300

            if is_deltascf and "gmx" not in [k.lower() for k in self.parsed_input.get("keywords", [])]:
                keywords_to_add.append("GMX")
                logger.info("Adding GMX keyword for deltaSCF calculation")

            return {
                "use_moread": True,
                "keywords_add": keywords_to_add,
                "scf_block": {
                    "MaxIter": maxiter,
                },
            }
        elif self.attempt == 2:
            # Attempt 2: VerySlowConv + KDIIS + high damping
            keywords_to_add = ["VerySlowConv", "KDIIS"]
            maxiter = 600 if is_deltascf else 500

            if is_deltascf and "gmx" not in [k.lower() for k in self.parsed_input.get("keywords", [])]:
                keywords_to_add.append("GMX")

            return {
                "use_moread": True,
                "keywords_add": keywords_to_add,
                "scf_block": {
                    "MaxIter": maxiter,
                    "DampFac": 0.9,  # High damping factor (default: 0.7)
                    "DampErr": 0.02,  # Keep damping longer
                },
            }
        else:
            # Attempt 3+: SOSCF (second-order) with very high damping
            keywords_to_add = ["VerySlowConv"]
            maxiter = 1000 if is_deltascf else 800

            if is_deltascf and "gmx" not in [k.lower() for k in self.parsed_input.get("keywords", [])]:
                keywords_to_add.append("GMX")

            return {
                "use_moread": True,
                "keywords_add": keywords_to_add,
                "scf_block": {
                    "MaxIter": maxiter,
                    "SOSCF": True,  # %scf block parameter
                    "DampFac": 0.95,  # Very high damping for pathological cases
                    "DampErr": 0.001,  # Damp until very converged
                },
            }

    def _leanscf_convergence_fixes(self) -> Dict:
        """Progressive fixes for LEANSCF (coupled-perturbed SCF) convergence failures.

        LEANSCF is used for analytical frequencies and can fail to converge for:
        - Highly excited deltaSCF states
        - Difficult electronic structures
        - Unstable wavefunctions

        Strategy:
        - Attempt 1: Increase main SCF convergence (tighter threshold helps LEANSCF)
        - Attempt 2: Add GMX for deltaSCF + tighten SCF further
        - Attempt 3: Skip FREQ (use MOREAD to keep geometry)
        """
        # Check if this is a deltaSCF calculation
        is_deltascf = False
        if self.parsed_input and "keywords" in self.parsed_input:
            keywords_lower = [k.lower() for k in self.parsed_input["keywords"]]
            is_deltascf = "deltascf" in keywords_lower

        if self.attempt == 1:
            # Attempt 1: Tighter SCF convergence + increased MaxIter
            # This helps LEANSCF by providing a better starting wavefunction
            keywords_to_add = ["TightSCF", "SlowConv"]
            maxiter = 500 if is_deltascf else 400

            if is_deltascf and "gmx" not in [k.lower() for k in self.parsed_input.get("keywords", [])]:
                keywords_to_add.append("GMX")
                logger.info("Adding GMX keyword for deltaSCF LEANSCF convergence")

            return {
                "use_moread": True,
                "keywords_add": keywords_to_add,
                "scf_block": {
                    "MaxIter": maxiter,
                    "ConvForced": "true",  # Force convergence for better LEANSCF starting point
                },
            }
        elif self.attempt == 2:
            # Attempt 2: VeryTightSCF + high damping + GMX for deltaSCF
            keywords_to_add = ["VeryTightSCF", "VerySlowConv"]
            maxiter = 700 if is_deltascf else 600

            if is_deltascf and "gmx" not in [k.lower() for k in self.parsed_input.get("keywords", [])]:
                keywords_to_add.append("GMX")

            return {
                "use_moread": True,
                "keywords_add": keywords_to_add,
                "scf_block": {
                    "MaxIter": maxiter,
                    "ConvForced": "true",
                    "DampFac": 0.9,  # High damping
                    "DampErr": 0.01,
                },
            }
        else:
            # Attempt 3: Skip FREQ entirely, just keep the optimized geometry
            logger.info("LEANSCF failed repeatedly - skipping FREQ calculation")
            return {
                "use_moread": True,
                "skip_freq": True,
            }

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
                    "SOSCF": True,  # %scf block parameter
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
                    "SOSCF": True,
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

    def _mpi_crash_fixes(self) -> Dict:
        """Fixes for MPI crashes.

        MPI crashes often happen early before orbitals are written.
        DO NOT use MOREAD - it will fail with "No orbitals found".
        Just reduce cores and retry from scratch.
        """
        return {
            "reduce_pal": 0.5,  # Reduce cores by half
            "env_vars": {
                "OMPI_MCA_btl": "^vader",  # Disable vader transport
                "OMPI_MCA_btl_base_warn_component_unused": "0",
            },
        }

    def _frequency_failure_fixes(self) -> Dict:
        """Fixes for frequency calculation failures."""
        if self.attempt == 1:
            return {
                "use_moread": True,
                "freq_method": "NumFreq",  # Try different method
                "scf_block": {
                    "Convergence": "VeryTightConv",  # Tighter for numerical derivatives
                },
            }
        else:
            return {
                "use_moread": True,
                "skip_freq": True,  # Skip frequency calculation
            }

    def _memory_error_fixes(self) -> Dict:
        """Fixes for memory errors.

        Memory errors often occur with SOSCF (needs Hessian) or high parallelization.
        Strategy:
        1. Disable SOSCF (most memory-intensive)
        2. Reduce cores significantly (more memory per core)
        3. Optionally increase MaxCore
        """
        if self.attempt == 1:
            return {
                "keywords_remove": ["SOSCF"],  # Disable SOSCF keyword if present
                "scf_block_remove": ["SOSCF", "SOSCFStart", "SOSCFConvFactor", "SOSCFMaxStep", "SOSCFHESSUP"],
                "reduce_pal": 0.5,  # Half the cores = double memory per core
            }
        else:
            # Attempt 2: Even fewer cores
            return {
                "keywords_remove": ["SOSCF"],
                "scf_block_remove": ["SOSCF", "SOSCFStart", "SOSCFConvFactor", "SOSCFMaxStep", "SOSCFHESSUP"],
                "reduce_pal": 0.25,  # Quarter cores = 4x memory per core
            }




class OrcaInputModifier:
    """Modifies ORCA input files by applying recovery strategies."""

    def __init__(self, inp_file: Path, config: Dict):
        """Initialize input modifier.

        Args:
            inp_file: Path to original ORCA .inp file
            config: DELFIN configuration dict
        """
        self.inp_file = inp_file
        self.config = config
        self.work_dir = inp_file.parent

    def apply_recovery(self, strategy: RecoveryStrategy) -> Path:
        """Apply recovery strategy and create modified input file.

        Args:
            strategy: Recovery strategy to apply

        Returns:
            Path to modified input file
        """
        try:
            # Parse original input first so strategy can access it
            parsed = self._parse_input()

            # Provide parsed input to strategy for context-aware modifications
            strategy.parsed_input = parsed

            mods = strategy.get_modifications()

            if not mods:
                logger.warning(f"No modifications defined for {strategy.error_type.value}")
                return self.inp_file

            # Apply modifications
            if mods.get("use_moread"):
                parsed = self._add_moread(parsed)
                # CRITICAL: Update coordinates from .xyz to match .gbw geometry
                parsed = self._update_geometry_from_xyz(parsed)

            if "scf_block" in mods:
                parsed = self._modify_scf_block(parsed, mods["scf_block"])

            if "geom_block" in mods:
                parsed = self._modify_geom_block(parsed, mods["geom_block"])

            if mods.get("reduce_pal"):
                parsed = self._reduce_parallelism(parsed, mods["reduce_pal"])

            if mods.get("reduce_maxcore"):
                parsed = self._reduce_maxcore(parsed, mods["reduce_maxcore"])

            if mods.get("skip_freq"):
                parsed = self._remove_freq(parsed)

            if "keywords_add" in mods:
                parsed = self._add_keywords(parsed, mods["keywords_add"])

            if "keywords_remove" in mods:
                parsed = self._remove_keywords(parsed, mods["keywords_remove"])

            if "scf_block_remove" in mods:
                parsed = self._remove_scf_params(parsed, mods["scf_block_remove"])

            # Create modified input file
            new_inp = self.inp_file.with_suffix(f'.retry{strategy.attempt}.inp')
            self._write_input(new_inp, parsed)

            logger.info(
                f"Created recovery input: {new_inp.name} "
                f"(error={strategy.error_type.value}, attempt={strategy.attempt})"
            )

            return new_inp

        except Exception as e:
            logger.error(f"Failed to apply recovery strategy: {e}", exc_info=True)
            return self.inp_file

    def _parse_input(self) -> Dict:
        """Parse ORCA input file into structured format.

        Returns:
            Dictionary with parsed components
        """
        with open(self.inp_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        parsed = {
            "keywords": [],
            "blocks": {},
            "blocks_before_geom": [],  # Blocks that appear before geometry
            "blocks_after_geom": [],   # Blocks that appear after geometry
            "geometry": [],
            "charge_mult": None,
            "raw_lines": lines,
            "comments": [],  # Store comment lines
        }

        # Parse comment lines and keyword line
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("#"):
                parsed["comments"].append(line.rstrip())
            elif stripped.startswith("!"):
                keywords = stripped[1:].split()
                parsed["keywords"] = keywords
                break

        # First, find geometry position
        geom_line_num = None
        for i, line in enumerate(lines):
            if line.strip().startswith("* xyz") or line.strip().startswith("* xyzfile"):
                geom_line_num = i
                break

        # Parse blocks (%pal, %scf, %geom, etc.)
        # Need to handle nested blocks (like Constraints { ... } end in %geom)
        current_block = None
        block_lines = []
        in_subblock = False  # Track if we're in a nested sub-block (like Constraints)
        block_start_line = None

        for i, line in enumerate(lines):
            stripped = line.strip()

            if stripped.startswith("%"):
                # Save previous block
                if current_block:
                    block_content = "\n".join(block_lines)
                    parsed["blocks"][current_block] = block_content
                    # Categorize: before or after geometry?
                    if geom_line_num is not None:
                        if block_start_line < geom_line_num:
                            parsed["blocks_before_geom"].append(current_block)
                        else:
                            parsed["blocks_after_geom"].append(current_block)

                # Start new block
                block_name = stripped[1:].split()[0]
                current_block = block_name
                block_lines = [stripped]
                in_subblock = False
                block_start_line = i

            elif current_block:
                block_lines.append(line.rstrip())

                # Check for start of sub-blocks (Constraints, MDCI, etc.)
                if any(keyword in stripped for keyword in ['Constraints', 'MDCI', 'ModifyInternal']):
                    in_subblock = True

                # Check for "end"
                if stripped.lower() == "end":
                    if in_subblock:
                        # This "end" closes the sub-block, not the main block
                        in_subblock = False
                    else:
                        # This "end" closes the main block
                        block_content = "\n".join(block_lines)
                        parsed["blocks"][current_block] = block_content
                        # Categorize: before or after geometry?
                        if geom_line_num is not None:
                            if block_start_line < geom_line_num:
                                parsed["blocks_before_geom"].append(current_block)
                            else:
                                parsed["blocks_after_geom"].append(current_block)
                        current_block = None
                        block_lines = []
                        in_subblock = False
                        block_start_line = None

        # Parse geometry
        in_geom = False
        for i, line in enumerate(lines):
            stripped = line.strip()

            if stripped.startswith("* xyz") or stripped.startswith("* xyzfile"):
                parsed["charge_mult"] = stripped
                in_geom = True
                continue

            if in_geom:
                if stripped == "*":
                    break
                parsed["geometry"].append(line.rstrip())

        return parsed

    def _update_geometry_from_xyz(self, parsed: Dict) -> Dict:
        """Update geometry coordinates from latest .xyz file while preserving inline basis sets.

        This is CRITICAL for geometry optimizations:
        - .gbw contains wavefunction for OPTIMIZED geometry
        - .xyz contains the OPTIMIZED coordinates
        - We must use these new coordinates with MOREAD, not old ones!

        Args:
            parsed: Parsed input dict

        Returns:
            Modified parsed dict with updated coordinates
        """
        # Find latest .xyz file
        xyz_candidates = [
            self.inp_file.with_suffix('.xyz'),  # Same name as input
            *sorted(
                self.work_dir.glob(f"{self.inp_file.stem}*.xyz"),
                key=lambda p: p.stat().st_mtime,
                reverse=True
            )[:3]
        ]

        xyz_file = None
        for candidate in xyz_candidates:
            if candidate.exists():
                xyz_file = candidate
                break

        if not xyz_file:
            logger.info("No .xyz file found, keeping original coordinates")
            return parsed

        logger.info(f"Updating coordinates from: {xyz_file.name}")

        # Read .xyz file
        try:
            with open(xyz_file, 'r') as f:
                xyz_lines = f.readlines()
        except Exception as e:
            logger.warning(f"Could not read .xyz file: {e}")
            return parsed

        # Parse .xyz format:
        # Line 0: number of atoms
        # Line 1: comment
        # Line 2+: Element X Y Z
        try:
            n_atoms = int(xyz_lines[0].strip())
            xyz_coords = []
            for i in range(2, 2 + n_atoms):
                xyz_coords.append(xyz_lines[i].strip())
        except Exception as e:
            logger.warning(f"Could not parse .xyz file: {e}")
            return parsed

        # Extract inline basis sets from original geometry
        # Map: atom_index -> list of basis set lines following that atom
        atom_basis_map = {}
        current_atom_idx = 0
        basis_lines_for_current = []

        for line in parsed["geometry"]:
            stripped = line.strip()
            # Check if this is an atom line (starts with element symbol)
            if stripped and not stripped.startswith('newgto') and not stripped.startswith('newaux'):
                # Save basis lines for previous atom
                if current_atom_idx > 0 and basis_lines_for_current:
                    atom_basis_map[current_atom_idx - 1] = basis_lines_for_current
                    basis_lines_for_current = []
                current_atom_idx += 1
            elif stripped.startswith('newgto') or stripped.startswith('newaux'):
                # This is a basis set line
                basis_lines_for_current.append(line)

        # Don't forget last atom's basis sets
        if basis_lines_for_current:
            atom_basis_map[current_atom_idx - 1] = basis_lines_for_current

        # Build new geometry with updated coordinates + preserved basis sets
        new_geometry = []
        for i, xyz_line in enumerate(xyz_coords):
            # Add coordinate line from .xyz
            new_geometry.append(xyz_line)
            # Add inline basis sets if this atom had them
            if i in atom_basis_map:
                new_geometry.extend(atom_basis_map[i])

        # Update parsed geometry
        parsed["geometry"] = new_geometry

        # CRITICAL: Change "* xyzfile ..." to "* xyz ..." since we're using inline coords
        if parsed["charge_mult"] and "xyzfile" in parsed["charge_mult"]:
            # Extract charge and multiplicity, replace xyzfile with xyz
            parts = parsed["charge_mult"].split()
            if len(parts) >= 4 and parts[0] == "*" and parts[1] == "xyzfile":
                charge = parts[2]
                mult = parts[3]
                parsed["charge_mult"] = f"* xyz {charge} {mult}"
                logger.info(f"Changed geometry format from xyzfile to inline xyz")

        logger.info(f"Updated {len(xyz_coords)} atom coordinates from .xyz, preserved inline basis sets")

        return parsed

    def _add_moread(self, parsed: Dict) -> Dict:
        """Add MOREAD keyword and %moinp block to continue from last state.

        Strategy for RETRY (this is only called during recovery, not initial run):
        - Always check if a valid GBW exists
        - If GBW exists: Add/update MOREAD to use it (continue from partial run)
        - If no GBW: Don't add MOREAD (start from scratch)

        Args:
            parsed: Parsed input dict

        Returns:
            Modified parsed dict
        """
        # Find valid GBW file to use for continuation
        gbw_file = self._find_valid_gbw()

        if not gbw_file:
            # No valid GBW found -> Remove MOREAD if present and start from scratch
            logger.info("No valid GBW found - removing MOREAD (will start from scratch)")
            parsed["keywords"] = [k for k in parsed["keywords"] if k != "MOREAD"]
            if "moinp" in parsed["blocks"]:
                del parsed["blocks"]["moinp"]
            return parsed

        # Valid GBW found -> Use it for continuation
        logger.info(f"Found valid GBW - using {gbw_file.name} for continuation")

        # Add MOREAD keyword if not present
        if "MOREAD" not in parsed["keywords"]:
            parsed["keywords"].append("MOREAD")
            logger.info("Added MOREAD keyword")

        # Add/update %moinp block
        parsed["blocks"]["moinp"] = f'%moinp "{gbw_file.name}"'
        logger.info(f"Set %moinp to {gbw_file.name}")

        return parsed

    def _find_valid_gbw(self, prefer_job_gbw: bool = True) -> Optional[Path]:
        """Find a valid GBW file to use for MOREAD.

        Priority order when prefer_job_gbw=True:
        1. GBW file with same name as input (e.g., T1.gbw for T1.inp) - highest priority
        2. Other GBW files in directory (sorted by modification time)

        Args:
            prefer_job_gbw: If True, prioritize job's own GBW over others

        Returns:
            Path to valid GBW file (backed up to _old.gbw), or None if not found
        """
        if prefer_job_gbw:
            # First check job's own GBW (from partial run)
            job_gbw = self.inp_file.with_suffix('.gbw')
            if job_gbw.exists() and job_gbw.stat().st_size > 100_000:
                gbw_backup = self.work_dir / f"{job_gbw.stem}_old.gbw"
                try:
                    shutil.copy2(job_gbw, gbw_backup)
                    logger.info(f"Found job's own GBW (partial run): {job_gbw.name} -> {gbw_backup.name}")
                    return gbw_backup
                except Exception as e:
                    logger.warning(f"Failed to backup job GBW {job_gbw.name}: {e}")

        # Check other GBW files (e.g., from prerequisite jobs like S0)
        other_gbw_candidates = sorted(
            self.work_dir.glob("*.gbw"),
            key=lambda p: p.stat().st_mtime,
            reverse=True
        )[:10]

        for gbw in other_gbw_candidates:
            if not gbw.exists():
                continue

            # Skip job's own GBW if we already checked it
            if prefer_job_gbw and gbw == self.inp_file.with_suffix('.gbw'):
                continue

            # Basic sanity checks
            size = gbw.stat().st_size
            if size < 100_000:  # Less than 100 KB is suspicious
                logger.debug(f"Skipping {gbw.name} - too small ({size} bytes)")
                continue

            # Create backup to prevent ORCA from deleting it on crash
            gbw_backup = self.work_dir / f"{gbw.stem}_old.gbw"
            try:
                shutil.copy2(gbw, gbw_backup)
                logger.info(f"Found valid GBW: {gbw.name} -> {gbw_backup.name}")
                return gbw_backup
            except Exception as e:
                logger.warning(f"Failed to backup GBW {gbw.name}: {e}")
                continue

        return None

    def _modify_scf_block(self, parsed: Dict, scf_params: Dict) -> Dict:
        """Modify or create %scf block with new parameters.

        Args:
            parsed: Parsed input dict
            scf_params: SCF parameters to set

        Returns:
            Modified parsed dict
        """
        # Parse existing %scf block or create new one
        if "scf" in parsed["blocks"]:
            scf_lines = parsed["blocks"]["scf"].split("\n")
            # Remove %scf header and end
            scf_lines = [l for l in scf_lines if l.strip() and l.strip().lower() not in ["%scf", "end"]]
        else:
            scf_lines = []

        # Build new SCF block
        new_scf = ["%scf"]

        # Track which parameters we've added
        added_params = set()

        # Update existing parameters or add new ones
        for line in scf_lines:
            # Check if this line sets a parameter we want to modify
            param_updated = False
            for param, value in scf_params.items():
                if param.lower() in line.lower():
                    param_updated = True
                    added_params.add(param)
                    break

            if not param_updated:
                new_scf.append("  " + line.strip())

        # Add new parameters that weren't in original block
        for param, value in scf_params.items():
            if param not in added_params:
                if value is True:
                    new_scf.append(f"  {param}")
                elif value is not False:
                    new_scf.append(f"  {param} {value}")

        new_scf.append("end")

        parsed["blocks"]["scf"] = "\n".join(new_scf)

        return parsed

    def _modify_geom_block(self, parsed: Dict, geom_params: Dict) -> Dict:
        """Modify or create %geom block.

        Args:
            parsed: Parsed input dict
            geom_params: Geometry parameters to set

        Returns:
            Modified parsed dict
        """
        if "geom" in parsed["blocks"]:
            geom_lines = parsed["blocks"]["geom"].split("\n")
            geom_lines = [l for l in geom_lines if l.strip() and l.strip().lower() not in ["%geom", "end"]]
        else:
            geom_lines = []

        new_geom = ["%geom"]

        added_params = set()
        for line in geom_lines:
            param_updated = False
            for param in geom_params:
                if param.lower() in line.lower():
                    param_updated = True
                    added_params.add(param)
                    break
            if not param_updated:
                new_geom.append("  " + line.strip())

        for param, value in geom_params.items():
            if param not in added_params:
                if value is True:
                    new_geom.append(f"  {param}")
                else:
                    new_geom.append(f"  {param} {value}")

        new_geom.append("end")
        parsed["blocks"]["geom"] = "\n".join(new_geom)

        return parsed

    def _reduce_parallelism(self, parsed: Dict, factor: float) -> Dict:
        """Reduce number of cores in %pal block.

        Args:
            parsed: Parsed input dict
            factor: Multiplication factor (e.g., 0.5 = half cores)

        Returns:
            Modified parsed dict
        """
        if "pal" not in parsed["blocks"]:
            return parsed

        pal_block = parsed["blocks"]["pal"]
        match = re.search(r'nprocs\s+(\d+)', pal_block, re.IGNORECASE)

        if match:
            current_nprocs = int(match.group(1))
            new_nprocs = max(1, int(current_nprocs * factor))

            new_pal = re.sub(
                r'nprocs\s+\d+',
                f'nprocs {new_nprocs}',
                pal_block,
                flags=re.IGNORECASE
            )

            parsed["blocks"]["pal"] = new_pal
            logger.info(f"Reduced PAL nprocs: {current_nprocs} → {new_nprocs}")

        return parsed

    def _reduce_maxcore(self, parsed: Dict, factor: float) -> Dict:
        """Reduce maxcore memory allocation.

        Args:
            parsed: Parsed input dict
            factor: Multiplication factor (e.g., 0.7 = 70% of original)

        Returns:
            Modified parsed dict
        """
        if "maxcore" not in parsed["blocks"]:
            return parsed

        maxcore_block = parsed["blocks"]["maxcore"]
        match = re.search(r'(\d+)', maxcore_block)

        if match:
            current_maxcore = int(match.group(1))
            new_maxcore = max(1000, int(current_maxcore * factor))

            parsed["blocks"]["maxcore"] = f"%maxcore {new_maxcore}"
            logger.info(f"Reduced maxcore: {current_maxcore} → {new_maxcore}")

        return parsed

    def _remove_freq(self, parsed: Dict) -> Dict:
        """Remove frequency calculation keywords.

        Args:
            parsed: Parsed input dict

        Returns:
            Modified parsed dict
        """
        freq_keywords = ["FREQ", "NUMFREQ", "ANFREQ"]

        original_count = len(parsed["keywords"])
        parsed["keywords"] = [
            kw for kw in parsed["keywords"]
            if kw.upper() not in freq_keywords
        ]

        if len(parsed["keywords"]) < original_count:
            logger.info("Removed frequency keywords from input")

        return parsed

    def _add_keywords(self, parsed: Dict, keywords_to_add: List[str]) -> Dict:
        """Add keywords to the input file.

        Args:
            parsed: Parsed input dict
            keywords_to_add: List of keywords to add (e.g., ["LooseOpt", "Grid5"])

        Returns:
            Modified parsed dict
        """
        for keyword in keywords_to_add:
            # Only add if not already present (case-insensitive check)
            if keyword.upper() not in [kw.upper() for kw in parsed["keywords"]]:
                parsed["keywords"].append(keyword)
                logger.info(f"Added keyword '{keyword}' to input")

        return parsed

    def _remove_keywords(self, parsed: Dict, keywords_to_remove: List[str]) -> Dict:
        """Remove keywords from the input file.

        Args:
            parsed: Parsed input dict
            keywords_to_remove: List of keywords to remove (e.g., ["SOSCF", "AutoAux"])

        Returns:
            Modified parsed dict
        """
        original_count = len(parsed["keywords"])
        keywords_upper = [kw.upper() for kw in keywords_to_remove]
        parsed["keywords"] = [
            kw for kw in parsed["keywords"]
            if kw.upper() not in keywords_upper
        ]
        removed_count = original_count - len(parsed["keywords"])
        if removed_count > 0:
            logger.info(f"Removed {removed_count} keyword(s): {keywords_to_remove}")

        return parsed

    def _remove_scf_params(self, parsed: Dict, params_to_remove: List[str]) -> Dict:
        """Remove parameters from %scf block.

        Args:
            parsed: Parsed input dict
            params_to_remove: List of parameter names to remove (e.g., ["SOSCF", "MaxDIIS"])

        Returns:
            Modified parsed dict
        """
        if "scf" not in parsed["blocks"]:
            return parsed

        scf_lines = parsed["blocks"]["scf"].split("\n")
        # Filter out block delimiters and empty lines
        scf_lines = [l for l in scf_lines if l.strip() and l.strip().lower() not in ["%scf", "end"]]

        params_upper = [p.upper() for p in params_to_remove]
        new_scf = ["%scf"]
        removed_params = []

        for line in scf_lines:
            keep_line = True
            line_stripped = line.strip()
            # Check if this line contains any of the parameters to remove
            for param in params_upper:
                # Check if parameter name appears at start of line (before space or =)
                if line_stripped.upper().startswith(param) or \
                   line_stripped.upper().startswith(param + " ") or \
                   line_stripped.upper().startswith(param + "="):
                    keep_line = False
                    removed_params.append(param)
                    break
            if keep_line:
                new_scf.append("  " + line_stripped)

        new_scf.append("end")
        parsed["blocks"]["scf"] = "\n".join(new_scf)

        if removed_params:
            logger.info(f"Removed {len(removed_params)} SCF parameter(s): {list(set(removed_params))}")

        return parsed

    def _write_input(self, path: Path, parsed: Dict):
        """Write modified input file preserving original format.

        Args:
            path: Path to write modified input
            parsed: Parsed input dict
        """
        with open(path, 'w', encoding='utf-8') as f:
            # Write comment lines (if any)
            if parsed.get("comments"):
                for comment in parsed["comments"]:
                    f.write(comment + "\n")

            # Write keywords
            if parsed["keywords"]:
                f.write("! " + " ".join(parsed["keywords"]) + "\n\n")

            # Write blocks BEFORE geometry (preserving original order)
            blocks_before = parsed.get("blocks_before_geom", [])
            all_blocks = set(parsed["blocks"].keys())

            # Priority for new blocks to insert before geometry
            new_blocks_before = ["moinp", "pal", "maxcore", "scf", "geom"]

            # Write original blocks that were before geometry
            for block_name in blocks_before:
                if block_name in parsed["blocks"]:
                    f.write(parsed["blocks"][block_name] + "\n")
                    all_blocks.discard(block_name)

            # Add new blocks (moinp, scf, etc.) if they don't exist yet
            for block_name in new_blocks_before:
                if block_name in all_blocks and block_name not in blocks_before:
                    f.write(parsed["blocks"][block_name] + "\n")
                    all_blocks.discard(block_name)

            f.write("\n")

            # Write geometry
            if parsed["charge_mult"]:
                f.write(parsed["charge_mult"] + "\n")
                for line in parsed["geometry"]:
                    f.write(line + "\n")
                # Only write closing * for inline coordinates, not for xyzfile
                if "xyzfile" not in parsed["charge_mult"]:
                    f.write("*\n")

            # Write blocks AFTER geometry (preserving original order)
            blocks_after = parsed.get("blocks_after_geom", [])
            if blocks_after:
                f.write("\n")
                for block_name in blocks_after:
                    if block_name in parsed["blocks"]:
                        f.write(parsed["blocks"][block_name] + "\n")
                        all_blocks.discard(block_name)

            # Write any remaining blocks at the end
            if all_blocks:
                f.write("\n")
                for block_name in sorted(all_blocks):
                    f.write(parsed["blocks"][block_name] + "\n")


def prepare_input_for_continuation(
    inp_file: Path,
    *,
    always_use_moread: bool = True,
    backup_gbw: bool = True,
    update_geometry: bool = True,
) -> bool:
    """Prepare ORCA input file for continuation using latest xyz and gbw.

    This function should be called BEFORE every ORCA run (not just retries!)
    to ensure the input uses the most recent geometry and wavefunction.

    Critical for:
    - RECALC mode: Always use latest state
    - Geometry optimizations: Use optimized coordinates with corresponding wavefunction
    - Any calculation that might be restarted

    Args:
        inp_file: Path to ORCA .inp file to prepare
        always_use_moread: If True, add MOREAD if xyz/gbw exist (recommended)
        backup_gbw: If True, create a backup copy of gbw (prevents ORCA deletion on failure)
        update_geometry: If True, update coordinates from latest .xyz file

    Returns:
        bool: True if preparation was successful (or no changes needed), False on error

    Example:
        >>> # Before running ORCA, always prepare the input:
        >>> prepare_input_for_continuation(Path("opt.inp"))
        >>> run_orca("opt.inp", "opt.out")
    """
    if not inp_file.exists():
        logger.debug(f"Input file {inp_file} does not exist, skipping preparation")
        return False

    work_dir = inp_file.parent
    inp_stem = inp_file.stem

    try:
        # Find latest xyz and gbw files
        xyz_candidates = [
            inp_file.with_suffix('.xyz'),
            *sorted(work_dir.glob(f"{inp_stem}*.xyz"), key=lambda p: p.stat().st_mtime, reverse=True)[:3]
        ]
        xyz_file = next((f for f in xyz_candidates if f.exists() and f.stat().st_size > 0), None)

        gbw_candidates = [
            inp_file.with_suffix('.gbw'),
            *sorted(work_dir.glob(f"{inp_stem}*.gbw"), key=lambda p: p.stat().st_mtime, reverse=True)[:3]
        ]

        # CRITICAL: Only consider GBW if it's large enough (basic sanity check)
        gbw_file = None
        for candidate in gbw_candidates:
            if candidate.exists() and candidate.stat().st_size > 100_000:
                gbw_file = candidate
                break

        # If no xyz/gbw found, nothing to prepare
        if not xyz_file and not gbw_file:
            logger.debug(f"No xyz/gbw files found for {inp_file.name}, no preparation needed")
            return True

        # Read and parse input file
        with inp_file.open('r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()

        modified = False

        # Check if original input already had MOREAD/moinp
        bang_line_idx = next((i for i, line in enumerate(lines) if line.strip().startswith('!')), None)
        original_has_moread = bang_line_idx is not None and 'MOREAD' in lines[bang_line_idx].upper()
        original_has_moinp = any(re.match(r'^\s*%moinp\s+', line, re.IGNORECASE) for line in lines)

        # Strategy: If valid GBW exists, add/update MOREAD; otherwise remove MOREAD
        if not gbw_file:
            logger.info(f"{inp_file.name}: No valid GBW found - removing MOREAD if present")
            # Remove MOREAD keyword if present
            if bang_line_idx is not None and original_has_moread:
                lines[bang_line_idx] = lines[bang_line_idx].replace(' MOREAD', '').replace('MOREAD', '')
                modified = True
            # Remove %moinp lines
            moinp_pattern = re.compile(r'^\s*%moinp\s+', re.IGNORECASE)
            new_lines = [line for line in lines if not moinp_pattern.match(line)]
            if len(new_lines) != len(lines):
                lines = new_lines
                modified = True
            always_use_moread = False
        else:
            logger.info(f"{inp_file.name}: Found valid GBW - ensuring MOREAD is set")
            always_use_moread = True

        # Step 1: Add MOREAD keyword if needed and not already present
        if always_use_moread and gbw_file and not original_has_moread:
            if bang_line_idx is not None:
                lines[bang_line_idx] = lines[bang_line_idx].rstrip() + ' MOREAD\n'
                modified = True
                logger.info(f"Added MOREAD to {inp_file.name}")

        # Step 2: Ensure %moinp block exists with correct gbw path
        if always_use_moread and gbw_file:
            # Backup GBW before modifying input
            if backup_gbw:
                gbw_backup = work_dir / f"{gbw_file.stem}_old.gbw"
                try:
                    shutil.copy2(gbw_file, gbw_backup)
                    gbw_file = gbw_backup
                    logger.info(f"Created GBW backup: {gbw_backup.name}")
                except Exception as e:
                    logger.warning(f"Failed to backup GBW {gbw_file.name}: {e}")

            # Remove any existing %moinp lines
            moinp_pattern = re.compile(r'^\s*%moinp\s+', re.IGNORECASE)
            lines = [line for line in lines if not moinp_pattern.match(line)]

            # Find bang line to insert %moinp right after it
            if bang_line_idx is not None:
                insert_idx = bang_line_idx + 1
                lines.insert(insert_idx, f'%moinp "{gbw_file.name}"\n')
                modified = True
                logger.info(f"Added %moinp pointing to {gbw_file.name}")

        # Step 3: Update geometry from xyz file while preserving inline basis sets
        if update_geometry and xyz_file:
            # Find geometry block
            geom_start = next((i for i, line in enumerate(lines) if line.strip().lower().startswith('* xyz')), None)
            if geom_start is not None:
                geom_end = next(
                    (i for i in range(geom_start + 1, len(lines)) if lines[i].strip() == '*'),
                    None
                )
                if geom_end is not None:
                    # Read xyz coordinates
                    try:
                        xyz_lines = xyz_file.read_text(encoding='utf-8', errors='replace').splitlines()
                        # Skip count/comment if present
                        xyz_start = 2 if xyz_lines and xyz_lines[0].strip().isdigit() else 0
                        xyz_coords = []
                        for line in xyz_lines[xyz_start:]:
                            parts = line.split()
                            if len(parts) >= 4:
                                try:
                                    x, y, z = map(float, parts[1:4])
                                    xyz_coords.append(f"{parts[0]:<2} {x:>16.8f} {y:>16.8f} {z:>16.8f}\n")
                                except (ValueError, IndexError):
                                    continue

                        if xyz_coords:
                            # Extract inline basis sets from current geometry
                            old_geom_lines = lines[geom_start + 1:geom_end]
                            atom_basis_map = {}  # {atom_index: [basis_lines]}
                            current_atom = -1
                            current_basis = []

                            for line in old_geom_lines:
                                stripped = line.strip().lower()
                                if stripped and not stripped.startswith('newgto') and not stripped.startswith('newaux'):
                                    # This is an atom line - save previous basis if any
                                    if current_atom >= 0 and current_basis:
                                        atom_basis_map[current_atom] = current_basis
                                    current_atom += 1
                                    current_basis = []
                                elif stripped.startswith('newgto') or stripped.startswith('newaux'):
                                    # This is a basis set line
                                    current_basis.append(line)

                            # Don't forget last atom
                            if current_atom >= 0 and current_basis:
                                atom_basis_map[current_atom] = current_basis

                            # Build new geometry with updated coords + preserved basis
                            new_geom_lines = []
                            for i, coord_line in enumerate(xyz_coords):
                                new_geom_lines.append(coord_line)
                                if i in atom_basis_map:
                                    new_geom_lines.extend(atom_basis_map[i])

                            # Replace geometry block
                            lines[geom_start + 1:geom_end] = new_geom_lines
                            modified = True

                            # CRITICAL: Change "* xyzfile ..." to "* xyz ..." since we're using inline coords
                            geom_header = lines[geom_start].strip()
                            if "xyzfile" in geom_header:
                                parts = geom_header.split()
                                if len(parts) >= 4 and parts[0] == "*" and parts[1] == "xyzfile":
                                    charge = parts[2]
                                    mult = parts[3]
                                    lines[geom_start] = f"* xyz {charge} {mult}\n"
                                    logger.info(f"Changed geometry format from xyzfile to inline xyz")

                            logger.info(f"Updated {len(xyz_coords)} coordinates from {xyz_file.name}, preserved inline basis sets")

                    except Exception as e:
                        logger.warning(f"Failed to update geometry from {xyz_file.name}: {e}")

        # Write modified input if changes were made
        if modified:
            try:
                inp_file.write_text(''.join(lines), encoding='utf-8')
                logger.info(f"Successfully prepared {inp_file.name} for continuation")
                return True
            except Exception as e:
                logger.error(f"Failed to write prepared input {inp_file.name}: {e}")
                return False
        else:
            logger.debug(f"No modifications needed for {inp_file.name}")
            return True

    except Exception as e:
        logger.error(f"Error preparing input {inp_file.name}: {e}", exc_info=True)
        return False


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
