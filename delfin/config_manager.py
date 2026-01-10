# config_manager.py
# Centralized configuration management for DELFIN

from typing import Dict, Any, Optional, Union, List
from dataclasses import dataclass, field

from delfin.common.paths import resolve_path

from .config import read_control_file, get_E_ref


@dataclass
class DelfinConfig:
    """Centralized configuration management for DELFIN workflows.

    This class provides:
    - Type-safe access to configuration parameters
    - Default values and validation
    - Centralized parameter documentation
    - Easy parameter access methods
    """

    # Core settings
    NAME: str = ""
    charge: int = 0
    multiplicity: int = 1
    input_file: str = "input.txt"
    solvent: Optional[str] = None

    # Computational settings
    functional: str = "PBE0"
    main_basisset: str = "def2-SVP"
    metal_basisset: str = ""
    PAL: int = 1
    relativity: str = "none"

    # Workflow control
    method: str = "full"
    run_XTB_before_ORCA: str = "auto"
    frequency_calculation: str = "no"
    frequency_calculation_OCCUPIER: str = "no"

    # OCCUPIER selection parameters
    occupier_selection: str = "tolerance"
    occupier_precision: int = 3
    occupier_epsilon: float = 5e-4

    # Clean preference parameters (new)
    clean_override_window_h: float = 0.002
    clean_quality_improvement: float = 0.05
    clean_quality_good: float = 0.05
    clean_bias_window_h: float = 0.001
    quality_bias_window: float = 0.05

    # Energy bias parameters
    energy_bias_window_h: float = 0.001
    mismatch_bias_window: float = 0.05
    dev_similarity: float = 0.15
    bs_override_window_h: float = 0.001

    # Spin contamination bias parameters
    spin_bias_energy_window_h: float = 0.0
    spin_bias_min_gain: float = 0.003
    spin_bias_trigger_dev: float = 0.05
    spin_pair_bias_window_h: float = 0.0
    spin_pair_bias_dev_window: float = 0.20
    spin_pair_bias_min_gain: float = 0.10

    # Advanced settings
    dev_max: Optional[float] = None
    approximate_spin_projection_APMethod: Optional[Union[str, int]] = 2

    # File configuration
    out_files: Dict[int, str] = field(default_factory=dict)
    inp_files: Dict[int, str] = field(default_factory=dict)

    # Internal data
    _raw_config: Dict[str, Any] = field(default_factory=dict)
    _file_path: Optional[str] = None

    @classmethod
    def from_control_file(cls, file_path: str = "CONTROL.txt") -> 'DelfinConfig':
        """Create DelfinConfig from CONTROL.txt file.

        Args:
            file_path: Path to CONTROL.txt file

        Returns:
            DelfinConfig instance with loaded parameters

        Raises:
            FileNotFoundError: If CONTROL.txt doesn't exist
            ValueError: If configuration is invalid
        """
        control_path = resolve_path(file_path)
        if not control_path.exists():
            raise FileNotFoundError(f"CONTROL.txt not found at {file_path}")

        raw_config = read_control_file(str(control_path))
        instance = cls._from_dict(raw_config)
        instance._raw_config = raw_config
        instance._file_path = str(control_path)

        return instance

    @classmethod
    def _from_dict(cls, config_dict: Dict[str, Any]) -> 'DelfinConfig':
        """Create DelfinConfig from dictionary with type conversion."""
        kwargs = {}

        # Map config keys to class fields with type conversion
        for field_name, field_type in cls.__annotations__.items():
            if field_name.startswith('_'):
                continue

            value = config_dict.get(field_name)
            if value is not None:
                kwargs[field_name] = cls._convert_type(value, field_type)

        return cls(**kwargs)

    @staticmethod
    def _convert_type(value: Any, target_type: type) -> Any:
        """Convert value to target type with proper handling."""
        if target_type == Optional[str] or target_type == str:
            return str(value).strip() if value else None
        elif target_type == int:
            return int(value) if value is not None else 0
        elif target_type == float:
            return float(value) if value is not None else 0.0
        elif target_type == Optional[float]:
            return float(value) if value is not None else None
        else:
            return value

    def get_E_ref(self) -> float:
        """Get reference energy using the configured method."""
        return get_E_ref(self._raw_config)

    def get_file_for_index(self, idx: int, file_type: str = "output") -> str:
        """Get filename for given index and type.

        Args:
            idx: File index (1, 2, 3, ...)
            file_type: 'output' or 'input'

        Returns:
            Filename for the given index
        """
        if file_type == "output":
            files_dict = self.out_files
            default = "output.out" if idx == 1 else f"output{idx}.out"
        else:
            files_dict = self.inp_files
            default = "input.inp" if idx == 1 else f"input{idx}.inp"

        return files_dict.get(idx, default)

    def is_frequency_calculation_enabled(self, workflow: str = "main") -> bool:
        """Check if frequency calculations are enabled for a workflow."""
        if workflow == "OCCUPIER":
            return str(self.frequency_calculation_OCCUPIER).lower() == 'yes'
        else:
            return str(self.frequency_calculation).lower() == 'yes'

    def get_selection_method(self) -> tuple[str, int, float]:
        """Get OCCUPIER selection method, precision, and epsilon.

        Returns:
            Tuple of (method, precision, epsilon)
        """
        raw_sel = str(self.occupier_selection).lower()
        method = raw_sel.split('|')[0].strip()

        if method in {'rounding', 'round', 'gerundet', 'runden'}:
            method = 'rounding'
        elif method in {'tolerance', 'toleranz', 'toleranzband', 'epsilon'}:
            method = 'tolerance'
        else:
            method = 'truncation'

        precision = max(0, min(10, self.occupier_precision))
        epsilon = self.occupier_epsilon if self.occupier_epsilon > 0 else 10.0**(-precision)

        return method, precision, epsilon

    def should_run_xtb_before_orca(self, output_exists: bool = False) -> bool:
        """Determine if XTB should run before ORCA.

        Args:
            output_exists: Whether the output file already exists

        Returns:
            True if XTB should run
        """
        config_value = str(self.run_XTB_before_ORCA).lower()

        if config_value in {"yes", "true", "1", "on"}:
            return True
        elif config_value in {"no", "false", "0", "off"}:
            return False
        else:
            # Auto mode: "prep" or full automation
            return self.method.lower() in {"prep", "full"} and not output_exists

    def validate(self) -> List[str]:
        """Validate configuration and return list of issues.

        Returns:
            List of validation error messages (empty if valid)
        """
        issues = []

        # Check required files
        if not resolve_path(self.input_file).exists():
            issues.append(f"Input file '{self.input_file}' not found")

        # Validate charge and multiplicity
        if self.charge < -10 or self.charge > 10:
            issues.append(f"Charge {self.charge} seems unrealistic")

        if self.multiplicity < 1:
            issues.append(f"Multiplicity must be >= 1, got {self.multiplicity}")

        # Validate PAL
        if self.PAL < 1:
            issues.append(f"PAL must be >= 1, got {self.PAL}")

        return issues

    def __getitem__(self, key: str) -> Any:
        """Allow dict-style access for backward compatibility."""
        if hasattr(self, key):
            return getattr(self, key)
        return self._raw_config.get(key)

    def get(self, key: str, default: Any = None) -> Any:
        """Get configuration value with default (dict-style interface)."""
        if hasattr(self, key):
            value = getattr(self, key)
            return value if value is not None else default
        return self._raw_config.get(key, default)
