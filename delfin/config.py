import ast
import re
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from delfin.common.control_validator import validate_control_config
from delfin.common.tddft_settings import maxdim_hint
from delfin.occupier_sequences import remove_existing_sequence_blocks
from delfin.define import TEMPLATE as CONTROL_TEMPLATE

from delfin.common.logging import get_logger

logger = get_logger(__name__)


_SEQUENCE_BLOCK_HEADER = re.compile(r"^\s*([+\-0-9,\s]+)=\[\s*$")
_TEMPLATE_DEFAULTS_CACHE: Optional[Dict[str, Any]] = None
_TEMPLATE_REQUIRED_KEYS: Set[str] = set()
_OPTION_PLACEHOLDER_DEFAULTS: Dict[str, Any] = {
    "stability_constant_mode": "auto",
    "thdy_smiles_converter": "NORMAL",
    "thdy_preopt": "xtb",
    # Left as the shipped "[GOAT|CREST]" this means no global optimiser at all,
    # which is what the two keys it replaces (XTB_GOAT=no, CREST=no) used to say.
    "global_optimizer": "",
    # The ESD lists are opt-in: shipped bracketed so the template shows the
    # shape, but an untouched template must not switch ESD on.
    "states": "",
    "ISCs": "",
    "ICs": "",
    "emission_rates": "",
}
_PLACEHOLDER_VALIDATION_VALUES: Dict[str, Any] = {
    "charge": 0,
    "solvent": "DMF",
    "method": "classic",
    "smiles_converter": "NORMAL",
    "ESD_modus": "tddft",
    "ESD_T1_opt": "uks",
}
_PLACEHOLDER_MESSAGES: Dict[str, str] = {
    # Say what goes in, not only that something has to. A message that names
    # the accepted values answers the question on the spot; one that says
    # "replace this" sends the reader looking for a manual.
    # One line each, and no semicolons: the errors are joined with "; " for
    # display, so a semicolon inside a message splits it into two entries that
    # each read as half a sentence.
    "charge": (
        "Placeholder [CHARGE] must be replaced with the total charge as a whole "
        "number, e.g. 0, -1, 2 — the charge before any redox step"
    ),
    "solvent": (
        "Placeholder [SOLVENT] must be replaced with a solvent name, e.g. water, "
        "acetonitrile, DMF, DMSO, THF, toluene — any solvent ORCA knows for "
        "CPCM/SMD"
    ),
    "method": "Placeholder [METHOD] must be set to one of: classic, manually, OCCUPIER",
    "smiles_converter": "Placeholder [SMILES_CONVERTER] must be set to one of: QUICK, NORMAL, MANTA, ARCHITECTOR",
    "stability_constant_mode": "Placeholder [STABILITY_CONSTANT_MODE] must be set to one of: auto, reaction",
    "thdy_smiles_converter": "Placeholder [THDY_SMILES_CONVERTER] must be set to one of: QUICK, NORMAL, MANTA, ARCHITECTOR",
    "thdy_preopt": "Placeholder [THDY_PREOPT] must be set to one of: none, xtb, crest, goat",
    "ESD_modus": "Placeholder [ESD_MODUS] must be set to one of: TDDFT, deltaSCF, hybrid1",
    "ESD_T1_opt": "Placeholder [ESD_T1_OPT] must be set to one of: uks, tddft",
}
_MISSING_KEY_MESSAGES: Dict[str, str] = {
    # smiles_converter is only demanded when the run actually builds a
    # structure from a SMILES, so say that rather than naming a bare key.
    "smiles_converter": (
        "This run builds its structure from a SMILES, so smiles_converter must "
        "be set to one of: QUICK, NORMAL, MANTA, ARCHITECTOR"
    ),
}
_CONTROL_KEY_ALIASES: Dict[str, str] = {
    "speciesdelta": "co2_species_delta",
    "co2speciesdelta": "co2_species_delta",
    "co2delta": "co2_species_delta",
    "co2coordination": "co2_coordination",
    "co2coordinationonoff": "co2_coordination",
    "thermodynamics": "stability_constant",
    "thermodynamicconstant": "stability_constant",
    "thermodynamicsconstant": "stability_constant",
    "thermodynconstant": "stability_constant",
    "thermodynamicsmode": "stability_constant_mode",
    "thermodynamicmode": "stability_constant_mode",
    "thermodynmode": "stability_constant_mode",
    "thermodynamicsreaction": "stability_reaction",
    "thermodynamicreaction": "stability_reaction",
    "thermodynreaction": "stability_reaction",
    "thdysmilesconverter": "thdy_smiles_converter",
    "thdypreopt": "thdy_preopt",
    "scsmilesconverter": "thdy_smiles_converter",
    "scpreopt": "thdy_preopt",
    # XTB_preOPT is the name the CONTROL file shows; XTB_OPT is what the
    # workflow code has always called it, so the new spelling maps onto the old.
    "xtbpreopt": "XTB_OPT",
    # TD-DFT has one spelling per setting.  TDDFT_TDDFT_maxiter is what the
    # template shipped for years and no job ever read; the ESD_ names are
    # older than the TDDFT section.  The canonical names are listed too, so
    # that tddft_nroots and TDDFT_NROOTS mean TDDFT_nroots instead of nothing.
    "tddftnroots": "TDDFT_nroots",
    "tddftmaxdim": "TDDFT_maxdim",
    "tddftmaxiter": "TDDFT_maxiter",
    "tddfttda": "TDDFT_TDA",
    "tddftfollowiroot": "TDDFT_followiroot",
    "tddftsoc": "TDDFT_SOC",
    "tddftadditions": "TDDFT_additions",
    "tddfttddftmaxiter": "TDDFT_maxiter",
    "esdtddftmaxiter": "TDDFT_maxiter",
    "esdnroots": "TDDFT_nroots",
    "esdmaxdim": "TDDFT_maxdim",
    "esdtda": "TDDFT_TDA",
    "esdfollowiroot": "TDDFT_followiroot",
    "esdsoc": "TDDFT_SOC",
}
_COLON_ASSIGNMENT_KEYS: Set[str] = {"co2_coordination", "co2_species_delta"}
_STRING_ONLY_KEYS: Set[str] = {"smiles"}


def _is_placeholder_value(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    stripped = value.strip()
    return bool(stripped.startswith("[") and stripped.endswith("]") and len(stripped) > 2)


def _normalize_control_key_for_alias_lookup(raw_key: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(raw_key or "").strip().lower())


def _canonicalize_control_key(raw_key: str) -> str:
    key = str(raw_key or "").strip()
    if not key:
        return key
    alias = _normalize_control_key_for_alias_lookup(key)
    return _CONTROL_KEY_ALIASES.get(alias, key)


def _strip_wrapping_quotes(value: str) -> str:
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        return value[1:-1]
    return value


def _validate_orca_override_entries(control_text: str) -> List[str]:
    """Errors in the keyword:/additions: entries, read as the run reads them."""
    from delfin.common.orca_overrides import findings

    return findings(control_text)[0]


def get_orca_override_hints(control_text: str) -> List[str]:
    """Non-blocking notes on the keyword:/additions: entries of a CONTROL text."""
    from delfin.common.orca_overrides import findings

    try:
        return findings(control_text)[1]
    except Exception:
        return []


def _collect_literal_list(lines: List[str], start_idx: int, initial_value: str) -> Tuple[List[Any], int]:
    """Collect a multi-line Python literal (list/dict) starting at start_idx."""
    buffer = initial_value + "\n"
    depth = initial_value.count('[') - initial_value.count(']')
    idx = start_idx + 1

    while idx < len(lines) and depth > 0:
        line = lines[idx]
        buffer += line
        depth += line.count('[') - line.count(']')
        idx += 1

    try:
        parsed = ast.literal_eval(buffer)
    except Exception as exc:  # noqa: BLE001
        logger.error(f"Failed to parse list literal near line {start_idx + 1}: {exc}")
        parsed = []

    return parsed if isinstance(parsed, list) else [], idx


def _parse_delta_tokens(raw: str) -> List[int]:
    deltas: List[int] = []
    seen: set[int] = set()
    for token in raw.split(','):
        text = token.strip()
        if not text:
            continue
        try:
            value = int(text)
        except ValueError:
            continue
        if value in seen:
            continue
        seen.add(value)
        deltas.append(value)
    return deltas


def _parse_sequence_block(lines: List[str], start_idx: int) -> Tuple[Dict[str, Any], int]:
    block: Dict[str, Any] = {"even_seq": [], "odd_seq": []}
    idx = start_idx

    while idx < len(lines):
        line = lines[idx]
        stripped = line.strip()

        if not stripped or stripped.startswith("#") or stripped.startswith("---") or stripped.startswith("***"):
            idx += 1
            continue
        if stripped.endswith(":"):
            idx += 1
            continue
        if stripped == "]":
            return block, idx + 1

        if stripped.startswith("even_seq"):
            parts = line.split("=", 1)
            if len(parts) != 2:
                idx += 1
                continue
            value = parts[1].strip()
            if value.startswith('[') and not value.endswith(']'):
                parsed, idx = _collect_literal_list(lines, idx, value)
            else:
                try:
                    parsed = ast.literal_eval(value)
                    idx += 1
                except Exception as exc:  # noqa: BLE001
                    logger.error(f"Failed to parse even_seq block near line {idx + 1}: {exc}")
                    parsed = []
                    idx += 1
            block["even_seq"] = parsed if isinstance(parsed, list) else []
            continue

        if stripped.startswith("odd_seq"):
            parts = line.split("=", 1)
            if len(parts) != 2:
                idx += 1
                continue
            value = parts[1].strip()
            if value.startswith('[') and not value.endswith(']'):
                parsed, idx = _collect_literal_list(lines, idx, value)
            else:
                try:
                    parsed = ast.literal_eval(value)
                    idx += 1
                except Exception as exc:  # noqa: BLE001
                    logger.error(f"Failed to parse odd_seq block near line {idx + 1}: {exc}")
                    parsed = []
                    idx += 1
            block["odd_seq"] = parsed if isinstance(parsed, list) else []
            continue

        idx += 1

    logger.warning("Unterminated OCCUPIER sequence block; reached EOF without closing ']'.")
    return block, idx


def _load_template_defaults() -> Dict[str, Any]:
    """Parse define.py TEMPLATE and cache defaults for missing CONTROL keys."""
    global _TEMPLATE_DEFAULTS_CACHE
    if _TEMPLATE_DEFAULTS_CACHE is None:
        parsed = _parse_control_file("<template>", keep_steps_literal=True, content=CONTROL_TEMPLATE)
        sanitized: Dict[str, Any] = {}
        for key, value in parsed.items():
            if isinstance(value, str) and "|" in value and not _is_placeholder_value(value):
                value = value.split("|", 1)[0].strip()

            if _is_placeholder_value(value):
                if key in _OPTION_PLACEHOLDER_DEFAULTS:
                    sanitized[key] = _OPTION_PLACEHOLDER_DEFAULTS[key]
                    continue
                _TEMPLATE_REQUIRED_KEYS.add(key)
                continue

            sanitized[key] = value

        validation_seed = dict(sanitized)
        for key in _TEMPLATE_REQUIRED_KEYS:
            if key in validation_seed:
                continue
            if key not in _PLACEHOLDER_VALIDATION_VALUES:
                raise ValueError(f"No validation fallback configured for placeholder key '{key}'")
            validation_seed[key] = _PLACEHOLDER_VALIDATION_VALUES[key]

        defaults = validate_control_config(validation_seed)
        if "_occupier_sequence_blocks" in parsed:
            defaults["_occupier_sequence_blocks"] = parsed["_occupier_sequence_blocks"]

        # Override ESD-specific defaults to empty (opt-in only)
        # Users must explicitly set these in CONTROL.txt to enable ESD calculations
        esd_opt_in_keys = ['states', 'ISCs', 'ICs', 'emission_rates']
        for key in esd_opt_in_keys:
            if key in defaults:
                defaults[key] = ''

        _TEMPLATE_DEFAULTS_CACHE = defaults

    return deepcopy(_TEMPLATE_DEFAULTS_CACHE)


_GLOBAL_OPTIMIZER_CHOICES: Set[str] = {"GOAT", "CREST"}
_GLOBAL_OPTIMIZER_OFF: Set[str] = {"NONE", "NO", "OFF", "FALSE", "0"}


def _is_yes_token(value: Any) -> bool:
    return str(value).strip().lower() in {"yes", "true", "1", "on"}


def _normalize_global_optimizer(value: Any) -> str:
    """Return ``GOAT``, ``CREST`` or ``""`` for a ``global_optimizer`` value.

    Anything unrecognised comes back as ``""``: parsing stays lenient and the
    validator is the one that tells the user what belongs there.
    """
    text = str(value or "").strip()
    if not text or _is_placeholder_value(text):
        return ""
    upper = text.upper()
    if upper in _GLOBAL_OPTIMIZER_OFF:
        return ""
    return upper if upper in _GLOBAL_OPTIMIZER_CHOICES else ""


def _apply_global_optimizer_compat(config: Dict[str, Any]) -> None:
    """Keep ``global_optimizer`` and the older ``XTB_GOAT``/``CREST`` pair in step.

    A CONTROL file now names one global optimiser instead of carrying a yes/no
    switch per program. The workflow code still reads ``XTB_GOAT`` and ``CREST``
    — and reads them by bare index — so after a parse both spellings are always
    present, whichever one the file was written with.
    """
    chosen = _normalize_global_optimizer(config.get("global_optimizer"))
    legacy_goat = _is_yes_token(config.get("XTB_GOAT"))
    legacy_crest = _is_yes_token(config.get("CREST"))

    if chosen:
        conflicting = "CREST" if chosen == "GOAT" else "XTB_GOAT"
        if (chosen == "GOAT" and legacy_crest) or (chosen == "CREST" and legacy_goat):
            logger.warning(
                "CONTROL sets global_optimizer=%s while the older %s=yes is also "
                "present; global_optimizer wins.", chosen, conflicting,
            )
        config["XTB_GOAT"] = "yes" if chosen == "GOAT" else "no"
        config["CREST"] = "yes" if chosen == "CREST" else "no"
    else:
        # No choice made (absent, emptied, or the untouched "[GOAT|CREST]"):
        # whatever the old keys say stands, and both are guaranteed to exist
        # because the workflow indexes them without a default.
        config.setdefault("XTB_GOAT", "no")
        config.setdefault("CREST", "no")
        if "global_optimizer" not in config:
            config["global_optimizer"] = (
                "GOAT" if legacy_goat else "CREST" if legacy_crest else ""
            )

    # XTB_preOPT is resolved to XTB_OPT by the alias table; guarantee it exists
    # for the same bare-index reason.
    config.setdefault("XTB_OPT", "no")


def _apply_guppy_legacy(config: Dict[str, Any]) -> None:
    """Translate legacy ``GUPPY=yes`` to ``smiles_converter=MANTA``.

    Two generations of CONTROL file said the same thing three ways: the oldest
    with a bare ``GUPPY=yes``, the next with ``smiles_converter=GUPPY``, and the
    current one with ``smiles_converter=MANTA``.  They all select the same
    builder -- since MANTA v1 the GUPPY path has called MANTA's own entry point
    and added an energy ranking on top -- so all three are read as MANTA.

    This bridge is load-bearing, not decoration: 62 of the 126 archived GUPPY
    runs are configured through the bare ``GUPPY=yes`` and would stop parsing
    without it.
    """
    guppy_val = str(config.get("GUPPY", "no")).strip().lower()
    if guppy_val != "yes":
        return
    sc = config.get("smiles_converter", "")
    if not sc or _is_placeholder_value(sc):
        config["smiles_converter"] = "MANTA"


def _control_names_a_smiles(config: Dict[str, Any]) -> bool:
    """The CONTROL file itself carries a SMILES that will have to be built."""
    value = config.get("SMILES", "")
    return bool(str(value or "").strip()) and not _is_placeholder_value(value)


def _input_file_holds_a_smiles(config: Dict[str, Any], control_path: Optional[str]) -> bool:
    """Whether the geometry source named by the CONTROL file is a SMILES.

    The SMILES module is imported lazily: it pulls in openbabel and costs about
    0.4 s, and the question only has to be answered when smiles_converter is
    missing in the first place.
    """
    if not control_path:
        return False
    try:
        entry = str(config.get("input_file") or "input.txt").strip() or "input.txt"
        source = Path(control_path).parent / entry
        if not source.is_file():
            return False
        content = source.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return False
    if not content.strip():
        return False
    try:
        from delfin.smiles_converter import is_smiles_string
    except Exception:  # noqa: BLE001
        return False
    return bool(is_smiles_string(content))


def _run_converts_smiles(config: Dict[str, Any], control_path: Optional[str] = None) -> bool:
    """True when this run turns a SMILES into coordinates.

    Only then does the converter have to be chosen. A run whose geometry comes
    from an XYZ block has nothing to convert, and demanding a converter for it
    made every CONTROL file written for such a run unusable.
    """
    return _control_names_a_smiles(config) or _input_file_holds_a_smiles(config, control_path)


def _collect_missing_required_keys(user_keys: Set[str], placeholder_keys: Set[str],
                                    esd_enabled: bool = True,
                                    converts_smiles: Any = None) -> List[str]:
    """Required keys the file does not answer.

    `converts_smiles` may be a bool or a zero-argument callable, and is
    consulted only when smiles_converter would otherwise be reported — the
    answer can cost an expensive import.
    """
    required_missing = {key for key in _TEMPLATE_REQUIRED_KEYS if key not in user_keys}
    placeholder_required = _TEMPLATE_REQUIRED_KEYS & placeholder_keys
    missing = sorted(required_missing | placeholder_required)
    # ESD_T1_opt is optional (default: uks); keep it as non-blocking hint only.
    missing = [k for k in missing if k != "ESD_T1_opt"]
    if not esd_enabled:
        missing = [k for k in missing if k != "ESD_modus"]
    if "smiles_converter" in missing and converts_smiles is not None:
        answer = converts_smiles() if callable(converts_smiles) else converts_smiles
        if not answer:
            missing = [k for k in missing if k != "smiles_converter"]
    return missing


def _raise_if_missing_required(user_keys: Set[str], placeholder_keys: Set[str],
                                esd_enabled: bool = True,
                                converts_smiles: Any = None) -> None:
    missing_all = _collect_missing_required_keys(
        user_keys, placeholder_keys, esd_enabled=esd_enabled, converts_smiles=converts_smiles,
    )
    if missing_all:
        joined = ", ".join(missing_all)
        raise ValueError(
            f"Missing required CONTROL values for: {joined}. "
            "Replace template placeholders (e.g. [CHARGE], [SOLVENT]) with actual values."
        )


def _parse_control_file(file_path: str, *, keep_steps_literal: bool, content: Optional[str] = None) -> Dict[str, Any]:
    config: Dict[str, Any] = {}
    sequence_blocks: List[Dict[str, Any]] = []

    if content is None:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    else:
        lines = content.splitlines(keepends=True)

    idx = 0
    total_lines = len(lines)

    while idx < total_lines:
        line = lines[idx]
        stripped = line.strip()

        upper = stripped.upper()
        # Stop parsing before the documentation block at the end of CONTROL.txt.
        # Require the exact "ESD MODULE:" heading so we don't accidentally stop
        # at the real ESD configuration section ("ESD module (excited state dynamics):").
        if upper.startswith("INFOS:") or upper.startswith("ESD MODULE:"):
            # Stop parsing before template/documentation blocks to avoid overriding user values.
            break

        # Skip comments / separators / blanks
        if not stripped or stripped.startswith('#') or stripped.startswith('---') or stripped.startswith('***'):
            idx += 1
            continue

        # Sequence block definition (e.g. "-1,0,+1=[")
        block_match = _SEQUENCE_BLOCK_HEADER.match(stripped)
        if block_match:
            deltas = _parse_delta_tokens(block_match.group(1))
            block, next_idx = _parse_sequence_block(lines, idx + 1)
            if deltas and (block.get("even_seq") or block.get("odd_seq")):
                block["deltas"] = deltas
                sequence_blocks.append(block)
            idx = next_idx
            continue

        key = ""
        value = ""

        if '=' in line:
            key_raw, value_raw = line.split('=', 1)
            key = _canonicalize_control_key(key_raw)
            value = value_raw.strip()
        elif ':' in stripped:
            key_raw, value_raw = line.split(':', 1)
            key_candidate = _canonicalize_control_key(key_raw)
            value_candidate = value_raw.strip()

            # Most ":" lines are section headings. Only parse explicit aliases
            # we support as data keys (e.g. "Species delta: -2").
            if key_candidate not in _COLON_ASSIGNMENT_KEYS or not value_candidate:
                idx += 1
                continue
            key = key_candidate
            value = value_candidate
        else:
            idx += 1
            continue

        if keep_steps_literal and key in ('oxidation_steps', 'reduction_steps'):
            config[key] = value
            idx += 1
            continue

        if _normalize_control_key_for_alias_lookup(key) in _STRING_ONLY_KEYS:
            # Keep string-like keys as raw text so values like SMILES=[Ni](N)(N)
            # are not mistaken for multi-line list literals.
            config[key] = _strip_wrapping_quotes(value)
            idx += 1
            continue

        if value.startswith('[') and not value.endswith(']'):
            parsed, next_idx = _collect_literal_list(lines, idx, value)
            config[key] = parsed
            idx = next_idx
            continue

        if ',' in value and not value.startswith('{') and not value.startswith('['):
            config[key] = [v.strip() for v in value.split(',') if v.strip()]
            idx += 1
            continue

        if key == "charge":
            config[key] = value
            idx += 1
            continue

        try:
            config[key] = ast.literal_eval(value)
        except Exception:
            config[key] = value
        idx += 1

    if sequence_blocks:
        config["_occupier_sequence_blocks"] = sequence_blocks
    _apply_global_optimizer_compat(config)
    return config


def parse_control_text(control_text: str, *, keep_steps_literal: bool = True) -> Dict[str, Any]:
    """Parse CONTROL.txt content from a string without full validation."""
    return _parse_control_file("<in-memory>", keep_steps_literal=keep_steps_literal, content=control_text)


def set_control_value(control_text: str, key: str, value: Any) -> str:
    """Return `control_text` with `key` set to `value`, appending it if absent.

    The one place that rewrites a CONTROL line in text form, so the dashboard
    editor and the derived templates cannot disagree about what counts as a
    key line.
    """
    pattern = rf"(?m)^{re.escape(str(key))}\s*=.*$"
    replacement = f"{key}={value}"
    if re.search(pattern, control_text):
        return re.sub(pattern, lambda _m: replacement, control_text)
    return control_text.rstrip() + f"\n{replacement}\n"


def read_control_file(file_path: str) -> Dict[str, Any]:
    """Parse CONTROL.txt file and return configuration dictionary.

    Supports:
    - Key=value pairs with type inference
    - Multi-line lists in [...] format
    - Comma-separated values converted to lists
    - Comments starting with # or --- or ***

    Args:
        file_path: Path to CONTROL.txt file

    Returns:
        Dictionary containing parsed configuration parameters
    """
    sanitized = remove_existing_sequence_blocks(Path(file_path), persist=False)
    text = sanitized or Path(file_path).read_text(encoding="utf-8")
    config = _parse_control_file(file_path, keep_steps_literal=True, content=text)
    _apply_guppy_legacy(config)
    override_errors = _validate_orca_override_entries(text)
    if override_errors:
        raise ValueError("; ".join(override_errors))
    user_keys = set(config.keys())
    placeholder_keys = {key for key, value in config.items() if _is_placeholder_value(value)}
    defaults = _load_template_defaults()
    esd_enabled = str(config.get("ESD_modul", "no")).strip().lower() == "yes"
    _raise_if_missing_required(
        user_keys, placeholder_keys, esd_enabled=esd_enabled,
        converts_smiles=lambda: _run_converts_smiles(config, file_path),
    )
    validated = validate_control_config(config)
    if "_occupier_sequence_blocks" in config:
        validated["_occupier_sequence_blocks"] = config["_occupier_sequence_blocks"]
        user_keys.add("_occupier_sequence_blocks")

    merged = dict(validated)
    for key, value in defaults.items():
        if key not in user_keys:
            merged[key] = value

    return merged


def validate_control_text(control_text: str, *, converts_smiles: Optional[bool] = None) -> List[str]:
    """Validate CONTROL.txt content and return a list of errors (empty if valid).

    Collects ALL validation errors instead of stopping at the first one.
    This includes placeholder errors AND field validation errors.

    `converts_smiles` says whether this run will build a structure from a
    SMILES, which decides whether smiles_converter has to be chosen. There is
    no file to look at here, so a caller that knows — the Submit tab knows what
    is in its input box — should say; otherwise only the CONTROL file's own
    SMILES key can answer.
    """
    all_errors: List[str] = []

    try:
        _load_template_defaults()
    except Exception as exc:
        return [f"Failed to load template defaults: {exc}"]

    try:
        config = _parse_control_file("<in-memory>", keep_steps_literal=True, content=control_text)
    except Exception as exc:
        return [f"Failed to parse CONTROL content: {exc}"]

    _apply_guppy_legacy(config)
    all_errors.extend(_validate_orca_override_entries(control_text))

    user_keys = set(config.keys())
    placeholder_keys = {key for key, value in config.items() if _is_placeholder_value(value)}

    # Collect placeholder errors (don't raise, just collect)
    esd_enabled = str(config.get('ESD_modul', 'no')).strip().lower() == 'yes'
    missing_required = _collect_missing_required_keys(
        user_keys, placeholder_keys,
        converts_smiles=(
            converts_smiles if converts_smiles is not None
            else lambda: _run_converts_smiles(config)
        ),
    )
    for key in missing_required:
        # ESD settings only matter when ESD_modul=yes
        if key in {'ESD_modus', 'ESD_T1_opt'} and not esd_enabled:
            continue
        if key in placeholder_keys:
            msg = _PLACEHOLDER_MESSAGES.get(
                key, f"Placeholder [{key.upper()}] must be replaced with an actual value"
            )
            all_errors.append(msg)
        else:
            all_errors.append(_MISSING_KEY_MESSAGES.get(key, f"Missing required key: {key}"))

    # Replace placeholders with valid dummy values for further validation
    config_for_validation = dict(config)
    for key in placeholder_keys:
        if key in _PLACEHOLDER_VALIDATION_VALUES:
            config_for_validation[key] = _PLACEHOLDER_VALIDATION_VALUES[key]

    # Run full field validation to catch other errors (invalid functionals, solvents, etc.)
    try:
        validate_control_config(config_for_validation)
    except ValueError as exc:
        # Split by semicolon and add each error
        for err in str(exc).split(";"):
            err = err.strip()
            if err:
                all_errors.append(err)

    return all_errors


def _esd_parse_transitions(val: Any) -> List[str]:
    """The entries of an ESD list, read as the ESD module reads them (brackets or not)."""
    from delfin.common.control_validator import esd_list_items

    return esd_list_items(val)


def get_esd_hints(control_text: str) -> List[str]:
    """Return advisory hints when ESD_modul=yes but ESD fields need attention.

    These are non-blocking — submission proceeds normally.
    """
    try:
        config = parse_control_text(control_text)
    except Exception:
        return []
    if str(config.get('ESD_modul', 'no')).lower() != 'yes':
        return []
    hints = []
    method_val = str(config.get('method', '')).strip().lower()
    if method_val in ('occupier', 'manually'):
        hints.append(
            f"ESD_modul=yes with method={config.get('method')} — "
            f"ESD is only supported with method=classic"
        )
    esd_fields = {
        'states': '[S1,T1,S2,T2]',
        'ISCs': '[S1>T1,T1>S1]',
        'ICs': '[S1>S0]',
        'emission_rates': '[f,p]',
    }
    for field, example in esd_fields.items():
        val = config.get(field)
        if val is None:
            hints.append(
                f"ESD_modul=yes: {field} not set — "
                f"define if needed (e.g. {field}={example[1:-1]}), or clear with {field}= or {field}=[]"
            )

    # Say what will be computed.  The template's example lists are read like
    # any others -- brackets included -- so a file that only switched ESD on
    # computes them.
    planned = []
    for field, label in (('states', 'states'), ('ISCs', 'ISC'), ('ICs', 'IC'), ('emission_rates', 'emission')):
        items = _esd_parse_transitions(config.get(field))
        if field == 'states':
            items = ['S0'] + [i.upper() for i in items if i.upper() != 'S0']
        if items:
            planned.append(f"{label} {','.join(items)}")
    hints.append("ESD_modul=yes computes: " + "; ".join(planned))

    states_upper = {s.upper() for s in _esd_parse_transitions(config.get('states'))}
    if "T1" in states_upper:
        t1_opt = config.get('ESD_T1_opt')
        if t1_opt is None or _is_placeholder_value(t1_opt):
            hints.append(
                "ESD_modul=yes with states containing T1: ESD_T1_opt not set — "
                "set ESD_T1_opt=uks or ESD_T1_opt=tddft."
            )

    # Warn if ISCs contain same-spin transitions (those belong in ICs)
    for trans in _esd_parse_transitions(config.get('ISCs')):
        if '>' in trans:
            a, b = trans.split('>', 1)
            a_spin = a.strip()[:1].upper()
            b_spin = b.strip()[:1].upper()
            if a_spin and b_spin and a_spin == b_spin:
                hints.append(
                    f"ISCs={trans!r}: same-spin transition (S→S or T→T) — "
                    f"did you mean ICs={trans}? ISCs connect different spin states (S↔T)"
                )

    # Warn if ICs contain different-spin transitions (those belong in ISCs)
    for trans in _esd_parse_transitions(config.get('ICs')):
        if '>' in trans:
            a, b = trans.split('>', 1)
            a_spin = a.strip()[:1].upper()
            b_spin = b.strip()[:1].upper()
            if a_spin and b_spin and a_spin != b_spin:
                hints.append(
                    f"ICs={trans!r}: different-spin transition (S→T or T→S) — "
                    f"did you mean ISCs={trans}? ICs connect same spin states (S→S, T→T)"
                )

    # An IC ORCA cannot compute is skipped at run time; say so before it is.
    # Read the way the ESD module reads it: brackets or not.
    from delfin.common.control_validator import unsupported_ic_reason
    raw_ics = config.get('ICs')
    ic_items = raw_ics if isinstance(raw_ics, list) else str(raw_ics or '').strip().strip('[]').split(',')
    for trans in [str(t).strip() for t in ic_items if str(t).strip()]:
        reason = unsupported_ic_reason(trans)
        if reason and 'ISC' not in reason and 'not a transition' not in reason:
            hints.append(f"ICs={trans!r}: {reason}; it will be skipped")

    # An explicit MaxDim outside ORCA's range usually means the old reading
    # of it as an absolute size (the template shipped 30).
    maxdim_note = maxdim_hint(config)
    if maxdim_note:
        hints.append(maxdim_note)

    return hints


def validate_control_file(file_path: str) -> List[str]:
    """Validate CONTROL.txt file and return a list of errors (empty if valid)."""
    try:
        read_control_file(file_path)
    except FileNotFoundError:
        return [f"CONTROL file not found: {file_path}"]
    except ValueError as exc:
        return [err.strip() for err in str(exc).split(";") if err.strip()]
    return []

def OCCUPIER_parser(path: str) -> Dict[str, Any]:
    """Parse OCCUPIER-specific configuration file.

    Similar to read_control_file but with specialized handling for OCCUPIER workflow.

    Args:
        path: Path to configuration file

    Returns:
        Dictionary containing parsed OCCUPIER configuration
    """
    sanitized = remove_existing_sequence_blocks(Path(path), persist=False)
    text = sanitized or Path(path).read_text(encoding="utf-8")
    config = _parse_control_file(path, keep_steps_literal=False, content=text)
    _apply_guppy_legacy(config)
    override_errors = _validate_orca_override_entries(text)
    if override_errors:
        raise ValueError("; ".join(override_errors))
    user_keys = set(config.keys())
    placeholder_keys = {key for key, value in config.items() if _is_placeholder_value(value)}
    defaults = _load_template_defaults()
    esd_enabled = str(config.get("ESD_modul", "no")).strip().lower() == "yes"
    _raise_if_missing_required(
        user_keys, placeholder_keys, esd_enabled=esd_enabled,
        converts_smiles=lambda: _run_converts_smiles(config, path),
    )
    validated = validate_control_config(config)
    if "_occupier_sequence_blocks" in config:
        validated["_occupier_sequence_blocks"] = config["_occupier_sequence_blocks"]
        user_keys.add("_occupier_sequence_blocks")

    merged = dict(validated)
    for key, value in defaults.items():
        if key not in user_keys:
            merged[key] = value

    return merged


def _coerce_float(val: Any) -> Optional[float]:
    """Convert various types to float with robust error handling.

    Handles:
    - Integers and floats
    - String representations (including comma as decimal separator)
    - Boolean values (returns None)
    - Infinity and NaN checks

    Args:
        val: Value to convert to float

    Returns:
        Float value or None if conversion fails
    """
    if val is None:
        return None
    if isinstance(val, bool):
        return None
    if isinstance(val, (int, float)):
        try:
            from math import isfinite
            f = float(val)
            return f if isfinite(f) else None
        except Exception:
            return None
    if isinstance(val, str):
        s = val.strip()
        if not s:
            return None
        s = s.replace(",", ".")
        try:
            return float(s)
        except ValueError:
            return None
    return None


def _is_oniom_calculation(config: Dict[str, Any]) -> bool:
    """Check if this is an ONIOM (QM/QM2) calculation.

    Detects ONIOM by checking:
    1. If input_file exists and contains ONIOM markers
    2. If any generated ORCA input files contain QM/QM2 keywords

    Args:
        config: Configuration dictionary

    Returns:
        True if ONIOM calculation is detected, False otherwise
    """
    import os
    from pathlib import Path

    # QM/MM method patterns to detect
    qmmm_patterns = ['QM/XTB', 'QM/MM', 'QM/QM2', 'QM/PBEH-3C', 'QM/HF-3C', 'QM/r2SCAN-3C']

    # Check input_file specified in config
    input_file = config.get('input_file', 'input.txt')
    if os.path.exists(input_file):
        try:
            with open(input_file, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read(2000)  # Read first 2000 chars
                # Check for ONIOM markers
                if any(pattern in content for pattern in qmmm_patterns):
                    return True
        except Exception:
            pass

    # Check generated ORCA input files (initial.inp, etc.)
    input_files_to_check = [
        'initial.inp',
        'ox_step_1.inp',
        'red_step_1.inp',
    ]

    for fname in input_files_to_check:
        if os.path.exists(fname):
            try:
                with open(fname, 'r', encoding='utf-8', errors='ignore') as f:
                    first_line = f.readline()
                    # ORCA input files start with "! keywords"
                    if any(pattern in first_line for pattern in qmmm_patterns):
                        return True
            except Exception:
                pass

    return False


def get_E_ref(config: Dict[str, Any]) -> float:
    """Get reference electrode potential for redox calculations.

    Returns user-specified E_ref if available, otherwise looks up
    solvent-specific reference potentials vs. SHE.

    For ONIOM calculations (QM/QM2), automatically uses adjusted E_ref values
    to account for systematic DFT errors in electrostatic interactions.

    Args:
        config: Configuration dictionary containing 'E_ref' and 'solvent'

    Returns:
        Reference electrode potential in V vs. SHE (default: 4.345 V)
    """
    e_ref = _coerce_float(config.get('E_ref', None))
    if e_ref is not None:
        return e_ref

    solvent_raw = config.get('solvent', '')
    solvent_key = solvent_raw.strip().lower() if isinstance(solvent_raw, str) else ''

    # Check if this is an ONIOM calculation
    is_oniom = _is_oniom_calculation(config)

    if is_oniom:
        # ONIOM-specific E_ref values (can be customized per solvent)
        solvent_E_ref_oniom = {
            "dmf": -3.31, "n,n-dimethylformamide": -3.31,
            "dcm": -3.31, "ch2cl2": -3.31, "dichloromethane": -3.31,
            "acetonitrile": -3.31, "mecn": -3.31,
            "thf": -3.31, "tetrahydrofuran": -3.31,
            "dmso": -3.31, "dimethylsulfoxide": -3.31,
            "dme": -3.31, "dimethoxyethane": -3.31,
            "acetone": -3.31, "propanone": -3.31,
        }

        e_ref_value = solvent_E_ref_oniom.get(solvent_key, -3.31)

        # Log the automatic adjustment
        logger.info(f"ONIOM calculation detected: Using E_ref = {e_ref_value:.3f} V for {solvent_raw or 'default solvent'}")
        logger.info(f"  → This accounts for systematic DFT errors in QM/QM2 electrostatic interactions")
        logger.info(f"  → To override, set E_ref manually in CONTROL.txt")

        return e_ref_value
    else:
        # Standard E_ref values for non-ONIOM calculations
        solvent_E_ref = {
            "dmf": 4.795, "n,n-dimethylformamide": 4.795,
            "dcm": 4.805, "ch2cl2": 4.805, "dichloromethane": 4.805,
            "acetonitrile": 4.745, "mecn": 4.745,
            "thf": 4.905, "tetrahydrofuran": 4.905,
            "dmso": 4.780, "dimethylsulfoxide": 4.780,
            "dme": 4.855, "dimethoxyethane": 4.855,
            "acetone": 4.825, "propanone": 4.825,
        }

        return solvent_E_ref.get(solvent_key, 4.345)
