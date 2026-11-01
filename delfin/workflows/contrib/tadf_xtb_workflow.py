"""Workflow wrapper for TADF-xTB screening."""

from __future__ import annotations

from typing import Any, Dict, List

from delfin.workflows.registry import register


# Settings the workflow needs before it can run; without one the physics of
# the screening is undefined, so it is refused with a message, never filled
# with a silent fallback (values from the CLI parser are explicit kwargs and
# are not affected).  Multiplicity is not on the list: the ground state of
# a molecule entered for TADF screening is a singlet, so it defaults to 1.
_REQUIRED_KEYS: tuple[str, ...] = (
    "charge",
    "xtb_method",
    "excited_method",
    "energy_window",
    "cores",
    "maxcore",
)

# CONTROL.txt keys that map onto run_single_tadf_xtb parameters, with the
# defaults the shipped CONTROL template carries when CONTROL does not
# name them (delfin/define.py:269-275: tadf_xTB_preopt=none means
# optimize_s0 False; tadf_xTB_run_t1_opt=yes is the template's default
# and is looked up by its own name below, not hard-coded).
_CONFIG_DEFAULTS: dict[str, Any] = {
    "optimize_s0": False,
    "use_bfw": False,
    "use_crest": False,
    "use_goat": False,
    "t1_multiplicity": 3,
}


def _yes_no(value: Any) -> bool:
    """CONTROL yes/no strings and real booleans, both accepted."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("yes", "true", "1", "on")


class TadfXtbWorkflow:
    """TADF screening via xTB excited-state calculations."""

    name = "tadf_xtb"
    description = "TADF screening via xTB excited-state calculations"

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> Any:
        from delfin.tadf_xtb import run_single_tadf_xtb

        def cfg(key: str) -> Any:
            """CONTROL.txt lookups are case-insensitive (PAL vs pal)."""
            for variant in (key, key.lower(), key.upper()):
                if variant in config:
                    return config[variant]
            return None

        # CONTROL.txt values are the defaults; an explicit kwarg wins.
        params: Dict[str, Any] = {}
        for key, value in _CONFIG_DEFAULTS.items():
            params[key] = kwargs.get(key, cfg(key) if cfg(key) is not None else value)

        params["charge"] = kwargs.get("charge", cfg("charge"))
        params["multiplicity"] = kwargs.get("multiplicity", cfg("multiplicity"))
        if params["multiplicity"] is None:
            params["multiplicity"] = 1  # S0 ground state of the entry
        params["xtb_method"] = kwargs.get(
            "xtb_method",
            cfg("xTB_method")
            or cfg("xtb_method")
            or cfg("xtb-method"),
        )
        params["excited_method"] = kwargs.get(
            "excited_method",
            cfg("tadf_xTB_excited_method")
            or cfg("excited_method")
            or cfg("excited-method"),
        )
        params["energy_window"] = kwargs.get(
            "energy_window",
            cfg("tadf_xTB_energy_window")
            or cfg("energy_window"),
        )
        params["cores"] = kwargs.get("cores", cfg("PAL"))
        params["maxcore"] = kwargs.get("maxcore", cfg("maxcore"))
        # Template-prefixed keys (delfin/define.py:271-275) win over the
        # hard-coded defaults; an explicit kwarg wins over both.  The
        # preopt vocabulary (none|xtb|crest|goat) is the validator's own
        # (_as_preopt_mode, control_validator.py:1055-1063); the mapping
        # of a mode to the geometry-preparation switches is defined here:
        # crest -> use_crest, goat -> use_goat, xtb -> optimize_s0,
        # none -> none of them.
        preopt = kwargs.get("preopt", cfg("tadf_xTB_preopt") or cfg("preopt"))
        preopt = str(preopt).strip().lower() if preopt is not None else "none"
        params["run_t1_opt"] = _yes_no(
            kwargs.get(
                "run_t1_opt",
                cfg("tadf_xTB_run_t1_opt") or cfg("run_t1_opt") or False,
            )
        )
        params["use_bfw"] = _yes_no(
            kwargs.get("use_bfw", cfg("tadf_xTB_bfw") or cfg("use_bfw") or False)
        )
        params["use_crest"] = _yes_no(
            kwargs.get("use_crest", preopt == "crest" or params["use_crest"])
        )
        params["use_goat"] = _yes_no(
            kwargs.get("use_goat", preopt == "goat" or params["use_goat"])
        )
        params["optimize_s0"] = _yes_no(
            kwargs.get("optimize_s0", preopt == "xtb" or params["optimize_s0"])
        )

        # Entry / workdir are not CONTROL settings; they come from the caller.
        params["entry"] = kwargs.get("entry")
        if "workdir" in kwargs:
            params["workdir"] = kwargs["workdir"]

        # The refusal is per key: a required setting must come from an
        # explicit kwarg or from CONTROL.txt.  Plumbing kwargs (entry,
        # workdir, …) prove nothing about the physics settings, and one
        # explicit key does not vouch for the others.  An explicit None
        # counts as absent.  A call whose required settings are all given
        # (the CLI path passes every one) has an empty `missing` and runs;
        # anything else is refused with a message naming the settings,
        # never TypeError downstream.
        missing = [key for key in _REQUIRED_KEYS if params.get(key) is None]
        if missing:
            raise ValueError(
                f"tadf_xtb workflow is missing required setting(s): "
                f"{', '.join(missing)}. Set them in CONTROL.txt or pass them "
                "as arguments."
            )

        return run_single_tadf_xtb(**params)

    def run_cli(self, argv: List[str]) -> int:
        from delfin.tadf_xtb import run_cli

        return run_cli(argv)


register(TadfXtbWorkflow())
