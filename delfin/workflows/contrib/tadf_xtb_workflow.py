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
# defaults ORCA's own tools use when CONTROL does not name them.
_CONFIG_DEFAULTS: dict[str, Any] = {
    "optimize_s0": True,
    "use_bfw": False,
    "use_crest": False,
    "use_goat": False,
    "run_t1_opt": False,
    "t1_multiplicity": 3,
}


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
            "xtb_method", cfg("xtb_method") or cfg("xtb-method")
        )
        params["excited_method"] = kwargs.get(
            "excited_method", cfg("excited_method") or cfg("excited-method")
        )
        params["energy_window"] = kwargs.get("energy_window", cfg("energy_window"))
        params["cores"] = kwargs.get("cores", cfg("PAL"))
        params["maxcore"] = kwargs.get("maxcore", cfg("maxcore"))

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
