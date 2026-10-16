"""ORCA's structured output, read through OPI, as a second witness.

DELFIN reads ORCA's energies from the text output (:mod:`delfin.energies`).
ORCA also writes the same numbers in a structured form: ``<base>.property.txt``,
which ``orca_2json <base> -property`` turns into ``<base>.property.json``,
which ORCA's Python interface OPI (``pip install "delfin[orca]"``) reads.
When both are there, this module hands over what the JSON says, so a reader
can compare it with what it took from the text -- that comparison is what
found the 1.18 eV ESD error, where the text parser took another job's
energy and the property file had the right one.

It is a witness, never the reader of record, and absence means nothing:

* without OPI, or without the JSON beside the output, :func:`structured`
  returns ``None`` and the caller carries on with the text exactly as before;
* a ``None`` inside the result means the JSON does not carry that value --
  NOT that it does not exist. Measured 2026-09-25: ``orca_2json`` from 6.1.1
  turns a 6.0.1 property file into a JSON without its energies, silently,
  and OPI's ``get_final_energy()`` then answers ``None``;
* OPI's own version gate (ORCA >= 6.1.1) is switched off: the archive holds
  6.0.1 and 6.1.0 runs, and for those a missing value reads as "no witness",
  which is the honest answer.

Measured 2026-10-15 on six HF-3c runs (ORCA 6.1.1, OPI 2.0.0): G, E and the
charge agree with the text parser to its printed precision; one parse takes
0.04 s.
"""

from __future__ import annotations

from pathlib import Path

#: The text prints G to 8 decimals; anything closer is the same number.
AGREE_EH = 5e-7

_SUFFIX = ".property.json"


def available() -> bool:
    """Whether OPI is installed (the ``[orca]`` extra)."""
    try:
        import opi.output.core  # noqa: F401
    except Exception:
        return False
    return True


def json_beside(output: Path) -> Path | None:
    """``<base>.property.json`` next to an ORCA output, if it is there."""
    output = Path(output)
    candidate = output.with_name(output.name.split(".")[0] + _SUFFIX)
    return candidate if candidate.is_file() else None


def structured(json_path: Path) -> dict | None:
    """``{"gibbs", "electronic", "charge"}`` from ORCA's property JSON.

    ``None`` when OPI is missing or cannot read the file; a ``None`` value
    when the JSON does not carry it (see the module docstring -- that is no
    evidence the value does not exist).
    """
    try:
        from opi.output.core import Output
    except Exception:
        return None
    json_path = Path(json_path)
    if not json_path.name.endswith(_SUFFIX):
        return None
    base = json_path.name[: -len(_SUFFIX)]
    try:
        out = Output(base, working_dir=json_path.parent, version_check=False)
        out.parse()
    except Exception:
        return None

    def value(getter, kind=float):
        try:
            v = getter()
        except Exception:
            return None
        return None if v is None else kind(v)

    return {"gibbs": value(out.get_free_energy),
            "electronic": value(out.get_final_energy),
            "charge": value(out.get_charge, int)}


def compare(text_value: float, witness: dict | None, energy: str) -> str:
    """What the witness says about a value read from the text, as a note.

    ``""`` when there is no witness for it; ``"confirmed by OPI"`` when the
    two agree; a WARNING naming both numbers when they do not.
    """
    if not witness or witness.get(energy) is None:
        return ""
    other = witness[energy]
    if abs(other - text_value) <= AGREE_EH:
        return "confirmed by OPI"
    return (f"WARNING: ORCA's property.json (OPI) says {other:.8f} Eh, the "
            f"output's text {text_value:.8f} Eh")
