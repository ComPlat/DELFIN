"""What CONTROL.txt says about TD-DFT, as the ``%tddft`` block ORCA reads.

DELFIN writes a ``%tddft`` block into every ORCA job with a TD-DFT part: the
S0 absorption check, each excited-state optimisation, the check jobs of the
deltaSCF and hybrid1 modes, and the ESD(ISC/IC/FLUOR/PHOSP) rate jobs.  Each
of those writers used to read CONTROL on its own, under five spellings
(``TDDFT_*``, ``ESD_*``, bare ``TDA``, per-job ``ESD_ISC_NROOTS`` ...), so a
``TDDFT_*`` key reached some jobs and not others.  In deltaSCF mode and in
every rate job ``TDDFT_nroots``, ``TDDFT_TDA`` and ``TDDFT_SOC`` did nothing,
and ``TDDFT_TDDFT_maxiter`` -- the spelling the template shipped -- was read by
no job at all.  This module is the one place that reads them, and
:func:`tddft_block` the one place that writes the block.

**MaxDim is a multiplier.**  ORCA sizes the Davidson expansion space as
``MaxDim x NRoots`` (manual 6.1.1, Table 5.8).  Measured on 6.0.1 and 6.1.1:
``maxdim 3`` with 5 roots prints "Maximum size of the expansion space ... 15",
and without a MaxDim line ORCA uses 10 for 3, 5, 15 and 30 roots alike.  Below
2 it raises the value to 2 itself and rebuilds the space every few iterations.
So ``TDDFT_maxdim=auto`` writes ORCA's own 10, the upper end of the 5-10 range
the manual recommends.  The old template value 30 asked for three times that
(450 vectors at 15 roots, same iteration count on benzene), and the old unset
fallback ``max(5, nroots/2)`` was written as if MaxDim were an absolute size,
which made the space grow with the square of NRoots.

**What DELFIN owns and what passes through.**  NRoots, MaxDim, MaxIter, TDA,
FollowIRoot and DoSOC have keys of their own because the same value has to
reach every job and DELFIN reads it back.  IRoot, IRootMult, Triplets,
SRoot/TRoot/TRootSSL, NACME and ETF are decided per job by the ESD module.
Every other ``%tddft`` keyword ORCA knows -- ETol, RTol, NGuessMat, DoNTO,
NTOStates, EWin, OrbWin, TPrint, CPCMEQ, SOCGrad ... -- goes through
``TDDFT_additions`` verbatim into every block.  That is how CONTROL reaches all
of ORCA's TD-DFT without DELFIN having to know each keyword, and without two
spellings for one setting.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, List, Mapping, Optional, Sequence, Tuple, Union

__all__ = [
    'ORCA_DEFAULT_MAXDIM', 'ORCA_RECOMMENDED_MAXDIM', 'FROM_CONTROL',
    'TddftSettings', 'read_settings', 'tddft_block', 'job_nroots',
    'coerce_bool', 'coerce_maxdim', 'coerce_maxiter', 'coerce_additions',
    'maxdim_hint',
]

#: MaxDim ORCA uses when the input names none; measured on 6.0.1 and 6.1.1.
ORCA_DEFAULT_MAXDIM = 10

#: The range the ORCA manual recommends for MaxDim, in units of NRoots.
ORCA_RECOMMENDED_MAXDIM = (5, 10)

#: Marker for "take this from CONTROL", as opposed to None ("leave it to ORCA").
FROM_CONTROL: Any = object()

_TRUE = {'yes', 'y', 'true', '1', 'on'}
_FALSE = {'no', 'n', 'false', '0', 'off'}

#: Keywords with a CONTROL key of their own.  Naming one in TDDFT_additions
#: would give a setting two spellings, and the later one would silently win.
_OWN_KEYS = {
    'nroots': 'TDDFT_nroots',
    'maxdim': 'TDDFT_maxdim',
    'maxiter': 'TDDFT_maxiter',
    'tda': 'TDDFT_TDA',
    'followiroot': 'TDDFT_followiroot',
    'dosoc': 'TDDFT_SOC',
}

#: Keywords the ESD module decides per job.  One global value would point
#: every job at the same root, or give a singlet job triplets.
_PER_JOB_KEYS = frozenset({
    'iroot', 'irootmult', 'irootlist', 'triplets',
    'sroot', 'troot', 'trootssl', 'nacme', 'etf',
})

_KEYWORD_RE = re.compile(r'^[A-Za-z][A-Za-z0-9_]*(\[\d+\])?$')


def coerce_bool(value: Any, key: str = 'value') -> bool:
    """true/false in any of the spellings CONTROL uses; anything else is an error.

    Stricter than the validator's yes/no reading on purpose: there a typo
    means "no", and ``TDDFT_TDA=TURE`` would switch on full TD-DFT silently.
    """
    if isinstance(value, bool):
        return value
    text = str(value if value is not None else '').strip().lower()
    if text in _TRUE:
        return True
    if text in _FALSE:
        return False
    raise ValueError(f"{key} must be true or false (also yes/no, on/off, 1/0), got {value!r}")


def coerce_maxdim(value: Any) -> Union[str, int]:
    """``auto`` or a positive integer -- the multiplier ORCA reads."""
    text = str(value if value is not None else '').strip().lower()
    if text in ('', 'auto'):
        return 'auto'
    try:
        number = int(text)
    except ValueError:
        raise ValueError(
            "TDDFT_maxdim must be auto or a positive whole number; ORCA multiplies "
            f"it by nroots to size the Davidson space (auto = {ORCA_DEFAULT_MAXDIM}), "
            f"got {value!r}") from None
    if number < 1:
        raise ValueError(f"TDDFT_maxdim must be auto or a positive whole number, got {value!r}")
    return number


def coerce_maxiter(value: Any) -> Union[str, int]:
    """A positive integer, or ``auto``/empty for ORCA's own limit."""
    text = str(value if value is not None else '').strip().lower()
    if text in ('', 'auto'):
        return 'auto'
    try:
        number = int(text)
    except ValueError:
        raise ValueError(
            "TDDFT_maxiter must be auto or a positive whole number "
            f"(auto leaves ORCA's own limit), got {value!r}") from None
    if number < 1:
        raise ValueError(f"TDDFT_maxiter must be auto or a positive whole number, got {value!r}")
    return number


def _addition_lines(value: Any) -> List[str]:
    """Split TDDFT_additions into ORCA lines.

    The CONTROL reader turns a value with a comma into a list at the commas,
    and ORCA values have commas of their own (``OrbWin[0] 2,-1,-1,14``), so a
    list is glued back together at the commas before splitting at ``;``.
    """
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        text = ','.join(str(part).strip() for part in value)
    else:
        text = str(value)
    return [' '.join(line.split()) for line in text.split(';') if line.strip()]


def coerce_additions(value: Any) -> str:
    """Validate TDDFT_additions and return it as ``line; line; ...``."""
    lines = _addition_lines(value)
    for line in lines:
        keyword = line.split()[0].split('=')[0]
        if keyword.lower() in ('end', '$new_job') or keyword[:1] in ('%', '*', '!', ''):
            raise ValueError(
                f"TDDFT_additions takes %tddft keywords, not blocks or input lines: {line!r}. "
                "Separate keywords with ';', e.g. TDDFT_additions=DoNTO true; ETol 1e-7")
        base = re.sub(r'\[\d+\]$', '', keyword.lower())
        if base in _OWN_KEYS:
            raise ValueError(
                f"TDDFT_additions must not set {keyword}; use {_OWN_KEYS[base]}, "
                "so the value reaches every TD-DFT job exactly once")
        if base in _PER_JOB_KEYS:
            raise ValueError(
                f"TDDFT_additions must not set {keyword}; DELFIN sets it per job "
                "from states, ISCs, ICs and emission_rates")
        if not _KEYWORD_RE.match(keyword):
            raise ValueError(f"TDDFT_additions: {keyword!r} is not an ORCA keyword name")
    return '; '.join(lines)


def _first(config: Mapping[str, Any], *keys: str) -> Any:
    """The first key that is set, skipping empties, so a blank never shadows a legacy value."""
    for key in keys:
        value = config.get(key)
        if value is not None and str(value).strip() != '':
            return value
    return None


@dataclass(frozen=True)
class TddftSettings:
    nroots: int
    maxdim: int
    maxdim_auto: bool
    maxiter: Optional[int]
    tda: bool
    followiroot: bool
    soc: bool
    additions: Tuple[str, ...]

    @property
    def expansion_space(self) -> int:
        """Davidson vectors ORCA may hold, before its own cap at the problem size."""
        return self.maxdim * self.nroots


def read_settings(config: Mapping[str, Any]) -> TddftSettings:
    """TD-DFT settings from a CONTROL dictionary.

    ``TDDFT_*`` first, then the legacy ``ESD_*`` names, then bare ``TDA``.
    A file read by :func:`delfin.config.read_control_file` always carries the
    ``TDDFT_*`` keys (the template supplies them and the old names are aliased
    onto them), so the fallbacks only matter for dictionaries built by hand.
    """
    nroots_raw = _first(config, 'TDDFT_nroots', 'ESD_nroots')
    nroots = int(nroots_raw) if nroots_raw is not None else 15
    if nroots < 1:
        raise ValueError(f"TDDFT_nroots must be a positive whole number, got {nroots_raw!r}")

    maxdim = coerce_maxdim(_first(config, 'TDDFT_maxdim', 'ESD_maxdim'))
    maxiter = coerce_maxiter(_first(config, 'TDDFT_maxiter', 'ESD_TDDFT_maxiter'))
    tda_raw = _first(config, 'TDDFT_TDA', 'ESD_TDA', 'TDA')
    follow_raw = _first(config, 'TDDFT_followiroot', 'ESD_followiroot')
    soc_raw = _first(config, 'TDDFT_SOC', 'ESD_SOC')

    return TddftSettings(
        nroots=nroots,
        maxdim=ORCA_DEFAULT_MAXDIM if maxdim == 'auto' else int(maxdim),
        maxdim_auto=maxdim == 'auto',
        maxiter=None if maxiter == 'auto' else int(maxiter),
        tda=coerce_bool(tda_raw, 'TDDFT_TDA') if tda_raw is not None else True,
        followiroot=coerce_bool(follow_raw, 'TDDFT_followiroot') if follow_raw is not None else True,
        soc=coerce_bool(soc_raw, 'TDDFT_SOC') if soc_raw is not None else False,
        additions=tuple(_addition_lines(coerce_additions(config.get('TDDFT_additions')))),
    )


def job_nroots(config: Mapping[str, Any], override_key: str) -> int:
    """NRoots for one kind of rate job: its own override key if set, else TDDFT_nroots."""
    raw = _first(config, override_key)
    return int(raw) if raw is not None else read_settings(config).nroots


def _fmt(value: Any) -> str:
    if isinstance(value, bool):
        return 'true' if value else 'false'
    return str(value)


def tddft_block(
    config: Mapping[str, Any],
    *,
    nroots: Optional[int] = None,
    iroot: Optional[int] = None,
    irootmult: Optional[str] = None,
    triplets: bool = False,
    follow: bool = False,
    dosoc: Any = FROM_CONTROL,
    job_lines: Sequence[Tuple[str, Any]] = (),
) -> str:
    """The ``%tddft ... end`` block for one ORCA job, without a trailing newline.

    ``follow`` marks a job that optimises ``iroot``: only there does
    FollowIRoot mean anything, so only there is TDDFT_followiroot written.
    ``dosoc`` is taken from TDDFT_SOC unless the job needs a fixed value
    (ISC and PHOSP need SOC on) or none at all (None).  ``job_lines`` carries
    the ESD keywords a rate job adds (SRoot, NACME ...).
    """
    settings = read_settings(config)
    lines: List[Tuple[str, Any]] = [
        ('nroots', settings.nroots if nroots is None else int(nroots)),
        ('maxdim', settings.maxdim),
        ('tda', settings.tda),
    ]
    if settings.maxiter is not None:
        lines.append(('maxiter', settings.maxiter))
    if triplets:
        lines.append(('triplets', True))
    if iroot is not None:
        lines.append(('iroot', int(iroot)))
    if irootmult:
        lines.append(('irootmult', irootmult))
    if follow and iroot is not None and settings.followiroot:
        lines.append(('followiroot', True))
    soc = settings.soc if dosoc is FROM_CONTROL else dosoc
    if soc is not None:
        lines.append(('dosoc', bool(soc)))
    lines.extend(job_lines)

    body = [f"  {key} {_fmt(value)}" for key, value in lines]
    body.extend(f"  {line}" for line in settings.additions)
    return "\n".join(["%tddft", *body, "end"])


def maxdim_hint(config: Mapping[str, Any]) -> Optional[str]:
    """A sentence for the dashboard when an explicit MaxDim is outside ORCA's range."""
    try:
        settings = read_settings(config)
    except ValueError:
        return None  # the validator reports it as an error
    low, high = ORCA_RECOMMENDED_MAXDIM
    if settings.maxdim_auto or low <= settings.maxdim <= high:
        return None
    return (
        f"TDDFT_maxdim={settings.maxdim} gives a Davidson space of {settings.maxdim} x "
        f"{settings.nroots} roots = {settings.expansion_space} vectors: ORCA multiplies "
        f"MaxDim by NRoots. ORCA recommends {low}-{high} (its own default is "
        f"{ORCA_DEFAULT_MAXDIM}); TDDFT_maxdim=auto uses that."
    )
