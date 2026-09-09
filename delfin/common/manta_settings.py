"""What CONTROL.txt says about the structure builder, in the shapes it is needed in.

MANTA is configured three different ways at once -- keyword arguments to
:func:`delfin.smiles_converter.smiles_to_xyz_isomers`, ``DELFIN_FFFREE_*``
environment flags that must be set *before* that module is imported, and a
selection stage that runs afterwards.  This module is the one place that reads
CONTROL and hands each of the three what it needs, so a key is spelled once.

**Why the environment half exists at all.**  The construction presets are read
at import time by :mod:`delfin.smiles_converter`, so setting them inside the
conversion call is too late -- by then the module has already decided which
builder it is.  :func:`apply_construction_env` therefore has to run early, from
``cli.py``, right after the CONTROL file is read.  ``delfin-manta`` does the
same thing for the same reason (``cli_manta.main`` sets the environment before
its own import), and this reuses that function rather than restating its
34-flag list.

**Legacy names are accepted.**  Every ``MANTA_*`` key falls back to its
``GUPPY_*`` predecessor, so a CONTROL file written before the rename keeps
working and means the same thing.
"""

from __future__ import annotations

import logging
import os
from typing import Any, Dict, List, Mapping, Optional, Tuple

logger = logging.getLogger(__name__)

__all__ = [
    'QUALITY_MODES', 'CONSTRUCTION_MODES', 'RANK_METHODS',
    'builder_options', 'apply_construction_env', 'selection_options',
    'SCREEN_METHODS', 'OPT_MODES', 'REFINE_MODES',
]

#: Seed counts behind the presets, from ``smiles_converter._DELFIN_PROFILES``:
#: fast 12, normal 20, max 40, extreme 60.  ``extreme`` is the command line's
#: own default and the only one the convergence study finds reliable on
#: multi-isomer systems -- fast and normal miss the GFN2 global minimum by
#: about 2.5 kcal/mol.  The pipeline used to pass nothing at all, which is the
#: library default of 20, i.e. weaker than what a user gets by typing
#: ``delfin-manta``.
QUALITY_MODES: Tuple[str, ...] = ('fast', 'normal', 'max', 'extreme')

CONSTRUCTION_MODES: Tuple[str, ...] = ('champion', 'builder', 'default')

#: What the selection stage may rank with.  ``clash`` is the builder's own
#: steric ordering and costs nothing; the rest need a binary.  ``gxtb`` is
#: separate from the ``gfn*`` family because it is a different program -- an
#: ordinary xtb accepts ``--gxtb`` and silently runs GFN2 instead, so it can
#: never be reached by passing a flag to the xtb beside it.
SCREEN_METHODS: Tuple[str, ...] = (
    'none', 'clash', 'gfnff', 'gfn0', 'gfn1', 'gfn2', 'gxtb')

#: Kept under its old name so anything importing it still resolves.
RANK_METHODS: Tuple[str, ...] = SCREEN_METHODS

#: What may happen to a frame that survived the screen.
OPT_MODES: Tuple[str, ...] = ('none', 'xtb')

#: And to the best of the optimised.  GOAT and CREST are both already in the
#: tree and both already take a directory and a file, so neither needs the
#: other's machinery rebuilt.
REFINE_MODES: Tuple[str, ...] = ('none', 'goat', 'crest')

#: Mirrors control_validator._GOAT_TOPK_CEILING; imported rather than
#: restated so the two cannot drift apart.
try:
    from delfin.common.control_validator import _GOAT_TOPK_CEILING as GOAT_TOPK_CEILING
except Exception:  # noqa: BLE001 - the reader must not need the validator
    GOAT_TOPK_CEILING = 10

REFINE_TOPK_CEILING = GOAT_TOPK_CEILING

_TRUE = {'yes', 'y', 'true', '1', 'on'}
_FALSE = {'no', 'n', 'false', '0', 'off'}


def _raw(config: Mapping[str, Any], *names: str) -> Optional[str]:
    """The first of ``names`` that CONTROL actually carries a value for."""
    for name in names:
        value = config.get(name)
        if value is None:
            continue
        text = str(value).strip()
        if text and not (text.startswith('[') and text.endswith(']')):
            return text
    return None


def _boolean(config: Mapping[str, Any], *names: str,
             default: bool) -> bool:
    text = _raw(config, *names)
    if text is None:
        return default
    lowered = text.lower()
    if lowered in _TRUE:
        return True
    if lowered in _FALSE:
        return False
    return default


def _integer(config: Mapping[str, Any], *names: str,
             default: Optional[int]) -> Optional[int]:
    text = _raw(config, *names)
    if text is None:
        return default
    try:
        return int(float(text))
    except (TypeError, ValueError):
        return default


def _number(config: Mapping[str, Any], *names: str,
            default: Optional[float]) -> Optional[float]:
    text = _raw(config, *names)
    if text is None:
        return default
    try:
        return float(text)
    except (TypeError, ValueError):
        return default


def _choice(config: Mapping[str, Any], *names: str,
            allowed: Tuple[str, ...], default: str) -> str:
    text = _raw(config, *names)
    if text is None:
        return default
    lowered = text.lower()
    if lowered in allowed:
        return lowered
    logger.warning('%s=%s is not one of %s; using %s',
                   names[0], text, ', '.join(allowed), default)
    return default


def builder_options(config: Mapping[str, Any]) -> Dict[str, Any]:
    """Keyword arguments for ``smiles_to_xyz_isomers``, from CONTROL.

    ``max_isomers`` is deliberately **not** here: it travels its own way to the
    sampler, which needs it for the enumeration cap as well.  Note what it does
    when it gets there -- the builder's pre-UFF candidate budget is
    ``max_isomers * cap_mult``, so a small number shrinks the search rather than
    shortening the answer, and ``0`` means the complete manifold.
    """
    options: Dict[str, Any] = {
        'quality_mode': _choice(config, 'MANTA_QUALITY',
                                allowed=QUALITY_MODES, default='extreme'),
        'apply_uff': _boolean(config, 'MANTA_UFF', default=True),
        'deterministic': _boolean(config, 'MANTA_DETERMINISTIC', default=True),
        'include_binding_mode_isomers': _boolean(
            config, 'MANTA_BINDING_MODES', default=True),
        'collapse_label_variants': _boolean(
            config, 'MANTA_COLLAPSE_VARIANTS', default=False),
    }

    seeds = _integer(config, 'MANTA_SEEDS', default=None)
    if seeds is not None and seeds > 0:
        # The builder clamps to the 1024-entry seed schedule itself; clamping
        # here too would hide a typo rather than let it be seen in the log.
        options['seeds_override'] = seeds

    confs = _integer(config, 'MANTA_NUM_CONFS', default=None)
    if confs is not None and confs > 0:
        options['num_confs'] = confs

    hapto = _choice(config, 'MANTA_HAPTO',
                    allowed=('auto', 'on', 'off'), default='auto')
    if hapto != 'auto':
        # None means "decide from the molecule"; True/False force it.  False
        # makes hapto SMILES fail fast rather than be approximated.
        options['hapto_approx'] = (hapto == 'on')

    return options


#: What the last :func:`apply_construction_env` call actually applied.  The
#: pipeline records this rather than asking the config a second time, so the
#: provenance is the applied set by construction and cannot drift from it.
#: Mutated in place so that an early ``from ... import LAST_CONSTRUCTION``
#: still sees it.
LAST_CONSTRUCTION: Dict[str, Any] = {'config': None, 'flags': {}}


def apply_construction_env(config: Mapping[str, Any],
                           environ: Optional[Dict[str, str]] = None
                           ) -> Dict[str, str]:
    """Set the builder's environment flags.  Must run before the import.

    Returns what it set, so a caller can log it or a test can assert on it
    without reading ``os.environ`` back.

    The 34 flags behind ``champion`` are applied as a set and are not exposed
    one by one: the file that defines them records that an earlier 29-flag
    stack scored 13.7 % topology-correct against 33.6 % for no flags at all, so
    a hand-picked subset is a way to make the builder worse.  ``MANTA_ENV`` is
    the escape hatch for the one flag somebody genuinely needs.
    """
    target = os.environ if environ is None else environ
    applied: Dict[str, str] = {}

    construction = _choice(config, 'MANTA_CONSTRUCTION',
                           allowed=CONSTRUCTION_MODES, default='champion')
    before = dict(target)
    try:
        from delfin.cli_manta import _apply_construction_env
    except Exception as exc:                      # noqa: BLE001
        logger.warning('MANTA construction preset unavailable: %s', exc)
    else:
        if environ is None:
            _apply_construction_env(construction)
        else:
            saved = os.environ
            try:
                os.environ = target                # type: ignore[assignment]
                _apply_construction_env(construction)
            finally:
                os.environ = saved                 # type: ignore[assignment]
        for key, value in target.items():
            if before.get(key) != value:
                applied[key] = value

    # Keep the builder inside the allocation it was given.  The batch-UFF pool
    # is bounded by ``min(len(batch), os.cpu_count(), DELFIN_MAX_PROCESS_WORKERS)``
    # -- by the whole machine and by 64, never by PAL.  On a 384-core node a run
    # allocated PAL=8 would spawn up to 64 UFF processes, eight times its share,
    # and on a shared node that is somebody else's job it is taking.
    #
    # Not pinned to 1: this pipeline makes one builder call and parallelises the
    # frame optimisations afterwards, so the build and the optimisations never
    # overlap and there is nothing to oversubscribe against.  A harness that
    # runs many builds at once does need 1, for the opposite reason.
    if not os.environ.get('DELFIN_MAX_PROCESS_WORKERS'):
        try:
            pal = max(1, int(float(str(config.get('PAL') or 1).strip())))
        except (TypeError, ValueError):
            pal = 1
        target['DELFIN_MAX_PROCESS_WORKERS'] = str(pal)
        applied['DELFIN_MAX_PROCESS_WORKERS'] = str(pal)

    # The build budget.  This is the only wall-clock limit MANTA has, and it is
    # enforced where the build actually is: ``smiles_to_xyz_isomers`` runs in an
    # isolated subprocess by default, and that subprocess is what gets killed.
    # A cut build returns nothing rather than a smaller answer, so the number
    # is a real decision -- 1800 s cut 39.6 % of complexes above 80 atoms.
    budget = _number(config, 'MANTA_TIME_BUDGET', default=None)
    if budget is not None and budget >= 0:
        target['DELFIN_UI_ISOLATE_TIMEOUT'] = str(int(budget))
        applied['DELFIN_UI_ISOLATE_TIMEOUT'] = str(int(budget))

    #: What was applied, so the run can record it instead of recomputing it.
    #: Recomputing is how a provenance record comes to disagree with the run it
    #: describes: the same environment lookup written twice with two different
    #: fallbacks reports one thing and does another, and both look right.
    #: (Seen in the MANTA harness: builds ran at quality ``extreme`` for weeks
    #: while every results file recorded the shipped default, because the call
    #: site and the recording site each supplied their own default.)
    #: The drop-capable filters are all off by default, including under
    #: ``champion``.  For a pipeline that then optimises every surviving frame,
    #: a torn frame that leads the ordering costs a whole DFT chain -- and the
    #: builder's own ranker is documented as scoring a decoordinated frame
    #: perfectly, because it has no overlap to penalise.  Each gate is
    #: never-worse by construction: asked to empty the list, it returns the
    #: list unchanged.
    gates = (
        ('MANTA_CLEAN_GATE', 'DELFIN_FFFREE_CLEAN_GATE', True),
        ('MANTA_TOPOLOGY_GATE', 'DELFIN_FFFREE_TOPOLOGY_GATE', True),
        ('MANTA_DEDUP', 'DELFIN_FFFREE_PERMUTE_DEDUP', True),
        ('MANTA_COORD_INTEGRITY', 'DELFIN_FFFREE_COORD_INTEGRITY', False),
        ('MANTA_CONF_COMPLETE', 'DELFIN_FFFREE_CONF_COMPLETE', False),
    )
    for key, flag, default in gates:
        value = '1' if _boolean(config, key, default=default) else '0'
        target[flag] = value
        applied[flag] = value

    extra = _raw(config, 'MANTA_ENV')
    if extra:
        for piece in extra.replace(';', ',').split(','):
            if '=' not in piece:
                continue
            name, _, value = piece.partition('=')
            name, value = name.strip(), value.strip()
            if name:
                target[name] = value
                applied[name] = value

    # Mutated in place, never rebound.  ``from manta_settings import
    # LAST_CONSTRUCTION`` binds the object, so rebinding the module global
    # leaves every existing import pointing at the old empty dict -- a name
    # bound at one moment and read at another, which is the third time that
    # shape has cost something in this integration.
    LAST_CONSTRUCTION.clear()
    LAST_CONSTRUCTION.update({'config': construction, 'flags': dict(applied)})
    return applied


#: Frame count above which "optimise everything" stops being the cheap option
#: and a single-point screen is run first.  Measured over 5810 systems built at
#: champion/extreme with ``max_isomers=0``, the manifold size per system is
#:
#:     mean 26.1 - p10 3 - p25 4 - p50 14 - p75 33 - p90 64 - p95 90
#:     p99 190 - max 399
#:
#: so a fixed "optimise all" is right for the median system and wrong for its
#: tail: 4.1 % of systems return more than 100 frames, and each one of those is
#: an ORCA optimisation.  Screening above this leaves roughly three quarters of
#: systems on the path that needs no screen at all, and caps the worst case at
#: this many optimisations instead of 399.
#:
#: The number will bite slightly more often over time: the builder's
#: PUCKER_SYMM_ADD, added after that measurement, can take the pucker
#: combination union to twice the budget cap on systems with two or more rings
#: in one automorphism orbit.
SCREEN_ABOVE_FRAMES = 30


def selection_options(config: Mapping[str, Any]) -> Dict[str, Any]:
    """How the frames are funnelled down to the one geometry the pipeline takes.

    Three stages, each switchable on its own, because they cost three different
    amounts and a run has different reasons to want each:

    1. **screen** -- one single point per frame, no geometry change.  Cheap, and
       the only thing that gives the frames an energy order at all: the builder
       returns them ordered by least steric clash, and its own docstring records
       that on a clean ensemble 86 % of them tie at the top score, so the order
       delivered is the enumeration order.
    2. **optimise** -- a real geometry optimisation per survivor.  This is where
       the cost is.  ``MANTA_SCREEN_KEEP`` is what decides how much of it there
       is, which is why the screen exists: it picks *which* frames deserve it
       instead of paying for all of them or guessing from clash.
    3. **refine** -- GOAT or CREST on the best of the optimised.

    Setting the screen to ``none`` and the keep to ``all`` optimises everything,
    which is what the pipeline did before there was a choice.  Setting
    ``MANTA_OPT=none`` optimises nothing and ranks on the screen alone.
    """
    optimise = _choice(config, 'MANTA_OPT',
                       allowed=OPT_MODES, default='xtb')
    if _raw(config, 'MANTA_OPT') is None and _raw(config, 'MANTA_RANK_OPT') is not None:
        # the older two-stage spelling: a yes/no on whether to optimise
        optimise = 'xtb' if _boolean(config, 'MANTA_RANK_OPT', default=True) else 'none'
    refine = _choice(config, 'MANTA_REFINE',
                     allowed=REFINE_MODES, default='goat')
    screen_keep = _count(config, 'MANTA_SCREEN_KEEP', 'MANTA_KEEP')

    # What ranks the frames, when CONTROL does not say.  There are exactly two
    # things that can, and which one is wanted follows from how much is being
    # optimised:
    #
    #   optimise everything  -> the optimisation *is* the ranking; a single
    #                           point on every frame first would be measured
    #                           and then thrown away.
    #   optimise a subset    -> something has to choose the subset, and only a
    #                           single point on every frame can.
    #   optimise nothing     -> the single point is the whole ranking.
    #
    # An explicit MANTA_SCREEN always wins over this; it is a default, not a
    # rule.  What it prevents is the two silent wastes: paying for single
    # points nobody reads, and cutting an unranked list at N.
    optimises_everything = (optimise == 'xtb' and screen_keep is None)
    screen = _choice(config, 'MANTA_SCREEN', 'MANTA_RANK', 'GUPPY_RANK',
                     allowed=SCREEN_METHODS,
                     default='none' if optimises_everything else 'gfn2')
    screen_explicit = _raw(config, 'MANTA_SCREEN', 'MANTA_RANK',
                           'GUPPY_RANK') is not None

    return {
        'screen': screen,
        'screen_explicit': screen_explicit,
        'screen_keep': screen_keep,
        'screen_above': _integer(config, 'MANTA_SCREEN_ABOVE',
                                 default=SCREEN_ABOVE_FRAMES),
        'optimise': optimise,
        'optimise_method': _raw(config, 'MANTA_OPT_METHOD') or None,
        'multiplicity': _count(config, 'MANTA_OPT_MULTIPLICITY',
                               'MANTA_RANK_MULTIPLICITY', word='auto'),
        'multiplicities': _multiplicity_list(config, 'MANTA_MULTIPLICITIES'),
        'refine': refine,
        'refine_topk': max(0, min(REFINE_TOPK_CEILING,
                                  _integer(config, 'MANTA_REFINE_TOPK',
                                           'MANTA_GOAT', 'GUPPY_GOAT',
                                           default=0) or 0)),
        'rmsd_cutoff': _number(config, 'MANTA_RMSD_CUTOFF',
                               'GUPPY_RMSD_CUTOFF', default=0.3),
        'energy_window_kcal': _number(config, 'MANTA_ENERGY_WINDOW',
                                      'GUPPY_ENERGY_WINDOW_KCAL', default=25.0),
        'parallel_jobs': _parallel_jobs(config),
        'time_budget_s': _number(config, 'MANTA_TIME_BUDGET', default=1800.0),
    }


#: Cores one frame optimisation can actually use.  An xtb geometry optimisation
#: of a 40-100 atom complex stops scaling well before this; the number is a
#: compromise between that and not spawning more ORCA processes than the node
#: wants to schedule.
CORES_PER_FRAME_JOB = 4


def _parallel_jobs(config: Mapping[str, Any]) -> int:
    """How many frames are worked on at once.

    The old fixed 4 wasted the machine in both directions: on PAL=8 it split
    two cores per job, and on PAL=450 it left the node 96 % idle while frames
    queued.  Unset -- or ``auto`` -- this now follows PAL, at roughly
    ``CORES_PER_FRAME_JOB`` cores per frame.

    Two things bound it further at the point of use and are deliberately not
    duplicated here: the sampler clamps to the number of frames that actually
    exist (``min(parallel_jobs, total_jobs, pal)``), and the total memory is
    ``pal * maxcore`` however the split falls, so a wider split does not ask
    for more RAM.
    """
    text = (_raw(config, 'MANTA_PARALLEL_JOBS', 'GUPPY_PARALLEL_JOBS') or '').strip().lower()
    if text and text not in ('auto', 'default'):
        try:
            return max(1, int(float(text)))
        except (TypeError, ValueError):
            pass
    try:
        pal = int(float(str(config.get('PAL') or 1).strip()))
    except (TypeError, ValueError):
        pal = 1
    return max(1, pal // CORES_PER_FRAME_JOB)


def _multiplicity_list(config: Mapping[str, Any], *names: str) -> List[int]:
    """The spin states to test, or an empty list meaning "decide from parity".

    A metal complex is not one molecule per topology: which coordination isomer
    lies lowest depends on the spin state and vice versa, so the frames are
    offered at every multiplicity named here and compete in one list.

    Left unset -- or set to ``auto`` -- this stays empty and the sampler falls
    back to what DELFIN already does everywhere else: an even electron count is
    a closed-shell singlet, an odd one is a doublet.  That is not a spin-state
    prediction and is not meant to be; it is the default that does not quietly
    assume a d7 complex is closed-shell.  ``OCCUPIER`` is where the real spin
    question is asked, later and properly.

    Accepts ``1,3,5`` or ``1 3 5``; duplicates and values below 1 are dropped.
    """
    text = (_raw(config, *names) or '').strip().lower()
    if not text or text in ('auto', 'none', 'default', '0'):
        return []
    found: List[int] = []
    for token in text.replace(',', ' ').split():
        try:
            value = int(float(token))
        except (TypeError, ValueError):
            continue
        if value >= 1 and value not in found:
            found.append(value)
    return sorted(found)


def _count(config: Mapping[str, Any], *names: str,
           word: str = 'all') -> Optional[int]:
    """A positive integer, or ``None`` for the word that means "no limit".

    ``all`` and ``auto`` are the two such words in this block; both mean "do not
    impose a number here", which is not the same as zero.
    """
    text = (_raw(config, *names) or word).strip().lower()
    if text in (word, 'all', 'auto', '0', ''):
        return None
    try:
        return max(1, int(float(text)))
    except (TypeError, ValueError):
        return None
