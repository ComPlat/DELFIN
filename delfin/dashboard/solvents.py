"""What a solvent is, for both of the engines that can be told about one.

The viewer offers two families of method -- xtb's GFN and MOPAC's PM -- and
until now only the first of them was given the solvent the user chose.  Setting
Water and choosing PM7 produced a gas-phase answer, silently.  This module is
the one place that knows what a solvent is, so that both engines are told about
the same liquid and neither can quietly ignore it.

**The models, and which of them can carry an optimisation.**  xtb advertises
five.  Measured on a glycine, eight optimisation cycles each, against the same
run in the gas phase:

===========  =========================================================
ALPB         works for every method and every solvent; costs nothing
             (GFN-FF 50 ms against 44 in vacuum, GFN2 91 against 92)
GBSA         works, and is the older model: it is parametrised for 14 of
             the 25 solvents under GFN-FF and GFN2 and for 12 under GFN1
ddCOSMO      only GFN1 and GFN2, and five to six times the cost (GFN2
             499 ms against 92).  **Diverges under GFN-FF**: the energy
             climbed 892 kcal/mol in a single step and the structure came
             apart, so GFN-FF is not offered it
CPCM-X       refuses outright -- "CPCM-X not implemented for geometry
             optimization" in xtb's own words.  Not offered anywhere
COSMO        MOPAC's, and the only one it has.  Water shifted the heat of
             formation by -10.1 kcal/mol under PM7 and -9.5 under
             PM6-D3H4, against ALPB's -9.6 under GFN2 -- the two engines
             agree about water to within half a kcal/mol
===========  =========================================================

**Why the dielectric constants are not xtb's.**  MOPAC's continuum takes a
number, not a name, so each solvent needs one.  The obvious move is to read
them out of xtb -- it prints the constant it used -- and the obvious move is
wrong.  Asked, xtb 6.7.1 gives benzene **and** toluene as 7.0, where the
measured constants are 2.27 and 2.38, and it gives hexadecane as 1.88, which is
hexane's.  Those are placeholders, not liquids.  Using them would have had
MOPAC describe the nonpolar solvents as moderately polar: on the same glycine,
benzene at 7.0 shifted the heat of formation by -7.2 kcal/mol and benzene at
2.27 by -2.3.  The constants below are the ordinary measured ones at 25 C, and
where xtb's own table disagrees it is noted on the line.
"""

from typing import Any, Dict, List, Optional, Tuple

__all__ = ['SOLVENTS', 'MODELS', 'models_for', 'solvents_for', 'xtb_flags',
           'mopac_words', 'dielectric', 'label_of', 'model_label', 'note',
           'refusal']


#: Every solvent xtb is parametrised for, with the dielectric constant that
#: describes it.  ``xtb`` is what xtb's own ddCOSMO table says where that
#: differs, kept as a comment on the line rather than used.
SOLVENTS: Dict[str, Dict[str, Any]] = {
    'acetone': {'label': 'acetone', 'eps': 20.49},
    'acetonitrile': {'label': 'acetonitrile', 'eps': 35.69},
    'aniline': {'label': 'aniline', 'eps': 6.89},
    'benzaldehyde': {'label': 'benzaldehyde', 'eps': 17.85},
    'benzene': {'label': 'benzene', 'eps': 2.27},        # xtb says 7.0
    'ch2cl2': {'label': 'dichloromethane', 'eps': 8.93},
    'chcl3': {'label': 'chloroform', 'eps': 4.71},
    'cs2': {'label': 'carbon disulfide', 'eps': 2.61},
    'dioxane': {'label': 'dioxane', 'eps': 2.21},
    'dmf': {'label': 'DMF', 'eps': 37.22},
    'dmso': {'label': 'DMSO', 'eps': 46.83},
    'ether': {'label': 'diethyl ether', 'eps': 4.24},
    'ethanol': {'label': 'ethanol', 'eps': 24.85},
    'ethylacetate': {'label': 'ethyl acetate', 'eps': 5.99},
    'furane': {'label': 'furan', 'eps': 2.94},
    'hexadecane': {'label': 'hexadecane', 'eps': 2.04},  # xtb says 1.88
    'hexane': {'label': 'hexane', 'eps': 1.88},
    'methanol': {'label': 'methanol', 'eps': 32.61},
    'nitromethane': {'label': 'nitromethane', 'eps': 36.56},
    'octanol': {'label': 'octanol', 'eps': 9.86},
    'woctanol': {'label': 'octanol (wet)', 'eps': 8.1},  # xtb says 9.86, dry's
    'phenol': {'label': 'phenol', 'eps': 12.40},
    'thf': {'label': 'THF', 'eps': 7.43},
    'toluene': {'label': 'toluene', 'eps': 2.38},        # xtb says 7.0
    'water': {'label': 'water', 'eps': 78.36},
}


#: The models, and what each one is.  ``needs`` is the sentence a user should
#: read before choosing it, and every clause of it was measured.
MODELS: Dict[str, Dict[str, Any]] = {
    'alpb': {
        'label': 'ALPB', 'engine': 'xtb', 'flag': '--alpb',
        'needs': ('xtb\'s own recommendation, and the only model here that '
                  'covers all 25 solvents. It costs nothing measurable: 50 ms '
                  'against 44 in vacuum.'),
    },
    'gbsa': {
        'label': 'GBSA', 'engine': 'xtb', 'flag': '--gbsa',
        'needs': ('The older generalised-Born model, kept because published '
                  'numbers were computed with it. It knows fewer solvents '
                  'than ALPB, and which ones depends on the method.'),
    },
    'ddcosmo': {
        'label': 'ddCOSMO', 'engine': 'xtb', 'flag': '--cosmo',
        'needs': ('A conductor-like continuum, and the model to use when the '
                  'result is going on to COSMO-RS. Five to six times the cost '
                  'of ALPB, and not available for GFN-FF, where it diverges.'),
    },
    'cosmo': {
        'label': 'COSMO', 'engine': 'mopac',
        'needs': ('MOPAC\'s continuum, and the only one it has. It is given '
                  'the dielectric constant of the solvent you chose, so a PM '
                  'run and a GFN run are asked about the same liquid.'),
    },
}


#: Which solvents GBSA refuses, per method -- asked of xtb 6.7.1 one solvent at
#: a time rather than read off the help text, which is close but not complete.
#: GFN1 additionally has no DMF and no hexane, which the help text does say.
GBSA_REFUSES: Dict[str, frozenset] = {
    'gfnff': frozenset({
        'aniline', 'benzaldehyde', 'dioxane', 'ethanol', 'ethylacetate',
        'furane', 'hexadecane', 'nitromethane', 'octanol', 'woctanol',
        'phenol'}),
    'gfn1': frozenset({
        'aniline', 'benzaldehyde', 'dioxane', 'dmf', 'ethanol', 'ethylacetate',
        'furane', 'hexadecane', 'hexane', 'nitromethane', 'octanol',
        'woctanol', 'phenol'}),
    'gfn2': frozenset({
        'aniline', 'benzaldehyde', 'dioxane', 'ethanol', 'ethylacetate',
        'furane', 'hexadecane', 'nitromethane', 'octanol', 'woctanol',
        'phenol'}),
}

#: ddCOSMO under GFN-FF does not fail -- it runs, and destroys the structure.
#: Eight cycles on a glycine: the energy climbed 892 kcal/mol in one step,
#: after four ordinary downhill ones.  A model that is wrong quietly is worse
#: than one that refuses, so it is refused here.
DDCOSMO_REFUSES: frozenset = frozenset({'gfnff'})

#: What each method can be asked about a solvent, in the order it is offered.
#: g-xTB is absent on purpose: its build takes no solvation at all.
_BY_METHOD: Dict[str, Tuple[str, ...]] = {
    'gfnff': ('alpb', 'gbsa'),
    'gfn1': ('alpb', 'gbsa', 'ddcosmo'),
    'gfn2': ('alpb', 'gbsa', 'ddcosmo'),
    'gxtb': (),
    'pm6': ('cosmo',),
    'pm7': ('cosmo',),
    'pm6d3h4': ('cosmo',),
}


def _key(value: Any) -> str:
    return str(value or '').strip().lower()


def label_of(solvent: Any) -> str:
    """The solvent's name as a chemist writes it, or '' for the gas phase."""
    found = SOLVENTS.get(_key(solvent))
    return found['label'] if found else ''


def dielectric(solvent: Any) -> Optional[float]:
    """The dielectric constant, or None for a solvent that is not known."""
    found = SOLVENTS.get(_key(solvent))
    return found['eps'] if found else None


def model_label(model: Any) -> str:
    found = MODELS.get(_key(model))
    return found['label'] if found else ''


def models_for(method: Any) -> Tuple[str, ...]:
    """Which solvation models *method* can actually be run with."""
    return _BY_METHOD.get(_key(method), ())


def solvents_for(model: Any, method: Any) -> List[str]:
    """Which solvents *model* is parametrised for under *method*."""
    which, who = _key(model), _key(method)
    if which not in MODELS:
        return []
    refused = GBSA_REFUSES.get(who, frozenset()) if which == 'gbsa' \
        else frozenset()
    return [name for name in SOLVENTS if name not in refused]


def refusal(model: Any, solvent: Any, method: Any) -> str:
    """Why this combination cannot be run, or '' if it can.

    Said before the run rather than after it, because three of the four ways
    this can go wrong do not produce an error message from the program itself.
    """
    which, wet, who = _key(model), _key(solvent), _key(method)
    if not wet:
        return ''
    if wet not in SOLVENTS:
        return (f'{solvent!r} is not a solvent that is known here. '
                'It knows: ' + ', '.join(SOLVENTS) + '.')
    if which not in MODELS:
        return f'{model!r} is not a solvation model.'
    offered = models_for(who)
    if not offered:
        return ('This method has no implicit solvation in this build. '
                'Optimise it in the gas phase, or choose a method that can '
                'be given a solvent.')
    if which not in offered:
        if which == 'ddcosmo' and who in DDCOSMO_REFUSES:
            return ('ddCOSMO is not offered for GFN-FF: it runs, and it '
                    'destroys the structure -- measured on a glycine, the '
                    'energy climbed 892 kcal/mol in one optimisation step. '
                    'Use ALPB or GBSA.')
        return (f'{MODELS[which]["label"]} is not available for this method. '
                'It offers: '
                + ', '.join(MODELS[name]['label'] for name in offered) + '.')
    if which == 'gbsa' and wet in GBSA_REFUSES.get(who, frozenset()):
        return (f'GBSA is not parametrised for {SOLVENTS[wet]["label"]} under '
                'this method -- ALPB is, and covers every solvent here.')
    return ''


def xtb_flags(model: Any, solvent: Any) -> List[str]:
    """What to add to an xtb command line, or [] for the gas phase."""
    which, wet = _key(model), _key(solvent)
    found = MODELS.get(which)
    if not wet or found is None or found['engine'] != 'xtb':
        return []
    return [found['flag'], wet]


#: Where ORCA writes a liquid's name differently from xtb.
#:
#: The names in :data:`SOLVENTS` are xtb's, because xtb is what most of this
#: runs on.  ORCA takes the same continuum by keyword and agrees on 24 of the
#: 25 -- measured with ORCA 6.1.1 on this box, one ``! XTB2 SP ALPB(name)``
#: per name: 24 ran and ``ALPB(ether)`` stopped with "UNRECOGNIZED OR
#: DUPLICATED KEYWORD(S) IN SIMPLE INPUT LINE" before doing any work.  So
#: every press that goes through ORCA -- climb from here, the chain's second
#: half, the band -- died outright on the one solvent the dropdown calls
#: diethyl ether, while the same choice ran perfectly well on the xtb halves.
#:
#: ``diethylether`` is the same liquid and not merely a name ORCA accepts: on
#: one water, ORCA's ``ALPB(diethylether)`` gives -5.078841108810 Eh against
#: xtb's own ``--alpb ether`` -5.078841108805, which is the agreement the two
#: also have in water (-5.084878982390 against -5.084878982393).
_ORCA_NAMES: Dict[str, str] = {'ether': 'diethylether'}


def orca_keyword(solvent: Any) -> str:
    """ORCA's continuum keyword for this liquid, or '' for the gas phase.

    ALPB and nothing else, because ORCA's xtb driver has nothing else to
    offer: asked for ``CPCM(water)`` or ``SMD(water)`` beside ``! XTB2`` it
    stops in ``main_input_check`` before doing any arithmetic.  So a press
    that runs through ORCA runs ALPB whatever the model box says, and what it
    owes the user is to name ALPB rather than what was asked for.
    """
    wet = _key(solvent)
    if not wet:
        return ''
    # A name this table does not know is handed on rather than dropped: ORCA
    # refuses one it does not recognise and says so, and a press that stops
    # with a reason is better than one that quietly answers the gas-phase
    # question instead.
    return f'ALPB({_ORCA_NAMES.get(wet, wet)})'


def mopac_words(solvent: Any) -> List[str]:
    """What to add to a MOPAC keyword line, or [] for the gas phase.

    MOPAC switches its continuum on by being given a dielectric constant; there
    is no solvent name to give it.  ``RSOLV`` is left at MOPAC's default.
    """
    eps = dielectric(solvent)
    return [] if eps is None else [f'EPS={eps:g}']


def note(model: Any, solvent: Any) -> str:
    """Which solvent a result is about, for the status line.

    A geometry optimised in the gas phase and one optimised in water are two
    different answers to two different questions, and a result that does not
    say which it is invites them to be compared.  The model is named too: ALPB
    and ddCOSMO do not agree to better than a few kcal/mol.
    """
    wet = _key(solvent)
    if not wet or wet not in SOLVENTS:
        return ''
    which = model_label(model)
    return f' In {SOLVENTS[wet]["label"]}' + (f' ({which}).' if which else '.')
