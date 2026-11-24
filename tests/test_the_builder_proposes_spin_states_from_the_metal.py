"""OCCUPIER spin states follow the metal and its d-count, not the parity.

Before this package, ``workflow_builder.build_workflow`` chose the OCCUPIER
sequence from the parity of the electron count alone, plus preferences it
had learned. It never looked at the metal, its oxidation state or its
d-electron count, so an Fe(III) d5 centre and a closed-shell organic were
indistinguishable to it: neither the sextet a chemist scans for d5 nor the
diamagnetic state of a square-planar d8 centre could be proposed.

Red on the baseline: ``result["spin_states"]`` raises KeyError for every
request below, because the key does not exist. The control
(``test_a_closed_shell_organic_gets_no_proposal``) asserts ``is None`` from
``.get`` and therefore passes on the baseline AND after the change -- it is
the control that would catch the proposal leaking into the non-metal path.

Every assertion here names the input that would make it fail.
"""

from __future__ import annotations

from delfin.agent import workflow_builder as wb


# ---------------------------------------------------------------------
# The three cases the package names, through the builder entry.
# ---------------------------------------------------------------------

def test_an_fe_iii_d5_complex_gets_a_sequence_that_includes_m6():
    """d5 Fe(III): the high-spin sextet must be one of the states computed.

    Fails if: the proposal is absent (KeyError on baseline), or d5 yields
    a set without M=6. M=6 is the only way to put 5 unpaired electrons in
    five d orbitals, so its absence means the high-spin state is never
    computed and the scan cannot find it.
    """
    request = ("Spin states of an Fe(III) d5 complex, SMILES "
               "[Fe+3]N1C=CC=C1, opt+freq at B3LYP/def2-SVP.")
    result = wb.build_workflow(request)
    spin = result["spin_states"]
    assert spin is not None, "an Fe(III) centre must get a proposal"
    assert (spin["metal"], spin["oxidation_state"], spin["dcount"]) == ("Fe", 3, 5)
    assert 6 in spin["multiplicities"], (
        f"d5 must include the sextet M=6, got {spin['multiplicities']}")
    # The low-spin d5 state is t2g^5: one unpaired electron, M=2 -- not M=1,
    # which is impossible for an odd d-count.
    assert 2 in spin["multiplicities"], spin["multiplicities"]
    assert 1 not in spin["multiplicities"], (
        "M=1 needs an even d-count; d5 cannot have it")
    assert spin["first"] == 6, "a first-row centre tries high spin first"


def test_a_d6_low_spin_4d_metal_does_not_get_a_high_spin_only_sequence():
    """Ru(II) d6 on 4d: the diamagnetic state is the reference, not a guess.

    Fails if: the multiplicity set omits M=1 (the t2g^6 closed shell), or
    ``first`` is not 1 -- either would make the answer a high-spin-only
    sequence for a centre whose ligand field is large.
    """
    request = ("Set up a spin-state comparison for a Ru(II) d6 polypyridyl "
               "complex with OCCUPIER.")
    result = wb.build_workflow(request)
    spin = result["spin_states"]
    assert spin is not None
    assert (spin["metal"], spin["dcount"], spin["row"]) == ("Ru", 6, "4d/5d")
    assert 1 in spin["multiplicities"], (
        f"low-spin d6 is the closed shell M=1, got {spin['multiplicities']}")
    assert spin["multiplicities"] != [5], "a high-spin-only sequence is the regression"
    assert spin["first"] == 1, (
        "for a 4d centre the low-spin occupancy is the reference, so it is tried first")
    assert "4d/5d" in spin["reason"], spin["reason"]


def test_a_square_planar_d8_centre_includes_the_diamagnetic_state():
    """Ni(II) d8: square planar is diamagnetic, and it must come first.

    Fails if: the d8 row returns [3] alone (octahedral only), or names the
    triplet as the state to try first. The archive's nickel steps whose
    winner is d8-CONSISTENT -- i.e. won at an ODD multiplicity, which an
    even d-count requires -- number 148, and M=1 won 111 of them against
    M=3's 35. A scan that omits M=1 could not have found it, and one that
    tries M=3 first starts from the state that usually loses. The 301-step
    whole-element row is NOT the population: 153 of its winners sit at an
    even multiplicity, which no d8 centre can have.
    """
    result = wb.build_workflow("Set up a Ni(II) d8 square planar complex")
    spin = result["spin_states"]
    assert spin is not None
    assert spin["dcount"] == 8
    assert 1 in spin["multiplicities"], (
        f"a square-planar d8 centre is diamagnetic; got {spin['multiplicities']}")
    assert 3 in spin["multiplicities"], "the octahedral triplet is the comparison"
    assert spin["first"] == 1, (
        "the diamagnetic state is the one that usually wins for d8, so it is "
        f"tried first; got first={spin['first']}")
    assert f"try M={spin['first']} first" in spin["reason"], (
        "the reason must nominate the same state `first` names, otherwise the "
        f"proposal contradicts itself: first={spin['first']} "
        f"reason={spin['reason']!r}")


# ---------------------------------------------------------------------
# The control: nothing changes where no metal is named.
# ---------------------------------------------------------------------

def test_a_closed_shell_organic_gets_no_proposal():
    """CONTROL: no metal named -> parity rule stays in charge, unchanged.

    Green on the baseline AND after the change, by design: ``.get`` yields
    None either way. This is the test that fails if the proposal ever
    leaks into a non-metal workflow.
    """
    result = wb.build_workflow(
        "Run an opt+freq workflow on ethanol (SMILES CCO) at B3LYP/def2-SVP.")
    assert result.get("spin_states") is None
    assert result["spec"]["steps"], "the draft itself is unchanged"
    assert result["checks"]["ok"], result["checks"]["errors"]


def test_a_metal_without_an_oxidation_state_falls_back_to_parity():
    """A named but unquantified metal is NOT guessed at.

    Fails if: the parser invents an oxidation state for "an iron porphyrin".
    Green on baseline too (the key is absent there -> .get is None); it is
    a guard against over-reach, not the red test.
    """
    result = wb.build_workflow(
        "Run an opt+freq workflow on an iron porphyrin, SMILES C1=CC2=CC=C1, "
        "at B3LYP/def2-SVP.")
    assert result.get("spin_states") is None


# ---------------------------------------------------------------------
# The proposal is a proposal: it never rewrites what the caller gave.
# ---------------------------------------------------------------------

def test_an_explicitly_given_sequence_is_not_overridden():
    """A caller's own spec keeps its own steps; the proposal only informs."""
    provided = {
        "name": "mine",
        "steps": [{"step": "orca_opt", "label": "opt",
                   "kwargs": {"charge": 3, "mult": 6}}],
    }
    result = wb.build_workflow("Spin states of Fe(III) d5", base_spec=provided)
    steps = result["spec"]["steps"]
    assert steps[0]["kwargs"]["mult"] == 6, "the given mult is not replaced"
    assert len(steps) == 1, "no step is added to a spec the caller gave"
    assert result["spin_states"]["multiplicities"] == [2, 6]


# ---------------------------------------------------------------------
# The parser: explicit centres only, and a symbol must stand alone.
# ---------------------------------------------------------------------

def test_the_proposal_is_a_function_of_the_metal_not_the_wording():
    """Fe(III), iron(III), Fe(3+) and Fe3+ agree; nothing else is parsed."""
    from delfin.occupier_auto import propose_sequence_for_request as prop

    for wording in ("Fe(III)", "iron(III)", "Fe(3+)", "Fe3+"):
        got = prop(f"compute {wording} d5 spin states")
        assert got is not None, wording
        assert (got["metal"], got["oxidation_state"]) == ("Fe", 3), wording
        assert 6 in got["multiplicities"], wording
    for nothing in ("compute the spin states of an iron complex",
                    "compute an opt+freq for ethanol"):
        assert prop(nothing) is None, nothing


def test_a_word_containing_a_metal_symbol_is_not_a_metal_centre():
    """A substring match would read prose as chemistry.

    Fails if: the metal is found by ``symbol in text``. Measured on the
    alternative implementation, "a co-solvent study of a Ni(II) salt" and
    "a Ti(IV) centre in a co-factor" were read as complexes because "co"
    and "Ti" occur inside other words. A symbol must stand alone.
    """
    from delfin.occupier_auto import parse_metal_request as parse

    # These DO name a real centre after the noise word, so they parse to
    # the centre that is explicitly written with a charge.
    assert parse("a co-solvent study of a Ni(II) salt")["metal"] == "Ni"
    assert parse("a Ti(IV) centre in a co-factor")["metal"] == "Ti"
    # These name no centre at all: no metal, or no oxidation state.
    assert parse("compute the optics of a co-complex") is None
    assert parse("a co-factor study of the enzyme") is None
    assert parse("the Ni centre of the enzyme") is None, (
        "a metal symbol without an oxidation state is not a centre")
    assert parse("an Fe site") is None


def test_every_row_obeys_the_dcount_parity():
    """M=1 is impossible for an odd d-count; the table must not offer it.

    Fails if: any odd d-count yields an even multiplicity, which would be
    a state the centre cannot have. This is the assertion that caught the
    rejected implementation's d5 -> [1, 2].
    """
    from delfin.occupier_auto import _multiplicities_for_dcount as table

    for n in range(11):
        mults = table(n)
        assert mults, n
        for m in mults:
            assert (m - 1) % 2 == n % 2, (
                f"d{n} cannot have M={m}: {m - 1} unpaired electrons has the "
                f"wrong parity")
            assert 1 <= m <= 11, m
    # Both ligand-field extremes are present for the open-shell counts.
    assert table(5) == [2, 6], table(5)
    assert table(6) == [1, 5], table(6)
    assert table(8) == [1, 3], table(8)
