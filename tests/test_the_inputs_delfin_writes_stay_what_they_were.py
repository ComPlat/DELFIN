"""What DELFIN writes for a fixture is what it wrote before, or the change is deliberate.

An input that comes out a little different is invisible until somebody runs
again: a recalc then computes a finished job a second time, and two states of
one molecule can end up computed with different settings.  Both were found on
the archive months after the change -- the first coordination sphere's basis
(a ligand crossing the cutoff got another basis than the same ligand in the
next state), and the resolution of the automatic OCCUPIER sequences (an edit
that reached OCCUPIER made 434 configuration inputs of 186 archived jobs of
one user come out different).

Here every writer runs over the same two made-up systems and the result is
compared with the input beside this test.  A deliberate change is made by
running the suite once with DELFIN_UPDATE_GOLDEN=1 and committing what
changed, so it is read in the diff of the commit that causes it.
"""

from __future__ import annotations

import difflib
import os
from pathlib import Path

import pytest

from delfin.common import orca_input as oi
from orca_input_corpus import build_corpus

GOLDEN = Path(__file__).parent / "fixtures" / "orca_inputs"

#: Every input of the corpus, named so a missing one is a failure and not a
#: test that quietly checks nothing.
NAMES = (
    "initial.inp", "red_step_1.inp", "initial_with_properties.inp",
    "red_step_1_OCCUPIER/input.inp", "red_step_1_OCCUPIER/input3.inp",
    "ESD/S0.inp", "ESD/S1.inp", "ESD/T1.inp", "ESD/S1_T1_ISC_ms0.inp", "ESD/S1_S0_IC.inp",
    "ESD/S1_S0_FLUOR.inp", "ESD/T1_S0_PHOSP.inp",
    "overridden/initial.inp", "overridden/S1_S0_IC.inp",
    "ox_step_1_with_reorganisation.inp", "solv_complex.inp",
    "nmr_of_a_structure.inp", "nmr_reference.inp",
)


@pytest.fixture(scope="module")
def corpus(tmp_path_factory) -> dict:
    return build_corpus(tmp_path_factory.mktemp("corpus"))


def _golden_path(name: str) -> Path:
    return GOLDEN / name


def test_every_writer_of_the_corpus_is_covered(corpus):
    # a writer that stops being called here stops being watched
    assert sorted(corpus) == sorted(NAMES)


@pytest.mark.parametrize("name", NAMES)
def test_the_input_is_the_one_beside_this_test(corpus, name):
    written = corpus[name]
    golden = _golden_path(name)
    if os.environ.get("DELFIN_UPDATE_GOLDEN") == "1":
        golden.parent.mkdir(parents=True, exist_ok=True)
        golden.write_text(written, encoding="utf-8")
    assert golden.exists(), f"no recorded input for {name}; write it with DELFIN_UPDATE_GOLDEN=1"
    recorded = golden.read_text(encoding="utf-8")
    if written != recorded:
        diff = "".join(difflib.unified_diff(recorded.splitlines(True), written.splitlines(True),
                                            fromfile=f"recorded/{name}", tofile=f"written/{name}"))
        pytest.fail(f"{name} is written differently than it was.  If that is the intention, run the "
                    f"suite once with DELFIN_UPDATE_GOLDEN=1 and commit the change with its reason.\n{diff}")


def test_every_input_of_the_corpus_reads_as_orca_input(corpus):
    for name, text in corpus.items():
        for job in oi.split_jobs(text):
            oi.parse_job(job)                     # raises on anything ORCA could not read


def test_reading_an_input_and_writing_it_again_changes_nothing(corpus):
    # the recovery rewrites inputs through this model; what it does not
    # understand it must at least hand back unchanged
    for name, text in corpus.items():
        for job in oi.split_jobs(text):
            assert oi.render_job(oi.parse_job(job)) == job, name


def test_the_metal_keeps_its_own_basis_in_every_job_that_has_it(corpus):
    for name in ("red_step_1.inp", "red_step_1_OCCUPIER/input.inp", "red_step_1_OCCUPIER/input3.inp"):
        metal_line = [l for l in corpus[name].splitlines() if l.strip().startswith("Co ")]
        assert metal_line and 'NewGTO "def2-TZVP" end' in metal_line[0], name


def _atoms_with_their_own_basis(text: str) -> set:
    return {" ".join(line.split()[:4]) for line in text.splitlines() if "NewGTO" in line}


def test_every_writer_draws_the_first_coordination_sphere_the_same(corpus):
    """The same complex, the same CONTROL: the sphere cannot depend on the writer.

    It did.  The main workflow reads the covalent radii only when CONTROL
    gives no scale of its own, so with the shipped ``scale=1.3`` the built-in
    table decided; the stability-constant job read the radii regardless and
    drew its sphere with the other table (Co-O at 2.50 A against 2.26 A).  On
    the corpus complex that is the oxygen at 2.40 A: metal basis in the run's
    own jobs, main basis in the stability job of the same run -- and a
    stability constant is a difference of energies, so the step is in it.
    """
    spheres = {name: _atoms_with_their_own_basis(corpus[name])
               for name in ("red_step_1.inp", "red_step_1_OCCUPIER/input.inp", "solv_complex.inp")}
    main = spheres["red_step_1.inp"]
    assert len(main) == 4, main           # the metal, both oxygens, the nitrogen
    for name, sphere in spheres.items():
        assert sphere == main, f"{name} gives another set of atoms their own basis than red_step_1.inp"
