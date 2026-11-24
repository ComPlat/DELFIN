"""67 skills are installed and the agent is told about none of them.

MEASURED. The ``skill`` tool's whole schema is 107 tokens -- "Skill
name, no slash", no enum -- and of the 67 skill names, 6 appear in the
51k-character system prompt, all as coincidental matches of common words
(``loops``, ``merge``, ``pr``, ``push``). ``pipeline-build``, which says
"work only through the MCP server delfin-tools; never hand-roll pipeline
logic", is not in it. On a real 33-minute run the agent called
``write_file`` 20 times and ``bash`` 46 times to assemble and start ORCA
by hand, abandoned 12 of its 14 calculations, and called none of the 26
delfin-tools tools it was offered.

A POINTER, NOT A DIRECTORY. Measured: names plus the descriptions the
registry resolves cost 1,461 tokens, 12% of a prompt head whose budget
has no slack; names alone about 370. A pointer costs about 43 tokens
when it fires and nothing otherwise, and it answers the question a
directory answers badly -- not "what exists" but "what applies to
this". ``repo_error_hints`` already does this for one skill, section by
section; this is that mechanism generalised rather than a second one
beside it.

THE TRIGGERS ARE THE METADATA, so there is no second list to drift from
the first, and a skill nobody can find has a DESCRIPTION problem.

AND THE ROOT WAS THE METADATA, WHICH IS THE POINT OF THIS FILE. The
matcher first missed the very case it was built for. Three faults, each
upstream of it:

  * three skills carried ``>-`` as their whole description -- a YAML
    folded block the line-based front-matter reader took literally --
    and the text it hid is the best metadata in the set: "Use when the
    user asks for an absorption spectrum, excitation energies, or an
    excited state of a small molecule";
  * 28 skills carry no ``description:`` and were described by their
    title alone, so ``casscf-setup``'s whole retrievable vocabulary was
    {casscf, setup} while the ``>`` line beneath its heading said
    "CASSCF/NEVPT2 ... justified active space";
  * and my own workaround for both: treating a term unique to one
    skill's metadata as decisive on its own. It bought the TD-DFT case
    and cost two false positives no frequency trick removes -- "the CI
    pipeline is red" selected ``pipeline-build``, "thanks, that looks
    good" selected ``provenance``.

Fixing the two data faults took the distinctive vocabulary from 375
terms to 706 and the unfindable skills from 10 to 2, after which the
workaround was unnecessary: the same sentence matches ``energies``,
``excitation`` and ``state``. It was removed, and the false positives
went with it. A plaster that becomes unnecessary when the root is fixed
is how you know it was a plaster.
"""
from __future__ import annotations

import pytest

from delfin.agent import skill_pointers as sp
from delfin.agent import skills as sk

#: The task the real run failed at, verbatim in substance.
BENCHMARK_TASK = (
    "Compute the TD-DFT vertical excitation energies of p-nitroaniline "
    "with PBE0 and wB97X-D4, def2-TZVP, and compare the first bright state."
)


@pytest.fixture(scope="module")
def installed():
    return list(sk.discover_skills(None))


def _names(text, skills):
    return [n for n, _terms in sp.match(text, skills)]


# --- the case it exists for -------------------------------------------

def test_the_measured_failure_case_now_gets_a_pointer(installed):
    """The sentence the agent answered by hand-writing 20 inputs."""
    names = _names(BENCHMARK_TASK, installed)
    assert names, "the case this mechanism was built for still matches nothing"
    assert "orca-spectra-uvvis-esd" in names


@pytest.mark.parametrize("text,expected", [
    ("set up a CASSCF/NEVPT2 active space for this complex",
     "casscf-setup"),
    ("run a GOAT conformer search with CREST for the ensemble",
     "conformer-search"),
    ("I need an IR spectrum with scaling factors",
     "orca-spectra-ir-raman"),
    ("assemble a computational chemistry pipeline from the building blocks",
     "pipeline-build"),
])
def test_a_request_reaches_its_playbook(installed, text, expected):
    assert expected in _names(text, installed)


# --- and the silence, which is most text ------------------------------

@pytest.mark.parametrize("text", [
    "please rename this variable and add a unit test",
    "the CI pipeline is red again on main",
    "thanks, that looks good to me",
    "",
    "   ",
])
def test_text_that_calls_for_nothing_gets_nothing(installed, text):
    """Two of these were false positives of a workaround that is gone:
    `pipeline` and `good` are each unique to one skill's metadata and
    each ordinary English."""
    assert _names(text, installed) == []


def test_at_most_two_pointers(installed):
    """A turn handed five playbooks reads none."""
    greedy = ("conformer ensemble CREST GOAT casscf nevpt2 excitation "
              "spectrum absorption oscillator pipeline building blocks")
    assert len(sp.match(greedy, installed)) <= sp.MAX_POINTERS


def test_the_pointer_is_cheap_and_says_why(installed):
    p = sp.pointer(sp.match(BENCHMARK_TASK, installed))
    assert p
    assert len(p) <= 400, f"{len(p)} characters rides on a matched turn"
    assert "skill(name=" in p
    # the words that called for it, so a wrong pointer can be argued with
    assert "excitation" in p


def test_each_skill_is_named_once_per_session(installed):
    shown: set = set()
    first = sp.hint(BENCHMARK_TASK, installed, shown=shown)
    assert first
    assert sp.hint(BENCHMARK_TASK, installed, shown=shown) == ""


def test_the_module_keeps_no_state(installed):
    """``shown`` is passed in, so one test cannot decide another's
    answer -- which is how a shared registry made an earlier test in
    this repository pass alone and fail in company."""
    assert sp.hint(BENCHMARK_TASK, installed) == sp.hint(
        BENCHMARK_TASK, installed)


# --- the metadata is the trigger, so its faults are the subject -------

def test_a_folded_block_scalar_is_its_text_and_not_its_marker(installed):
    """`description: >-` with the text beneath. Three shipped skills
    carried ">-" as their whole description, and the text it hid is the
    best metadata in the set."""
    s = next(x for x in installed if x.name == "orca-spectra-uvvis-esd")
    assert s.description not in (">", ">-", "|", "|-")
    assert "excitation energies" in s.description
    assert "absorption" in s.description


def test_the_cue_line_is_read_where_no_description_key_exists(installed):
    """28 skills open with a heading and a `>` line and carry no
    `description:`. The title alone is not a trigger."""
    s = next(x for x in installed if x.name == "casscf-setup")
    assert "CASSCF setup" in s.description      # the heading
    assert "active space" in s.description      # and the line beneath it
    terms = sp.distinctive_terms(installed)["casscf-setup"]
    assert "nevpt2" in terms and "active" in terms


def test_a_block_scalar_is_folded_with_spaces():
    """Unit, so the rule is pinned without the shipped files."""
    meta, body = sk._parse_frontmatter(
        "---\nname: x\ndescription: >-\n  one line\n  and the next\n"
        "domains: code\n---\nbody\n")
    assert meta["description"] == "one line and the next"
    assert meta["domains"] == "code", "the key after the block must survive"
    assert body.strip() == "body"


def test_a_plain_value_is_unchanged_by_the_block_reader():
    meta, _body = sk._parse_frontmatter(
        "---\nname: x\ndescription: one line\n---\nb\n")
    assert meta["description"] == "one line"


# --- what cannot be found is reported, not patched around -------------

def test_the_unfindable_are_named_rather_than_given_a_trigger(installed):
    """A skill here needs a better DESCRIPTION. Adding a pattern for it
    would start the second list this module exists to avoid, so the
    count is asserted and may only go down."""
    gaps = sp.unfindable(installed)
    assert len(gaps) <= 2, (
        f"{len(gaps)} skills cannot be selected by their own metadata: "
        f"{gaps}. Write a description that says WHEN to reach for it.")


def test_a_term_shared_by_everything_is_not_a_trigger(installed):
    """"run", "method", "report", "delfin" cannot select between
    skills, and the cut is computed from the registry rather than kept
    as a stopword list."""
    terms = sp.distinctive_terms(installed)
    everywhere = {"delfin", "run", "report", "method", "code"}
    for name, ts in terms.items():
        assert not (ts & everywhere), f"{name} kept a useless trigger"
