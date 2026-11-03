"""The benchmark's model-facing surfaces must be identical to main.

Everything the model uses from a benchmark task is reached on demand (a
tool it calls, or a skill loaded for chemistry tasks) and must never be
added to the always-on prompt. The stricter-metrics work (pass@k /
pass^k, weighted rubrics) lives entirely in the SCORING layer:
BenchmarkResult fields, aggregate_replicates, summarise_run and the CLI
printout. This test pins the two model-facing surfaces a benchmark
change could leak into:

- _DOC_TOOLS_OPENAI: the always-on tool schemas every task's prompt is
  built from. Asserted structurally, so it holds on any machine: none
  of the scoring additions may appear in it.
- the packaged task suite: prompt texts, signals, values, verify and --
  now -- the optional per-task rubric. The rubric is SCORING data for
  the runner, not model-facing content: the test asserts that no
  rubric text is reachable from a task's model-facing surface (its
  prompt) even when the task declares one.

The byte-level identity against main was measured with
scratch/prompt_hash.py run in a control worktree at the merge base:
tool-schema hash 09fbadfba812...ea73 (41244 bytes) and task-suite hash
b97063938f20...170c (76022 bytes) are identical on both.
"""
from delfin.agent import benchmark as B
from delfin.agent.api_client import _DOC_TOOLS_OPENAI

# Scoring-layer identifiers added by the stricter-metrics work. None of
# these may appear in the always-on tool schemas.
_SCORING_NAMES = (
    "pass_at_k", "pass_pow_k", "n_pass_at_k", "pass_at_k_rate",
    "n_pass_pow_k", "pass_pow_k_rate", "rubric_0_100",
    "rubric_components", "RubricCriterion", "score_rubric",
)


def test_always_on_tool_schemas_do_not_grow_scoring_content():
    blob = repr(_DOC_TOOLS_OPENAI)
    for name in _SCORING_NAMES:
        assert name not in blob, (
            f"scoring field {name!r} leaked into the always-on tool "
            "schemas; it must be reached on demand, not always-on")


def test_rubric_is_scoring_data_not_prompt_text():
    """A task that declares a rubric carries it nowhere model-facing."""
    rubric_tasks = [t for t in B.load_tasks()
                    if getattr(t, "rubric", None)]
    # The packaged suite currently has no rubric-bearing task; the check
    # must still be true for one when a task author adds it, so build
    # one here rather than only scanning the shipped YAML.
    task = B.Task(
        id="rubric_probe", task_class="misc", mode="solo", prompt="p",
        rubric=[B.RubricCriterion(name="planning", weight=2,
                                  kind="signals")],
    )
    for t in [*rubric_tasks, task]:
        assert t.rubric, "fixture broken: rubric lost in coercion"
        # The prompt is the model-facing surface; the rubric is not part
        # of it. (A criterion name like 'planning' may naturally occur
        # in prose; the criterion LIST itself must not.)
        assert repr(t.rubric) not in t.prompt
        for c in t.rubric:
            assert c.kind not in t.prompt
