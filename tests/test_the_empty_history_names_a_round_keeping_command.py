"""The empty-history message may only name a command that keeps a round.

`bench history` with nothing stored tells the reader how to get a round
stored. Before this test the message named `delfin-agent bench nightly`,
and that verb does not keep anything: `cli._cmd_bench_nightly` runs
`bench_watch.nightly`, which ends at `benchmark.write_run`
(run-file + report under the runs dir). The only writer of the kept
history is `bench_nightly._persist`, whose sole caller of
`bench_history.write_round` sits behind `bench_nightly.run`, reached by
`python -m delfin.agent.bench_nightly`. The scheduler's bench hook
(`scheduler add-bench`) calls `bench_watch.nightly` as well, so the
recurring path never stores a round either.

A message that promises a stored round from a command that stores none
is the silent absence this store exists to end: a person runs it nightly
for a week and `bench history` is still empty.

What would make this fail: a message naming `bench nightly` again.
"""

from __future__ import annotations

from delfin.agent import bench_history


def test_the_empty_history_message_names_only_a_round_keeping_command():
    text = bench_history.format_history([])

    # The documented working route: the module's own entry point, whose
    # run() persists through _persist().
    assert "python -m delfin.agent.bench_nightly" in text, text

    # The route that keeps nothing must not be offered as one that does.
    assert "bench nightly" not in text, (
        "the empty-history hint promises a kept round from `bench nightly`, "
        "which runs bench_watch.nightly and never calls write_round: "
        f"{text!r}"
    )


def test_a_stored_round_still_reads_back_through_the_formatter():
    """Control case: passes before and after the message fix.

    It pins that the message change did not touch the populated path, so a
    failure of the test above cannot be blamed on the formatter as a whole.
    """
    from delfin.agent import bench_history as bh

    assert bh.format_history([]).startswith("no benchmark rounds are kept yet")
    assert "benchmark history" in bh.format_history([
        bh.Round(
            round_id="ts-abc",
            commit="0" * 40,
            model="m",
            provider="",
            seed_mode="practice",
            task_set_version="v1",
            started_at="2026-01-01T00:00:00Z",
            ended_at="2026-01-01T00:01:00Z",
            instances=[],
        )
    ])
