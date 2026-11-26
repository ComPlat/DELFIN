"""Who a commit of this session is attributed to, and whether at all.

The product writes two trailer lines into every ``git commit`` a model runs:
``Co-Authored-By: DELFIN-Agent (<model>) <email>`` and, when the session was
started by a wave or a swarm job, ``Agent: <run>/<id>``. The executor splices
them in through git's own ``--trailer`` mechanism
(``api_client._commit_trailer_rewrite``, called from ``_bash_isolation_argv``,
the one choke point every bash command passes through), and
``prompt_loader._commit_attribution_block`` tells the model to add the first one
itself on the paths that do not go through the executor.

One half reads the decision computed here today: the prompt block. The executor
does NOT yet -- ``api_client._commit_trailer_enabled`` still reads
``agent.commit_coauthor`` directly, so a wave member is told nothing about the
trailer while its ``git commit`` still gets the line spliced in. The wiring is
``tests/fixtures/h1-trailer-off-for-members.patch`` (a unified diff against
``api_client.py`` and ``user_settings.py``, which sit behind the
self-modification guard and so cannot be edited in the session that writes the
patch) and the contract for it is
``tests/test_the_executor_asks_the_same_attribution_question.py``, whose case is
``xfail(strict=True)`` until that patch is applied. Read this paragraph before
trusting the module's name: two callers, one answer is the design, and the
second caller is not wired yet.

Why a wave member answers "no" by default: its commits are read later as the
wave's work, on a branch someone else merges, and the line is public text — the
product name, the model alias, and a second trailer naming the run and the
member. A merged history cannot be reworded, so the direction here is off, and
a deployment that wants the line sets ``agent.swarm_commit_coauthor`` to true.
An ordinary session is unchanged: ``agent.commit_coauthor``, default true.

Rejected: one key for both kinds of session (a person who wants the line in
their own commits would have switched it on for the wave too, and the two
questions are different); deciding inside each caller (the prompt block and the
executor disagreed for a whole release because each read the settings itself).
"""
from __future__ import annotations

#: The run marker a wave member carries. Written for a member by
#: ``swarm/wave._launch`` and ``swarm/job.run_in_job`` and inherited by
#: everything the member starts, its subagents included. Same name as
#: ``swarm_memory.ENV_RUN`` — that module owns it; the value is repeated rather
#: than imported because importing swarm_memory pulls the store's dependencies
#: into every prompt build, and a name shared by two modules is exactly what
#: :func:`test_the_run_marker_is_the_one_the_swarm_writes` is for.
ENV_RUN = "DELFIN_SWARM_RUN"


def in_swarm_run() -> bool:
    """Was this process started as a member of a wave or a swarm job?

    Input: ``DELFIN_SWARM_RUN``. Output: whether it names a run.

    Deliberately the environment ALONE and not
    ``swarm_memory.current_run_id()``, which falls back to
    ``perms.swarm_run_id``: a plain session sets that when it delegates a
    subagent (``subagents._derive_perms``), and an ordinary session's subagent
    is not a swarm member — deciding on the fallback would quietly strip the
    attribution from every delegation in the product.
    """
    import os

    try:
        return bool(str(os.environ.get(ENV_RUN, "") or "").strip())
    except Exception:
        return False


def commit_coauthor_for(agent: dict | None) -> bool:
    """Whether THIS session ends its commits with the attribution trailer.

    Input: the ``agent`` block of the user settings (``None``/empty = the
    defaults). Output: true to attribute, false to stay silent.

    A member of a wave answers on ``agent.swarm_commit_coauthor`` (default
    false); any other session on ``agent.commit_coauthor`` (default true).
    """
    settings = agent if isinstance(agent, dict) else {}
    if in_swarm_run():
        return bool(settings.get("swarm_commit_coauthor", False))
    return bool(settings.get("commit_coauthor", True))
