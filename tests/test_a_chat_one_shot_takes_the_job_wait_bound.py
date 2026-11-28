"""``--job-wait`` reaches the one-shot wake loop through BOTH doors.

``delfin-agent chat -p ...`` ends in ``cmd_run`` (cli.py, ``_cmd_chat_one_shot``
returns ``cmd_run(args, prepare=_prepare)``), and ``cmd_run`` reads
``args.job_wait`` to bound the wake loop. The flag was declared only on the
``run`` subparser, so the door a script or supervisor drives -- ``chat -p ...
--output-format stream-json`` -- refused the flag at argparse (rc 2) and the
wake loop always ran at its 3600 s default: the caller could not bound it.

The control case parses ``run --job-wait`` and passes on the base commit too:
it binds the flag where it already existed, so a green here alone proves
nothing about the fix.

A parser test is the whole door for the refusal: argparse answers before the
engine exists, so the front door's answer IS the parser's answer.
"""

import pytest

from delfin.agent import cli


def test_run_door_takes_job_wait_baseline():
    """Control: the flag exists on ``run`` on the base commit too."""
    parser = cli.build_parser()
    assert parser.parse_args(["run", "hi"]).job_wait == cli._JOB_WAIT_DEFAULT_S
    assert parser.parse_args(["run", "--job-wait", "0", "hi"]).job_wait == 0


def test_chat_door_takes_job_wait():
    """The fix: ``chat -p --job-wait N`` parses and binds N.

    Red on base: argparse answers ``unrecognized arguments: --job-wait``
    because the flag was declared on ``run`` only.
    """
    parser = cli.build_parser()
    argv = ["chat", "-p", "hi", "--output-format", "stream-json",
            "--job-wait", "45"]
    args = parser.parse_args(argv)
    assert args.job_wait == 45
    assert args.print_prompt == "hi"


def test_chat_door_job_wait_default_is_the_shared_default():
    """No ``--job-wait`` given: the chat door defaults to the same value
    ``cmd_run`` reads, so the two doors wake for the same length of time
    unless the caller bounds one of them."""
    parser = cli.build_parser()
    args = parser.parse_args(["chat", "-p", "hi"])
    assert args.job_wait == cli._JOB_WAIT_DEFAULT_S


def test_chat_door_job_wait_zero_reaches_cmd_run():
    """``--job-wait 0`` means 'end with the answer and stop the shells' --
    the bound must arrive at ``cmd_run`` as 0, not as the default, because
    ``cmd_run`` distinguishes them (cli.py: ``wait_s <= 0`` branch)."""
    parser = cli.build_parser()
    args = parser.parse_args(["chat", "-p", "hi", "--job-wait", "0"])
    assert args.job_wait == 0


def test_run_door_still_has_the_flag_after_the_change():
    """Guard: adding the argument to ``chat`` must not disturb ``run``."""
    parser = cli.build_parser()
    assert parser.parse_args(["run", "--job-wait", "45", "hi"]).job_wait == 45
