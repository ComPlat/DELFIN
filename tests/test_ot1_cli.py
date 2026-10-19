"""OT1 phase 3: the CLI `trace` subcommand (real CLI path).

Committed cases go through delfin.agent.cli (the real parser the product
uses). The px-copy probe cases stay in tests/test_zz_trace_cli.py, which
is never committed (wave-17 lesson 4).
"""


def test_trace_parser_present_on_real_cli():
    """The real CLI parses `trace <session>` and dispatches it.

    History: this case began as the red control for the protected-file
    patch (SystemExit on 'trace' in the unpatched cli.py). The operator
    applied .gate/ot1_trace_cli.patch to the worktree (verified:
    cmd_trace at cli.py:4153, parser registration at cli.py:5520), which
    superseded the red direction; per the QS ruling (nacht-s29) the
    real-CLI case is the one that matters and lives committed here.
    """
    from delfin.agent import cli as real_cli
    parser = real_cli.build_parser()
    args = parser.parse_args(["trace", "whatever"])
    assert args.session_id == "whatever"
    assert args.otlp == ""
