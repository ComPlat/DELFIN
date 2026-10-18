"""Every bare-slash form the dashboard prompt promises is one the parser accepts.

``dashboard_agent.md`` (the dashboard agent's own prompt) teaches that a
bare slash command on its own line — ``/tab calc`` without the ``ACTION:``
prefix — works, and lists the prefixes it works for. The dashboard parser
only accepts a bare-slash line when the command starts with a prefix in
``_KNOWN_SLASH_PREFIXES`` (tab_agent.py, inside ``create_tab``); anything
else on a bare line is silently dropped.

F1 (wave 17 inventory): the doc listed ``/effort``, ``/mode``,
``/provider`` and ``/model`` as bare-slash forms, but the parser's
whitelist did not contain them — a model that wrote the bare form the
prompt showed it had the line dropped in silence, while the same line
with the ``ACTION:`` prefix executed. The benchmark mirror
(``benchmark_runner._BARE_SLASH_PREFIXES``) accepted all of them, so the
suite never saw the drop either.
"""

import re
from pathlib import Path

_PACK = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack"
_DOC = _PACK / "agents" / "dashboard_agent.md"


def _dashboard_whitelist() -> set[str]:
    """The bare-slash prefixes the dashboard parser actually accepts."""
    import inspect

    from delfin.dashboard import tab_agent

    src = inspect.getsource(tab_agent.create_tab)
    m = re.search(r"_KNOWN_SLASH_PREFIXES\s*=\s*frozenset\(\((.*?)\)\)", src, re.S)
    assert m, "dashboard parser whitelist (_KNOWN_SLASH_PREFIXES) not found"
    return set(re.findall(r'"(/\w+)"', m.group(1)))


def _bare_prefixes_promised() -> set[str]:
    """The bare-slash prefixes the doc promises, from its own list."""
    body = _DOC.read_text(encoding="utf-8")
    # The doc names the bare-slash prefixes in backticks in the sentence
    # about the bare slash on its own line (item 5 of the accepted forms).
    m = re.search(
        r"bare slash on its own line.*?(`/[a-z]+`(?:,\s*`/[a-z]+`)*)",
        body, re.S | re.I)
    assert m, "dashboard doc's bare-slash prefix list not found"
    return {tok for tok, _ in re.findall(r"`(/(\w+))`", m.group(1))}


def test_every_promised_bare_slash_prefix_is_accepted_by_the_dashboard():
    promised = _bare_prefixes_promised()
    assert promised, "no bare-slash prefixes found in the doc's list"
    accepted = _dashboard_whitelist()
    assert accepted, "dashboard parser whitelist came back empty"
    missing = promised - accepted
    assert not missing, (
        f"dashboard_agent.md promises bare-slash forms {sorted(missing)} "
        f"that the dashboard parser's whitelist rejects; a model writing "
        f"the bare form has the line dropped in silence")


def test_the_benchmark_mirror_and_the_dashboard_agree():
    """The benchmark sees the routing a user would get; the two
    bare-slash whitelists must not drift apart (the mirror accepted
    four prefixes the dashboard dropped)."""
    from delfin.agent import benchmark_runner

    mirror = set(benchmark_runner._BARE_SLASH_PREFIXES)
    accepted = _dashboard_whitelist()
    # The mirror holds "prefix " with a trailing space (it requires the
    # argument); the dashboard holds bare "/tab". Compare stems.
    stems = {p.strip() for p in mirror}
    for stem in stems:
        assert stem in accepted, (
            f"benchmark mirror accepts bare {stem!r} but the dashboard "
            f"parser rejects it — the suite no longer measures the product")
