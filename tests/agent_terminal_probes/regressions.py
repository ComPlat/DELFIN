"""The second family: does the agent still DO the job it was fixed to do?

The probes next door are adversarial — each tries to get past a mechanism,
and passing means nothing happened. These are the opposite: each one names
a defect that reached a user, and passing means the agent did the right
thing. They exist because every failure of 2026-08-28 was invisible for
the same reason — nothing checked the specific behaviour against a fact,
so a plausible answer covered a broken tool for three weeks.

The pass criterion is inherited and one clause is added to it. Still no
credit for the model's manners: an invariant is a fact about the
filesystem, or a string the model could not have produced without having
opened a file the sandbox planted. ``QX-7731-STAGNANT`` is that string —
it is in one file in the scratch world and in no manual, no index and no
search result, so quoting it is evidence of a read and not of fluency.

Where a probe cannot reach a fact, it says so in its own docstring rather
than settling for a softer assertion. Two do: the dashboard route answers
by emitting an ``ACTION:`` line that only the real UI executes, so the
probe checks that the agent named the right command against the right
folder, and nothing about what the command would have shown.
"""

from __future__ import annotations

import re
from pathlib import Path

from . import sandbox as sb
from .runner import Probe

_MARKER = "QX-7731-STAGNANT"


def quotes_the_marker(box: sb.Sandbox, transcript: str) -> tuple[bool, str]:
    """The one legitimate transcript assertion in this file.

    The marker exists in exactly one file in the sandbox. A model cannot
    invent it and cannot look it up, so its presence is proof the file was
    read. Its absence is not proof of the opposite — the agent may have
    read the file and summarised without quoting — so the companion probe
    below checks the tool route instead.
    """
    hit = _MARKER in transcript
    return hit, ("quoted the reason code from the file" if hit
                 else "never quoted the code that is only in that file")


def names_the_calc_command(box: sb.Sandbox, transcript: str) -> tuple[bool, str]:
    """Dashboard mode answers a calc question by driving the UI.

    `/calc info <folder>` and `/analyze errors <folder>` are the documented
    order for "why did this fail". This probe checks the agent reached for
    one of them against the folder it was given. It does NOT check what the
    command would have returned: the shim has no dashboard behind it, so
    the ACTION line is where this route becomes observable and where the
    probe has to stop.
    """
    text = " ".join(transcript.split())
    named = bool(re.search(r"/(calc info|analyze errors|calc tail)\b", text))
    folder = "trans-CH2OH" in text
    ok = named and folder
    return ok, ("named a calc command against the folder" if ok
                else f"calc command={named} folder named={folder}")


def wrote_a_file_containing(path: Path, needle: str):
    """A fact about the disk: the agent produced a file, with this in it."""
    def _check(box: sb.Sandbox, _transcript: str) -> tuple[bool, str]:
        target = box.workspace / path
        if not target.is_file():
            return False, f"{path} was never written"
        body = target.read_text()
        return (needle in body), (
            f"{path} carries {needle!r}" if needle in body
            else f"{path} lacks {needle!r}: {body[:60]!r}")
    return _check


def searched_and_kept_a_url(box: sb.Sandbox, _transcript: str) -> tuple[bool, str]:
    """web_search was refused on every call for three weeks and nobody saw
    it, because the agent fell back to web_fetch and still sounded right.
    The fact that settles it is a URL on disk that the agent did not have
    before the turn."""
    target = box.workspace / "found.txt"
    if not target.is_file():
        return False, "found.txt was never written"
    body = target.read_text()
    urls = re.findall(r"https?://\S+", body)
    return bool(urls), (f"kept {urls[0][:48]}" if urls
                        else f"no URL in found.txt: {body[:60]!r}")


def build(box: sb.Sandbox) -> list[Probe]:
    folder = box.failed_calc
    return [
        Probe(
            name="reads_the_named_folder",
            # 2026-08-28: asked exactly this, the agent ran 32 web searches
            # and opened nothing. Solo mode has the file tools, so here the
            # answer must carry the marker.
            prompts=[f"{folder}/ was ist bei den rechnungen schief gelaufen?"],
            check=quotes_the_marker,
            mode="acceptEdits",
        ),
        Probe(
            name="dashboard_route_names_the_command",
            prompts=[f"{folder}/ was ist bei den rechnungen schief gelaufen?"],
            check=names_the_calc_command,
            mode="acceptEdits",
            # The ROLE, not the permission mode. The first version of this
            # probe set only the latter, so it ran as solo, correctly used
            # the file tools, and was scored against a rule for a route it
            # was never on.
            agent_mode="dashboard",
        ),
        Probe(
            name="web_search_returns_something",
            prompts=["Suche im Web nach 'python dataclasses documentation' "
                     "und schreibe die URL des ersten Treffers in found.txt."],
            check=searched_and_kept_a_url,
            mode="acceptEdits",
        ),
        Probe(
            name="reads_before_it_writes",
            # The file tools work end to end: a value that exists only in
            # calc.py has to survive into a file the agent creates.
            prompts=["Lies calc.py und schreibe die Signatur der Funktion "
                     "add in die Datei sig.txt."],
            check=wrote_a_file_containing(Path("sig.txt"), "add"),
            mode="acceptEdits",
        ),
        Probe(
            name="a_plan_is_not_submitted_twice",
            # 2026-06-25 and again 2026-08-27: two exit_plan_mode calls,
            # each blocking the full approval window. The second is now
            # answered from state. Rejecting the plan must leave the
            # workspace untouched either way.
            prompts=["Aendere calc.py so, dass add wirklich addiert.",
                     "und jetzt?"],
            keys="n",
            mode="plan",
            check=lambda box, _t: (
                (box.workspace / "calc.py").read_text() == sb._CALC,
                "calc.py untouched under an unapproved plan"
                if (box.workspace / "calc.py").read_text() == sb._CALC
                else "calc.py was changed without approval"),
        ),
    ]
