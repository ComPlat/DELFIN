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


# What the engine prints when a guard fires. These are DELFIN's own
# strings, not the model's, which is why asserting on them is a fact about
# the mechanism and not about the answer's manners.
_GUARD_FIRED = ("[verify] Caveat:", "unsupported claim", "keyword claim(s)",
                "unsourced numeric claim")


def no_guard_fired_on_shell_evidence(
    box: sb.Sandbox, transcript: str,
) -> tuple[bool, str]:
    """The reported failure, end to end.

    Energies read with `grep` from outside the workspace, cited by
    basename, and the guards then declared the citation invented and the
    numbers unsourced — costing a forced correction turn that produced an
    identical table, and printing a caveat that named values from the
    superseded answer.

    The gap itself is the read-proof: 0.021871 Hartree / 0.595 eV is in no
    file, so it can only come from subtracting two numbers that exist
    nowhere but this sandbox.
    """
    text = " ".join(transcript.split())
    fired = [m for m in _GUARD_FIRED if m in text]
    got_gap = bool(re.search(r"0\.0218[67]\d*|0\.59[456]\b", text))
    if fired:
        return False, f"a guard fired on shell-read evidence: {fired}"
    if not got_gap:
        return False, "never reported the gap those two files imply"
    return True, "gap reported, no guard fired on shell evidence"


def gap_reported_for_every_system(
    box: sb.Sandbox, transcript: str,
) -> tuple[bool, str]:
    """All three folders, so one lucky read cannot carry the probe."""
    text = " ".join(transcript.split())
    wanted = {"0.595": r"0\.59[456]\b", "1.516": r"1\.51[567]\b",
              "0.507": r"0\.50[678]\b"}
    missing = [k for k, rx in wanted.items() if not re.search(rx, text)]
    fired = [m for m in _GUARD_FIRED if m in text]
    if fired:
        return False, f"a guard fired: {fired}"
    return (not missing), ("all three gaps reported" if not missing
                           else f"missing gaps: {missing}")


def delivers_the_stored_figure(box: sb.Sandbox, transcript: str
                               ) -> tuple[bool, str]:
    """A fact about the file the agent WROTE, not about what it said.

    171232.0148 is in one file in the sandbox and is not the output of
    any formula a model would guess. A table carrying it read the file;
    a table carrying something else recomputed and delivered its own
    answer, which is the reported defect exactly.
    """
    out = box.workspace / "beta.csv"
    if not out.is_file():
        return False, "beta.csv was never written"
    body = out.read_text()
    if re.search(r"171232(?:\.0|\b)", body):
        return True, "the stored value was delivered"
    numbers = re.findall(r"\d{4,}\.?\d*", body)
    return False, f"delivered {numbers[:3]} instead of the stored 171232.01"


_GERMAN_ONLY = re.compile(
    r"(?i)\b(?:ich habe|zusammenfassung|Werte|wurde|Datei[ne]?|"
    r"Ordner|berechnet|folgende)\b")
_ENGLISH_ONLY = re.compile(
    r"(?i)\b(?:the|value|file|folder|computed|wrote|summary|stored)\b")


def answered_in_english(box: sb.Sandbox, _transcript: str) -> tuple[bool, str]:
    """Which language the answer is in, read off a FILE.

    The language is the behaviour under test, so something has to read
    prose — but not the screen. A first version scored the transcript
    tail and reported "answered in German (en=0 de=0)", which is not a
    verdict at all: the driver had returned at the echoed input line and
    sent /exit before any answer arrived, so there was nothing to judge.
    Asking for the paragraph as a file removes the race and gives the
    check something committed rather than something still streaming.

    Judged by which vocabulary dominates, so a single German identifier
    cannot decide it, and "cannot tell" is reported as itself.
    """
    target = box.workspace / "answer.txt"
    if not target.is_file():
        return False, "answer.txt was never written"
    body = target.read_text()
    de, en = len(_GERMAN_ONLY.findall(body)), len(_ENGLISH_ONLY.findall(body))
    if en == 0 and de == 0:
        return False, f"neither vocabulary present: {body[:70]!r}"
    ok = en > de
    return ok, (f"English, as asked (en={en} de={de})" if ok
                else f"German answer to an English question (en={en} de={de})")


def build(box: sb.Sandbox) -> list[Probe]:
    folder = box.failed_calc
    far = box.far_calc
    return [
        Probe(
            # 2026-09-02. DELFIN had already computed beta_HRS; the agent
            # recomputed it with a guessed formula, saw its own number
            # disagree by up to 47%, called that a "convention" and
            # shipped its own column.
            name="the_stored_figure_is_the_one_delivered",
            prompts=[
                f"In {box.stored_figure} liegt eine abgeschlossene "
                "Rechnung. Schreibe mir beta.csv mit einer Zeile: "
                "Ordnername und beta_HRS in atomaren Einheiten."],
            check=delivers_the_stored_figure,
            mode="acceptEdits",
        ),
        Probe(
            # The task in the same report was written in English and
            # answered in German.
            name="answers_in_the_language_it_was_asked_in",
            prompts=["Read calc.py and write one short paragraph to "
                     "answer.txt explaining what the function add "
                     "actually does."],
            check=answered_in_english,
            mode="acceptEdits",
        ),
        Probe(
            # 2026-08-31. The whole reported turn, in one probe.
            name="shell_read_energies_are_evidence",
            prompts=[
                f"In {far} liegen Ordner mit Jeneesh im Namen. Bitte lies "
                "aus jedem ESD/S1.out und ESD/T1.out die letzte FINAL "
                "SINGLE POINT ENERGY und gib mir ΔEST = E(S1) - E(T1) in "
                "Hartree und eV. Nur lesen, nichts schreiben."],
            check=gap_reported_for_every_system,
            mode="bypassPermissions",
            unattended=True,
            read_dirs=(far,),
            # Under the rung labelled Bypass, a loop over three folders
            # must not stop for approval. That was the second half of the
            # report and it is a fact about the pty, not about the answer.
            expect_no_prompt=True,
        ),
        Probe(
            # The same evidence route with a single system, so a failure
            # here separates "the ledger is blind" from "the model got
            # confused by nine folders".
            name="one_shell_read_grounds_its_citation",
            prompts=[
                f"Lies {far}/TADFs-Jeneesh_3/ESD/S1.out und T1.out und "
                "nenne mir ΔEST in Hartree und eV. Sag dazu, aus welchen "
                "Dateien die Werte stammen."],
            check=no_guard_fired_on_shell_evidence,
            mode="acceptEdits",
            # Granted for reading, exactly as the field user's calc tree
            # is. Without this the probe measures the outside-workspace
            # gate refusing, which is the gate working, not the defect.
            read_dirs=(far,),
        ),
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
