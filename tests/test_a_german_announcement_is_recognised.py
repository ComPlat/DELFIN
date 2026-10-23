"""A German first-person announcement draws the follow-up note; a German
hand-back, completion, report, question or wait does not.

Input: the closing text of a turn, as ``announced_action`` receives it.
Output: whether the turn ended announcing an action it did not take.

Before this file the detector knew German only as a hand-back ("sobald",
"wenn", "ich warte"); "Ich prüfe jetzt die Logs." drew no note while
"Let me check the logs now." did. Replayed on the 434 unique closing
texts of the shared bug archive (245 of them German), the German count
was 0 of 245.

The cost asymmetry is the one the English table is built on: a false
note costs a model round and makes the model invent an action, a missed
one leaves the session idle. So the German forms are first person plus
an immediate or future marker, and every closing that addresses the
user (imperative, request, question, offer, approval) is a hand-back.

Fixtures are written for this test or taken from closing sentences of
archived sessions with paths and names replaced.
"""

from __future__ import annotations

import pytest

from delfin.agent import turn_continuation as TC


ANNOUNCES = [
    "Ich prüfe jetzt die Logs.",
    "Lass mich das kurz testen.",
    "Ich schaue mir das an.",
    "Als Nächstes starte ich die Rechnung.",
    "Ich werde den Branch pushen.",
    "Der Import ist an zwei Stellen falsch. Ich korrigiere ihn jetzt.",
    "Lass mich zuerst prüfen, was die Orakel tatsächlich liefern:",
    "Ich analysiere zuerst, was in `/trans-CH2OH` schiefgelaufen ist.",
    "Ich analysiere die Rechnung jetzt mit den internen Tools. Zuerst "
    "prüfe ich den Status der Dateien im Ordner.",
    "Ich extrahiere die Tensoren aus allen Ordnern. Zuerst verschaffe ich "
    "mir einen Überblick über die Struktur.",
    "Ich korrigiere die dokumentierten Fehler. Zuerst lese ich die "
    "betroffenen Zeilen, um sauber zu editieren:",
    "Ich beginne mit der Implementierung:",
    "Die Aufgabenliste ist erstellt. Ich beginne jetzt mit Task 1.",
    "Jetzt führe ich die Tests noch einmal aus.",
    "Ich kümmere mich gleich darum.",
    "Ich mache mich an den Fix für den Parser.",
    "Als Nächstes werde ich die Ergebnistabelle erstellen.",
    "Lassen Sie mich die Konfiguration kurz prüfen.",
    "Ich sehe mir die Ausgabe von ORCA genauer an.",
    "Der Test schlägt fehl. Ich baue jetzt einen Minimalfall.",
    "Ich schreibe jetzt den roten Test.",
    "Ich lese gleich die CONTROL.txt nach.",
    "Ich teste das sofort mit xtb.",
    "Nun pushe ich den Branch auf origin.",
]

DOES_NOT_ANNOUNCE = [
    # completion
    "Ich habe die Logs geprüft.",
    "Ich habe den Branch gepusht und die Tests sind grün.",
    "Ich habe die beiden Fakten im Memory gespeichert. Kein Code, kein "
    "Commit.",
    # a report in the present tense
    "Die Rechnung läuft jetzt.",
    "Der Job ist jetzt in der Warteschlange.",
    "Ich sehe jetzt, dass der Fehler im Parser liegt.",
    "Ich lese aus dem Log, dass die Rechnung jetzt konvergiert.",
    "Ich ignoriere das Ereignis weiter — an meiner Arbeit ist nichts offen.",
    "Ich lasse den Test bewusst rot stehen und ändere nichts an fremdem "
    "Code.",
    # a question to the user
    "Soll ich jetzt die Rechnung starten?",
    "Soll ich den Push jetzt ausführen?",
    "Möchtest du, dass ich direkt den PR eröffne?",
    "Oder soll ich zuerst die Logs prüfen?",
    "Was soll ich als Nächstes tun?",
    # a request, an imperative, an offer
    "Lass mich wissen, ob ich den PR eröffnen soll.",
    "Lassen Sie mich wissen, welche Variante Sie bevorzugen.",
    "Sag Bescheid, dann lege ich den Worktree an.",
    "Sag kurz Go, dann pushe ich alle drei Commits.",
    "Bestätige mit einem Wort, dann führe ich ihn aus.",
    "Bitte klicke im Dashboard auf Plan akzeptieren, damit ich jetzt "
    "starten kann.",
    "Danach prüfe ich gern per git ls-remote, dass der Branch steht.",
    "Wenn du willst, fixe ich das jetzt als eigenen Commit.",
    "Falls du zustimmst, starte ich jetzt die Rechnung.",
    # waiting, or a step that depends on someone else
    "Ich warte auf deine Freigabe.",
    "Ich melde mich, sobald der Job fertig ist.",
    "Ich beende diese Runde und melde mich mit dem Ergebnis.",
    "Sobald du den Plan akzeptierst, lege ich sofort los.",
    "Ich werde auf deine Rückmeldung warten.",
    "Der Abschluss meldet sich automatisch — danach erstelle ich die "
    "Ergebnistabelle.",
    "Ich fahre fort mit Abschnitt 7.3, sofern du nichts anderes vorgibst.",
    "Plan ist eingereicht. Nach deiner Freigabe setze ich ihn jetzt um.",
    # a negation, a habit
    "Ich werde nichts pushen, solange CI rot ist.",
    "Ich werde das künftig beachten.",
    "Das werde ich mir merken.",
    "Ich prüfe das jetzt nicht weiter, der Befund ist eindeutig.",
]


@pytest.mark.parametrize("text", ANNOUNCES)
def test_a_german_announcement_is_followed_up(text):
    assert TC.announced_action(text), text
    assert TC.should_continue(text, 0)[0], text


@pytest.mark.parametrize("text", DOES_NOT_ANNOUNCE)
def test_a_german_closing_that_announces_nothing_is_not(text):
    assert not TC.announced_action(text), text
    assert TC.should_continue(text, 2) == (False, ""), text


def test_a_german_announcement_in_an_earlier_paragraph_does_not_count():
    text = ("Ich lese zuerst die Datei.\n\n"
            "Die Datei definiert drei Gates; Abschnitt 7 fehlt.")
    assert not TC.announced_action(text)


def test_an_english_announcement_beside_german_text_still_counts():
    assert TC.announced_action(
        "Die Tests schlagen am Import fehl. Let me fix the import now.")
    assert not TC.announced_action(
        "Let me fix the import. Sag Bescheid, ob das passt.")


def test_the_plan_mode_sentence_sees_a_german_announcement():
    turn = "Ich pushe jetzt den Branch.Der Aufruf wurde abgelehnt."
    assert TC.blocked_by_plan_mode(turn, refusals=1)
