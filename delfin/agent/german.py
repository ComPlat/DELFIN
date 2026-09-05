"""What the framework has to recognise, and write back, in German.

This framework's users write German and work on German files. The code
and its comments are English; the *vocabulary* below is German because it
is what a user types and what a customer reads. Collecting it in one
module exists for three reasons, each of them a defect that was shipped:

* **Two guards disagreed about one sentence.** "Die Summe liegt bei rund
  45.000 EUR." was an honest hedge to the office guard and a confident
  claim to the answer guard, in the same answer, because each carried its
  own uncertainty list. There is now one list and both import it.
* **Two slugifiers produced three answers for one title.** ``Müller GmbH``
  and ``Möller GmbH`` both slugged to ``mller-gmbh`` — one customer's note
  overwrote the other's — and a second implementation elsewhere returned
  ``m-ller-gmbh`` for the same string. Umlauts are transliterated before
  anything is stripped, once, here.
* **German verbs split.** German puts the prefix of a separable verb at
  the end of the clause, so a pattern for ``anpass\\w*`` cannot see
  "Passe die Tabelle an" — and a write task went out "verified" with
  nothing written. The verb tables below match the stem wherever the
  prefix ended up.

Everything here is data plus small pure helpers. Nothing reads a file,
nothing raises.
"""

from __future__ import annotations

import re
from datetime import date, datetime, time
from typing import Any, Optional

# ---------------------------------------------------------------------------
# Umlauts: transliterate before stripping
# ---------------------------------------------------------------------------

# ä -> ae, not ä -> "". Dropping the diaeresis collapses distinct names
# ("Müller"/"Möller", "Größe"/"Grüße") onto one key, and a slug is a key:
# it names the file, it is what a merge and a delete match on, and it is
# the link text the user reads.
_TRANSLITERATION = (
    ("ä", "ae"), ("ö", "oe"), ("ü", "ue"),
    ("Ä", "Ae"), ("Ö", "Oe"), ("Ü", "Ue"),
    ("ß", "ss"),
    # Neighbouring languages appear in German address data often enough
    # that losing them the same way is the same defect.
    ("é", "e"), ("è", "e"), ("ê", "e"), ("à", "a"), ("á", "a"), ("â", "a"),
    ("ç", "c"), ("í", "i"), ("ï", "i"), ("ñ", "n"), ("ó", "o"), ("ô", "o"),
    ("ø", "oe"), ("å", "aa"), ("ú", "u"), ("û", "u"), ("ý", "y"),
)


def transliterate(text: Any) -> str:
    """German (and neighbouring) letters as their ASCII spelling."""
    out = str(text or "")
    for source, target in _TRANSLITERATION:
        if source in out:
            out = out.replace(source, target)
    return out


_SLUG_STRIP_RE = re.compile(r"[^a-z0-9]+")


def slugify(text: Any, max_len: int = 60, fallback: str = "eintrag") -> str:
    """One slug rule for the whole framework.

    Transliterate, lowercase, collapse everything else to a single dash,
    trim. ``Müller GmbH`` -> ``mueller-gmbh`` and ``Möller GmbH`` ->
    ``moeller-gmbh``: two customers, two keys.
    """
    lowered = transliterate(text).lower()
    slug = _SLUG_STRIP_RE.sub("-", lowered).strip("-")
    slug = slug[: max(1, int(max_len or 1))].rstrip("-")
    return slug or fallback


# ---------------------------------------------------------------------------
# Uncertainty: one vocabulary, imported by every guard
# ---------------------------------------------------------------------------
#
# A guard exists to stop CONFIDENT fabrication. An answer that says it is
# guessing has already disclosed its grounding, so every scanner skips a
# claim whose sentence carries a marker from this list.
#
# The German half used to be adverbs only, while the English half had a
# full first-person vocabulary ("i think", "from memory", "iirc"). So
# "Ich schätze die Summe auf 45.000 EUR." was read as an assertion and
# "I'd guess 45,000" was not — in a framework whose users write German.
# The first-person forms below are the fix, and they are the reason this
# list is shared rather than copied.
#
# Abbreviated hedges END in a period. A trailing ``\b`` after "." never
# matches (both sides are non-word), which is how "ca." — the most
# ordinary hedge there is — sat in a list and could not fire. They stand
# OUTSIDE the ``\b``-wrapped alternation.
HEDGE_RE = re.compile(
    r"(?i)\bca\.|\bz\.\s?T\.|\bz\.\s?B\.|\bggf\.|\bevtl\.|\bs\.\s?o\.|"
    r"\b(?:"
    # -- German, adverbial ------------------------------------------------
    r"vermutlich|wahrscheinlich|m(?:ö|oe)glicherweise|vielleicht|"
    r"ungef(?:ä|ae)hr|etwa|circa|zirka|sch(?:ä|ae)tzungsweise|"
    r"gesch(?:ä|ae)tzt|grob|rund|knapp|(?:ü|ue)berschl(?:ä|ae)gig|"
    r"n(?:ä|ae)herungsweise|gr(?:ö|oe)ssenordnung|gr(?:ö|oe)ßenordnung|"
    r"unsicher|ungepr(?:ü|ue)ft|unverifiziert|womöglich|womoeglich|"
    r"angenommen|annahme|hypothetisch|fiktiv|platzhalter|"
    r"beispielsweise|beispiel|ohne\s+gew(?:ä|ae)hr|"
    # -- German, explicit non-measurement ---------------------------------
    r"nicht\s+(?:gepr(?:ü|ue)ft|verifiziert|(?:ü|ue)berpr(?:ü|ue)ft|sicher|"
    r"nachgesehen|nachgeschaut|berechnet|gerechnet|nachgerechnet|"
    r"kontrolliert)|"
    # -- German, Konjunktiv II: how German says "not measured" ------------
    r"w(?:ä|ae)re|w(?:ä|ae)ren|w(?:ü|ue)rde|w(?:ü|ue)rden|l(?:ä|ae)ge|"
    r"erg(?:ä|ae)be|betr(?:ü|ue)ge|k(?:ö|oe)nnte|d(?:ü|ue)rfte|"
    r"scheint|offenbar|"
    # -- German, first person: the half that was missing -------------------
    r"ich\s+(?:denke|glaube|vermute|sch(?:ä|ae)tze|meine|nehme\s+an|"
    r"gehe\s+davon\s+aus|erinnere\s+mich|bin\s+mir\s+nicht\s+sicher|"
    r"wei(?:ß|ss)\s+(?:es\s+)?nicht)|"
    r"so(?:weit|viel)\s+ich\s+wei(?:ß|ss)|meines\s+wissens|"
    r"wenn\s+ich\s+mich\s+(?:recht|richtig)\s+erinnere|"
    r"nach\s+meiner\s+erinnerung|"
    r"aus\s+dem\s+(?:ged(?:ä|ae)chtnis|kopf)|"
    # -- English ----------------------------------------------------------
    r"probably|likely|perhaps|maybe|possibly|approximately|roughly|around|"
    r"unsure|uncertain|unverified|presumably|assumed|example|placeholder|"
    r"not\s+(?:sure|certain|verified|checked|confirmed)|"
    r"i\s+(?:think|believe|guess|assume|recall|suspect|expect)|"
    r"have\s*n(?:o|')t\s+(?:checked|verified|looked|read|confirmed)|"
    r"did\s*n(?:o|')t\s+(?:check|verify|look|read|confirm)|"
    r"without\s+(?:checking|looking|verifying)|"
    r"from\s+memory|iirc|if\s+i\s+recall|"
    r"might|may\s+be|should\s+be|seems?|estimated?"
    r")\b"
)


# ---------------------------------------------------------------------------
# Separable verbs: the stem stays, the prefix moves
# ---------------------------------------------------------------------------
#
# "Trage die Werte ein" and "eintragen" are the same verb. In verb-second
# order the prefix lands at the end of the clause, so a pattern written
# for the infinitive sees nothing. Each entry below is
# ``(prefix, attached stems, finite forms)`` and produces two alternatives:
#
#   attached   ``ein(?:ge)?trag<suffix>``  — eintragen, eingetragen, Eintragung
#   split      ``trage ... ein``           — with the prefix at a clause end
#
# The split form REQUIRES the prefix to close its clause. Without that
# "ich trage ein Kleid" would read as a booking, because "ein" is also
# the article.

# Verb endings, as a closed set rather than ``\w*``: "angelegenheit"
# starts with "angeleg" and is not a booking.
_VERB_TAIL = r"(?:e|en|et|st|t|te|ten|test|tet|end|ung|ungen)?"

# Words that continue a clause without ending the sentence. "Trage die
# Werte ein und speichere die Datei" keeps the prefix at a clause end.
_CLAUSE_END = (
    r"(?:\s*[.,;:!?)]|\s+(?:und|oder|sowie|damit|dann|danach|"
    r"anschlie(?:ß|ss)end|bevor|nachdem|bitte)\b|\s*$)"
)

# How far the prefix may sit from its stem. One clause of an ordinary
# office sentence; beyond that the two words are not the same verb.
_PREFIX_WINDOW = 60


def _separable(prefix: str, attached: tuple[str, ...],
               finite: tuple[str, ...]) -> str:
    """Regex source matching one separable verb, joined or split."""
    joined = "|".join(attached)
    forms = "|".join(sorted(finite, key=len, reverse=True))
    return (
        rf"\b{prefix}(?:ge)?(?:{joined}){_VERB_TAIL}\b"
        rf"|\b(?:{forms})\b(?=[^;:!?\n]{{0,{_PREFIX_WINDOW}}}?"
        rf"\b{prefix}\b{_CLAUSE_END})"
    )


def _plain(*stems: str) -> str:
    """Regex source for verbs whose prefix never moves (übertragen,
    verbuchen) and for plain stems."""
    return r"\b(?:" + "|".join(stems) + rf"){_VERB_TAIL}\b"


def _plain_unless(stem: str, *prefixes: str) -> str:
    """A plain stem that is NOT the finite half of a separable verb.

    "Setze das Datum auf den 31.07." changes a file; "Setze die Recherche
    fort" is the same finite form of a different verb, and reading it as a
    write accuses a finished reading task of having written nothing. The
    prefix is what tells them apart, and in verb-second order it sits at
    the end of the clause -- so this is :func:`_separable` inverted: match
    the stem only when none of these prefixes closes the clause after it.
    """
    escaped = "|".join(prefixes)
    return (rf"\b{stem}{_VERB_TAIL}\b"
            rf"(?![^;:!?\n]{{0,{_PREFIX_WINDOW}}}?\b(?:{escaped})"
            rf"{_CLAUSE_END})")


# The office work a German user asks for. These CHANGE a file, which is
# what separates them from the reading half below.
_SEPARABLE_WRITE: tuple[tuple[str, tuple[str, ...], tuple[str, ...]], ...] = (
    ("ein", ("trag", "träg"), ("trage", "trägst", "trägt", "tragt", "trag")),
    ("ein", ("füg",), ("füge", "fügst", "fügt", "füg")),
    ("ein", ("geb", "gib"), ("gebe", "gibst", "gibt", "gebt", "gib")),
    ("ein", ("pfleg",), ("pflege", "pflegst", "pflegt", "pfleg")),
    ("ein", ("setz",), ("setze", "setzt", "setz")),
    ("ein", ("buch",), ("buche", "buchst", "bucht", "buch")),
    ("an", ("leg",), ("lege", "legst", "legt", "leg")),
    ("an", ("pass",), ("passe", "passt", "pass")),
    ("aus", ("füll",), ("fülle", "füllst", "füllt", "füll")),
    ("aus", ("fertig",), ("fertige", "fertigst", "fertigt")),
    ("hinzu", ("füg",), ("füge", "fügst", "fügt", "füg")),
    ("nach", ("trag", "träg"), ("trage", "trägst", "trägt", "trag")),
    ("nach", ("pfleg",), ("pflege", "pflegst", "pflegt")),
    ("ab", ("leg",), ("lege", "legst", "legt", "leg")),
    ("ab", ("speicher",), ("speichere", "speicherst", "speichert")),
    ("ab", ("hak",), ("hake", "hakst", "hakt", "hak")),
    ("um", ("buch",), ("buche", "buchst", "bucht", "buch")),
    ("auf", ("nehm", "nimm"), ("nehme", "nimmst", "nimmt", "nimm")),
    ("auf", ("teil",), ("teile", "teilst", "teilt", "teil")),
    ("zu", ("ordn",), ("ordne", "ordnest", "ordnet")),
    ("zusammen", ("führ",), ("führe", "führst", "führt")),
    # A rename IS a write, and only the joined form was ever matched --
    # "Benenne die Datei um", which is how somebody actually asks for it,
    # matched nothing at all.
    ("um", ("benenn",), ("benenne", "benennst", "benennt", "benenn")),
)

# Verbs that only LOOK at a file. "Rechne die Summe aus" is answered by
# reading and arithmetic; the answer is a number in the chat, not a
# changed file. Calling it a write would accuse an honest, finished task
# of having written nothing — which is how a check stops being read.
_SEPARABLE_READ: tuple[tuple[str, tuple[str, ...], tuple[str, ...]], ...] = (
    ("aus", ("rechn",), ("rechne", "rechnest", "rechnet")),
    ("aus", ("wert",), ("werte", "wertest", "wertet")),
    ("aus", ("les",), ("lese", "liest", "lest", "lies")),
    ("durch", ("rechn",), ("rechne", "rechnest", "rechnet")),
    ("durch", ("seh",), ("sehe", "siehst", "sieht", "sieh", "seh")),
    ("nach", ("rechn",), ("rechne", "rechnest", "rechnet")),
    ("nach", ("seh", "schau"), ("sehe", "siehst", "sieht", "sieh", "schau",
                                "schaue", "schaust", "schaut")),
    ("an", ("schau", "seh"), ("schaue", "schaust", "schaut", "schau",
                              "sehe", "siehst", "sieht", "sieh")),
    ("gegen", ("rechn",), ("rechne", "rechnest", "rechnet")),
    ("ab", ("gleich",), ("gleiche", "gleichst", "gleicht")),
    ("zusammen", ("fass",), ("fasse", "fasst", "fass")),
)

# Office verbs whose prefix never moves, plus plain stems. Kept apart from
# the coding vocabulary the write check already had (erstell, refaktor,
# implementier): a bookkeeping task uses none of those words.
_PLAIN_WRITE = (
    "übertrag", "überträg", "uebertrag", "ueberträg",
    "verbuch", "verrechn", "berichtig", "vervollständig", "vervollstaendig",
    "sortier", "gruppier", "kategorisier", "beschrift", "markier",
    "quittier", "stornier", "fakturier", "archivier", "duplizier",
    "protokollier", "nummerier", "durchnummerier", "formatier",
    "eintrag", "austrag", "umbuch",
    "buche", "buchen", "gebucht",
    "fülle", "füllen", "gefüllt", "fuelle", "fuellen",
    # Measured against the shipped matcher, not guessed: "Setze das Datum
    # auf den 31.07.2026", "Hinterlege die Bankverbindung" and "Streiche
    # die Position 7" were all read as tasks that want no write, so a
    # session that had only READ the file reported them verified.
    "hinterleg", "streich",
)

# "setzen" needs the inverted-prefix guard: see _plain_unless.
_PLAIN_WRITE_GUARDED = (
    _plain_unless("setz", "fort", "voraus", "durch", "aus", "ab", "um"),
)

_PLAIN_READ = (
    "auswert", "ermittl", "ermittel", "berechn", "errechn", "ausrechn",
    "nachrechn", "abgleich", "gegenrechn", "zusammenfass", "aufliste",
    "kontrollier", "durchseh", "sichte", "zähl", "zaehl",
    # Same measurement, other direction. A read task that matches no read
    # verb is not merely unrecognised: the artefact branch runs on it, and
    # "Öffne den Bericht" was reported unmet for having produced no PDF.
    # "benennen" and "umbenennen" are writes and the word boundary keeps
    # them out of "nenn".
    "zeig", "öffn", "oeffn", "nenn",
)

# A question is a read intent even when it carries no read verb, and the
# most ordinary German way to ask for a figure is exactly that: "Wie hoch
# ist die Gesamtsumme?", "Was steht in Spalte D?". Left out of the verb
# tables because they are not verbs -- naming them separately keeps
# _PLAIN_READ honest about what it holds.
GERMAN_ASKING_SOURCE = (
    r"\bwie\s+hoch\b|\bwie\s+viele?\b|\bwieviel\w*\b|\bwas\s+steht\b"
    r"|\bwas\s+ergibt\b|\bwelche\s+(?:summe|beträge|betraege|posten|"
    r"zeilen|werte)\b|\bberichte?\s+mir\b"
)

#: Everything a German office task can say that means "change the file".
GERMAN_WRITE_VERB_SOURCE = "|".join(
    [_separable(*entry) for entry in _SEPARABLE_WRITE]
    + [_plain(*_PLAIN_WRITE)] + list(_PLAIN_WRITE_GUARDED)
)

#: Everything that means "look at the file and tell me" — the verbs, and
#: the question forms that say the same thing without one.
GERMAN_READ_VERB_SOURCE = "|".join(
    [_separable(*entry) for entry in _SEPARABLE_READ]
    + [_plain(*_PLAIN_READ), GERMAN_ASKING_SOURCE]
)

GERMAN_WRITE_VERB_RE = re.compile("(?i)" + GERMAN_WRITE_VERB_SOURCE)
GERMAN_READ_VERB_RE = re.compile("(?i)" + GERMAN_READ_VERB_SOURCE)


# ---------------------------------------------------------------------------
# What a delegate says it did
# ---------------------------------------------------------------------------
#
# A report claims a change in the perfect: "Ich habe die Werte
# eingetragen." Listed rather than derived, because German participles are
# not derivable from the stem (tragen -> getragen, legen -> gelegt).
GERMAN_WRITE_PARTICIPLES: tuple[str, ...] = (
    "geschrieben", "erstellt", "eingetragen", "nachgetragen", "ausgetragen",
    "übertragen", "uebertragen", "angelegt", "abgelegt", "eingefügt",
    "eingefuegt", "eingegeben", "eingepflegt", "ergänzt", "ergaenzt",
    "gespeichert", "abgespeichert", "verbucht", "umgebucht", "eingebucht",
    "gebucht", "ausgefüllt", "ausgefuellt", "sortiert", "gruppiert",
    "angepasst", "geändert", "geaendert", "überarbeitet", "ueberarbeitet",
    "gelöscht", "geloescht", "gestrichen", "korrigiert", "berichtigt",
    "markiert", "beschriftet", "archiviert", "exportiert", "storniert",
    "vervollständigt", "vervollstaendigt", "formatiert", "nummeriert",
)

#: English participles the mutation family was missing. It carried
#: "rewrote" and not "wrote", and no "saved" at all.
ENGLISH_WRITE_PARTICIPLES: tuple[str, ...] = (
    "wrote", "written", "saved", "stored", "appended", "inserted",
    "filled", "filled in", "filled out", "populated", "exported",
    "generated", "recorded", "entered", "sorted", "replaced",
)


# ---------------------------------------------------------------------------
# The label on a closing row
# ---------------------------------------------------------------------------
#
# A bookkeeping sheet ends in a total, and an append lands BELOW it —
# outside the range the total sums, so the workbook's own number stops
# matching its rows. Finding that row means reading its label.
#
# The label is matched by TOKENS, not by a prefix plus a length clamp. The
# clamp ("summe" plus at most two characters) was silent on every label a
# German sheet actually carries: "Summe EUR", "Endsumme", "Gesamtbetrag",
# "Summe gesamt", "Gesamt netto".

# A token containing one of these is a total word, compounded or not:
# Endsumme, Zwischensumme, Gesamtbetrag, Jahressumme.
_TOTAL_TOKEN_PARTS = ("summe", "summen", "gesamt", "total", "saldo",
                      "endbetrag", "endergebnis")

# A token that is exactly one of these is a total word too — short ones
# that must not match as substrings ("sum" inside "consumption").
_TOTAL_TOKEN_WORDS = frozenset({
    "sum", "sums", "subtotal", "subtotals", "summe", "summen", "gesamt",
    "total", "totals", "saldo", "ergebnis", "endstand", "bilanz",
})

# Words that may stand NEXT to a total word without turning the row into
# prose: a currency, a tax qualifier, a scope. "Summe der Belege folgt" is
# a sentence and not a closing row, and "der/folgt" is what says so.
_TOTAL_QUALIFIERS = frozenset({
    "eur", "euro", "chf", "usd", "€", "$", "netto", "brutto", "inkl",
    "inkl.", "exkl", "exkl.", "zzgl", "zzgl.", "mwst", "mwst.", "ust",
    "ust.", "insgesamt", "gesamt", "jahr", "monat", "quartal", "periode",
    "grand", "net", "gross", "incl", "excl", "vat", "amount", "sum",
    "betrag", "summe",
})

_LABEL_TOKEN_RE = re.compile(r"[^\W\d_]+|[€$]", re.UNICODE)

# A closing-row label is short. Four tokens covers "Summe gesamt netto
# EUR"; past that it is a sentence.
_TOTAL_LABEL_MAX_TOKENS = 4


def is_total_row_label(label: Any) -> bool:
    """Whether *label* is what a sheet writes on its closing row.

    True for ``Summe``, ``Summe EUR``, ``Endsumme``, ``Gesamtbetrag``,
    ``Summe gesamt``, ``Gesamt netto``, ``Total``. False for prose that
    merely starts with one of those words — ``Summe der Belege folgt``.
    """
    tokens = _LABEL_TOKEN_RE.findall(str(label or "").lower())
    if not tokens or len(tokens) > _TOTAL_LABEL_MAX_TOKENS:
        return False
    saw_total = False
    for token in tokens:
        is_total = (token in _TOTAL_TOKEN_WORDS
                    or any(part in token for part in _TOTAL_TOKEN_PARTS))
        if is_total:
            saw_total = True
            continue
        if token not in _TOTAL_QUALIFIERS:
            return False
    return saw_total


# ---------------------------------------------------------------------------
# Writing German back out
# ---------------------------------------------------------------------------
#
# The read side of this framework detects decimal commas and day-first
# dates. The write side was ``str(value)``, so a Serienbrief went to a
# customer reading "1234.5 EUR ... fällig am 2026-07-31 00:00:00" and the
# run reported ok. These are the two primitives that stop that.

_THOUSANDS = "."
_DECIMAL = ","


def format_number(value: Any, decimals: Optional[int] = None,
                  grouping: bool = False, percent: bool = False) -> str:
    """A number as German writes it: ``1.234,50``.

    *decimals* fixes the number of places (from the cell's own format);
    ``None`` keeps whatever the value carries. *grouping* inserts the
    thousands separator, which the cell's format decides — Excel's
    ``#,##0.00`` means "grouped, two places" in every locale, and the
    locale only decides which characters those are.
    """
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if percent:
        number *= 100.0
    if decimals is None:
        text = repr(number)
        if text.endswith(".0"):
            text = text[:-2]
        if "e" in text or "E" in text:          # keep exponent form intact
            return text.replace(".", _DECIMAL)
        whole, _, frac = text.partition(".")
    else:
        text = f"{number:.{max(0, int(decimals))}f}"
        whole, _, frac = text.partition(".")
    negative = whole.startswith("-")
    digits = whole.lstrip("-")
    if grouping and len(digits) > 3:
        parts = []
        while len(digits) > 3:
            parts.insert(0, digits[-3:])
            digits = digits[:-3]
        parts.insert(0, digits)
        digits = _THOUSANDS.join(parts)
    out = ("-" if negative else "") + digits
    if frac:
        out += _DECIMAL + frac
    return out + (" %" if percent else "")


# Excel's date tokens, longest first, as strftime directives.
_DATE_TOKENS = (
    ("yyyy", "%Y"), ("yy", "%y"),
    ("mmmm", "%B"), ("mmm", "%b"),
    ("dddd", "%A"), ("ddd", "%a"),
    ("mm", "%m"), ("m", "%m"),
    ("dd", "%d"), ("d", "%d"),
    ("hh", "%H"), ("h", "%H"),
    ("ss", "%S"),
)

#: Day first, dots, four-digit year — what a German document shows when
#: the cell says nothing more specific.
DEFAULT_DATE_PATTERN = "%d.%m.%Y"
DEFAULT_DATETIME_PATTERN = "%d.%m.%Y %H:%M"


def excel_date_pattern(number_format: Any) -> str:
    """An Excel date format as a strftime pattern ("" when it is not one).

    ``DD.MM.YYYY`` -> ``%d.%m.%Y``. Minutes are not distinguished from
    months (Excel decides that by position); a bookkeeping sheet formats
    dates, not durations, so the month reading is taken.
    """
    text = str(number_format or "").strip()
    if not text or text.lower() == "general":
        return ""
    text = text.split(";")[0]
    text = re.sub(r"\[[^\]]*\]", "", text)
    text = re.sub(r'"[^"]*"', "", text)
    lowered = text.lower()
    if not any(ch in lowered for ch in "ymdh"):
        return ""
    if any(ch in lowered for ch in "#0?"):        # a number format, not a date
        return ""
    out = []
    i = 0
    while i < len(lowered):
        for token, directive in _DATE_TOKENS:
            if lowered.startswith(token, i):
                out.append(directive)
                i += len(token)
                break
        else:
            out.append(text[i])
            i += 1
    pattern = "".join(out).strip()
    return pattern if "%" in pattern else ""


def format_date(value: Any, number_format: Any = "") -> str:
    """A date/datetime as the cell says, or day-first German if it does not."""
    pattern = excel_date_pattern(number_format)
    if isinstance(value, datetime):
        if not pattern:
            pattern = (DEFAULT_DATE_PATTERN
                       if (value.hour, value.minute, value.second) == (0, 0, 0)
                       else DEFAULT_DATETIME_PATTERN)
    elif isinstance(value, date):
        if not pattern or "%H" in pattern:
            pattern = DEFAULT_DATE_PATTERN
    elif isinstance(value, time):
        return value.strftime("%H:%M" if not pattern else pattern)
    else:
        return str(value)
    try:
        return value.strftime(pattern)
    except (ValueError, TypeError):
        return value.strftime(DEFAULT_DATE_PATTERN)
