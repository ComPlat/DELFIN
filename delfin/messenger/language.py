"""Which language a person's phone is answered in: English unless they write German.

Input: the text of each message a verified phone sends. Output: the
language of that contact's phone texts (``en`` or ``de``), kept in the
gateway's store per contact.

Rules (deterministic, no model):

* a new contact starts in DELFIN_MESSENGER_LANG (default ``en``);
* :func:`detect` reads one whole message: German when its words are mostly
  German (the German command words of lang.COMMAND_WORDS and common German
  words, umlauts and ß count), English when mostly English, None when it
  cannot tell (a number, a name, a single unknown word);
* hysteresis: the language switches after SWITCH_AFTER messages in a row
  detected as the other language; one stray word or message never flips it;
* ``language en|de|auto`` / ``sprache en|de|auto`` sets it explicitly; once
  set, automatic switching stops until ``auto``.

The agent's own answers are the model's: it answers in the language the
user writes. Only the gateway's texts (help, notices, approvals, the command
menu) follow this setting.
"""
from __future__ import annotations

import os
import re

#: Messages in a row in the other language before the phone switches.
SWITCH_AFTER = 2
LANGS = ("en", "de")

#: Common German words (function words and the words people write to an
#: agent), and their English counterparts. Matchers for what the user types.
_DE = frozenset("""
und oder aber der die das den dem des ein eine einen einem einer ist sind war
nicht kein keine bitte danke mit für fuer auf aus bei von vom zum zur nach
ich du er sie wir ihr mir mich dir dich uns euch mein dein sein ihr unser
was wie wo wann warum welche welcher welches noch schon jetzt dann auch nur
mal hier dort gut ja nein doch sehr mehr alle alles etwas hallo
schreib schreibe zeig zeige mach mache starte rechne berechne prüfe pruefe
hol holen neu sitzungen sitzung hilfe weg trennen entkoppeln handys geräte
datei dateien meldungen zeitpläne pläne ausschalten einschalten sprache
""".split())
_EN = frozenset("""
and or but the a an is are was not no please thanks thank with for on of from
to at by in i you he she we they me my your our what how where when why which
still already now then also only here there good yes very more all something
hello write show make start run compute check get new sessions session help
off disconnect unpair phones devices file files notify schedules language
""".split())
_WORD = re.compile(r"[a-zA-ZäöüÄÖÜß]+")


def default() -> str:
    lang = os.environ.get("DELFIN_MESSENGER_LANG", "").strip().lower()
    return lang if lang in LANGS else "en"


#: Phone command words (either language): a command says nothing about the
#: language its sender reads ("sessions" twice is no English).
_COMMANDS = frozenset("""
status pending stop update details more mehr show zeige zeigen neu new name forum meldungen
notify dateien datei files file austausch exchange profil pläne plaene zeitpläne schedules plans
plan ok nein sync language sprache lang approvals freigaben sessions sitzungen hol holen get
fetch weg trennen off disconnect entkoppeln abmelden revoke unpair handys geräte geraete phones
devices hilfe help chat
""".split())


def detect(text: str) -> str | None:
    """The language of one whole message, or None when it cannot tell (a
    command, a number, a name, a single unknown word)."""
    words = [w.lower() for w in _WORD.findall(str(text or ""))]
    if not words or words[0] in _COMMANDS:
        return None
    de = sum(1 for w in words if w in _DE or any(ch in w for ch in "äöüß"))
    en = sum(1 for w in words if w in _EN and w not in _DE)
    if de == en:
        return None
    if de > en and de >= max(1, len(words) // 3):
        return "de"
    if en > de and en >= max(1, len(words) // 3):
        return "en"
    return None


def explicit(text: str) -> str | None:
    """``language en|de|auto`` / ``sprache en|de|auto`` as the setting it
    names (``auto`` included), or None."""
    words = str(text or "").strip().lower().lstrip("/").split()
    if len(words) == 2 and words[0] in {"language", "sprache", "lang"} and \
            words[1] in {"en", "de", "auto", "english", "deutsch", "englisch"}:
        return {"english": "en", "englisch": "en", "deutsch": "de"}.get(words[1], words[1])
    return None


def next_state(state: dict, text: str) -> dict:
    """The contact's language state after one message: {lang, mode, streak}."""
    lang = state.get("lang") if state.get("lang") in LANGS else default()
    mode = state.get("mode") if state.get("mode") in {"auto", "fixed"} else "auto"
    streak = int(state.get("streak") or 0)
    if mode == "fixed":
        return {"lang": lang, "mode": mode, "streak": 0}
    seen = detect(text)
    if seen is None or seen == lang:
        return {"lang": lang, "mode": mode, "streak": 0}
    streak += 1
    if streak >= SWITCH_AFTER:
        return {"lang": seen, "mode": mode, "streak": 0}
    return {"lang": lang, "mode": mode, "streak": streak}


# -- the stored state: one entry per contact, in the gateway's state dir --------

def _file():
    from . import store
    return store.state_dir() / "languages.json"


def _load() -> dict:
    import json
    try:
        value = json.loads(_file().read_text())
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def _save(data: dict) -> None:
    import json
    from delfin.agent.state_paths import write_text_atomic
    write_text_atomic(_file(), json.dumps(data, sort_keys=True))


def of_contact(cid: str) -> str:
    """The language contact ``cid`` is answered in."""
    lang = (_load().get(str(cid)) or {}).get("lang")
    return lang if lang in LANGS else default()


def set_contact(cid: str, lang: str, *, mode: str = "fixed") -> None:
    """Set contact ``cid``'s language (``mode`` fixed: no automatic
    switching until ``auto``)."""
    import time
    data = _load()
    data[str(cid)] = {"lang": lang if lang in LANGS else default(), "mode": mode,
                      "streak": 0, "at": time.time()}
    _save(data)


def choose(cid: str, setting: str) -> str:
    """``language en|de|auto``: fix the language, or switch automatic
    choosing back on (keeping the current language). Output: the
    language now in force."""
    if setting == "auto":
        set_contact(cid, of_contact(cid), mode="auto")
    else:
        set_contact(cid, setting)
    return of_contact(cid)


def observe(cid: str, text: str) -> str:
    """One message of contact ``cid``: the state after it (written only
    when it changed), and the language now in force."""
    import time
    data = _load()
    old = data.get(str(cid)) or {}
    new = next_state(old, text)
    if {k: old.get(k) for k in new} != new or time.time() - float(old.get("at") or 0) > 60:
        data[str(cid)] = {**new, "at": time.time()}
        _save(data)
    return new["lang"]


def of_person(person: str) -> str:
    """The language a person's chats get: the language of the person's
    phone that wrote last."""
    from . import store
    data = _load()
    seen = [(float((data.get(cid) or {}).get("at") or 0), cid) for cid, _ in store.phones(person)]
    seen = [s for s in seen if s[0]]
    return of_contact(max(seen)[1]) if seen else default()
