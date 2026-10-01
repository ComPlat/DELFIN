"""Questions in somebody else's words, so retrieval can be judged on them.

Input: an index and a callable that answers a prompt with one line.
Output: a fixture of (question, the section it came from), saved so it is
built once and re-used.

Why this exists. ``retrieval_eval`` derives its queries from the corpus
verbatim, which measures whether a document can be reached at all -- the
question that caught a chunker indexing a quarter of each paper. It
cannot judge a change that trades exact matching for meaning: measured on
that set, adding character n-grams or LSA to the lexical scorer moves
recall@1 from 0.975 to 0.959 and 0.943. That is the set answering the
question it was built for, not evidence against the method.

A question a researcher actually types shares few words with the passage
that answers it. So the fixture has to contain questions that do NOT
reuse the passage's language, and this module refuses the ones that do --
mechanically, with the same token measure the memory store dedups by,
rather than by trusting the writer of the question to have obeyed.

Model-free by construction. The caller passes ``ask``; this module never
chooses a backend and never reaches a network. That is not tidiness: the
passages are the user's literature, and sending them anywhere is a
decision for whoever runs this, made with the backend in front of them.
"""

from __future__ import annotations

import json
import random
import time
from pathlib import Path
from typing import Callable

from delfin.agent.memory_store import _tokenize

#: Above this share of its OWN words taken from the passage, a question is
#: a quotation. Containment, not Jaccard: Jaccard is symmetric, so a short
#: quotation out of a long passage scores low because the passage's other
#: words fill the union -- the filter then failed exactly where passages
#: are long, which is every real one. Measured on a 500-character passage,
#: a verbatim 12-word fragment scores Jaccard 0.19 and containment 1.00.
#:
#: 0.60 is strict on purpose. A real question shares its topic words with
#: the passage and little else; the whole value of this fixture is that a
#: lexical scorer cannot win it by matching words it was shown.
_MAX_OVERLAP = 0.60

#: A question shorter than this says nothing; longer than this is a summary.
_MIN_WORDS = 5
_MAX_WORDS = 30

_DEFAULT_PATH = Path.home() / ".delfin" / "retrieval_paraphrase.json"


def prompt_for(passage: str, title: str = "") -> str:
    """What to ask a model for one passage."""
    return (
        "Below is a passage from a scientific document.\n\n"
        "Write ONE question that a researcher might type into a search box "
        "and that this passage answers.\n\n"
        "Rules:\n"
        "- Use the researcher's own words, not the passage's. Avoid its "
        "distinctive terms, names and numbers.\n"
        "- Ask about the subject, not about the text ('what basis set ...', "
        "not 'what does this passage say').\n"
        "- One line, no quotes, no preamble.\n\n"
        f"Passage{(' from ' + title) if title else ''}:\n{passage[:1800]}"
    )


def accept(question: str, passage: str) -> tuple[bool, str]:
    """Whether *question* may enter the fixture, and why not when it may not.

    The overlap test is the one that matters. A model asked not to reuse
    the passage's words will often reuse them anyway, and a fixture full
    of quotations measures the same thing the verbatim set already does --
    twice, at the cost of a model call.
    """
    text = " ".join(str(question or "").split())
    if not text:
        return False, "empty"
    words = text.split()
    if len(words) < _MIN_WORDS:
        return False, f"too short ({len(words)} words)"
    if len(words) > _MAX_WORDS:
        return False, f"too long ({len(words)} words)"
    if "\n" in str(question or "").strip():
        return False, "more than one line"
    asked = _tokenize(text)
    if not asked:
        return False, "no content words"
    taken = len(asked & _tokenize(passage)) / len(asked)
    if taken > _MAX_OVERLAP:
        return False, f"quotes the passage (overlap {taken:.2f})"
    return True, ""


def build(index: dict, ask: Callable[[str], str], *,
          per_document: int = 6, seed: int = 0,
          min_chars: int = 400) -> dict:
    """Generate the fixture. Returns it; the caller decides to save it.

    *ask* takes a prompt and returns one line. Every failure -- a refusal,
    an exception, a rejected question -- is recorded and skipped, never
    retried with a softer rule: a fixture that argues with its own filter
    is not a fixture.
    """
    rng = random.Random(seed)
    cases: list[dict] = []
    rejected: list[dict] = []
    for doc_id, doc in (index.get("documents") or {}).items():
        sections = [(sid, s) for sid, s in (doc.get("sections") or {}).items()
                    if len(s.get("text", "")) >= min_chars]
        rng.shuffle(sections)
        for section_id, section in sections[:per_document]:
            passage = section.get("text", "")
            try:
                question = ask(prompt_for(passage, doc.get("title", "")))
            except Exception as exc:
                rejected.append({"doc_id": doc_id, "section_id": section_id,
                                 "reason": f"ask failed: {exc}"})
                continue
            ok, why = accept(question, passage)
            if not ok:
                rejected.append({"doc_id": doc_id, "section_id": section_id,
                                 "reason": why, "question": question})
                continue
            cases.append({"question": " ".join(str(question).split()),
                          "doc_id": doc_id, "section_id": section_id})
    return {
        "version": 1,
        "built_at": time.time(),
        # Pinned to the index it was made from: a fixture whose documents
        # have changed measures something nobody can name.
        "index_built_at": index.get("built_at", ""),
        "document_ids": sorted((index.get("documents") or {})),
        "cases": cases,
        # Kept, not counted. The rejected questions say what the model did
        # with the instruction, which is the only way to tell "the model
        # quotes" from "the passage has no question in it".
        "rejected": rejected,
    }


def save(fixture: dict, path: Path | None = None) -> Path:
    target = Path(path or _DEFAULT_PATH)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(fixture, indent=2, ensure_ascii=False),
                      encoding="utf-8")
    return target


def load(path: Path | None = None) -> dict:
    try:
        return json.loads(Path(path or _DEFAULT_PATH).read_text(
            encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def cases(fixture: dict, index: dict | None = None) -> list[tuple[str, str]]:
    """(question, doc_id) pairs, or [] when the fixture does not fit *index*.

    Refusing a mismatched fixture rather than scoring against it: a set
    built for other documents produces a number, and the number is about
    nothing.
    """
    if index is not None:
        here = sorted((index.get("documents") or {}))
        if fixture.get("document_ids") and fixture["document_ids"] != here:
            return []
    return [(c["question"], c["doc_id"]) for c in fixture.get("cases", [])
            if c.get("question") and c.get("doc_id")]


def summary(fixture: dict) -> dict:
    """What the fixture is, in numbers a reader can check."""
    rejected = fixture.get("rejected") or []
    reasons: dict[str, int] = {}
    for item in rejected:
        key = str(item.get("reason", "")).split(" (")[0]
        reasons[key] = reasons.get(key, 0) + 1
    return {"accepted": len(fixture.get("cases") or []),
            "rejected": len(rejected), "reasons": reasons}
