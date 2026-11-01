"""A streaming index of the sections ORCA prints in its output file.

This complements DELFIN's typed parsers (``extract_*`` / ``parse_orca_output``);
it replaces none. It answers a different question: *where in the output file is
a given section?* — the section banners ORCA emits, as ``(title, start, end)``
ranges — so an agent can page straight to, say, the ``FINAL SINGLE POINT
ENERGY`` block of a 100 MB output without any parser having to understand the
lines between.

The index is built from the banner patterns ORCA prints — two kinds, both
derived from the four ``tests/fixtures/orca611`` outputs and the manual's
descriptions (``quickstartguide/hellowater.txt``, ``structurereactivity/
optimizations.txt``):

* framed ``* Title *`` banner blocks: a content line ``<spaces>* Title *``
  between two asterisk-only lines, e.g. ``* Geometry Optimization Run *`` and
  ``* GEOMETRY OPTIMIZATION CYCLE   1 *``;
* plain heading lines notched into the text, e.g. ``ORBITAL ENERGIES``,
  ``VIBRATIONAL FREQUENCIES``, ``FINAL SINGLE POINT ENERGY  -76.32 ...``,
  ``TD-DFT/TDA EXCITED STATES (SINGLETS)``;
* the SCF iteration table header ``----D-I-I-S----``;
* end-of-run markers such as ``****ORCA TERMINATED NORMALLY****`` and the
  ``* ERROR *`` banner of an aborted run.

``orca_out_index`` streams the file line by line (it never loads the whole
file, so a 100 MB output is cheap, and ``orca_out_section`` likewise reads only
the line range it returns). Both keep to the marker catalog below; sections
outside it are not reported.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# --- Marker catalog, derived from tests/fixtures/orca611 + the manual ---------

# A framed banner content line: optional spaces, asterisk, the title, asterisk.
_BANNER_LINE = re.compile(r"^\s*\*\s*(?P<title>.+?)\s*\*\s*$")


def _is_aster_only(text: str) -> bool:
    s = text.strip()
    return len(s) > 3 and set(s) == {"*"}


# Banner titles (the ``* Title *`` inner text, upper-cased for matching) that
# open a section. ``GEOMETRY OPTIMIZATION CYCLE n`` is matched by prefix.
_BANNER_TITLES = {
    "SINGLE POINT CALCULATION",
    "GEOMETRY OPTIMIZATION RUN",
    "REAL PERTURBATIONS",
    "MULLIKEN POPULATION ANALYSIS",
    "LOEWDIN POPULATION ANALYSIS",
    "MAYER POPULATION ANALYSIS",
    "ERROR",
}
_BANNER_PREFIX = ("GEOMETRY OPTIMIZATION CYCLE",)

# Plain headings matched on their leading line stub.
_HEADING_TITLES = (
    "FINAL SINGLE POINT ENERGY",
    "ORBITAL ENERGIES",
    "MULLIKEN ATOMIC CHARGES",
    "TOTAL SCF ENERGY",
    "VIBRATIONAL FREQUENCIES",
    "NORMAL MODES",
    "IR SPECTRUM",
    "TD-DFT/TDA EXCITED STATES",
    "CIS/TD-DFT TOTAL ENERGY",
)

# End-of-run markers: title -> substring that closes the run.
_END_MARKERS = {
    "ORCA TERMINATED NORMALLY": "****ORCA TERMINATED NORMALLY****",
}

# SCF iteration table header (ORCA 6, ``----D-I-I-S----``...).
_SCF_TABLE = re.compile(r"-{10,}D-I-I-S-{10,}")


@dataclass
class Section:
    """One indexed region of an ORCA output file.

    ``start_line``/``end_line`` are 1-based and inclusive; ``end_line`` is the
    last line of the section, i.e. the line before the next section's start
    (or the last line of the file for the final section).
    """

    title: str
    start_line: int
    end_line: int
    note: str = ""


def _scan_markers(path) -> list[tuple[int, str]]:
    """Return ``(line_no, title)`` markers in file order, streaming.

    Only lines that open a recognized section are reported; the line number is
    the banner/heading line itself (1-based), not its asterisk frame.
    """
    found: list[tuple[int, str]] = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh, 1):
            s = line.rstrip("\n")
            if _is_aster_only(s):
                continue
            m = _BANNER_LINE.match(s)
            if m:
                title = m.group("title").strip()
                up = " ".join(title.upper().split())
                if up in _BANNER_TITLES or up.startswith(_BANNER_PREFIX):
                    found.append((lineno, up))
                    continue
            for title, marker in _END_MARKERS.items():
                if marker in s:
                    found.append((lineno, title))
                    break
            else:
                for head in _HEADING_TITLES:
                    if s.lstrip().startswith(head):
                        found.append((lineno, head))
                        break
                else:
                    if _SCF_TABLE.match(s.strip()):
                        found.append((lineno, "SCF ITERATIONS"))
    return found


def orca_out_index(path) -> list[Section]:
    """Index an ORCA output file: ordered sections with title and line range.

    The result is entirely marker-driven: each recognized banner or heading
    opens a section that runs until the next section's marker (or end of
    file). Sections are returned in file order. Works on outputs of any size:
    the file is streamed, never loaded whole.
    """
    markers = _scan_markers(path)
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        total = sum(1 for _ in fh)
    sections: list[Section] = []
    for i, (line_no, title) in enumerate(markers):
        end_line = (markers[i + 1][0] - 1) if i + 1 < len(markers) else total
        note = _note_for(title)
        sections.append(Section(title=title, start_line=line_no,
                                end_line=end_line, note=note))
    return sections


def _note_for(title: str) -> str:
    up = title.upper()
    if up.startswith("GEOMETRY OPTIMIZATION CYCLE"):
        return "optimization cycle block"
    if up == "SCF ITERATIONS":
        return "SCF iteration table"
    if up == "ERROR":
        return "error/abort block"
    if up == "FINAL SINGLE POINT ENERGY":
        return "final electronic energy"
    if up == "TD-DFT/TDA EXCITED STATES":
        return "TD-DFT excited states"
    if up == "ORCA TERMINATED NORMALLY":
        return "run ended normally"
    if "POPULATION ANALYSIS" in up or up.endswith("ATOMIC CHARGES"):
        return "population analysis"
    return ""


def orca_out_section(path, title: str, n: int | None = None,
                     max_lines: int = 200) -> list[str]:
    """Return the lines of a section in an ORCA output, as a string list.

    ``title`` names the section as returned by :func:`orca_out_index`
    (case-insensitive on the canonical form). ``n`` selects the *n*-th
    occurrence (1-based); ``None`` means the first. At most ``max_lines``
    lines are returned, starting at the section's marker line. Raises
    ``KeyError`` when no section matches (or ``IndexError`` for an ``n``
    beyond the occurrences). Only the returned range is read.
    """
    want = title.strip().upper()
    sections = [s for s in orca_out_index(path) if s.title.upper() == want]
    if not sections:
        raise KeyError(f"no ORCA section named {title!r} in {path}")
    if n is None:
        n = 1
    if n < 1 or n > len(sections):
        raise IndexError(f"{title!r} has {len(sections)} occurrence(s); "
                         f"requested n={n}")
    sec = sections[n - 1]
    start = sec.start_line
    end = sec.end_line
    lines: list[str] = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh, 1):
            if lineno < start:
                continue
            if lineno > end or len(lines) >= max_lines:
                break
            lines.append(line.rstrip("\n"))
    return lines
