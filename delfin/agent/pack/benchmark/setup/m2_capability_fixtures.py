#!/usr/bin/env python3
"""Build the workspace for the m2 capability-measurement tasks.

Input: the workspace path. Output: one subdirectory per task, each holding
the fixture that task's question is asked about.

Everything here is written with the standard library only -- a ``.docx``
is a zip of XML parts, not an output of ``python-docx`` -- because the
benchmark grades a MODEL's answer, and a fixture that needs a package the
grading host may not have would turn "the model could not answer" into
"the fixture could not be built". The runner treats a non-zero exit from
a setup script as UNMEASURED rather than as a model failure, and that is
the correct stance: a question never put is not a wrong answer.

The content is chosen so that the facts a correct answer must state are
NOT recoverable from a plausible-sounding guess, and so that each fixture
holds the near-miss a shallow answer reaches for:

  g1   one definition, one real call, a prose mention, and a same-named
       local closure that must NOT be reported as a call site
  b2   three CSV runs with unequal row counts -> the pooled mean matches
       none of the per-run means, so averaging everything is detectable
  a3   two in-service spec sheets and one SPARE whose flow rate is the
       highest in the folder -> reporting the spare is detectable
  d4   a real .docx (its text is compressed, so a plain read cannot
       recover it) beside a superseded .txt that states a DIFFERENT limit
  sv1  three ethanol calculations at two levels of theory plus a much
       lower METHANOL folder -> any cross-method ranking is detectable

Exit non-zero on any failure.
"""
from __future__ import annotations

import json
import sys
import zipfile
from pathlib import Path

# --------------------------------------------------------------------------
# A5 read_document -- the .docx fixture (binary to a plain read)
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# A3 read_file -- one answer, split across three files, with a decoy
# --------------------------------------------------------------------------

_PUMP = {
    "specs/pump_a.txt": (
        "Specification sheet -- circulation pump A\n"
        "model: PA-200\n"
        "maximum flow rate: 18.4 L/min\n"
        "maximum head: 6.2 m\n"
    ),
    "specs/pump_b.txt": (
        "Specification sheet -- circulation pump B\n"
        "model: PB-310\n"
        "maximum flow rate: 11.9 L/min\n"
        "maximum head: 9.8 m\n"
    ),
    # The decoy: a pump that is NOT one of the two quoted, with a flow
    # rate that would win if the question were "which pump flows most".
    "specs/pump_spare.txt": (
        "Specification sheet -- spare unit, not in service\n"
        "model: PX-900\n"
        "maximum flow rate: 24.0 L/min\n"
        "maximum head: 4.1 m\n"
    ),
}


def _build_a3(root: Path) -> None:
    d = root / "a3"
    for rel, body in _PUMP.items():
        p = d / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(body, encoding="utf-8")


# --------------------------------------------------------------------------
# A1 grep_file -- a symbol defined once, referenced twice, shadowed once
# --------------------------------------------------------------------------

_GREP = {
    "src/core.py": (
        "# Not implemented, never executed: a CAD placeholder.\n"
        "def _settle(grid, dt):\n"
        '    """Propagate one time step of dt over grid."""\n'
        "    return grid\n"
        "\n"
        "\n"
        "def step(grid, n):\n"
        "    for _ in range(n):\n"
        "        grid = _settle(grid, 0.01)\n"
        "    return grid\n"
    ),
    "src/run.py": (
        "from core import step, _settle\n"
        "\n"
        "\n"
        "def main():\n"
        "    grid = [0.0] * 8\n"
        "    grid = step(grid, 10)\n"
        "    # Final tidy-up shares the name but is a LOCAL closure, not the\n"
        "    # helper in core.py -- reporting it would be a wrong answer.\n"
        "    def _settle():\n"
        "        return grid\n"
        "    return _settle()\n"
    ),
    "notes/README.md": (
        "The solver lives in src/. `_settle` is internal and may change.\n"
    ),
}


def _build_g1(root: Path) -> None:
    d = root / "g1"
    for rel, body in _GREP.items():
        p = d / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(body, encoding="utf-8")


# --------------------------------------------------------------------------
# A2 bash -- numbers in files that must be computed, not guessed
# --------------------------------------------------------------------------

_BASH = {
    "run_a/cycle_times.csv": (
        "run,seconds\n"
        "a1,120\na2,98\na3,143\na4,131\na5,109\n"
    ),
    "run_b/cycle_times.csv": (
        "run,seconds\n"
        "b1,87\nb2,92\nb3,101\nb4,79\nb5,95\n"
    ),
    # A third run whose header names a DIFFERENT column; a blind
    # concatenation of every seconds-like column is wrong here.
    "run_c/cycle_times.csv": (
        "run,seconds,warmup_seconds\n"
        "c1,150,20\nc2,141,18\nc3,162,25\n"
    ),
}


def _build_b2(root: Path) -> None:
    d = root / "b2"
    for rel, body in _BASH.items():
        p = d / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(body, encoding="utf-8")


# --------------------------------------------------------------------------
# A4 read_document, second shape -- a real .docx (zip of XML parts)
# --------------------------------------------------------------------------

_W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"

_DOCX_PARAS = [
    "Method note MN-14: trace-metal digestion",
    "Digest 0.25 g of sample in 5 mL of HNO3 at 95 C for 2 h.",
    "CAUTION: do not exceed 120 C -- selenium is lost above it.",
    "Cool, dilute to 25 mL, filter through 0.45 um.",
    "Report results as mg/kg dry mass, not as mg/L.",
]

_CONTENT_TYPES = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
    '<Default Extension="xml" ContentType="application/xml"/>'
    '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
    "</Types>"
)

_RELS = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>'
    "</Relationships>"
)


def _document_xml() -> str:
    paras = "".join(
        f"<w:p><w:r><w:t>{p}</w:t></w:r></w:p>" for p in _DOCX_PARAS)
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        f'<w:document xmlns:w="{_W}"><w:body>{paras}</w:body></w:document>')


def _build_d4(root: Path) -> None:
    d = root / "d4"
    d.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(d / "method_note.docx", "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", _CONTENT_TYPES)
        z.writestr("_rels/.rels", _RELS)
        z.writestr("word/document.xml", _document_xml())
    # A plain-text distractor that states a DIFFERENT temperature limit.
    # An answer that reads only the .txt gives 150 C and is wrong.
    (d / "method_note_old.txt").write_text(
        "Method note MN-14 (superseded 2019)\n"
        "Digest at 95 C. Do not exceed 150 C.\n", encoding="utf-8")


# --------------------------------------------------------------------------
# A6 search_calcs -- two methods, a decoy folder, and a truth file
# --------------------------------------------------------------------------

_SV_ENERGIES = {
    "ethanol_pbe0": ("PBE0", "def2-SVP", -154.981234),
    "ethanol_b3lyp": ("B3LYP", "def2-SVP", -154.907142),
    "ethanol_pbe0_tzvp": ("PBE0", "def2-TZVP", -155.041188),
}


def _build_sv1(root: Path) -> None:
    base = root / "sv1" / "calc"
    for name, (func, basis, energy) in _SV_ENERGIES.items():
        d = base / name
        d.mkdir(parents=True, exist_ok=True)
        (d / "CONTROL.txt").write_text(
            f"functional: {func}\nbasis_set: {basis}\n"
            "charge: 0\nmultiplicity: 1\n", encoding="utf-8")
        out = d / f"{name}.out"
        lines = [
            f"! {func} {basis} Opt Freq".replace("def2-SVP", "def2-SVP"),
            "",
            "FINAL SINGLE POINT ENERGY      "
            f"{energy:.6f}",
            "",
            "VIBRATIONAL FREQUENCIES",
            "",
            "   0:       0.00 cm**-1",
            "   1:       0.00 cm**-1",
            "   2:       0.00 cm**-1",
        ]
        for i in range(3, 21):
            lines.append(f"  {i}:     {1050 + 37 * i:.2f} cm**-1")
        lines += ["", "ORCA TERMINATED NORMALLY", ""]
        out.write_text("\n".join(lines), encoding="utf-8")
    # The decoy: a folder that looks like a run of the same molecule but
    # is a DIFFERENT species, so ranking it in would be wrong.
    decoy = base / "methanol_pbe0"
    decoy.mkdir(parents=True, exist_ok=True)
    (decoy / "CONTROL.txt").write_text(
        "functional: PBE0\nbasis_set: def2-SVP\n", encoding="utf-8")
    (decoy / "methanol_pbe0.out").write_text(
        "! PBE0 def2-SVP\n\nFINAL SINGLE POINT ENERGY      -115.501000\n"
        "\nORCA TERMINATED NORMALLY\n", encoding="utf-8")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    root = Path(argv[1])
    try:
        _build_b2(root)
        _build_g1(root)
        _build_a3(root)
        _build_d4(root)
        _build_sv1(root)
    except Exception as exc:                              # pragma: no cover
        print(f"setup failed: {exc}", file=sys.stderr)
        return 1
    # A manifest so an acceptance script can tell "setup ran" from "the
    # fixture was never built", instead of inferring it from a missing
    # directory and reporting a model failure.
    (root / "m2_fixture_manifest.json").write_text(
        json.dumps({"built": ["A1", "A2", "A3", "A4", "A5", "A6"]},
                   indent=1), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
