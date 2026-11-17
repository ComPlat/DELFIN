"""delfin-manta never leaves a frame file or the manifest half written.

A reader of the output directory (a second delfin-manta run into the same -o, the
parity side script, a dashboard tailing the directory) must see either the old
file or the new one, never a truncated one.  Writing with Path.write_text opens
the destination with truncation first, so a reader in between sees an empty
frame.  The build is replaced by a fixed manifold here; no structure is built.
"""

from __future__ import annotations

import io
import json
from pathlib import Path

from delfin import cli_manta
from delfin.common import manta_build

_FRAMES = [
    ("Pt 0.0 0.0 0.0\nCl 2.3 0.0 0.0\nCl -2.3 0.0 0.0", "trans"),
    ("Pt 0.0 0.0 0.0\nCl 2.3 0.0 0.0\nCl 0.0 2.3 0.0", "cis"),
]


def _is_whole_output_file(path: Path) -> bool:
    text = path.read_text()
    if path.name == "manifest.json":
        try:
            json.loads(text)
        except ValueError:
            return False
        return True
    lines = text.splitlines()
    return bool(lines) and lines[0].strip().isdigit() and len(lines) == int(lines[0]) + 2


def test_every_destination_file_is_whole_at_every_open(tmp_path, monkeypatch):
    monkeypatch.setattr(manta_build, "build_isomers_isolated",
                        lambda smiles, kwargs, timeout=None, env=None: (list(_FRAMES), None))
    out = tmp_path / "out"
    assert cli_manta.main(["[Pt](Cl)Cl", "-o", str(out), "-q"]) == 0
    names = sorted(p.name for p in out.iterdir())

    torn = []
    real_open = io.open
    looking = []

    def _watching_open(file, mode="r", *a, **k):
        fh = real_open(file, mode, *a, **k)
        if looking:                      # the look below reads files itself
            return fh
        looking.append(1)
        try:
            # What a concurrent reader would find right now in the destination files.
            for p in out.iterdir():
                if p.name in names and not _is_whole_output_file(p):
                    torn.append((p.name, mode))
        finally:
            looking.clear()
        return fh

    monkeypatch.setattr(io, "open", _watching_open)
    assert cli_manta.main(["[Pt](Cl)Cl", "-o", str(out), "-q"]) == 0
    monkeypatch.setattr(io, "open", real_open)

    assert not torn, torn
    assert sorted(p.name for p in out.iterdir()) == names      # no temporary file left behind
    assert all(_is_whole_output_file(out / n) for n in names)
