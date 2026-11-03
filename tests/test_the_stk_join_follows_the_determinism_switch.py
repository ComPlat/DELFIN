"""The two stk-conversion wall-clock joins follow the determinism master switch.

The MANTA build paths ship two stk-conversion daemon threads whose join is the
only wall-clock decision on the build path NOT gated on
``DELFIN_DETERMINISTIC``:

  * delfin/manta/single_structure.py:473  ``_stk_t.join(timeout=_EMBED_TIMEOUT)``
  * delfin/manta/embed_strategies.py:121  ``_st.join(timeout=_EMBED_TIMEOUT)``

Every sibling join on this path IS gated (embed_timeout.py:83,152;
single_structure.py:1081; binding_modes.py:390). When the join fires under CPU
load the daemon stk thread (``stk.BuildingBlock(...).to_rdkit_mol()``, which
internally runs RDKit ETKDG) is orphaned and keeps mutating the process's RDKit
state while the main thread embeds with RDKit in parallel. That unsynchronised
overlap is the load-dependent hazard that can tear an embedded conformer --
the reported empty ``tetrahedral trans`` frame (0 atom lines) in 1 of 32
parallel builds -- and it depends on how fast the machine was at that moment.

These tests pin the requirement that both joins are UNBOUNDED (wait for the
thread) under ``DELFIN_DETERMINISTIC=1``, exactly as the embed joins are
pinned in test_the_embed_join_follows_the_determinism_switch.py, so a loaded
machine can never again orphan a stk thread into the build.

The primary test asserts the source (both callsites gate the join on
``_deterministic_mode()``), the repo's established idiom for a determinism
switch that is read directly at the callsite (cf. test_cli_dashboard_parity.py
asserts source text). A behavioral test then builds the metal SMILES whose
frame is reported empty and asserts no emitted frame ever has zero atom lines
under the deterministic default with a forced-tiny timeout.
"""
from __future__ import annotations

import re
from pathlib import Path

_SINGLE = Path(__file__).resolve().parents[1] / "delfin/manta/single_structure.py"
_STRATS = Path(__file__).resolve().parents[1] / "delfin/manta/embed_strategies.py"

# The deterministic-branch expression every wall-clock join on this path must
# use, matching the existing pattern at single_structure.py:1081 and
# embed_timeout.py:83,152.
_DET_GUARD = re.compile(
    r"join\(\s*timeout\s*=\s*None\s+if\s+_deterministic_mode\(\)\s+else\s+(_EMBED_TIMEOUT|timeout)\s*\)"
)


def _stk_join_statements(source: str):
    """Yield the inner expressions of every ``.join(...)`` line in *source*."""
    for line in source.splitlines():
        if "join(" in line and ("_TIMEOUT" in line or "timeout" in line):
            yield line


# --- Primary: both stk joins must be gated on the master switch ------------

def test_the_two_stk_joins_gate_on_the_determinism_switch():
    """single_structure.py:473 and embed_strategies.py:121 must join the stk
    thread with ``timeout=None if _deterministic_mode() else ...`` so that no
    wall-clock decision survives under DELFIN_DETERMINISTIC=1. Before this fix
    both passes a bare ``timeout=_EMBED_TIMEOUT``, so the test is RED."""
    failing = []
    for path, name in ((_SINGLE, "single_structure.py"),
                       (_STRATS, "embed_strategies.py")):
        source = path.read_text()
        stk_joins = [ln for ln in _stk_join_statements(source)
                     if "stk" in ln.lower() or "_stk" in ln or "_st" in ln]
        if not stk_joins:
            # Fall back to any join line that is not already the gated pattern.
            stk_joins = [ln for ln in _stk_join_statements(source)
                         if not _DET_GUARD.search(ln)]
        for ln in stk_joins:
            if not _DET_GUARD.search(ln):
                failing.append(f"{name}: {ln.strip()}")
    assert not failing, (
        "A stk-conversion join on the build path still cuts on wall-clock even "
        "under DELFIN_DETERMINISTIC=1; under CPU load it orphans a daemon stk "
        "thread into the RDKit build (the reported empty frame). Fix these to "
        "join with timeout=None if _deterministic_mode() else ...:\n"
        + "\n".join(failing))


# --- Behavioral: never emit an empty frame under the deterministic default ---

def test_build_of_the_reported_smiles_emits_no_empty_frame(monkeypatch, tmp_path):
    """Building the metal SMILES whose frame was reported empty must never emit
    a 0-atom frame, even with the embed timeout forced to ~0 (the wall-clock
    cutoff a loaded machine hits). Asserts every emitted frame carries atom
    lines and that two builds agree byte-for-byte."""
    from delfin import cli_manta

    import os
    for k in [k for k in os.environ if k.startswith("DELFIN_")]:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("DELFIN_DETERMINISTIC", "1")
    monkeypatch.setenv("DELFIN_EMBED_TIMEOUT", "0")

    manifolds = []
    for tag in ("a", "b"):
        out = tmp_path / f"m_{tag}"
        code = cli_manta.main(["[Pt](Cl)(Cl)(N)N", "-o", str(out), "-q"])
        assert code == 0, code
        manifest = __import__("json").loads(
            (out / "manifest.json").read_text())
        manifolds.append([
            (item["label"], cli_manta._atom_lines((out / item["file"]).read_text()))
            for item in manifest["isomers"]
        ])
    for frames in manifolds:
        assert frames, "no frames emitted"
        empty = [lbl for lbl, lines in frames if not lines]
        assert not empty, f"empty (0-atom) frames emitted: {empty}"
    assert manifolds[0] == manifolds[1], (
        "two builds under DELFIN_DETERMINISTIC=1 with forced-0 embed timeout "
        "disagree — a build decision still depends on wall-clock timing")
