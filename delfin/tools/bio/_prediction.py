"""Structure-prediction engines: detected, described, never run or installed.

Input: nothing but a mapping to read ``PATH`` from (``prediction_engines``),
or a FASTA path plus the same mapping (the
``structure_prediction_engine`` block in
:mod:`delfin.tools.adapters.structure_prediction`).  Output: one entry per
engine — whether its command is on PATH, where, and the command line that
would run it on a FASTA.

Semantics, and the line the interface does not cross:

* **Detection only.**  ``shutil.which`` under the caller's ``PATH``, the same
  probe every external-program block in this layer makes.  The engine is
  never executed here, never installed, never fetched: DELFIN steers real
  calculations, and a prediction engine is exactly the kind of program a
  caller must decide to have.  A command is *described*, built to be read —
  the whole of what this interface does.
* **The command is described from the sequence the caller names.**  Each
  engine's invocation is what its own documentation publishes (see the
  table): input FASTA, output directory, the engine's own controls left at
  their defaults.  A described command is a starting point for a caller who
  has read the engine's docs, not a promise.
* **The environment is read where the caller points.**  ``os.environ`` is
  the default; a caller can hand a different mapping (tests do), so the
  interface reads no environment it was not given.

Licences of the programs named here (code licences; trained weights carry
their own terms in every case):

* ColabFold — MIT code; the AlphaFold weights it loads are licensed
  separately by their own terms.
* ESMFold (Meta ESM) — MIT.
* Boltz — MIT, code and weights.
* AlphaFold 3 (DeepMind) — inference code under Creative Commons
  BY-NC-SA 4.0 (non-commercial); the model terms must be accepted before
  use.

Where a prediction job would run on a GPU: DELFIN drives machine-learned
potentials through the ``orca_mlips`` layer (``delfin/mlp_tools/orca_mlips.py``,
``detect_orca_mlips``), which resolves the ``orca-mlips`` environment and its
ProgExt entry points and marks GPU-capable backends for Slurm GPU partitions
(docs/MLIP_GOAT_GPU.md).  A prediction engine would be driven the same way —
detection, an environment resolved through the builder registry, and a job
submitted to a GPU partition — but that path is deliberately left unchanged
here: this interface detects and describes, it does not submit.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import List, Mapping, Optional


class _Engine:
    """One predictor: the command that identifies it and how it is invoked."""

    __slots__ = ("name", "binary", "licence", "argv")

    def __init__(self, name: str, binary: str, licence: str,
                 argv) -> None:
        self.name = name
        self.binary = binary
        self.licence = licence
        self.argv = argv

    def entry(self, env: Mapping[str, str],
              fasta: Optional[str] = None) -> dict:
        """The engine's report entry: detected?, path, described command.

        The command is described whether or not a sequence was named: with
        *fasta* it names the file, without it it carries the ``{fasta}``
        placeholder for a caller to fill.
        """
        path = shutil.which(self.binary, path=env.get("PATH"))
        if fasta is not None:
            command = [part.replace("{fasta}", str(fasta)) for part in self.argv]
        else:
            command = list(self.argv)
        return {
            "name": self.name,
            "command_name": self.binary,
            "detected": path is not None,
            "path": path or "",
            "licence": self.licence,
            "command": command,
        }


#: The invocation each engine's own documentation publishes.  ``{fasta}`` is
#: the sequence file; every control is left at the engine's default.
_ENGINES = (
    _Engine(
        "colabfold", "colabfold_batch",
        "MIT (code); the AlphaFold weights it loads are licensed separately",
        ("colabfold_batch", "{fasta}", "predictions_out"),
    ),
    _Engine(
        "esmfold", "esmfold",
        "MIT (code and weights)",
        ("esmfold", "--input", "{fasta}", "--output", "predictions_out"),
    ),
    _Engine(
        "boltz", "boltz",
        "MIT (code and weights)",
        ("boltz", "predict", "{fasta}", "--outdir", "predictions_out"),
    ),
    _Engine(
        "alphafold3", "alphafold3",
        "CC BY-NC-SA 4.0 (code, non-commercial); model terms accepted before use",
        ("alphafold3", "--fasta_path={fasta}", "--output_dir=predictions_out"),
    ),
)

#: The output directory a described command names.  The caller replaces it
#: with their own when they run the engine.
DEFAULT_OUTPUT_DIR = "predictions_out"


def prediction_engines(env: Optional[Mapping[str, str]] = None,
                       fasta: Optional[str] = None) -> List[dict]:
    """One entry per engine: detected under *env*'s PATH, and its command.

    With *fasta*, the described command names the sequence file; without it,
    the command is empty and only detection is reported.  No engine is run,
    installed or downloaded by this call.
    """
    mapping = dict(os_environ()) if env is None else env
    return [engine.entry(mapping, fasta) for engine in _ENGINES]


def os_environ() -> Mapping[str, str]:
    """The real environment, named for callers who want the default."""
    import os
    return os.environ


def fasta_path_for(record, work_dir: Path) -> Path:
    """Write the first sequence record to a FASTA an engine can read."""
    path = work_dir / "sequence.fasta"
    path.write_text(f">{record.id}\n{record.sequence}\n", encoding="utf-8")
    return path


__all__ = [
    "DEFAULT_OUTPUT_DIR",
    "prediction_engines",
    "fasta_path_for",
]
