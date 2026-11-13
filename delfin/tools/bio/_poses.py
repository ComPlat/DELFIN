"""A pose set: one multi-record SDF, every record one placement, scores in SD tags.

Input: the path of an SDF holding one or more records of the same molecule,
each record one placement (pose).  Optional ``score_tag`` names the SD tag
whose value is the pose's score (a docking score, an affinity, an energy);
``lower_is_better`` says which direction ranks first (docking scores are
energies, so the default is ``True``).

Output: a :class:`PoseSet` whose ``poses`` are the records **in file order** —
the order a docking program or depositor wrote is itself information, and a
reordered pose file would silently put a different placement first — each with
its score (``None`` where the tag is absent).  ``ranking`` lists the poses
best-first as indices into ``poses``: lowest score first when
``lower_is_better``, ties in file order, unscored poses after every scored one
(in file order among themselves).

Semantics: a score is read only when it parses as a number.  A tag whose
value is not a number is refused with the record named
(:class:`~delfin.tools.bio.StructureFormatError`) — treating a broken score
as a good pose would put a wrong placement first, and the caller cannot see
the file at the point the number is used.  Without ``score_tag`` no scores
are read and the ranking is the file order.

The records are read through :func:`delfin.tools.bio.parse_sdf`, which
enforces the atom/bond counts and the V2000 layout; this module adds only the
pose semantics on top.

Rejected alternatives: sorting ``poses`` in place (loses the file order the
ranking is expressed against); a NaN or infinite score is a number only to
``float`` and not to a ranking — such a value is refused like any other
non-numeric one.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union

from delfin.tools.bio._errors import StructureFormatError
from delfin.tools.bio._model import SmallMolecule
from delfin.tools.bio._ligand import parse_sdf

_PathLike = Union[str, Path]


@dataclass(frozen=True)
class Pose:
    """One placement of the molecule, with its score when the file gave one."""

    molecule: SmallMolecule
    score: Optional[float] = None
    index: int = 0           # its position in the file (0-based)

    @property
    def name(self) -> str:
        return self.molecule.name


@dataclass(frozen=True)
class PoseSet:
    """Every record of a pose SDF, in file order, plus the score ranking."""

    poses: Tuple[Pose, ...] = ()
    lower_is_better: bool = True

    @property
    def ranking(self) -> Tuple[int, ...]:
        """Indices into :attr:`poses`, best first."""
        scored = [(i, p.score) for i, p in enumerate(self.poses)
                  if p.score is not None]
        scored.sort(key=lambda item: (item[1] if self.lower_is_better
                                      else -item[1], item[0]))
        unscored = [i for i, p in enumerate(self.poses) if p.score is None]
        return tuple(i for i, _ in scored) + tuple(unscored)


def read_pose_set(path: _PathLike, *, score_tag: Optional[str] = None,
                  lower_is_better: bool = True) -> PoseSet:
    """Read a pose SDF.  Raises :class:`StructureFormatError`."""
    p = Path(path)
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise StructureFormatError(
            f"{p.name}: cannot be read ({exc.strerror})"
        ) from None
    molecules = parse_sdf(text, source_name=p.name)

    poses: List[Pose] = []
    for i, molecule in enumerate(molecules):
        score: Optional[float] = None
        if score_tag is not None:
            raw = molecule.properties.get(score_tag)
            if raw is not None:
                try:
                    score = float(raw.strip())
                    if not math.isfinite(score):
                        raise ValueError
                except ValueError:
                    raise StructureFormatError(
                        f"{p.name}: record {i + 1} "
                        f"({molecule.name or 'unnamed'}): the {score_tag} "
                        f"tag is {raw.strip()!r}, not a number"
                    ) from None
        poses.append(Pose(molecule=molecule, score=score, index=i))

    return PoseSet(poses=tuple(poses), lower_is_better=lower_is_better)


__all__ = ["Pose", "PoseSet", "read_pose_set"]
