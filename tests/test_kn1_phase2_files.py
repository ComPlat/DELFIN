"""KN1 red control: the wave-22 distilled lessons must exist as shipped
operating-knowledge entries. RED on a4c9ebaf (files do not exist), GREEN
after phase 2 ships them. Slugs are fixed here so a renamed lesson is a
visible change, not a silent one."""
from pathlib import Path

from delfin.agent import knowledge as kn

OPERATING = kn.global_dir() / "operating"

WAVE22_SLUGS = [
    "count-and-loader-tests-pin-their-own-tasks",
    "a-fresh-folder-for-every-run",
]


def test_wave22_lessons_are_shipped():
    missing = [s for s in WAVE22_SLUGS
               if not (OPERATING / f"{s}.md").is_file()]
    assert missing == [], missing
