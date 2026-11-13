"""A pose set is one multi-record SDF where every record is one placement of
the same molecule, and its SD tags carry the scores.

The records come back in file order — a pose file's record order is meaningful
(depositors and docking programs rank by it).  With ``score_tag`` given, the
value under that SD tag is parsed as a number; a tag whose value is not a
number is refused with the record named, because silently treating a
missing/broken score as a good pose would put a wrong placement first.
Without ``score_tag`` every pose simply has no score and the ranking is the
file order.  The ranking lists the poses best-first; with
``lower_is_better=True`` (docking scores are energies) the lowest score leads.
"""

import pytest

from delfin.tools.bio import Pose, PoseSet, read_pose_set
from delfin.tools.bio._errors import StructureFormatError


def _record(name, score=None, extra_tag="CLIPBOARD"):
    lines = [
        name,
        "  pose",
        "",
        "  2  1  0  0  0  0  0  0  0  0999 V2000",
        "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0",
        "    1.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
        "  1  2  1  0  0  0  0",
        "M  END",
    ]
    if score is not None:
        lines.append(f">  <{extra_tag}>")
        lines.append(str(score))
        lines.append("")
    lines.append("$$$$")
    return lines


def _write(tmp_path, name, *records):
    path = tmp_path / name
    path.write_text("\n".join("\n".join(r) for r in records) + "\n",
                    encoding="utf-8")
    return path


def test_the_records_come_back_in_file_order(tmp_path):
    path = _write(tmp_path, "poses.sdf", _record("b"), _record("a"),
                  _record("c"))
    ps = read_pose_set(path)
    assert isinstance(ps, PoseSet)
    assert [p.molecule.name for p in ps.poses] == ["b", "a", "c"]


def test_a_score_tag_gives_every_pose_a_number(tmp_path):
    path = _write(
        tmp_path, "poses.sdf",
        _record("one", score="-7.5", extra_tag="SCORE"),
        _record("two", score="-9.1", extra_tag="SCORE"),
    )
    ps = read_pose_set(path, score_tag="SCORE")
    assert [p.score for p in ps.poses] == [-7.5, -9.1]
    assert all(isinstance(p, Pose) for p in ps.poses)


def test_with_lower_is_better_the_lowest_score_ranks_first(tmp_path):
    path = _write(
        tmp_path, "poses.sdf",
        _record("one", score="-7.5", extra_tag="SCORE"),
        _record("two", score="-9.1", extra_tag="SCORE"),
        _record("three", score="-8.0", extra_tag="SCORE"),
    )
    ps = read_pose_set(path, score_tag="SCORE")
    assert ps.ranking == (1, 2, 0)


def test_higher_is_better_reverses_the_ranking(tmp_path):
    path = _write(
        tmp_path, "poses.sdf",
        _record("one", score="3.0", extra_tag="AFFINITY"),
        _record("two", score="7.0", extra_tag="AFFINITY"),
    )
    ps = read_pose_set(path, score_tag="AFFINITY", lower_is_better=False)
    assert ps.ranking == (1, 0)


def test_a_non_numeric_score_is_refused_and_names_the_record(tmp_path):
    path = _write(
        tmp_path, "poses.sdf",
        _record("one", score="-7.5", extra_tag="SCORE"),
        _record("two", score="n/a", extra_tag="SCORE"),
    )
    with pytest.raises(StructureFormatError, match="two"):
        read_pose_set(path, score_tag="SCORE")


def test_a_pose_without_the_score_tag_has_no_score_and_ranks_last(tmp_path):
    path = _write(
        tmp_path, "poses.sdf",
        _record("scored", score="-7.5", extra_tag="SCORE"),
        _record("unscored", score=None, extra_tag="SCORE"),
    )
    ps = read_pose_set(path, score_tag="SCORE")
    assert ps.poses[1].score is None
    assert ps.ranking == (0, 1)


def test_without_a_score_tag_no_scores_are_read_and_the_ranking_is_file_order(tmp_path):
    path = _write(tmp_path, "poses.sdf", _record("b"), _record("a"))
    ps = read_pose_set(path)
    assert [p.score for p in ps.poses] == [None, None]
    assert ps.ranking == (0, 1)


def test_ties_rank_in_file_order(tmp_path):
    path = _write(
        tmp_path, "poses.sdf",
        _record("first", score="-7.0", extra_tag="SCORE"),
        _record("second", score="-7.0", extra_tag="SCORE"),
    )
    ps = read_pose_set(path, score_tag="SCORE")
    assert ps.ranking == (0, 1)


def test_a_file_that_is_not_an_sdf_at_all_is_refused(tmp_path):
    path = tmp_path / "poses.sdf"
    path.write_text("nothing here\n", encoding="utf-8")
    with pytest.raises(StructureFormatError):
        read_pose_set(path)
