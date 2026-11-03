"""A review set shared out to a team: disjoint shares, a shared overlap, the
tab showing only one's own share, the agreement, unpacking and sending back.

Synthetic data only: invented file names, a water molecule moved a little.
"""
from __future__ import annotations

import csv
import gzip
import json
import subprocess
from pathlib import Path

import pytest

from delfin import review as rv
from delfin import review_team as rt


def _water(stem: str, k: int = 0) -> str:
    return (f"3\n{stem}\n"
            f"O 0.000 0.000 {0.117 + 0.05 * k:.3f}\n"  # far enough apart: never identical
            "H 0.000 0.757 -0.467\nH 0.000 -0.757 -0.467\n")


def _set(folder: Path, n: int) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    for k in range(1, n + 1):
        (folder / f"T{k:04d}.xyz").write_text(_water(f"T{k:04d}", k))
    return folder


@pytest.fixture
def team_dir(tmp_path: Path) -> Path:
    return _set(tmp_path / "set" / "structures", 40)


# ------------------------------------------------------------------ assignment
def test_shares_are_disjoint_the_overlap_shared_and_the_draw_repeatable(team_dir):
    a, _ = rt.assign(team_dir, ["Ann", "Ben", "Cem"], overlap=4, seed=11, redundancy=0)
    b, _ = rt.assign(team_dir, ["Ann", "Ben", "Cem"], overlap=4, seed=11, redundancy=0)
    assert {k: v for k, v in a.items() if k not in ("created", "updated")} == \
        {k: v for k, v in b.items() if k not in ("created", "updated")}
    assert a["files_total"] == 40 and len(a["shared"]) == 4 and a["seed"] == 11
    shares = [set(v) for v in a["raters"].values()]
    assert all(len(s) == 12 for s in shares)                       # (40 - 4) // 3
    assert not (shares[0] & shares[1]) and not (shares[1] & shares[2])
    assert not set(a["shared"]) & set().union(*shares)
    for v in a["raters"].values():
        assert v == sorted(v)                                       # in file order


def test_overlap_zero_and_a_block_leave_the_rest_for_later(team_dir):
    a, _ = rt.assign(team_dir, ["Ann"], overlap=0, block=10, seed=3)
    assert a["shared"] == [] and len(a["raters"]["Ann"]) == 10


def test_a_reviewer_added_later_changes_nobody_elses_share(team_dir):
    a, _ = rt.assign(team_dir, ["Ann"], overlap=4, block=10, seed=5, redundancy=0)
    before = json.loads(json.dumps(a))
    b, msgs = rt.assign(team_dir, ["Ann", "Ben"], block=10, existing=a)
    assert b["raters"]["Ann"] == before["raters"]["Ann"] and b["shared"] == before["shared"]
    assert len(b["raters"]["Ben"]) == 10
    assert not set(b["raters"]["Ben"]) & (set(b["raters"]["Ann"]) | set(b["shared"]))
    assert any("Ann: already in the assignment" in m for m in msgs)
    with pytest.raises(ValueError):
        rt.assign(team_dir, ["Cem"], overlap=7, existing=b)        # overlap is fixed


def test_what_a_reviewer_already_rated_stays_theirs(team_dir):
    s = rv.ReviewSession.open(team_dir, name="Ben")
    s.go_to_file("T0007.xyz")
    s.rate("pass")
    s.go_to_file("T0031.xyz")
    s.rate("block", ["X"])
    a, _ = rt.assign(team_dir, ["Ann", "Ben"], overlap=0, block=5, seed=2, redundancy=0)
    assert {"T0007.xyz", "T0031.xyz"} <= set(a["raters"]["Ben"])
    assert len(a["raters"]["Ben"]) == 5
    assert not {"T0007.xyz", "T0031.xyz"} & set(a["raters"]["Ann"])


def test_the_test_split_is_fixed_and_drawn_apart_from_the_shares(team_dir):
    a, _ = rt.assign(team_dir, ["Ann", "Ben"], overlap=4, seed=8, test_fraction=0.2)
    b, _ = rt.assign(team_dir, ["Ann"], overlap=9, seed=99, test_fraction=0.2)
    assert len(a["split"]["test"]) == 8
    c, _ = rt.assign(team_dir, ["Ann"], overlap=4, seed=8, test_fraction=0.2)
    assert c["split"] == a["split"] and b["split"] != a["split"]


def test_the_order_mixes_the_shared_ones_early_among_the_own(team_dir):
    a, _ = rt.assign(team_dir, ["Ann", "Ben"], overlap=6, seed=4)
    order = rt.rater_order(a, "ann")                               # case is forgiven
    assert sorted(order) == sorted(a["shared"] + a["raters"]["Ann"])
    head = order[:12]
    assert set(a["shared"]) <= set(head)
    assert head[:6] != sorted(a["shared"]) and order[12:] == a["raters"]["Ann"][6:]
    with pytest.raises(KeyError):
        rt.rater_order(a, "Zoe")


def test_two_names_writing_one_review_file_are_refused(team_dir):
    with pytest.raises(ValueError):
        rt.assign(team_dir, ["Ann Lee", "Ann_Lee"], seed=1)


# ------------------------------------------------------------------ the tab
def _panel(where: Path, name: str):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard.tab_review import ReviewPanel
    panel = ReviewPanel(where)
    panel.name_input.value = name
    return panel, panel.open(str(where))


def test_the_tab_shows_only_the_reviewers_share(team_dir):
    a, _ = rt.assign(team_dir, ["Ann", "Ben"], overlap=4, block=8, seed=6, redundancy=0)
    rt.write_assignment(team_dir.parent / rt.ASSIGNMENT_NAME, a)  # next to kategorien.txt
    panel, s = _panel(team_dir, "ann")
    assert s is not None, panel.status_html.value
    assert s.review_path.name == "review_Ann.json"
    assert [f.file for f in s.frames] == rt.rater_order(a, "Ann")
    assert "0 / 12 rated" in panel.progress_html.value
    assert "Your share: 12 structures of 40" in panel.status_html.value
    names = [v for _, v in panel.file_list.options]
    assert names == rt.rater_order(a, "Ann")
    assert panel.go_to(a["raters"]["Ben"][0][:5]) is False
    assert panel.go_to(a["raters"]["Ann"][-1]) is True
    for _ in range(12):
        panel.rate("pass")
    assert "12 / 12 rated" in panel.progress_html.value
    assert set(rv.load_review(s.review_path)["reviews"]) == set(a["shared"] + a["raters"]["Ann"])


def test_a_new_name_draws_its_own_block_and_keeps_it(team_dir, monkeypatch):
    # No name and no GitHub account to take it from: the tab asks for one.
    monkeypatch.setattr(rt, "_gh_login", lambda run: "")
    a, _ = rt.assign(team_dir, ["Ann"], overlap=2, block=3, seed=6, self_block=5,
                     redundancy=0)
    rt.write_assignment(team_dir / rt.ASSIGNMENT_NAME, a)
    panel, s = _panel(team_dir, "")
    assert s is None and "enter your name first" in panel.status_html.value
    panel, s = _panel(team_dir, "Zoe")
    assert s is not None, panel.status_html.value
    assert "Your share: 7 structures" in panel.status_html.value
    assert "drawn for you" in panel.status_html.value
    stored = rv.load_review(team_dir / "review_Zoe.json")["assignment"]
    assert stored["self"] and len(stored["files"]) == 5 and stored["name"] == "Zoe"
    mine = set(stored["files"])
    assert not mine & (set(a["shared"]) | set(a["raters"]["Ann"]))
    assert {f.file for f in s.frames} == mine | set(a["shared"])
    assert "0 / 7 rated" in panel.progress_html.value
    # The same name draws the same block, also on another machine ...
    other = _set(team_dir.parent / "copy", 40)
    assert rt.self_share(other, a, "Zoe")["files"] == stored["files"]
    # ... and keeps it although the assignment changes later.
    b, _ = rt.assign(team_dir, ["Cem"], block=30, existing=a)
    rt.write_assignment(team_dir / rt.ASSIGNMENT_NAME, b)
    panel, s = _panel(team_dir, "Zoe")
    assert set(rv.load_review(team_dir / "review_Zoe.json")["assignment"]["files"]) == mine
    assert {f.file for f in s.frames} == mine | set(a["shared"])
    # The coordinator does not give a self-drawn block out again where it can be seen.
    assert not mine & set(b["raters"]["Cem"])


def test_without_an_assignment_the_tab_shows_the_whole_folder(team_dir):
    panel, s = _panel(team_dir, "Ann")
    assert s is not None and len(s) == 40 and s.assigned is None
    assert "Your share" not in panel.status_html.value


# ------------------------------------------------------------------ agreement
def test_cohen_and_fleiss_kappa_on_known_tables():
    a = [True] * 20 + [False] * 5 + [True] * 10 + [False] * 15
    b = [True] * 20 + [True] * 5 + [False] * 10 + [False] * 15
    assert rt.cohen_kappa(a, b) == pytest.approx(0.4)              # po .7, pe .5
    assert rt.cohen_kappa([True, True], [True, True]) is None
    # three raters, items: all yes, all no, 2 of 3 yes
    k = rt.fleiss_kappa([3, 0, 2], 3)
    p = 5 / 9
    pbar = (1 + 1 + 1 / 3) / 3
    assert k == pytest.approx((pbar - (p * p + (1 - p) ** 2)) / (1 - (p * p + (1 - p) ** 2)))


def _rate(folder, name, plan):
    s = rv.ReviewSession.open(folder, name=name, categories=["H", "TOPO"])
    for fname, verdict, cats in plan:
        s.go_to_file(fname)
        s.rate(verdict, cats)


def test_merge_writes_every_rating_and_measures_the_shared_ones(team_dir, tmp_path):
    a, _ = rt.assign(team_dir, ["Ann", "Ben", "Cem"], overlap=3, block=2, seed=1,
                     test_fraction=0.25)
    rt.write_assignment(team_dir / rt.ASSIGNMENT_NAME, a)
    sh = a["shared"]
    _rate(team_dir, "Ann", [(sh[0], "block", ["H"]), (sh[1], "pass", []), (sh[2], "block", ["TOPO"]),
                            (a["raters"]["Ann"][0], "pass", [])])
    _rate(team_dir, "Ben", [(sh[0], "block", ["H"]), (sh[1], "pass", []), (sh[2], "pass", [])])
    two = rt.agreement({n: rv.load_review(team_dir / f"review_{n}.json") for n in ("Ann", "Ben")}, sh)
    assert two["method"] == "Cohen" and two["frames"] == 3
    assert two["verdict"]["agreement"] == pytest.approx(2 / 3)
    _rate(team_dir, "Cem", [(sh[0], "block", ["H", "TOPO"]), (sh[1], "pass", []), (sh[2], "block", [])])
    out = tmp_path / "merged.csv"
    assert rv.main(["merge", str(team_dir), "--csv", str(out)]) == 0
    rows = list(csv.DictReader(out.open()))
    assert len(rows) == 4 + 3 + 3 and {r["rater"] for r in rows} == {"Ann", "Ben", "Cem"}
    assert sum(r["shared"] == "yes" for r in rows) == 9
    test = set(a["split"]["test"])
    assert all(r["split"] == ("test" if r["file"] in test else "calibration") for r in rows)
    res = rt.agreement(rt.collect_reviews(team_dir), sh)
    assert res["method"] == "Fleiss" and res["frames"] == 3 and res["raters"] == ["Ann", "Ben", "Cem"]
    assert res["categories"]["H"]["agreement"] == pytest.approx(1.0)
    assert res["categories"]["H"]["n"] == 3
    text = rt.format_agreement(res)
    assert "ok/defect" in text and "TOPO" in text
    calib = tmp_path / "calib.csv"
    rt.write_merged_csv(rt.collect_reviews(team_dir), calib, a, split="calibration")
    assert all(r["split"] == "calibration" for r in csv.DictReader(calib.open()))


# ------------------------------------------------------------------ unpack and submit
def _published(tmp_path: Path) -> Path:
    set_dir = tmp_path / "repo" / "review_sets" / "demo"
    set_dir.mkdir(parents=True)
    with gzip.open(set_dir / rt.SET_ARCHIVE_NAME, "wt", encoding="utf-8") as fh:
        for k in range(1, 9):
            fh.write(_water(f"T{k:04d}", k))
    (set_dir / "kategorien.txt").write_text("H\nTOPO\n")
    return set_dir


def test_unpack_writes_single_files_and_never_alters_a_review_in_progress(tmp_path):
    set_dir = _published(tmp_path)
    dest = tmp_path / "mine"
    assert rt.unpack(set_dir, dest) == (9, 0)
    assert (dest / "T0003.xyz").read_text() == _water("T0003", 3)
    assert rt.unpack(set_dir, dest) == (0, 9)
    (dest / "T0003.xyz").write_text("changed")
    with pytest.raises(ValueError):
        rt.unpack(set_dir, dest)


def _git(cwd, *args):
    return subprocess.run(["git", *args], cwd=cwd, check=True, capture_output=True, text=True)


def test_submit_commits_only_the_own_review_file_on_a_new_branch(tmp_path):
    origin = tmp_path / "origin.git"
    _git(tmp_path, "init", "--quiet", "--bare", "-b", "main", str(origin))
    repo = tmp_path / "clone"
    _git(tmp_path, "clone", "--quiet", str(origin), str(repo))
    for k, v in (("user.name", "T"), ("user.email", "t@example.org")):
        _git(repo, "config", k, v)
    (repo / "README").write_text("x\n")
    _git(repo, "add", "README")
    _git(repo, "commit", "--quiet", "-m", "init")
    _git(repo, "push", "--quiet", "origin", "main")
    folder = _set(tmp_path / "work", 6)
    a, _ = rt.assign(folder, ["Ann"], overlap=0, block=6, seed=1, set_name="demo")
    rt.write_assignment(folder / rt.ASSIGNMENT_NAME, a)
    with pytest.raises(FileNotFoundError):
        rt.submit(folder, "Ann", repo=repo, push=False)
    _rate(folder, "Ann", [("T0002.xyz", "pass", [])])
    res = rt.submit(folder, "ann", repo=repo, push=False)
    assert res["path"] == "ratings/demo/review_Ann.json" and res["rated"] == 1
    files = _git(repo, "show", "--name-only", "--format=", res["branch"]).stdout.split()
    assert files == ["ratings/demo/review_Ann.json"]
    sent = json.loads(_git(repo, "show", f"{res['branch']}:ratings/demo/review_Ann.json").stdout)
    assert sent["submitted"]["rater"] == "Ann" and sent["reviews"]
    assert _git(repo, "worktree", "list").stdout.count("\n") == 1     # temporary worktree gone
    # Someone else on another machine picked the same name: refused.
    _git(repo, "push", "--quiet", "origin", f"{res['branch']}:main")
    other = _set(tmp_path / "other", 6)
    rt.write_assignment(other / rt.ASSIGNMENT_NAME, a)
    _rate(other, "Ann", [("T0003.xyz", "block", ["H"])])
    with pytest.raises(ValueError, match="another reviewer"):
        rt.submit(other, "Ann", repo=repo, push=False)
    # The same reviewer may send again (the token matches).
    _rate(folder, "Ann", [("T0004.xyz", "pass", [])])
    assert rt.submit(folder, "Ann", repo=repo, push=False)["rated"] == 2


def test_submit_without_gh_pushes_the_branch_and_says_how_to_open_the_pull_request(
        tmp_path, monkeypatch):
    origin = tmp_path / "origin.git"
    _git(tmp_path, "init", "--quiet", "--bare", "-b", "main", str(origin))
    repo = tmp_path / "clone"
    _git(tmp_path, "clone", "--quiet", str(origin), str(repo))
    for k, v in (("user.name", "T"), ("user.email", "t@example.org")):
        _git(repo, "config", k, v)
    (repo / "README").write_text("x\n")
    _git(repo, "add", "README")
    _git(repo, "commit", "--quiet", "-m", "init")
    _git(repo, "push", "--quiet", "origin", "main")
    folder = _set(tmp_path / "work", 3)
    _rate(folder, "Ann", [("T0001.xyz", "pass", [])])
    monkeypatch.setattr(rt.shutil, "which", lambda _name: None)
    res = rt.submit(folder, "Ann", repo=repo, set_name="demo")
    assert res["pushed_to"] == "origin" and res["pr"] is None
    assert "open a pull request" in res["next"]
    assert res["branch"] in _git(origin, "branch", "--list").stdout


def test_a_published_set_is_assigned_from_its_archive(tmp_path):
    set_dir = _published(tmp_path)
    assert rt.structure_files(set_dir) == [f"T{k:04d}.xyz" for k in range(1, 9)]
    assert rv.main(["assign", str(set_dir), "Ann,Ben", "--overlap", "2", "--seed", "3", "--block", "3",
                    "--test-fraction", "0.25", "--set", "demo"]) == 0
    a = rt.load_assignment(set_dir / rt.ASSIGNMENT_NAME)
    assert a["set"] == "demo" and len(a["split"]["test"]) == 2
    assert len(a["raters"]["Ann"]) == len(a["raters"]["Ben"]) == 3
    dest = tmp_path / "mine"
    rt.unpack(set_dir, dest)
    assert rt.load_assignment(dest / rt.ASSIGNMENT_NAME) == a


def test_a_group_never_straddles_the_test_split(team_dir):
    groups = {f"T{k:04d}.xyz": f"g{(k - 1) // 2}" for k in range(1, 41)}   # pairs
    a, _ = rt.assign(team_dir, ["Ann"], overlap=0, seed=13, test_fraction=0.3, groups=groups)
    test = set(a["split"]["test"])
    assert a["split"]["grouped"] and 12 <= len(test) <= 13
    for k in range(1, 41, 2):
        pair = {f"T{k:04d}.xyz", f"T{k + 1:04d}.xyz"}
        assert pair <= test or not pair & test



def test_overlapping_self_drawn_blocks_count_for_the_agreement(team_dir, tmp_path):
    a, _ = rt.assign(team_dir, ["Ann"], overlap=0, block=4, seed=2, self_block=30)
    rt.write_assignment(team_dir / rt.ASSIGNMENT_NAME, a)
    zoe = rt.self_share(team_dir, a, "Zoe")["files"]
    yan = rt.self_share(team_dir, a, "Yan")["files"]
    both = sorted(set(zoe) & set(yan))
    assert len(both) >= 2                                           # 30 of 36 each
    only_zoe = sorted(set(zoe) - set(yan))[0]
    _rate(team_dir, "Zoe", [(both[0], "block", ["H"]), (both[1], "pass", []), (only_zoe, "pass", [])])
    _rate(team_dir, "Yan", [(both[0], "block", ["H"]), (both[1], "block", ["TOPO"])])
    for n, f in (("Zoe", zoe), ("Yan", yan)):
        path = team_dir / f"review_{n}.json"
        d = rv.load_review(path)
        d["assignment"] = rt.self_share(team_dir, a, n)
        rv.write_json_atomic(path, d)
    res = rt.agreement(rt.collect_reviews(team_dir))
    assert res["method"] == "Cohen" and res["frames"] == 2 and res["raters"] == ["Yan", "Zoe"]
    assert res["verdict"]["agreement"] == pytest.approx(0.5)
    out = tmp_path / "m.csv"
    rt.write_merged_csv(rt.collect_reviews(team_dir), out, a)
    rows = {(r["rater"], r["file"]): r["share"] for r in csv.DictReader(out.open())}
    assert rows[("Zoe", only_zoe)] == "own" and rows[("Yan", both[0])] == "own"


def test_fleiss_with_a_varying_number_of_raters():
    assert rt.fleiss_kappa([3, 0, 2], [3, 3, 3]) == pytest.approx(rt.fleiss_kappa([3, 0, 2], 3))
    # items: 2 of 2 yes, 0 of 3 yes, 1 of 2 yes
    p = 3 / 7
    pbar = (1 + 1 + 0) / 3
    pe = p * p + (1 - p) ** 2
    assert rt.fleiss_kappa([2, 0, 1], [2, 3, 2]) == pytest.approx((pbar - pe) / (1 - pe))
    assert rt.fleiss_kappa([1], [1]) is None


def test_unpack_with_a_new_name_draws_the_block_into_the_review_file(tmp_path):
    set_dir = _published(tmp_path)
    a, _ = rt.assign(set_dir, ["Ann"], overlap=2, block=2, seed=3, self_block=3)
    rt.write_assignment(set_dir / rt.ASSIGNMENT_NAME, a)
    dest = tmp_path / "mine"
    assert rv.main(["unpack", str(set_dir), str(dest), "--name", "Zoe"]) == 0
    own = rv.load_review(dest / "review_Zoe.json")["assignment"]
    assert len(own["files"]) == 3 and not set(own["files"]) & set(a["shared"] + a["raters"]["Ann"])
    assert rt.unpack_for(set_dir, dest, "Zoe")[2].startswith("Zoe: a block of 3")
    assert rt.unpack_for(set_dir, dest, "ann")[2] == "Ann: your share is 4 structures"
    assert not (dest / "review_ann.json").exists() and not (dest / "review_Ann.json").exists()
    s = rv.ReviewSession.open(dest, name="Zoe")
    assert rv.load_review(dest / "review_Zoe.json")["assignment"] == own and len(s) == 8


def test_when_nothing_is_free_a_new_block_holds_second_ratings(team_dir):
    a, _ = rt.assign(team_dir, ["Ann"], overlap=0, block=40, seed=6, self_block=10)
    rt.write_assignment(team_dir / rt.ASSIGNMENT_NAME, a)
    panel, s = _panel(team_dir, "Zoe")
    assert s is not None, panel.status_html.value
    own = rv.load_review(team_dir / "review_Zoe.json")["assignment"]["files"]
    assert len(own) == 10 and set(own) <= set(a["raters"]["Ann"])
    assert rv.main(["assign", str(team_dir), "Ben", "--block", "0"]) == 0
    assert rt.load_assignment(team_dir / rt.ASSIGNMENT_NAME)["self_block"] == 10  # kept
