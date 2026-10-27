"""A review set rated by a collective: the least rated first, overlapping
blocks, reservations that run out, reliability and consensus, the status of
a set, the check of a helper's pull request, and a new set.

Synthetic data only: invented names, a water molecule moved a little, and a
made-up chiral molecule; git runs through delfin.agent.safe_git in temporary
repositories, nothing goes to the network.
"""
from __future__ import annotations

import csv
import gzip
import json
import random
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pytest

from delfin import review as rv
from delfin import review_check as rc
from delfin import review_stats as rs
from delfin import review_team as rt
from delfin.agent import safe_git

NOW = datetime(2026, 10, 27, 12, 0, tzinfo=timezone.utc)


def _water(stem: str, k: int = 0) -> str:
    return (f"3\n{stem}\n"
            f"O 0.000 0.000 {0.117 + 0.05 * k:.3f}\n"  # far enough apart: never identical
            "H 0.000 0.757 -0.467\nH 0.000 -0.757 -0.467\n")


def _set(folder: Path, n: int) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    for k in range(1, n + 1):
        (folder / f"T{k:04d}.xyz").write_text(_water(f"T{k:04d}", k))
    return folder


def _rate(folder, name, plan, categories=("H", "TOPO")):
    s = rv.ReviewSession.open(folder, name=name, categories=list(categories))
    for fname, verdict, cats in plan:
        s.go_to_file(fname)
        s.rate(verdict, cats)


def _review(ratings: Path, name: str, files, verdict="pass", **extra) -> None:
    """A review file as it lands in ratings/<set>/ (written here, not by the tab)."""
    data = {"version": 1, "categories": ["H", "TOPO"], "files": {},
            "reviews": {f: {"0": {"verdict": verdict, "categories": [], "note": "",
                                  "file": f, "frame": 0}} for f in files}, **extra}
    rv.write_json_atomic(ratings / f"review_{name}.json", data)


@pytest.fixture
def team(tmp_path):
    folder = _set(tmp_path / "set", 60)
    ratings = tmp_path / "ratings"
    ratings.mkdir()
    return folder, ratings


# ------------------------------------------------------------------ the draw
def test_a_new_block_takes_the_unrated_first_then_second_ratings(team):
    folder, ratings = team
    a, _ = rt.assign(folder, [], overlap=0, seed=3, redundancy=0)
    files = rt.structure_files(folder)
    _review(ratings, "Ann", files[:50])
    own = rt.self_share(folder, a, "Zoe", 8, ratings=[ratings])["files"]
    assert set(own) <= set(files[50:]) and len(own) == 8          # nothing rated by Ann
    own = rt.self_share(folder, a, "Yan", 15, ratings=[ratings])["files"]
    assert set(files[50:]) <= set(own) and len(own) == 15         # the 10 free, then 5 twice
    # The same state draws the same block; another name another one.
    assert rt.self_share(folder, a, "Yan", 15, ratings=[ratings])["files"] == own
    assert rt.self_share(folder, a, "Xia", 15, ratings=[ratings])["files"] != own


def test_reservations_count_until_they_run_out(team):
    folder, ratings = team
    a, _ = rt.assign(folder, [], overlap=0, seed=4, redundancy=0)
    files = rt.structure_files(folder)
    (ratings / "claims").mkdir()
    claim = {"name": "Ann", "files": files[:55], "date": (NOW - timedelta(days=3)).isoformat()}
    (ratings / "claims" / "Ann.json").write_text(json.dumps(claim))
    own = rt.self_share(folder, a, "Zoe", 5, ratings=[ratings], now=NOW)["files"]
    assert set(own) == set(files[55:])
    # 29 days after the claim it has run out: its unrated structures are free again ...
    later = NOW + timedelta(days=26)
    res = rt.reservations(a, {}, rt.collect_claims(ratings), now=later)
    assert not res["ann"]["active"] and res["ann"]["kind"] == "claim"
    assert not any(rt.coverage(files, {}, res).values())
    # ... but what Ann has rated stays rated.
    _review(ratings, "Ann", files[:20])
    reviews = rt.collect_reviews(ratings)
    res = rt.reservations(a, reviews, rt.collect_claims(ratings), now=later)
    assert sum(rt.coverage(files, reviews, res).values()) == 20
    # A block drawn for oneself and kept in a review file reserves as well.
    _review(ratings, "Bob", [], assignment={"self": True, "name": "Bob", "files": files[40:50],
                                            "drawn": NOW.isoformat()})
    reviews = rt.collect_reviews(ratings)
    res = rt.reservations(a, reviews, {}, now=NOW)
    assert res["bob"]["active"] and res["bob"]["kind"] == "self"
    assert rt.coverage(files, reviews, res)[files[45]] == 1


def test_the_redundant_share_overlaps_many_reviewers_and_is_spread(tmp_path):
    folder = _set(tmp_path / "set", 300)
    ratings = tmp_path / "ratings"
    ratings.mkdir()
    a, _ = rt.assign(folder, [], overlap=10, seed=5, redundancy=0.3)
    wanted = rt.raters_wanted(a, rt.structure_files(folder))
    assert len(wanted) == 290 and sum(n > 1 for n in wanted.values()) == 87
    assert sum(n == 3 for n in wanted.values()) == 29
    assert not set(a["shared"]) & set(wanted)
    blocks = {}
    for name in ("Ann", "Ben", "Cem", "Dan"):
        own = rt.self_share(folder, a, name, 40, ratings=[ratings], now=NOW)
        assert not set(own["files"]) & set(a["shared"])
        blocks[name] = set(own["files"])
        _review(ratings, name, [], assignment=own)
    # Every later block shares structures with earlier ones (12 places of 40 are
    # redundant), only ones meant for more than one reviewer.
    for later, earlier in (("Ben", ["Ann"]), ("Cem", ["Ann", "Ben"]),
                           ("Dan", ["Ann", "Ben", "Cem"])):
        common = blocks[later] & set().union(*(blocks[e] for e in earlier))
        assert 1 <= len(common) <= 12
        assert all(wanted[f] > 1 for f in common)
    assert len([e for e in ("Ann", "Ben", "Cem") if blocks["Dan"] & blocks[e]]) >= 2
    # No structure gets more reviewers than it is meant for while others have none.
    reviews = rt.collect_reviews(ratings)
    cover = rt.coverage(rt.structure_files(folder), reviews,
                        rt.reservations(a, reviews, {}, now=NOW))
    assert all(cover[f] <= wanted[f] for f in wanted)


def test_assign_reads_the_ratings_and_the_claims_too(team):
    folder, ratings = team
    a, _ = rt.assign(folder, [], overlap=0, seed=6, redundancy=0)
    files = rt.structure_files(folder)
    _review(ratings, "Zoe", files[:30])
    b, msgs = rt.assign(folder, ["Ann"], block=30, existing=a, ratings=[ratings], now=NOW)
    assert set(b["raters"]["Ann"]) == set(files[30:])
    assert b["redundancy"] == 0 and any("Ann: 30 own" in m for m in msgs)
    with pytest.raises(ValueError):
        rt.assign(folder, ["Cem"], existing=b, redundancy=2)


def test_the_set_version_and_the_redundancy_are_kept_in_the_assignment(tmp_path):
    set_dir = tmp_path / "review_sets" / "demo"
    set_dir.mkdir(parents=True)
    with gzip.open(set_dir / rt.SET_ARCHIVE_NAME, "wt") as fh:
        for k in range(1, 9):
            fh.write(_water(f"T{k:04d}", k))
    a, _ = rt.assign(set_dir, [], overlap=2, seed=1)
    assert a["redundancy"] == rt.DEFAULT_REDUNDANCY
    assert a["set_version"] == rt.set_version(set_dir) and a["set_version"].startswith("sha256:")
    with gzip.open(set_dir / rt.SET_ARCHIVE_NAME, "wt") as fh:
        fh.write(_water("T0001", 99))
    with pytest.raises(ValueError):
        rt.assign(set_dir, ["Ann"], existing=a)


def test_unpack_points_the_folder_at_the_ratings_of_the_clone(tmp_path):
    repo = tmp_path / "clone"
    set_dir = repo / "review_sets" / "demo"
    set_dir.mkdir(parents=True)
    with gzip.open(set_dir / rt.SET_ARCHIVE_NAME, "wt") as fh:
        for k in range(1, 13):
            fh.write(_water(f"T{k:04d}", k))
    (set_dir / "kategorien.txt").write_text("H\nTOPO\n")
    a, _ = rt.assign(set_dir, [], overlap=0, seed=2, redundancy=0, self_block=4)
    rt.write_assignment(set_dir / rt.ASSIGNMENT_NAME, a)
    ratings = repo / "ratings" / "demo"
    ratings.mkdir(parents=True)
    _review(ratings, "Ann", [f"T{k:04d}.xyz" for k in range(1, 9)])
    dest = tmp_path / "mine"
    assert rt.unpack(set_dir, dest)[0] == 14                        # 12 + 2, not the marker
    assert rt.ratings_dirs(dest)[-1].resolve() == ratings.resolve()
    _, _, msg = rt.unpack_for(set_dir, dest, "Zoe")
    own = rv.load_review(dest / "review_Zoe.json")["assignment"]
    assert set(own["files"]) == {f"T{k:04d}.xyz" for k in range(9, 13)}
    assert own["drawn_from"]["reviews"] == 1 and "a block of 4" in msg
    # The coordinator's folder review_sets/<set> sees ratings/<set> by itself.
    assert ratings.resolve() in [d.resolve() for d in rt.ratings_dirs(set_dir)]


# ------------------------------------------------------------------ statistics
def test_krippendorffs_alpha_on_the_published_example():
    # Krippendorff (2011), "Computing Krippendorff's alpha-reliability": four
    # observers, twelve units, nominal alpha = 0.743.
    units = [[1, 1, 1], [2, 2, 3, 2], [3, 3, 3, 3], [3, 3, 3, 3], [2, 2, 2, 2],
             [1, 2, 3, 4], [4, 4, 4, 4], [1, 1, 2, 1], [2, 2, 2, 2], [5, 5, 5],
             [1, 1], [3]]
    assert rs.krippendorff_alpha(units) == pytest.approx(0.743, abs=5e-4)
    assert rs.krippendorff_alpha([[1], [2]]) is None                 # nothing pairable
    assert rs.krippendorff_alpha([[1, 1], [1, 1]]) is None           # one value only
    a = rs.alpha_with_ci(units, n_boot=300, seed=1)
    b = rs.alpha_with_ci(units, n_boot=300, seed=1)
    assert a == b and a["low"] <= a["alpha"] <= a["high"] and a["items"] == 11


def _simulated(seed=7, n=300):
    rng = random.Random(seed)
    truth = {f"I{k}": rng.random() < 0.4 for k in range(n)}
    skill = {"good1": (0.9, 0.95), "good2": (0.85, 0.9), "good3": (0.9, 0.85),
             "coin": (0.5, 0.5)}
    votes = {}
    for item, t in truth.items():
        votes[item] = {}
        for r, (sens, spec) in skill.items():
            if r != "coin" and rng.random() < 0.3:
                continue                                             # missing vote
            votes[item][r] = int(rng.random() < (sens if t else 1 - spec))
    return truth, votes


def test_dawid_skene_finds_the_reliable_reviewers_and_beats_the_majority():
    truth, votes = _simulated()
    fit = rs.dawid_skene(votes)
    assert fit["converged"]
    sens, spec = rs.sensitivity_specificity(fit["confusion"]["good1"])
    assert sens > 0.8 and spec > 0.85
    sens, spec = rs.sensitivity_specificity(fit["confusion"]["coin"])
    assert 0.3 < sens < 0.7 and 0.3 < spec < 0.7
    ds_right = sum((fit["posterior"][i][1] >= 0.5) == t for i, t in truth.items())
    major = sum((sum(votes[i].values()) * 2 > len(votes[i])) == t for i, t in truth.items())
    assert ds_right >= major and ds_right / len(truth) > 0.9
    assert rs.dawid_skene(votes)["posterior"] == fit["posterior"]    # deterministic


def _chiral():
    # C with four different substituents: its mirror image is another molecule.
    return ["C", "H", "F", "Cl", "Br"], np.array(
        [[0, 0, 0], [0.63, 0.63, 0.63], [-0.8, -0.8, 0.8], [-1.0, 1.0, -1.0], [1.1, -1.1, -1.1]])


def _turn(xyz, seed):
    q, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    q = q if np.linalg.det(q) > 0 else -q
    return np.round(xyz @ q.T + np.array([0.3, -0.2, 1.0]), 6)


def test_a_structure_shown_twice_is_found_but_not_its_mirror_image():
    elems, xyz = _chiral()
    mirror = xyz * np.array([-1, 1, 1])
    other = xyz.copy()
    other[2] *= 1.3
    perm = [0, 1, 2, 4, 3]                       # another atom order within the file
    frames = [("A.xyz", elems, xyz), ("B.xyz", elems, _turn(xyz, 1)),
              ("M.xyz", elems, _turn(mirror, 2)), ("O.xyz", elems, _turn(other, 3)),
              ("P.xyz", [elems[i] for i in perm], _turn(xyz, 4)[perm])]
    assert rs.find_repeats(frames) == [("A.xyz", "B.xyz"), ("A.xyz", "P.xyz"), ("B.xyz", "P.xyz")]
    verdicts = {"Ann": {"A.xyz": True, "B.xyz": True, "P.xyz": False}, "Ben": {"A.xyz": False}}
    res = rs.self_consistency(verdicts, [("A.xyz", "B.xyz"), ("A.xyz", "P.xyz")])
    assert res["pairs_rated"] == 2 and res["same"] == pytest.approx(0.5)
    assert res["raters"] == {"Ann": {"pairs": 2, "same": 1}}


def _rated_team(tmp_path):
    folder = _set(tmp_path / "set", 30)
    (folder / "kategorien.txt").write_text("H\nTOPO\n")
    a, _ = rt.assign(folder, [], overlap=6, seed=9, test_fraction=0.3, redundancy=0)
    rt.write_assignment(folder / rt.ASSIGNMENT_NAME, a)
    files = rt.structure_files(folder)
    rng = random.Random(2)
    for name in ("Ann", "Ben", "Cem"):
        plan = []
        for f in files:
            if rng.random() < 0.7:
                bad = (int(f[1:5]) % 3 == 0) != (rng.random() < 0.1)
                plan.append((f, "block" if bad else "pass", ["H"] if bad else []))
        _rate(folder, name, plan)
    return folder, a


def test_merge_reports_alpha_dawid_skene_and_keeps_the_test_part_out(tmp_path, capsys):
    folder, a = _rated_team(tmp_path)
    test = set(a["split"]["test"])
    reviews = rt.collect_reviews(folder)
    calib = rt.split_items(a, "calibration", rt.structure_files(folder))
    res = rt.reliability(reviews, calib, categories=["H", "TOPO"], n_boot=200)
    assert set(res["questions"]) == {"ok/defect", "H", "TOPO"}
    assert res["questions"]["ok/defect"]["alpha"]["alpha"] > 0.5
    assert not {f for f, _ in res["consensus"]} & test
    assert set(res["raters"]) == {"Ann", "Ben", "Cem"}
    key = tmp_path / "private_key.csv"
    key.write_text("id,source\nT0001,crystal_cod\nT0002,manta\nT0004,crystal\n")
    pairs = tmp_path / "pairs.json"
    pairs.write_text(json.dumps({"T0001.xyz": "g", "T0002.xyz": "g", "T0004.xyz": "h"}))
    out = tmp_path / "m.csv"
    assert rv.main(["merge", str(folder), "--csv", str(out), "--split", "calibration",
                    "--bootstrap", "100", "--controls", str(key),
                    "--repeats", str(pairs)]) == 0
    text = capsys.readouterr().out
    assert "Krippendorff" in text and "sensitivity" in text and "Controls" in text
    rows = list(csv.DictReader((tmp_path / "m_consensus.csv").open()))
    assert rows and not {r["file"] for r in rows} & test
    assert all(r["split"] == "calibration" for r in rows)
    assert {r["consensus"] for r in rows} <= {"block", "pass"} and "p_H" in rows[0]
    raters = list(csv.DictReader((tmp_path / "m_raters.csv").open()))
    assert {r["question"] for r in raters} == {"ok/defect", "H", "TOPO"}
    assert rv.main(["merge", str(folder), "--no-csv", "--split", "test",
                    "--bootstrap", "0"]) == 0
    assert rt.load_controls(key) == {"T0001.xyz", "T0004.xyz"}
    assert rt.load_repeats(pairs) == [("T0001.xyz", "T0002.xyz")]


def test_the_status_counts_ratings_reservations_and_contributors(tmp_path):
    folder = _set(tmp_path / "set", 20)
    ratings = tmp_path / "ratings"
    (ratings / "claims").mkdir(parents=True)
    a, _ = rt.assign(folder, ["Max"], overlap=2, block=4, seed=1, test_fraction=0.25,
                     redundancy=0)
    rt.write_assignment(folder / rt.ASSIGNMENT_NAME, a)
    files = [f for f in rt.structure_files(folder) if f not in a["shared"]]
    _review(ratings, "Ann", files[:6], credit="Ann Lee")
    _review(ratings, "Ben", files[3:6] + a["shared"])
    _review(ratings, "Cem", files[5:6])
    for name, days in (("Dan", 2), ("Eve", 40)):
        (ratings / "claims" / f"{name}.json").write_text(json.dumps(
            {"name": name, "files": files[10:12],
             "date": (NOW - timedelta(days=days)).isoformat()}))
    st = rt.status(folder, [ratings], now=NOW, n_boot=50)
    assert st["hist"] == {"0": 12, "1": 3, "2": 2, "3+": 1}
    assert st["shared"] == {"files": 2, "raters": [1, 1]}
    states = {r["name"]: r["active"] for r in st["reservations"]}
    assert states == {"Max": True, "Dan": True, "Eve": False}
    assert st["credits"] == ["Ann Lee"]
    assert sum(t for _, t in st["parts"].values()) == 20
    assert set(st["parts"]) == {"calibration", "test"}
    text = rt.format_status(st)
    assert "2 active, 1 expired" in text
    assert "contributors: Ann Lee" in text and "alpha" in text
    assert rv.main(["status", str(folder), "--ratings", str(ratings), "--bootstrap", "10"]) == 0


# ------------------------------------------------------------------ git: claim, submit, check
def _g(cwd, *args):
    res = safe_git.run(list(args), cwd)
    assert res.returncode == 0, res.stderr
    return res.stdout


def _published_repo(tmp_path, n=8):
    origin = tmp_path / "origin.git"
    _g(tmp_path, "init", "--quiet", "--bare", "-b", "main", str(origin))
    repo = tmp_path / "clone"
    _g(tmp_path, "clone", "--quiet", str(origin), str(repo))
    for k, v in (("user.name", "T"), ("user.email", "t@example.org")):
        _g(repo, "config", k, v)
    set_dir = repo / "review_sets" / "demo"
    set_dir.mkdir(parents=True)
    with gzip.open(set_dir / rt.SET_ARCHIVE_NAME, "wt", encoding="utf-8") as fh:
        for k in range(1, n + 1):
            fh.write(_water(f"T{k:04d}", k))
    (set_dir / "kategorien.txt").write_text("H\nTOPO\n")
    a, _ = rt.assign(set_dir, [], overlap=2, seed=3, self_block=3, set_name="demo")
    rt.write_assignment(set_dir / rt.ASSIGNMENT_NAME, a)
    (repo / "ratings").mkdir()
    (repo / "ratings" / "README.md").write_text("ratings\n")
    _g(repo, "add", "-A")
    _g(repo, "commit", "--quiet", "-m", "the set")
    _g(repo, "push", "--quiet", "origin", "main")
    return repo, set_dir, a


def _pr_check(repo, branch, base="origin/main", author="zoe-gh"):
    return rc.check_pull_request(repo, base, branch, author=author, now=NOW)


def test_claim_reserves_the_block_in_a_pull_request_of_its_own(tmp_path):
    repo, set_dir, a = _published_repo(tmp_path)
    dest = tmp_path / "mine"
    rt.unpack_for(set_dir, dest, None)
    res = rt.claim(dest, "Zoe", repo=repo, push=False, github="zoe-gh")
    assert res["path"] == "ratings/demo/claims/Zoe.json" and res["files"] == 3
    assert res["owner"] == "ratings/demo/owners/Zoe.json"
    assert sorted(_g(repo, "show", "--name-only", "--format=", res["branch"]).split()) == \
        sorted([res["path"], res["owner"]])
    sent = json.loads(_g(repo, "show", f"{res['branch']}:{res['path']}"))
    own = rv.load_review(dest / "review_Zoe.json")["assignment"]
    assert sent["files"] == own["files"] and sent["set"] == "demo" and sent["days"] == 28
    assert sent["set_version"] == a["set_version"] and sent["token"]
    bound = json.loads(_g(repo, "show", f"{res['branch']}:{res['owner']}"))
    assert bound["login"] == "zoe-gh" and bound["name"] == "Zoe"
    first = _pr_check(repo, res["branch"])
    assert first["problems"] == [] and first["first"] and not first["automerge"]
    # Claiming again keeps the block (a new date, the same files).
    again = rt.claim(dest, "Zoe", repo=repo, push=False)
    assert json.loads(_g(repo, "show", f"{again['branch']}:{again['path']}"))["files"] == \
        own["files"]
    assert again["owner"] is None and "--github" in again["next"]
    a2 = rt.load_assignment(dest / rt.ASSIGNMENT_NAME)
    a2["raters"]["Max"] = ["T0001.xyz"]
    rt.write_assignment(dest / rt.ASSIGNMENT_NAME, a2)
    with pytest.raises(ValueError, match="nothing to claim"):
        rt.claim(dest, "max", repo=repo, push=False)


def test_a_helpers_pull_request_is_checked(tmp_path, monkeypatch):
    repo, set_dir, a = _published_repo(tmp_path)
    dest = tmp_path / "mine"
    rt.unpack_for(set_dir, dest, "Zoe")
    own = rv.load_review(dest / "review_Zoe.json")["assignment"]["files"]
    _rate(dest, "Zoe", [(own[0], "block", ["H"]), (own[1], "pass", [])])
    res = rt.submit(dest, "Zoe", repo=repo, push=False, credit="Zoe Q.", github="zoe-gh")
    ok = _pr_check(repo, res["branch"])
    assert ok["helper"] and ok["kind"] == "review" and ok["problems"] == [], ok
    assert ok["first"] and not ok["known"] and not ok["automerge"]
    assert "maintainer review" in rc.format_check(ok)
    sent = json.loads(_g(repo, "show", f"{res['branch']}:{res['path']}"))
    assert sent["credit"] == "Zoe Q." and sent["submitted"]["set_version"] == a["set_version"]
    assert rv.load_review(dest / "review_Zoe.json")["credit"] == "Zoe Q."      # kept here too
    # The first pull request binds the name to the account it comes from.
    assert "comes from mallory" in " ".join(_pr_check(repo, res["branch"], author="mallory")
                                            ["problems"])

    def branch_with(name, files, delete=()):
        work = tmp_path / f"wt_{name}"
        _g(repo, "worktree", "add", "--quiet", "-b", name, str(work), "origin/main")
        for rel, body in files.items():
            (work / rel).parent.mkdir(parents=True, exist_ok=True)
            (work / rel).write_text(body if isinstance(body, str) else json.dumps(body))
        for rel in delete:
            (work / rel).unlink()
        _g(work, "add", "-A")
        _g(work, "commit", "--quiet", "-m", name)
        return name

    def problems(name, files, author="zoe-gh", **kw):
        return " | ".join(_pr_check(repo, branch_with(name, files, **kw), author=author)
                          ["problems"])

    rel = "ratings/demo/review_Zoe.json"
    orel = "ratings/demo/owners/Zoe.json"
    owner = {"version": 1, "set": "demo", "name": "Zoe", "login": "zoe-gh", "since": "x"}
    bad = json.loads(json.dumps(sent))
    bad["reviews"][own[0]]["0"]["categories"] = ["NOPE"]
    bad["reviews"][own[0]]["0"]["atoms"] = [0, 7]
    bad["reviews"]["T9999.xyz"] = {"0": {"verdict": "pass"}}
    text = problems("b1", {rel: bad, orel: owner})
    assert "NOPE" in text and "atoms outside 0..2" in text and "T9999.xyz is not in the set" in text
    assert "is not the name of the file" in problems("b2", {"ratings/demo/review_Yan.json": sent})
    assert "exactly one" in problems("b3", {rel: sent, orel: owner, "delfin/x.py": "x = 1\n"})
    assert "exactly one" in problems("b4", {rel: sent, "ratings/demo/review_Q.json": sent})
    old = dict(sent, submitted=dict(sent["submitted"], set_version="sha256:0"))
    assert "set version" in problems("b5", {rel: old, orel: owner})
    assert "credit" in problems("b6", {rel: dict(sent, credit="x" * 81), orel: owner})
    assert "not valid JSON" in problems("b7", {rel: "{"})
    assert "must add ratings/demo/owners/Zoe.json" in problems("b7b", {rel: sent})
    assert "not that of" in problems("b7c", {rel: sent, "ratings/demo/owners/Yan.json": owner})
    assert "no GitHub login" in problems("b7d", {rel: sent, orel: dict(owner, login="a b")})
    # Not a helper's pull request: no rating, claim or owner file.
    other = _pr_check(repo, branch_with("b8", {"ratings/README.md": "more\n"}))
    assert not other["helper"] and other["problems"] == []
    # Once Zoe's file and owner file are on main: only Zoe's token, only her account.
    _g(repo, "push", "--quiet", "origin", f"{res['branch']}:main")
    _g(repo, "fetch", "--quiet", "origin")
    stolen = dict(sent, submitted=dict(sent["submitted"], token="another"))
    assert "another reviewer" in problems("b9", {rel: stolen})
    again = branch_with("b10", {rel: dict(sent, credit="")})
    known = _pr_check(repo, again)
    assert known["problems"] == [] and known["known"] and known["automerge"]
    assert "belongs to the GitHub account zoe-gh" in \
        " ".join(_pr_check(repo, again, author="mallory")["problems"])
    assert "may not change it" in problems("b12", {rel: sent, orel: dict(owner, login="mallory")},
                                           author="mallory")
    assert "deleted" in problems("b11", {}, delete=[rel])
    # Claims: a valid first one; then too large, from the future, another length.
    claim = {"version": 1, "set": "demo", "set_version": a["set_version"], "name": "Yan",
             "files": own, "date": NOW.isoformat(), "days": 28, "token": "t"}
    crel = "ratings/demo/claims/Yan.json"
    yan = dict(owner, name="Yan", login="yan-gh")
    assert problems("c1", {crel: claim, "ratings/demo/owners/Yan.json": yan},
                    author="yan-gh") == ""
    files = rt.structure_files(set_dir)
    text = problems("c2", {crel: dict(claim, files=files, days=90,
                                       date=(NOW + timedelta(days=9)).isoformat())},
                    author="yan-gh")
    assert "more than one block" in text and "future" in text and "28 days" in text
    out = tmp_path / "gh_output"
    monkeypatch.setenv("GITHUB_OUTPUT", str(out))
    monkeypatch.delenv("GITHUB_STEP_SUMMARY", raising=False)
    assert rc.main(["--repo", str(repo), "--base", "origin/main", "--head", "b10",
                    "--author", "zoe-gh"]) == 0
    assert rc.main(["--repo", str(repo), "--base", "origin/main", "--head", "b1"]) == 1
    text = out.read_text()
    assert "automerge=true" in text and "valid=true" in text and "valid=false" in text


def test_the_workflow_checks_from_the_base_and_merges_only_known_accounts():
    import yaml
    root = Path(__file__).resolve().parents[1]
    wf = yaml.safe_load((root / ".github" / "workflows" / "ratings.yml").read_text())
    on = wf.get("on") or wf.get(True)
    assert list(on) == ["pull_request_target"]
    assert on["pull_request_target"]["paths"] == ["ratings/**"]
    assert wf["permissions"] == {}
    jobs = wf["jobs"]
    check = jobs["check"]
    assert check["permissions"] == {"contents": "read"}
    steps = check["steps"]
    assert steps[0]["with"]["ref"] == "${{ github.event.pull_request.base.sha }}"
    assert steps[0]["with"]["persist-credentials"] is False
    # Only one checkout, of the base; the pull request is fetched as objects.
    assert sum("actions/checkout" in str(s.get("uses", "")) for s in steps) == 1
    fetch = next(s for s in steps if "git fetch" in s.get("run", ""))
    assert "refs/pull/${PR_NUMBER}/head" in fetch["run"] and "checkout" not in fetch["run"]
    run = steps[-1]
    assert "-m delfin.review_check" in run["run"] and '--author "$AUTHOR"' in run["run"]
    assert run["env"]["AUTHOR"] == "${{ github.event.pull_request.user.login }}"
    for name in ("automerge", "first"):
        job = jobs[name]
        assert job["needs"] == "check"
        assert not any("uses" in s for s in job["steps"])           # nothing checked out
        assert all("${{" not in s["run"] for s in job["steps"])      # event data via env only
    assert jobs["automerge"]["if"] == "needs.check.outputs.automerge == 'true'"
    assert jobs["automerge"]["permissions"] == {"contents": "write", "pull-requests": "write"}
    assert "gh pr merge --auto --squash" in jobs["automerge"]["steps"][0]["run"]
    assert jobs["first"]["permissions"] == {"pull-requests": "write"}
    assert "first contribution, maintainer review" in jobs["first"]["steps"][0]["run"]


# ------------------------------------------------------------------ identical structures
def _chiral_xyz(ident, seed, mirror=False):
    elems, xyz = _chiral()
    if mirror:
        xyz = xyz * np.array([-1, 1, 1])
    xyz = _turn(xyz, seed)
    rows = "".join(f"{e} {x:.6f} {y:.6f} {z:.6f}\n" for e, (x, y, z) in zip(elems, xyz))
    return f"{len(elems)}\n{ident}\n{rows}"


def _with_identical(folder: Path) -> list[set[str]]:
    """Waters T0001..T0020 and three copies (turned) of one chiral molecule,
    plus its mirror image, which is another structure."""
    _set(folder, 20)
    for k, ident in enumerate(("T0021", "T0022", "T0023")):
        (folder / f"{ident}.xyz").write_text(_chiral_xyz(ident, k))
    (folder / "T0024.xyz").write_text(_chiral_xyz("T0024", 9, mirror=True))
    return [{"T0021.xyz", "T0022.xyz", "T0023.xyz"}]


def test_identical_structures_stay_in_one_part_of_the_split(tmp_path):
    folder = tmp_path / "set"
    groups = _with_identical(folder)
    group_of, joined = rt.identical_groups(folder)
    assert joined == 3 and len({group_of[f] for f in groups[0]}) == 1
    assert group_of["T0024.xyz"] == "T0024.xyz"                       # the mirror image is not
    for seed in range(12):
        a, _ = rt.assign(folder, [], overlap=0, seed=seed, test_fraction=0.3)
        test = set(a["split"]["test"])
        assert groups[0] <= test or not groups[0] & test
        assert a["split"]["identical_kept_together"] == 3


def test_a_split_that_straddles_is_corrected_with_a_dated_record(tmp_path):
    folder = tmp_path / "set"
    groups = _with_identical(folder)
    a, _ = rt.assign(folder, ["Max"], overlap=2, block=4, seed=3, test_fraction=0.25)
    # As drawn before identical structures were kept together: one in each part.
    test = (set(a["split"]["test"]) - groups[0]) | {"T0021.xyz"}
    a["split"]["test"] = sorted(test)
    before = json.loads(json.dumps(a))
    b, msgs = rt.assign(folder, [], existing=a, regroup=True)
    new = set(b["split"]["test"])
    assert groups[0] <= new or not groups[0] & new
    assert len(new ^ test) == 1                                     # the minority moved
    change = b["split"]["changes"][-1]
    assert change["groups"] == 1 and change["moved_to_test"] + change["moved_to_calibration"] == 1
    assert change["test_size"] == len(new) and change["date"] and change["reason"]
    assert "T0021" not in json.dumps(change)                         # counts, no names
    assert b["shared"] == before["shared"] and b["raters"] == before["raters"]
    assert any("1 group(s)" in m for m in msgs)
    c, _ = rt.assign(folder, [], existing=b, regroup=True)
    assert c["split"]["test"] == b["split"]["test"] and len(c["split"]["changes"]) == 1
    rt.write_assignment(folder / rt.ASSIGNMENT_NAME, a)
    assert rv.main(["assign", str(folder), "--regroup-split"]) == 0
    assert rt.load_assignment(folder / rt.ASSIGNMENT_NAME)["split"]["test"] == b["split"]["test"]


def test_a_new_set_keeps_identical_structures_together(tmp_path):
    folder = tmp_path / "frames"
    groups = _with_identical(folder)
    src = tmp_path / "all.xyz"
    src.write_text("".join((folder / f).read_text() for f in rt.structure_files(folder)))
    cats = tmp_path / "kategorien.txt"
    cats.write_text("H\n")
    out = tmp_path / "review_sets" / "twins"
    rt.new_set(src, out, cats, seed=4, test_fraction=0.3, overlap=0)
    a = rt.load_assignment(out / rt.ASSIGNMENT_NAME)
    test = set(a["split"]["test"])
    assert groups[0] <= test or not groups[0] & test
    assert a["split"]["identical_kept_together"] == 3 and a["split"]["grouped"]


# ------------------------------------------------------------------ a new set
def test_a_new_set_gets_its_test_part_before_any_rating(tmp_path):
    src = tmp_path / "all.xyz"
    src.write_text("".join(_water(f"mol {k}", k) for k in range(1, 21)))
    cats = tmp_path / "kategorien.txt"
    cats.write_text("# the categories\nH\nTOPO\n")
    out = tmp_path / "repo" / "review_sets" / "fresh"
    with pytest.raises(ValueError, match="renumber"):
        rt.new_set(src, out, cats)
    with pytest.raises(ValueError, match="test part"):
        rt.new_set(src, out, cats, renumber=True, test_fraction=0)
    assert rv.main(["new-set", str(src), str(out), "--categories", str(cats), "--renumber",
                    "--seed", "5", "--overlap", "4"]) == 0
    a = rt.load_assignment(out / rt.ASSIGNMENT_NAME)
    assert a["set"] == "fresh" and len(a["shared"]) == 4 and len(a["split"]["test"]) == 4
    assert a["raters"] == {} and a["redundancy"] == rt.DEFAULT_REDUNDANCY
    assert a["set_version"] == rt.set_version(out)
    assert rt.structure_files(out)[0] == "T0001.xyz" and len(rt.structure_files(out)) == 20
    assert (out / "kategorien.txt").read_text() == "H\nTOPO\n"
    readme = (out / "README.md").read_text()
    assert "TODO" in readme and "seed\n`5:split`" in readme and "Cambridge" in readme
    first = (out / rt.SET_ARCHIVE_NAME).read_bytes()
    with pytest.raises(ValueError, match="not empty"):
        rt.new_set(src, out, cats, renumber=True)
    again = tmp_path / "again"
    rt.new_set(src, again, cats, renumber=True, seed=5, overlap=4)
    assert (again / rt.SET_ARCHIVE_NAME).read_bytes() == first          # reproducible

    def read(p):
        local = out / p.split("/", 2)[2]
        return local.read_bytes() if local.is_file() else None

    the_set = rc.base_set(read, "fresh")
    assert the_set["atoms"]["T0020.xyz"] == 3 and the_set["version"] == a["set_version"]
