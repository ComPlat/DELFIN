"""Reviewing as a team: share a set of structures out, bring the ratings back together.

A set of single-structure ``*.xyz`` files is too large for one reviewer, so it
is shared out: every reviewer rates a block of their own, and all of them rate
a few structures in common, so that their agreement can be measured.  The
share lives in ``assignment.json`` next to the structures (or in the folder
above, like ``kategorien.txt``)::

    {"version": 1, "seed": 7, "files_total": 5000,
     "shared": ["T0042.xyz", ...],            # rated by everyone
     "raters": {"Max": ["T0003.xyz", ...]},   # disjoint, in file order
     "split": {"seed": ..., "test_fraction": 0.2, "test": [...]}}   # optional

* :func:`assign` draws the share: the files are shuffled with the seed, the
  first ``overlap`` become the shared ones, and every new reviewer takes
  ``block`` files (by default the free files split evenly) by
  :func:`draw_block`: the least rated first.  Every ``review_<name>.json``
  and every active reservation that can be seen (``ratings/<set>/`` and its
  ``claims/``) counts, so structures nobody has rated or reserved go out
  first.  A ``redundancy`` share of the set (default 15 %, in
  ``assignment.json``), drawn once with the seed and spread over the whole
  set, is meant for two or three reviewers: every block fills that share of
  its places with such structures that already have a reviewer, so every
  reviewer overlaps with many others (an incomplete block design).  Second
  ratings of the rest come only when everything has one.  Reviewers added
  later change nothing given out before; files a reviewer has already rated
  stay theirs.  An optional split marks a fixed fraction of the files, drawn
  once with its own seed, as *test*: a pre-registered hold-out whose ratings
  are never used to tune anything.
* :func:`claim` reserves a reviewer's block at once, as a small pull request
  that adds only ``ratings/<set>/claims/<name>.json``; a reservation runs
  out after ``CLAIM_DAYS`` (28) days from its date, and the structures its
  reviewer has not rated are then given out again.
* :func:`rater_order` is the order the Review tab shows a reviewer: the
  shared files mixed among as many of their own, then the rest of their own
  in file order, so the shared ones are rated early but do not stand out.
* :func:`self_share` lets anyone help without being added first: a name not
  in the assignment draws a block of its own (``self_block``, default 500)
  by :func:`draw_block`, from the seed, the name and the ratings and
  reservations it can see, kept in its own ``review_<name>.json`` rather
  than in ``assignment.json``.
* :func:`agreement` measures how well the reviewers agree on every frame two
  or more of them rated (the shared ones, and overlapping blocks): Cohen's
  kappa for two, Fleiss' kappa for more, for the verdict (pass/block) and for
  every category, with the share of frames on which all agree and the number
  of frames.  :func:`reliability` adds Krippendorff's alpha with a bootstrap
  interval and the Dawid-Skene consensus (:mod:`delfin.review_stats`), and
  :func:`status` counts how often every structure is rated so far.
  :func:`write_merged_csv` writes every rating of every reviewer in one CSV.
* :func:`unpack` writes the structures of a published set
  (``structures.xyz.gz``, one frame per structure) as single files, and
  :func:`submit` brings one reviewer's ``review_<name>.json`` back into the
  repository as a pull request that adds only that one file.  :func:`new_set`
  makes the public part of a new set from a multi-frame XYZ file.

``delfin review assign|merge|unpack|submit|claim|status|new-set`` are the
command line to all of it.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import random
import re
import secrets
import shutil
import subprocess
import tempfile
from collections.abc import Iterable, Sequence
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from delfin import review as rv
from delfin.agent import safe_git

ASSIGNMENT_NAME = 'assignment.json'
SET_ARCHIVE_NAME = 'structures.xyz.gz'
DEFAULT_OVERLAP = 50
DEFAULT_BLOCK = 500
SCHEMA_VERSION = 1
UPSTREAM = 'ComPlat/DELFIN'
RATINGS_DIR = 'ratings'
CLAIMS_DIR = 'claims'
OWNERS_DIR = 'owners'
CLAIM_DAYS = 28
DEFAULT_REDUNDANCY = 0.15
DEFAULT_TEST_FRACTION = 0.2
# Written by unpack next to the structures: where the set and its ratings are.
SOURCE_NAME = '.review_set.json'
VERDICT_QUESTION = 'ok/defect'
COMMANDS = ('assign', 'merge', 'unpack', 'submit', 'claim', 'status', 'new-set')


# ---------------------------------------------------------------------------
# The assignment
# ---------------------------------------------------------------------------

def structure_files(folder) -> list[str]:
    """The ``*.xyz`` file names of *folder*, in name order (the order of the Review tab);
    for a published set (only ``structures.xyz.gz``) the names it unpacks to."""
    folder = Path(folder).expanduser()
    names = sorted(p.name for p in folder.glob('*.xyz') if p.is_file())
    if not names and (folder / SET_ARCHIVE_NAME).is_file():
        names = sorted(f'{ident}.xyz' for ident, _ in archive_frames(folder))
    return names


def find_assignment(folder) -> Path | None:
    """``assignment.json`` in *folder* or the folder above, else None."""
    folder = Path(folder).expanduser()
    for where in (folder, folder.parent):
        path = where / ASSIGNMENT_NAME
        if path.is_file():
            return path
    return None


def load_assignment(path) -> dict[str, Any]:
    data = json.loads(Path(path).read_text(encoding='utf-8'))
    if not isinstance(data, dict) or not isinstance(data.get('raters'), dict):
        raise ValueError(f'{path} is not an assignment file')
    data.setdefault('shared', [])
    return data


def find_rater(assignment: dict[str, Any], name: str) -> str | None:
    """The reviewer's name as the assignment spells it (exact, else ignoring case)."""
    want = (name or '').strip()
    if not want:
        return None
    raters = assignment.get('raters') or {}
    if want in raters:
        return want
    return next((r for r in raters if r.casefold() == want.casefold()), None)


def _file_stem(name: str) -> str:
    return rv.default_review_path(Path('.'), name).name


def _shuffled(files: Sequence[str], seed) -> list[str]:
    order = sorted(files)
    random.Random(seed).shuffle(order)
    return order


# ---------------------------------------------------------------------------
# Who has what: ratings, reservations, how often each structure is covered
# ---------------------------------------------------------------------------

def _key(name: str) -> str:
    """One reviewer, however the name is written: the review file's stem, case folded."""
    return _file_stem(name)[len('review_'):-len('.json')].casefold()


def set_version_of(frames_text: bytes) -> str:
    """``sha256:<16 hex>`` of a set's frames (the decompressed ``structures.xyz.gz``)."""
    return 'sha256:' + hashlib.sha256(frames_text).hexdigest()[:16]


def set_version(set_dir) -> str | None:
    """The version of the published set in *set_dir*, None without ``structures.xyz.gz``."""
    archive = Path(set_dir).expanduser() / SET_ARCHIVE_NAME
    if not archive.is_file():
        return None
    with gzip.open(archive, 'rb') as handle:
        return set_version_of(handle.read())


def ratings_dirs(folder, extra: Sequence = ()) -> list[Path]:
    """The folders whose ``review_<name>.json`` and ``claims/`` count for *folder*:
    the folder itself, *extra*, the ratings folder an unpacked set names in
    ``.review_set.json``, and ``ratings/<set>`` of the clone for a folder
    ``review_sets/<set>``."""
    folder = Path(folder).expanduser()
    found = [folder, *(Path(e).expanduser() for e in extra)]
    marker = folder / SOURCE_NAME
    if marker.is_file():
        try:
            where = json.loads(marker.read_text(encoding='utf-8')).get('ratings')
        except (OSError, ValueError, AttributeError):
            where = None
        if where:
            found.append(Path(where))
    if folder.parent.name == 'review_sets':
        found.append(folder.parent.parent / RATINGS_DIR / folder.name)
    out: list[Path] = []
    seen: set[str] = set()
    for d in found:
        real = str(d.resolve())
        if real not in seen and d.is_dir():
            seen.add(real)
            out.append(d)
    return out


def collect_claims(*folders) -> dict[str, dict[str, Any]]:
    """``{name: claim}`` from every ``claims/<name>.json`` under *folders*."""
    out: dict[str, dict[str, Any]] = {}
    for folder in folders:
        for path in sorted((Path(folder).expanduser() / CLAIMS_DIR).glob('*.json')):
            try:
                data = json.loads(path.read_text(encoding='utf-8'))
            except (OSError, ValueError):
                continue
            if isinstance(data, dict) and isinstance(data.get('files'), list):
                out.setdefault(path.stem, data)
    return out


def _parse_time(text) -> datetime | None:
    try:
        when = datetime.fromisoformat(str(text))
    except (TypeError, ValueError):
        return None
    return when if when.tzinfo else when.replace(tzinfo=timezone.utc)


def reservations(assignment: dict[str, Any], reviews: dict[str, dict[str, Any]],
                 claims: dict[str, dict[str, Any]], *, now: datetime | None = None,
                 days: int = CLAIM_DAYS) -> dict[str, dict[str, Any]]:
    """Every reviewer's reserved block, ``{key: {name, kind, files, date, expires, active}}``.

    A block named in ``assignment.json`` is reserved for good (``assigned``).
    A block drawn for oneself (kept in ``review_<name>.json``, ``self``) and
    a ``claims/<name>.json`` (``claim``) are reserved from their date for
    *days* (a claim may give its own ``days``); one reviewer's blocks are
    one reservation, dated by the latest of them.  Once it has run out, the
    structures its reviewer has not rated count as free again."""
    now = now or datetime.now(timezone.utc)
    out: dict[str, dict[str, Any]] = {}

    def add(name, files, kind, date, length):
        k = _key(name)
        entry = out.setdefault(k, {'name': name, 'kind': kind, 'files': set(),
                                   'date': None, 'expires': None, 'active': False,
                                   'assigned': False})
        entry['files'] |= {f for f in files or [] if isinstance(f, str)}
        if kind == 'assigned':
            entry.update(kind='assigned', assigned=True)
        elif not entry['assigned']:
            entry['kind'] = 'claim' if 'claim' in (kind, entry['kind']) else kind
        when = _parse_time(date)
        if when is not None and (entry['date'] is None or when > entry['date']):
            entry['date'] = when
            entry['expires'] = when + timedelta(days=int(length or days))

    for name, files in (assignment.get('raters') or {}).items():
        add(name, files, 'assigned', None, None)
    for name, review in reviews.items():
        own = review.get('assignment') or {}
        if isinstance(own, dict) and own.get('files'):
            add(own.get('name') or name, own['files'], 'self', own.get('drawn'), None)
    for name, data in claims.items():
        add(data.get('name') or name, data.get('files'), 'claim', data.get('date'),
            data.get('days'))
    for entry in out.values():
        entry['active'] = entry['assigned'] or (entry['expires'] is not None
                                                and now < entry['expires'])
    return out


def _rated_by(reviews: dict[str, dict[str, Any]]) -> dict[str, set[str]]:
    """``{file: {reviewer key, ...}}`` of every current rating."""
    who: dict[str, set[str]] = {}
    for name, data in reviews.items():
        for f, _i in current_ratings(data):
            who.setdefault(f, set()).add(_key(name))
    return who


def coverage(files: Sequence[str], reviews: dict[str, dict[str, Any]],
             reserved: dict[str, dict[str, Any]], *, exclude: str | None = None) -> dict[str, int]:
    """How many reviewers (other than *exclude*) have rated or actively reserved each file."""
    who = {f: set(v) for f, v in _rated_by(reviews).items()}
    for k, entry in reserved.items():
        if entry['active']:
            for f in entry['files']:
                who.setdefault(f, set()).add(k)
    return {f: len(who.get(f, set()) - {exclude}) for f in files}


def redundancy_of(assignment: dict[str, Any]) -> float:
    try:
        return max(0.0, min(1.0, float(assignment.get('redundancy', DEFAULT_REDUNDANCY))))
    except (TypeError, ValueError):
        return DEFAULT_REDUNDANCY


def raters_wanted(assignment: dict[str, Any], files: Sequence[str]) -> dict[str, int]:
    """How many reviewers each structure is meant to get: one, or for the
    ``redundancy`` share two (two thirds of it) or three (one third), drawn
    once with the seed and spread over the whole set.  The shared structures
    are rated by everyone anyway and are left out."""
    shared = set(assignment.get('shared', []))
    pool = [f for f in files if f not in shared]
    n = round(redundancy_of(assignment) * len(pool))
    out = dict.fromkeys(pool, 1)
    for k, f in enumerate(_shuffled(pool, f'{assignment.get("seed", 0)}:redundant')[:n]):
        out[f] = 3 if k % 3 == 2 else 2
    return out


def draw_block(order: Sequence[str], size: int, cover: dict[str, int],
               wanted: dict[str, int], redundancy: float) -> list[str]:
    """*size* files from *order* (shuffled, without the shared ones and the
    drawer's own), the least covered first.

    First the redundant share of the places (``round(redundancy * size)``):
    structures meant for two or three reviewers that have fewer, those that
    already have a reviewer first, so blocks overlap.  Then structures nobody
    has rated or reserved.  Only when none is left, the rest, those below
    their target first, then by how often they are covered."""
    picked: list[str] = []
    seen: set[str] = set()

    def take(candidates, limit):
        for f in candidates:
            if len(picked) >= limit:
                return
            if f not in seen:
                seen.add(f)
                picked.append(f)

    size = max(0, int(size))
    under = [f for f in order if wanted.get(f, 1) > 1 and cover.get(f, 0) < wanted[f]]
    take([f for f in under if cover.get(f, 0) >= 1] + [f for f in under if not cover.get(f, 0)],
         min(size, round(redundancy * size)))
    take([f for f in order if not cover.get(f, 0)], size)
    take(sorted(order, key=lambda f: (cover.get(f, 0) >= wanted.get(f, 1), cover.get(f, 0))),
         size)
    return picked


def _coverage_mark(cover: dict[str, int]) -> str:
    """A short mark of the coverage a block was drawn from (to reproduce the draw)."""
    text = '\n'.join(f'{f}\t{n}' for f, n in sorted(cover.items()) if n)
    return 'sha256:' + hashlib.sha256(text.encode()).hexdigest()[:16]


def identical_groups(folder, groups: dict[str, str] | None = None) -> tuple[dict[str, str], int]:
    """``({file: group}, joined)``: the *groups* given (file -> group), joined
    with the structures that are geometrically the same
    (:func:`delfin.review_stats.find_repeats`: shown twice, or built alike);
    *joined* counts the pairs found.  A group is named by the smallest name
    among its members' own groups, so without identical structures the
    groups are those given."""
    from delfin import review_stats as rs
    files = structure_files(folder)
    label = {f: str((groups or {}).get(f, f)) for f in files}
    parent = {f: f for f in files}

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def join(a, b):
        ra, rb = root(a), root(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    first: dict[str, str] = {}
    for f in files:
        if label[f] in first:
            join(first[label[f]], f)
        else:
            first[label[f]] = f
    pairs = [(a, b) for a, b in rs.find_repeats(set_frames(folder)) if a in parent and b in parent]
    for a, b in pairs:
        join(a, b)
    name: dict[str, str] = {}
    for f in files:
        r = root(f)
        name[r] = min(name.get(r, label[f]), label[f])
    return {f: name[root(f)] for f in files}, len(pairs)


def regroup_split(folder, data: dict[str, Any], *, groups: dict[str, str] | None = None,
                  reason: str = 'structures that are geometrically the same kept in one part',
                  date: str | None = None) -> list[str]:
    """Move every group (:func:`identical_groups`) that straddles the test split
    wholly into one part, changing as little as possible: into the part that
    holds most of it, on a tie into the one that keeps the test part nearest
    its size (then by the seed).  Records the change with its date, reason
    and counts (not the files) in ``split.changes``.  For a correction made
    before any rating; returns messages."""
    split = data.get('split')
    if not split:
        raise ValueError('the assignment has no test split')
    files = structure_files(folder)
    place = {f: i for i, f in enumerate(files)}
    group_of, joined = identical_groups(folder, groups)
    members: dict[str, list[str]] = {}
    for f in files:
        members.setdefault(group_of[f], []).append(f)
    test = set(split.get('test', []))
    target = round(float(split['test_fraction']) * len(files))
    crossing = [m for m in members.values() if 0 < sum(f in test for f in m) < len(m)]
    crossing.sort(key=lambda m: m[0])
    random.Random(f'{split.get("seed", "")}:regroup').shuffle(crossing)
    coin = random.Random(f'{split.get("seed", "")}:regroup:tie')
    to_test = to_cal = 0
    for m in crossing:
        k_t = sum(f in test for f in m)
        k_c = len(m) - k_t
        if k_t != k_c:
            into_test = k_t > k_c
        elif len(test) != target:
            into_test = len(test) < target
        else:
            into_test = coin.random() < 0.5
        if into_test:
            to_test += k_c
            test |= set(m)
        else:
            to_cal += k_t
            test -= set(m)
    msgs = [f'{len(crossing)} group(s) of identical structures straddled the split; '
            f'{to_test} structure(s) moved to the test part, {to_cal} to calibration; '
            f'test part now {len(test)} of {len(files)}']
    if crossing:
        split['test'] = sorted(test, key=place.get)
        split['grouped'] = True
        split['identical_kept_together'] = joined
        split.setdefault('changes', []).append({
            'date': date or datetime.now(timezone.utc).date().isoformat(), 'reason': reason,
            'groups': len(crossing), 'moved_to_test': to_test, 'moved_to_calibration': to_cal,
            'test_size': len(test)})
    return msgs


def assign(folder, raters: Iterable[str], *, overlap: int | None = None,
           block: int | None = None, seed: int | None = None,
           test_fraction: float | None = None,
           existing: dict[str, Any] | None = None,
           set_name: str | None = None,
           groups: dict[str, str] | None = None,
           ratings: Sequence = (),
           self_block: int | None = None,
           redundancy: float | None = None,
           now: datetime | None = None,
           regroup: bool = False) -> tuple[dict[str, Any], list[str]]:
    """Share the structures of *folder* out to *raters*; returns (assignment, messages).

    With *existing*, the reviewers there and their files stay as they are
    and only the new names get files, the least rated and reserved first
    (:func:`draw_block`; *ratings* = further folders with review files and
    ``claims/``).  *groups* (file -> group) keeps the files of one group in
    one part of the test split; structures that are geometrically the same
    are kept together by themselves (:func:`identical_groups`); the groups
    are not stored.  *redundancy* is the share of the set meant for two or
    three reviewers (default 0.15).  *regroup* (an existing assignment, before
    any rating) moves groups that straddle its split into one part
    (:func:`regroup_split`).
    """
    folder = Path(folder).expanduser()
    files = structure_files(folder)
    if not files:
        raise ValueError(f'no .xyz files in {folder}')
    place = {f: i for i, f in enumerate(files)}
    msgs: list[str] = []
    names = [n.strip() for n in raters if n and n.strip()]
    if existing:
        data = json.loads(json.dumps(existing))
        if int(data.get('files_total', len(files))) != len(files):
            raise ValueError(f'the folder has {len(files)} structures, the assignment was '
                             f'made for {data.get("files_total")}; make a new one')
        given = set(data.get('shared', [])) | {f for v in data['raters'].values() for f in v}
        if given - set(files):
            raise ValueError(f'{len(given - set(files))} assigned file(s) are no longer in {folder}')
        if seed is not None and int(seed) != int(data['seed']):
            raise ValueError(f'the assignment has seed {data["seed"]}; it cannot change')
        if overlap is not None and int(overlap) != len(data.get('shared', [])):
            raise ValueError(f'the assignment shares {len(data.get("shared", []))} structures; '
                             'that cannot change once given out')
        if test_fraction is not None and 'split' not in data:
            raise ValueError('a test split can only be drawn with a new assignment')
        if regroup:
            msgs += regroup_split(folder, data, groups=groups)
    else:
        seed = int(seed) if seed is not None else secrets.randbelow(2**31 - 1) + 1
        overlap = DEFAULT_OVERLAP if overlap is None else int(overlap)
        if overlap < 0:
            raise ValueError('overlap must be 0 or more')
        perm = _shuffled(files, seed)
        data = {'version': SCHEMA_VERSION, 'seed': seed, 'files_total': len(files),
                'shared': sorted(perm[:overlap], key=place.get), 'raters': {},
                'created': _now()}
        if test_fraction:
            frac = float(test_fraction)
            if not 0 < frac < 1:
                raise ValueError('test fraction must lie between 0 and 1')
            split_seed = f'{seed}:split'
            n_test = round(frac * len(files))
            # Structures that are geometrically the same (shown twice, or built
            # alike) stay in one part, with the groups given.
            group_of, joined = identical_groups(folder, groups)
            members: dict[str, list[str]] = {}
            for f in files:
                members.setdefault(group_of[f], []).append(f)
            test: list[str] = []
            for g in _shuffled(list(members), split_seed):
                if len(test) >= n_test:
                    break
                test += members[g]
            data['split'] = {'seed': split_seed, 'test_fraction': frac,
                             'grouped': bool(groups) or bool(joined),
                             'identical_kept_together': joined,
                             'test': sorted(test, key=place.get)}
    if set_name:
        data['set'] = set_name
    if self_block is not None:
        data['self_block'] = int(self_block)
    else:
        data.setdefault('self_block', DEFAULT_BLOCK)
    if redundancy is not None:
        if not 0 <= float(redundancy) <= 1:
            raise ValueError('redundancy must lie between 0 and 1')
        if existing and float(redundancy) != redundancy_of(data):
            msgs.append(f'redundancy {redundancy_of(data):g} -> {float(redundancy):g} '
                        '(only blocks drawn from now on follow it)')
        data['redundancy'] = float(redundancy)
    else:
        data.setdefault('redundancy', DEFAULT_REDUNDANCY)
    version = set_version(folder)
    if version and data.get('set_version') not in (None, version):
        raise ValueError(f'the structures are not those of the assignment (set version '
                         f'{data.get("set_version")}, now {version}); make a new set')
    if version:
        data['set_version'] = version
    raters_now: dict[str, list[str]] = data['raters']
    stems = {_key(n): n for n in raters_now}
    new: list[str] = []
    for name in names:
        if find_rater(data, name):
            msgs.append(f'{name}: already in the assignment, kept as it is')
            continue
        stem = _key(name)
        if stem in stems:
            raise ValueError(f'"{name}" and "{stems[stem]}" would write the same review file')
        stems[stem] = name
        new.append(name)
    shared = set(data['shared'])
    # Everything that can be seen counts: ratings, blocks reviewers drew for
    # themselves (kept in their review files) and active claims.
    dirs = ratings_dirs(folder, ratings)
    reviews = collect_reviews(*dirs)
    reserved = reservations(data, reviews, collect_claims(*dirs), now=now)
    wanted = raters_wanted(data, files)
    order = [f for f in _shuffled(files, data['seed']) if f not in shared]
    if new:
        free = coverage(files, reviews, reserved)
        per = (int(block) if block is not None
               else sum(1 for f in order if not free[f]) // len(new))
        if per < 0:
            raise ValueError('block must be 0 or more')
        rated = _rated_by(reviews)
        for name in new:
            k = _key(name)
            cover = coverage(files, reviews, reserved, exclude=k)
            own = [f for f in order if k in rated.get(f, ())]
            fresh = draw_block([f for f in order if k not in rated.get(f, ())],
                               per - len(own), cover, wanted, redundancy_of(data))
            raters_now[name] = sorted(own + fresh, key=place.get)
            reserved[k] = {'name': name, 'kind': 'assigned', 'files': set(raters_now[name]),
                           'date': None, 'expires': None, 'active': True, 'assigned': True}
            twice = sum(1 for f in fresh if cover[f])
            note = (f' ({len(own)} already rated)' if own else '') + \
                   (f', {twice} also with another reviewer' if twice else '')
            msgs.append(f'{name}: {len(raters_now[name])} own + {len(shared)} shared{note}')
    data['updated'] = _now()
    cover = coverage(files, reviews, reserved)
    left = sum(1 for f in order if not cover[f])
    msgs.append(f'{len(files)} structures: {len(shared)} shared, '
                f'{len(order) - left} given out or rated, {left} free')
    return data, msgs


def _mixed_order(shared: Sequence[str], own: Sequence[str], seed) -> list[str]:
    own = list(own)
    mine = set(own)
    shared = [f for f in shared if f not in mine]
    head = shared + own[:len(shared)]
    random.Random(seed).shuffle(head)
    return head + own[len(shared):]


def rater_order(assignment: dict[str, Any], name: str) -> list[str]:
    """The files the reviewer *name* sees, in the order shown.

    The shared files come mixed among as many of the reviewer's own (shuffled
    with the seed and the name), then the rest of their own in file order: the
    shared ones are rated early, so the agreement can be measured from the
    start, but they do not stand out as a block.
    """
    rater = find_rater(assignment, name)
    if rater is None:
        raise KeyError(name)
    return _mixed_order(assignment.get('shared', []), assignment['raters'][rater],
                        f'{assignment.get("seed", 0)}:{rater}')


def self_share(folder, assignment: dict[str, Any], name: str,
               size: int | None = None, *, ratings: Sequence = (),
               now: datetime | None = None) -> dict[str, Any]:
    """The block a reviewer not in the assignment draws for themselves.

    Anyone can help without being added first: a new name gets ``self_block``
    (default 500) structures by :func:`draw_block`, the least rated and
    reserved first, with the redundant share overlapping other reviewers.
    Ratings and reservations count where they can be seen
    (:func:`ratings_dirs`: the folder, *ratings*, and the clone's
    ``ratings/<set>`` named by ``unpack``), so a fresh ``git pull`` before
    unpacking makes the block fit the newest state.  The draw depends only
    on the seed, the name and that state, so the same name always gets the
    same block from the same state.  The block is kept in the reviewer's own
    ``review_<name>.json`` (key ``assignment``), not in ``assignment.json``,
    so helpers never edit a common file; once drawn it is read from there (or
    from the reviewer's submitted file or claim) and stays, whatever the
    assignment does later.
    """
    name = (name or '').strip()
    if not name:
        raise ValueError('a name is needed')
    folder = Path(folder).expanduser()
    try:
        stored = rv.load_review(rv.default_review_path(folder, name)).get('assignment')
    except (OSError, ValueError):
        stored = None
    if isinstance(stored, dict) and stored.get('files'):
        return stored
    dirs = ratings_dirs(folder, ratings)
    reviews = collect_reviews(*dirs)
    claims = collect_claims(*dirs)
    k = _key(name)
    files = structure_files(folder)
    place = {f: i for i, f in enumerate(files)}
    shared = [f for f in assignment.get('shared', []) if f in place]
    # Drawn before, on another machine: the block sent back or claimed.
    for found in ([v.get('assignment') for n, v in reviews.items() if _key(n) == k]
                  + [c for n, c in claims.items() if _key(c.get('name') or n) == k]):
        if isinstance(found, dict) and found.get('files'):
            return {'self': True, 'name': name, 'seed': found.get('seed', ''),
                    'size': found.get('size', len(found['files'])),
                    'files': sorted((f for f in found['files'] if f in place), key=place.get),
                    'shared': shared, 'drawn': found.get('drawn') or found.get('date') or _now()}
    reserved = reservations(assignment, reviews, claims, now=now)
    cover = coverage(files, reviews, reserved, exclude=k)
    done = {f for f, who in _rated_by(reviews).items() if k in who and f in place}
    size = int(size if size is not None else assignment.get('self_block', DEFAULT_BLOCK))
    seed = f'{assignment.get("seed", 0)}:self:{name}'
    order = [f for f in _shuffled(files, seed) if f not in set(shared)]
    own = [f for f in order if f in done]
    block = own + draw_block([f for f in order if f not in done], size - len(own), cover,
                             raters_wanted(assignment, files), redundancy_of(assignment))
    return {'self': True, 'name': name, 'seed': seed, 'size': size,
            'files': sorted(block, key=place.get), 'shared': shared, 'drawn': _now(),
            'drawn_from': {'reviews': len([n for n in reviews if _key(n) != k]),
                           'claims': len([n for n in claims if _key(n) != k]),
                           'redundancy': redundancy_of(assignment),
                           'coverage': _coverage_mark(cover)}}


def share_of(folder, assignment: dict[str, Any], name: str):
    """(name, files in the order shown, self-drawn block or None) of a reviewer."""
    rater = find_rater(assignment, name)
    if rater is not None:
        return rater, rater_order(assignment, rater), None
    own = self_share(folder, assignment, name)
    order = _mixed_order(own.get('shared', assignment.get('shared', [])), own['files'],
                         f'{assignment.get("seed", 0)}:{own["name"]}')
    return own['name'], order, own


def write_assignment(path, data: dict[str, Any]) -> Path:
    path = Path(path)
    rv.write_json_atomic(path, data)
    return path


# ---------------------------------------------------------------------------
# Bringing the ratings together
# ---------------------------------------------------------------------------

def collect_reviews(*folders) -> dict[str, dict[str, Any]]:
    """``{name: review}`` from every ``review_<name>.json`` in *folders*."""
    out: dict[str, dict[str, Any]] = {}
    for folder in folders:
        for path in sorted(Path(folder).expanduser().glob('review_*.json')):
            name = path.stem[len('review_'):]
            try:
                data = rv.load_review(path)
            except (OSError, ValueError):
                continue
            if data and name not in out:
                out[name] = data
    return out


def current_ratings(data: dict[str, Any]) -> dict[tuple[str, int], dict[str, Any]]:
    return {(f, i): rec for f, i, rec, stale in rv.iter_records(data) if not stale}


MERGED_FIELDS = ('rater',) + rv.CSV_FIELDS + ('stale', 'shared', 'share', 'split')


def merged_rows(reviews: dict[str, dict[str, Any]],
                assignment: dict[str, Any] | None = None, split: str = 'all'):
    """Every row of every reviewer (:func:`delfin.review.csv_rows`), with the
    reviewer, whether the file is shared, whose share it is (``shared``,
    ``own`` -- assigned or self-drawn -- or ``extra``), and its part of the split."""
    shared = set((assignment or {}).get('shared', []))
    test = set(((assignment or {}).get('split') or {}).get('test', []))
    has_split = bool((assignment or {}).get('split'))
    for name in sorted(reviews):
        rater = find_rater(assignment or {}, name)
        own = set((assignment or {}).get('raters', {}).get(rater, []) if rater else
                  (reviews[name].get('assignment') or {}).get('files') or [])
        for row in rv.csv_rows(reviews[name]):
            part = ('test' if row['file'] in test else 'calibration') if has_split else ''
            if split != 'all' and part != split:
                continue
            share = ('shared' if row['file'] in shared else
                     'own' if row['file'] in own else 'extra')
            yield {'rater': name, **row, 'shared': 'yes' if row['file'] in shared else '',
                   'share': share, 'split': part}


def write_merged_csv(reviews, path, assignment=None, split: str = 'all') -> int:
    import csv
    n = 0
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=MERGED_FIELDS)
        writer.writeheader()
        for row in merged_rows(reviews, assignment, split):
            writer.writerow(row)
            n += 1
    return n


def cohen_kappa(a: Sequence[bool], b: Sequence[bool]) -> float | None:
    """Cohen's kappa of two raters on a yes/no question; None when undefined."""
    n = len(a)
    if not n or n != len(b):
        return None
    po = sum(x == y for x, y in zip(a, b)) / n
    pa, pb = sum(a) / n, sum(b) / n
    pe = pa * pb + (1 - pa) * (1 - pb)
    return None if pe >= 1 else (po - pe) / (1 - pe)


def fleiss_kappa(counts: Sequence[int], raters) -> float | None:
    """Fleiss' kappa on a yes/no question: *counts* = yes-votes per item out of
    *raters* (one number for all items, or one per item when the number of
    raters varies; every item needs two or more)."""
    N = len(counts)
    m = [int(raters)] * N if isinstance(raters, int) else [int(r) for r in raters]
    if not N or len(m) != N or min(m) < 2:
        return None
    p_items = [(k * k + (n - k) * (n - k) - n) / (n * (n - 1)) for k, n in zip(counts, m)]
    p_bar = sum(p_items) / N
    p = sum(counts) / sum(m)
    pe = p * p + (1 - p) * (1 - p)
    return None if pe >= 1 else (p_bar - pe) / (1 - pe)


def _measure(votes: list[list[bool]], cohen: bool) -> dict[str, Any]:
    """Agreement of one yes/no question; *votes* = one list per item, one vote per
    rater who rated it (in a fixed rater order for Cohen's two)."""
    n = len(votes)
    if cohen:
        kappa = cohen_kappa([v[0] for v in votes], [v[1] for v in votes])
    else:
        kappa = fleiss_kappa([sum(v) for v in votes], [len(v) for v in votes])
    return {'kappa': kappa, 'agreement': (sum(len(set(v)) == 1 for v in votes) / n) if n else None,
            'n': n, 'yes': sum(sum(v) for v in votes)}


def agreement(reviews: dict[str, dict[str, Any]],
              items: Iterable[str] | None = None) -> dict[str, Any] | None:
    """How well the reviewers agree, on every frame two or more of them rated
    (only frames of the files *items* when given): the shared structures, and
    any a reviewer rated besides, such as overlapping self-drawn blocks.
    Cohen's kappa when all these frames are rated by the same two reviewers,
    else Fleiss' kappa (the number of reviewers may differ per frame).  None
    when no frame was rated twice."""
    rated = {name: current_ratings(data) for name, data in reviews.items()}
    wanted = None if items is None else set(items)
    by_key: dict[tuple[str, int], list[str]] = {}
    for name in sorted(rated):
        for k in rated[name]:
            if wanted is None or k[0] in wanted:
                by_key.setdefault(k, []).append(name)
    keys = sorted(k for k, who in by_key.items() if len(who) >= 2)
    if not keys:
        return None
    names = sorted({n for k in keys for n in by_key[k]})
    cohen = len(names) == 2
    cats: list[str] = []
    for n in names:
        for c in reviews[n].get('categories') or []:
            if c not in cats:
                cats.append(c)
    for k in keys:
        for n in by_key[k]:
            for c in rated[n][k].get('categories') or []:
                if c not in cats:
                    cats.append(c)
    out = {'raters': names, 'method': 'Cohen' if cohen else 'Fleiss',
           'frames': len(keys), 'votes': sum(len(by_key[k]) for k in keys),
           'verdict': _measure([[rated[n][k]['verdict'] == 'block' for n in by_key[k]]
                                for k in keys], cohen),
           'unsure': sum(rv.is_unsure(rated[n][k].get('note', ''))
                         for k in keys for n in by_key[k]),
           'categories': {}}
    for c in cats:
        out['categories'][c] = _measure(
            [[c in (rated[n][k].get('categories') or []) for n in by_key[k]] for k in keys],
            cohen)
    return out


def format_agreement(result: dict[str, Any] | None) -> str:
    if not result:
        return 'Agreement: no frame has been rated by two reviewers yet.'

    def fmt(m):
        kappa = 'n/a' if m['kappa'] is None else f'{m["kappa"]:.3f}'
        agree = 'n/a' if m['agreement'] is None else f'{100 * m["agreement"]:.1f} %'
        return f'{kappa:>7}  {agree:>8}  {m["n"]:5d}  {m["yes"]:5d}'

    width = max([len('ok/defect')] + [len(c) for c in result['categories']])
    lines = [f'Agreement of {len(result["raters"])} reviewers ({result["method"]} kappa) '
             f'on {result["frames"]} frames rated by two or more ({result["votes"]} votes): '
             f'{", ".join(result["raters"])}',
             f'{"question".ljust(width)}    kappa  all agree      n    yes',
             f'{"ok/defect".ljust(width)}  {fmt(result["verdict"])}']
    for c, m in result['categories'].items():
        lines.append(f'{c.ljust(width)}  {fmt(m)}')
    lines.append('"yes" counts the votes for block / for the category over all reviewers; '
                 f'{result["unsure"]} of the votes are marked unsure.')
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Reliability: Krippendorff's alpha, Dawid-Skene, self-consistency, controls
# ---------------------------------------------------------------------------

def split_items(assignment: dict[str, Any] | None, split: str,
                files: Sequence[str]) -> set[str] | None:
    """The files of one part of the test split (None for ``all``)."""
    if split == 'all':
        return None
    part = (assignment or {}).get('split')
    if not part:
        raise ValueError('the assignment has no test split; use --split all')
    test = set(part.get('test', []))
    return {f for f in files if (f in test) == (split == 'test')}


def question_votes(reviews: dict[str, dict[str, Any]], items: Iterable[str] | None = None,
                   categories: Sequence[str] = ()) -> dict[str, dict[tuple[str, int], dict[str, int]]]:
    """``{question: {(file, frame): {rater: 0/1}}}`` for ok/defect (1 = block) and
    every category (1 = marked), over the frames of *items* (all when None)."""
    wanted = None if items is None else set(items)
    rated = {name: current_ratings(data) for name, data in reviews.items()}
    cats = list(categories)
    for n in sorted(rated):
        for c in list(reviews[n].get('categories') or []) + [
                c for rec in rated[n].values() for c in rec.get('categories') or []]:
            if c not in cats:
                cats.append(c)
    out: dict[str, dict] = {VERDICT_QUESTION: {}, **{c: {} for c in cats}}
    for n in sorted(rated):
        for key, rec in rated[n].items():
            if wanted is not None and key[0] not in wanted:
                continue
            out[VERDICT_QUESTION].setdefault(key, {})[n] = int(rec['verdict'] == 'block')
            got = set(rec.get('categories') or [])
            for c in cats:
                out[c].setdefault(key, {})[n] = int(c in got)
    return out


def reliability(reviews: dict[str, dict[str, Any]], items: Iterable[str] | None = None, *,
                categories: Sequence[str] = (), n_boot: int | None = None,
                seed: int = 0) -> dict[str, Any]:
    """Krippendorff's alpha (nominal, bootstrap interval) and the Dawid-Skene fit
    for ok/defect and every category.

    Returns ``questions`` ({question: {alpha, items, multi, prevalence,
    converged}}), ``raters`` ({rater: {question: {n, sensitivity,
    specificity}}}) and ``consensus`` ({(file, frame): {'n': raters,
    question: probability of yes}}).  Only the frames of *items* enter, so a
    fit on the calibration part never sees the test part."""
    from delfin import review_stats as rs
    n_boot = rs.DEFAULT_BOOTSTRAP if n_boot is None else int(n_boot)
    out: dict[str, Any] = {'questions': {}, 'raters': {}, 'consensus': {}}
    for q, votes in question_votes(reviews, items, categories).items():
        alpha = rs.alpha_with_ci([list(v.values()) for v in votes.values()],
                                 n_boot=n_boot, seed=seed)
        ds = rs.dawid_skene(votes, iterations=100, tol=1e-6)
        out['questions'][q] = {'alpha': alpha, 'items': len(votes),
                               'multi': sum(len(v) >= 2 for v in votes.values()),
                               'prevalence': ds['prevalence'][1], 'converged': ds['converged']}
        for r, m in ds['confusion'].items():
            sens, spec = rs.sensitivity_specificity(m)
            out['raters'].setdefault(r, {})[q] = {'n': ds['n'][r], 'sensitivity': sens,
                                                  'specificity': spec}
        for key, p in ds['posterior'].items():
            out['consensus'].setdefault(key, {'n': len(votes[key])})[q] = p[1]
    return out


def _fmt_ci(a: dict[str, Any]) -> tuple[str, str]:
    if a['alpha'] is None:
        return 'n/a', ''
    ci = '' if a['low'] is None else f'[{a["low"]:.2f}, {a["high"]:.2f}]'
    return f'{a["alpha"]:.3f}', ci


def format_reliability(res: dict[str, Any]) -> str:
    qs = res['questions']
    if not qs or not qs[VERDICT_QUESTION]['items']:
        return 'Reliability: no ratings.'
    width = max(len(q) for q in qs)
    lines = ["Krippendorff's alpha (nominal; 95 % bootstrap interval over the structures) "
             'and the Dawid-Skene prevalence:',
             f'{"question".ljust(width)}    alpha  interval        rated  twice+  prevalence']
    for q, m in qs.items():
        alpha, ci = _fmt_ci(m['alpha'])
        lines.append(f'{q.ljust(width)}  {alpha:>7}  {ci:<14}  {m["items"]:5d}  {m["multi"]:6d}'
                     f'  {100 * m["prevalence"]:9.1f} %')
    raters = {r: v[VERDICT_QUESTION] for r, v in res['raters'].items() if VERDICT_QUESTION in v}
    if raters:
        w = max(len('reviewer'), *(len(r) for r in raters))
        lines += ['', 'Dawid-Skene, ok/defect (defect = yes), per reviewer:',
                  f'{"reviewer".ljust(w)}  rated  sensitivity  specificity']
        for r, m in sorted(raters.items()):
            lines.append(f'{r.ljust(w)}  {m["n"]:5d}  {m["sensitivity"]:11.3f}  '
                         f'{m["specificity"]:11.3f}')
    return '\n'.join(lines)


def write_consensus_csv(res: dict[str, Any], path, assignment: dict[str, Any] | None = None) -> int:
    """One row per frame: the raters, the probability of a defect, the consensus
    and the probability of every category (Dawid-Skene)."""
    import csv
    test = set(((assignment or {}).get('split') or {}).get('test', []))
    has_split = bool((assignment or {}).get('split'))
    cats = [q for q in res['questions'] if q != VERDICT_QUESTION]
    fields = ['file', 'frame', 'split', 'raters', 'p_block', 'consensus'] + [f'p_{c}' for c in cats]
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for (f, i), row in sorted(res['consensus'].items()):
            p = row[VERDICT_QUESTION]
            writer.writerow({'file': f, 'frame': i,
                             'split': ('test' if f in test else 'calibration') if has_split else '',
                             'raters': row['n'], 'p_block': f'{p:.4f}',
                             'consensus': 'block' if p >= 0.5 else 'pass',
                             **{f'p_{c}': f'{row[c]:.4f}' for c in cats if c in row}})
    return len(res['consensus'])


def write_raters_csv(res: dict[str, Any], path) -> int:
    """Sensitivity and specificity of every reviewer for every question (Dawid-Skene)."""
    import csv
    n = 0
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['rater', 'question', 'rated', 'sensitivity', 'specificity'])
        for r in sorted(res['raters']):
            for q, m in res['raters'][r].items():
                writer.writerow([r, q, m['n'], f'{m["sensitivity"]:.4f}',
                                 f'{m["specificity"]:.4f}'])
                n += 1
    return n


def verdicts_by_file(reviews: dict[str, dict[str, Any]],
                     items: Iterable[str] | None = None) -> dict[str, dict[str, bool]]:
    """``{rater: {file: block?}}`` of the first frame of every rated file."""
    wanted = None if items is None else set(items)
    out: dict[str, dict[str, bool]] = {}
    for name, data in reviews.items():
        for (f, i), rec in current_ratings(data).items():
            if i == 0 and (wanted is None or f in wanted):
                out.setdefault(name, {})[f] = rec['verdict'] == 'block'
    return out


def set_frames(folder):
    """(file, elements, coordinates) of every structure of *folder* (or its archive)."""
    import numpy as np
    folder = Path(folder).expanduser()
    if (folder / SET_ARCHIVE_NAME).is_file():
        blocks = [(f'{ident}.xyz', text) for ident, text in archive_frames(folder)]
    else:
        blocks = [(f, (folder / f).read_text(encoding='utf-8')) for f in structure_files(folder)]
    for fname, text in blocks:
        frames = rv.parse_frames(text, fname)
        if not frames:
            continue
        rows = [line.split() for line in frames[0].coords.splitlines()]
        try:
            xyz = np.array([[float(x) for x in r[1:4]] for r in rows])
        except (ValueError, IndexError):
            continue
        yield fname, [r[0] for r in rows], xyz


def load_repeats(path) -> list[tuple[str, str]]:
    """Pairs of the same structure shown twice, from a private file: JSON
    ``{file: group}`` (e.g. the groups of the test split) or ``[[a, b], ...]``."""
    data = json.loads(Path(path).expanduser().read_text(encoding='utf-8'))
    if isinstance(data, dict):
        members: dict[str, list[str]] = {}
        for f, g in data.items():
            members.setdefault(str(g), []).append(f)
        return sorted((a, b) for m in members.values() for i, a in enumerate(sorted(m))
                      for b in sorted(m)[i + 1:])
    return sorted((min(a, b), max(a, b)) for a, b in data)


def load_controls(path) -> set[str]:
    """The control structures named in a private key: JSON ``{file: true}`` or a
    list of files, or a CSV with a column ``file`` or ``id`` and a column
    ``control`` (true/1/yes) or ``source`` (a control when it starts with
    ``crystal``)."""
    import csv
    path = Path(path).expanduser()

    def name(x):
        x = str(x).strip()
        return x if x.endswith('.xyz') else f'{x}.xyz'

    if path.suffix.lower() == '.json':
        data = json.loads(path.read_text(encoding='utf-8'))
        if isinstance(data, dict):
            return {name(f) for f, v in data.items()
                    if v is True or str(v).lower() in ('1', 'true', 'yes')
                    or str(v).lower().startswith('crystal')}
        return {name(f) for f in data}
    out: set[str] = set()
    with open(path, newline='', encoding='utf-8') as handle:
        for row in csv.DictReader(handle):
            ident = row.get('file') or row.get('id')
            flag = row.get('control')
            src = row.get('source') or ''
            if ident and ((flag is not None and flag.strip().lower() in ('1', 'true', 'yes'))
                          or (flag is None and src.lower().startswith('crystal'))):
                out.add(name(ident))
    return out


def format_self_consistency(res: dict[str, Any], source: str) -> str:
    if not res['pairs_known']:
        return f'Self-consistency: no structure is shown twice ({source}).'
    if not res['pairs_rated']:
        return (f'Self-consistency: {res["pairs_known"]} structures are shown twice ({source}); '
                'no reviewer has rated both of a pair yet.')
    kappa = 'n/a' if res['kappa'] is None else f'{res["kappa"]:.3f}'
    lines = [f'Self-consistency ({source}): {res["pairs_rated"]} pairs rated twice by the same '
             f'reviewer, same verdict {100 * res["same"]:.1f} %, kappa {kappa}']
    for r, m in res['raters'].items():
        lines.append(f'  {r}: {m["same"]} of {m["pairs"]} the same')
    return '\n'.join(lines)


def format_controls(res: dict[str, Any]) -> str:
    lines = [f'Controls (private key): {res["controls"]} structures, {res["votes"]} votes, '
             f'{res["passed"]} passed'
             + (f' ({100 * res["passed"] / res["votes"]:.1f} %)' if res['votes'] else '')]
    if res['consensus_block'] is not None:
        lines.append(f'  Dawid-Skene probability of a defect, mean over '
                     f'{res["consensus_items"]} rated controls: {res["consensus_block"]:.3f}')
    for r, m in res['raters'].items():
        lines.append(f'  {r}: {m["passed"]} of {m["rated"]} passed')
    lines.append('  (printed only; never write this into the repository)')
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Status of a set
# ---------------------------------------------------------------------------

def status(folder, ratings: Sequence = (), *, now: datetime | None = None,
           n_boot: int | None = None, seed: int = 0, repeats: bool = True) -> dict[str, Any]:
    """How far a set is: how often each structure is rated, the reservations,
    the reviewers, the agreement, the self-consistency on the structures shown
    twice (found from the geometry, :func:`delfin.review_stats.find_repeats`; only
    counts are reported) and the progress of each part."""
    from delfin import review_stats as rs
    folder = Path(folder).expanduser()
    files = structure_files(folder)
    found = find_assignment(folder)
    assignment = load_assignment(found) if found else {'raters': {}, 'shared': []}
    dirs = ratings_dirs(folder, ratings)
    reviews = collect_reviews(*dirs)
    claims = collect_claims(*dirs)
    now = now or datetime.now(timezone.utc)
    reserved = reservations(assignment, reviews, claims, now=now)
    who = _rated_by(reviews)
    shared = [f for f in assignment.get('shared', []) if f in set(files)]
    rest = [f for f in files if f not in set(shared)]
    hist = {'0': 0, '1': 0, '2': 0, '3+': 0}
    for f in rest:
        n = len(who.get(f, ()))
        hist['3+' if n >= 3 else str(n)] += 1
    wanted = raters_wanted(assignment, files)
    planned = [f for f, n in wanted.items() if n > 1]
    test = set((assignment.get('split') or {}).get('test', []))
    parts: dict[str, list[int]] = {}
    for f in files:
        part = ('test' if f in test else 'calibration') if assignment.get('split') else 'all'
        row = parts.setdefault(part, [0, 0])
        row[1] += 1
        row[0] += bool(who.get(f))
    res_rows = []
    for k, e in sorted(reserved.items()):
        if e['kind'] == 'assigned' and not e['files']:
            continue
        res_rows.append({'name': e['name'], 'kind': e['kind'], 'files': len(e['files']),
                         'rated': sum(1 for f in e['files'] if k in who.get(f, ())),
                         'active': e['active'],
                         'expires': e['expires'].date().isoformat() if e['expires'] else ''})
    alpha = rs.alpha_with_ci(
        [list(v.values()) for v in question_votes(reviews)[VERDICT_QUESTION].values()],
        n_boot=rs.DEFAULT_BOOTSTRAP if n_boot is None else n_boot, seed=seed)
    pairs = rs.find_repeats(set_frames(folder)) if repeats else []
    return {'set': assignment.get('set') or folder.name, 'files': len(files),
            'shared': {'files': len(shared),
                       'raters': sorted(len(who.get(f, ())) for f in shared)},
            'hist': hist, 'redundant': {'planned': len(planned),
                                        'reached': sum(len(who.get(f, ())) >= wanted[f]
                                                       for f in planned)},
            'parts': parts, 'reservations': res_rows,
            'raters': {n: {'rated': len(current_ratings(d)), 'credit': d.get('credit') or ''}
                       for n, d in sorted(reviews.items())},
            'credits': sorted({str(d['credit']) for d in reviews.values() if d.get('credit')}),
            'alpha': alpha,
            'self': rs.self_consistency(verdicts_by_file(reviews), pairs) if repeats else None}


def format_status(st: dict[str, Any]) -> str:
    h = st['hist']
    rest = sum(h.values())
    sh = st['shared']
    lines = [f'Review set {st["set"]}: {st["files"]} structures, {sh["files"]} rated by '
             f'everyone, {rest} shared out',
             '',
             'rated by     0      1      2     3+',
             f'structures {h["0"]:5d}  {h["1"]:5d}  {h["2"]:5d}  {h["3+"]:5d}']
    if sh['files']:
        lines.append(f'the {sh["files"]} rated by everyone: {sh["raters"][0]} to '
                     f'{sh["raters"][-1]} reviewers each')
    red = st['redundant']
    if red['planned']:
        lines.append(f'meant for 2-3 reviewers: {red["planned"]}, reached {red["reached"]}')
    lines.append('')
    lines.append('part          rated  of  total')
    for part, (n, total) in st['parts'].items():
        lines.append(f'{part:<12} {n:6d}     {total:6d}  ({100 * n / total:.1f} %)')
    act = [r for r in st['reservations'] if r['active']]
    gone = [r for r in st['reservations'] if not r['active']]
    lines += ['', f'reservations: {len(act)} active, {len(gone)} expired']
    for r in st['reservations']:
        state = 'active' if r['active'] else 'expired'
        until = f' until {r["expires"]}' if r['expires'] else ''
        lines.append(f'  {r["name"]:<20} {r["kind"]:<8} {r["rated"]:5d} of {r["files"]:5d} rated'
                     f'  {state}{until}')
    lines += ['', f'reviewers: {len(st["raters"])}']
    for n, r in st['raters'].items():
        lines.append(f'  {n:<20} {r["rated"]:5d} rated')
    alpha, ci = _fmt_ci(st['alpha'])
    lines += ['', f"Krippendorff's alpha, ok/defect: {alpha} {ci} "
                  f'({st["alpha"]["items"]} structures rated twice or more)']
    if st.get('self') is not None:
        lines += ['', format_self_consistency(st['self'], 'found in the set')]
    if st['credits']:
        lines += ['', 'contributors: ' + ', '.join(st['credits'])]
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# A new set
# ---------------------------------------------------------------------------

README_TEMPLATE = """# Review set {name}

{n} three-dimensional structures to be rated by eye in the Review tab of the
DELFIN dashboard. TODO: say what the structures are and what the ratings
are for. How to take part:
[docs/REVIEW_TEAM_GUIDE.md](../../docs/REVIEW_TEAM_GUIDE.md).

* `structures.xyz.gz` -- all structures, one XYZ frame each; the comment line
  is only the structure's id. `delfin review unpack` writes them as single
  files.
* `kategorien.txt` -- the {n_cats} defect categories the Review tab offers.
* `assignment.json` -- the {overlap} structures everyone rates, the test part
  and the parameters of the draw (seed {seed}, block {block}, redundancy
  {redundancy:g}).

**Just pick your name.** Anyone can help without being added:
`delfin review unpack review_sets/{name} <folder> --name <your name>` draws a
block of {block} structures, the least rated first, plus the {overlap} everyone
rates; `delfin review claim <folder> --name <your name>` reserves it for
{days} days.

**What is in it.** TODO: the sources and their licences (no data of the
Cambridge Structural Database), how the structures were made unrecognisable
(centred, rotated, atoms sorted), and what is kept back until the ratings
are in.

**Calibration and test.** Before any rating, a fixed fraction
({test_fraction:g}) of the structures (`split.test` in `assignment.json`, seed
`{split_seed}`) was set aside as the test part. Ratings of the test part are
never used to tune anything; they only measure the tuned result.

Ratings come back as `ratings/{name}/review_<name>.json`
(`delfin review submit`), reservations as `ratings/{name}/claims/<name>.json`.
"""


def new_set(source, out_dir, categories_file, *, set_name: str | None = None,
            seed: int | None = None, test_fraction: float = DEFAULT_TEST_FRACTION,
            overlap: int = DEFAULT_OVERLAP, self_block: int = DEFAULT_BLOCK,
            redundancy: float = DEFAULT_REDUNDANCY, groups: dict[str, str] | None = None,
            renumber: bool = False) -> tuple[Path, list[str]]:
    """The public part of a new review set, from a multi-frame XYZ file.

    Writes ``<out_dir>/structures.xyz.gz`` (one frame per structure, the
    comment line only its id), ``kategorien.txt``, ``assignment.json`` with
    the shared structures, the seed and the test part drawn now, before any
    rating, and a ``README.md`` to complete.  The ids are the comment lines
    (plain and unique), or ``T0001`` ... in file order with *renumber*.
    Blinding (which structure comes from where) is not done here; the frames
    are taken as they are.
    """
    out_dir = Path(out_dir).expanduser()
    if out_dir.exists() and any(out_dir.iterdir()):
        raise ValueError(f'{out_dir} is not empty')
    frames = rv.parse_frames(Path(source).expanduser().read_text(encoding='utf-8'))
    if not frames:
        raise ValueError(f'no XYZ frames in {source}')
    cats = rv.clean_categories(
        line for line in Path(categories_file).expanduser().read_text(encoding='utf-8').splitlines()
        if not line.strip().startswith('#'))
    if not cats:
        raise ValueError(f'{categories_file} lists no category')
    if not test_fraction or not 0 < float(test_fraction) < 1:
        raise ValueError('a new set needs a test part: --test-fraction between 0 and 1')
    if renumber:
        width = max(4, len(str(len(frames))))
        ids = [f'T{k:0{width}d}' for k in range(1, len(frames) + 1)]
    else:
        ids = [f.comment.strip() for f in frames]
        bad = [i for i in ids if not re.fullmatch(r'[A-Za-z0-9_.-]+', i)]
        if bad:
            raise ValueError(f'{len(bad)} comment line(s) are no plain id (first: {bad[0]!r}); '
                             'give --renumber')
        if len(set(ids)) != len(ids):
            raise ValueError('the ids (comment lines) are not unique; give --renumber')
    out_dir.mkdir(parents=True, exist_ok=True)
    text = ''.join(f'{f.n_atoms}\n{ident}\n{f.coords}\n' for f, ident in zip(frames, ids))
    with open(out_dir / SET_ARCHIVE_NAME, 'wb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, filename='') as handle:
            handle.write(text.encode('utf-8'))
    (out_dir / rv.CATEGORY_FILE_NAME).write_text('\n'.join(cats) + '\n', encoding='utf-8')
    name = set_name or out_dir.name
    data, msgs = assign(out_dir, [], overlap=overlap, seed=seed, test_fraction=test_fraction,
                        set_name=name, groups=groups, self_block=self_block,
                        redundancy=redundancy)
    write_assignment(out_dir / ASSIGNMENT_NAME, data)
    readme = out_dir / 'README.md'
    readme.write_text(README_TEMPLATE.format(
        name=name, n=len(frames), n_cats=len(cats), overlap=len(data['shared']),
        seed=data['seed'], block=data['self_block'], redundancy=data['redundancy'],
        days=CLAIM_DAYS, test_fraction=data['split']['test_fraction'],
        split_seed=data['split']['seed']), encoding='utf-8')
    msgs.append(f'test part: {len(data["split"]["test"])} structures, drawn now (seed '
                f'{data["split"]["seed"]})')
    msgs.append(f'written: {out_dir} (complete the TODOs in {readme.name})')
    return out_dir, msgs


# ---------------------------------------------------------------------------
# A published set: unpacking it, and sending the ratings back
# ---------------------------------------------------------------------------

def archive_frames(set_dir) -> list[tuple[str, str]]:
    """(id, frame text) of every frame in ``<set_dir>/structures.xyz.gz``; the id is the comment line."""
    archive = Path(set_dir).expanduser() / SET_ARCHIVE_NAME
    with gzip.open(archive, 'rt', encoding='utf-8') as handle:
        return frames_of_text(handle.read(), str(archive))


def frames_of_text(text: str, where: str = '') -> list[tuple[str, str]]:
    """(id, frame text) of every frame of a set's multi-frame XYZ text."""
    lines = text.splitlines(keepends=True)
    archive = where or SET_ARCHIVE_NAME
    blocks: list[tuple[str, str]] = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        ident = lines[i + 1].strip()
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', ident):
            raise ValueError(f'{archive}: frame {len(blocks) + 1} has no plain id ({ident!r})')
        blocks.append((ident, ''.join(lines[i:i + 2 + n])))
        i += 2 + n
    return blocks


def unpack_for(set_dir, dest, name: str | None = None) -> tuple[int, int, str]:
    """:func:`unpack`, then, for *name*, their share: an assigned one is only
    reported; a new name draws its block (:func:`self_share`) into
    ``review_<name>.json`` at once.  Returns (written, unchanged, message)."""
    written, same = unpack(set_dir, dest)
    msg = ''
    found = find_assignment(dest)
    if name and found:
        assignment = load_assignment(found)
        rater, order, own = share_of(dest, assignment, name)
        if own is not None:
            path = rv.default_review_path(Path(dest).expanduser(), rater)
            data = rv.load_review(path)
            if data.get('assignment') != own:
                data['assignment'] = own
                rv.write_json_atomic(path, data)
            msg = (f'{rater}: a block of {len(own["files"])} drawn for you, plus '
                   f'{len(assignment.get("shared", []))} everyone rates (kept in {path.name})')
        else:
            msg = f'{rater}: your share is {len(order)} structures'
    return written, same, msg


def unpack(set_dir, dest) -> tuple[int, int]:
    """Write the frames of ``<set_dir>/structures.xyz.gz`` as ``<dest>/<id>.xyz``
    (the id is the comment line), next to the set's ``kategorien.txt`` and
    ``assignment.json``.  Returns (written, unchanged); a file that exists
    with other content is an error, so a review in progress is never altered.
    """
    set_dir, dest = Path(set_dir).expanduser(), Path(dest).expanduser()
    blocks = archive_frames(set_dir)
    dest.mkdir(parents=True, exist_ok=True)
    wanted = [(f'{ident}.xyz', block) for ident, block in blocks]
    for extra in (rv.CATEGORY_FILE_NAME, ASSIGNMENT_NAME):
        if (set_dir / extra).is_file():
            wanted.append((extra, (set_dir / extra).read_text(encoding='utf-8')))
    clash = [name for name, body in wanted
             if (dest / name).exists() and (dest / name).read_text(encoding='utf-8') != body
             and name != ASSIGNMENT_NAME]
    if clash:
        raise ValueError(f'{len(clash)} file(s) in {dest} differ from the set '
                         f'(first: {clash[0]}); unpack into an empty folder')
    written = same = 0
    for name, body in wanted:
        path = dest / name
        if path.exists() and path.read_text(encoding='utf-8') == body:
            same += 1
            continue
        path.write_text(body, encoding='utf-8')
        written += 1
    # Where the set and its ratings are, so that a block drawn here sees the
    # ratings and claims of the others (after a git pull).
    source = {'set': set_dir.name, 'set_dir': str(set_dir.resolve())}
    if set_dir.parent.name == 'review_sets':
        source['ratings'] = str((set_dir.parent.parent / RATINGS_DIR / set_dir.name).resolve())
    rv.write_json_atomic(dest / SOURCE_NAME, source)
    return written, same


def _git(repo, *args, check=True, timeout=300.0) -> subprocess.CompletedProcess:
    res = safe_git.run(list(args), repo, timeout=timeout)
    if check and res.returncode:
        raise subprocess.CalledProcessError(res.returncode, res.args, res.stdout, res.stderr)
    return res


def _token(folder: Path, name: str) -> str:
    """A random mark of this reviewer on this machine (``.review_<name>.token``): a
    second person who picks the same name by mistake cannot overwrite the ratings."""
    path = folder / ('.' + _file_stem(name)[:-len('.json')] + '.token')
    if path.is_file():
        tok = path.read_text(encoding='utf-8').strip()
        if tok:
            return tok
    tok = secrets.token_hex(12)
    path.write_text(tok + '\n', encoding='utf-8')
    return tok


def _repo_top(repo) -> Path:
    if repo is None:
        repo = Path(__file__).resolve().parent.parent
    top = _git(repo, 'rev-parse', '--show-toplevel', check=False)
    if top.returncode:
        raise ValueError(f'{repo} is not a git checkout of DELFIN; give --repo <your clone>')
    return Path(top.stdout.strip())


def _gh_login(run) -> str:
    """The GitHub account gh is logged in with ('' when unknown)."""
    if shutil.which('gh') is None:
        return ''
    res = run(['gh', 'api', 'user', '-q', '.login'])
    return res.stdout.strip() if res.returncode == 0 else ''


def owner_path(set_name: str, stem: str) -> str:
    """Where the GitHub account a reviewer's name belongs to is kept."""
    return f'{RATINGS_DIR}/{set_name}/{OWNERS_DIR}/{stem}.json'


def _send(repo: Path, rel: str, data: dict[str, Any], *, token: str, name: str,
          branch: str, message: tuple[str, str], title: str, body: str,
          token_of, remote: str, base: str, upstream: str, push: bool, fetch: bool,
          run, set_name: str, stem: str, login: str | None = None) -> dict[str, Any]:
    """Commit *data* as *rel* on a new *branch* from ``<remote>/<base>``, push it
    and open a pull request (to a fork without write access).  A file of that
    name already on *base* is replaced only when ``token_of`` reads this
    reviewer's *token* from it.

    The first pull request of a name also adds ``ratings/<set>/owners/<name>.json``
    with the GitHub account it comes from (*login*, else the one ``gh`` is
    logged in with when pushing): from then on the check accepts pull
    requests for this name only from that account.  Nothing else changes."""
    if fetch:
        res = _git(repo, 'fetch', remote, base, check=False)
        if res.returncode:
            raise RuntimeError(f'git fetch {remote} {base} failed: {res.stderr.strip()}')
    start = f'{remote}/{base}'
    there = _git(repo, 'show', f'{start}:{rel}', check=False)
    if there.returncode == 0:
        try:
            theirs = token_of(json.loads(there.stdout))
        except (ValueError, AttributeError):
            theirs = ''
        if theirs and theirs != token:
            raise ValueError(f'{rel} on {base} was sent by another reviewer named "{name}"; '
                             'pick a name of your own')
    login = (login or '').strip() or (_gh_login(run) if push else '')
    owner_rel = owner_path(set_name, stem)
    files = {rel: data}
    bound = _git(repo, 'show', f'{start}:{owner_rel}', check=False)
    note = ''
    if bound.returncode == 0:
        try:
            who = str(json.loads(bound.stdout).get('login') or '')
        except (ValueError, AttributeError):
            who = ''
        if login and who and who.casefold() != login.casefold():
            raise ValueError(f'the name "{name}" belongs to the GitHub account {who}; '
                             'pick a name of your own')
    elif login:
        files[owner_rel] = {'version': SCHEMA_VERSION, 'set': set_name, 'name': name,
                            'login': login, 'since': _now()}
    else:
        note = (f' First contribution under this name: the check asks for {owner_rel} '
                'with your GitHub account; send again with --github <your login>.')
    stem_branch, k = branch, 1
    while _git(repo, 'rev-parse', '--verify', '--quiet', f'refs/heads/{branch}',
               check=False).returncode == 0:
        k += 1
        branch = f'{stem_branch}-{k}'
    tmp = Path(tempfile.mkdtemp(prefix='delfin-submit-'))
    work = tmp / 'wt'
    try:
        _git(repo, 'worktree', 'add', '--quiet', '-b', branch, str(work), start)
        for path, content in files.items():
            target = work / path
            target.parent.mkdir(parents=True, exist_ok=True)
            rv.write_json_atomic(target, content)
            _git(work, 'add', '--', path)
        changed = _git(work, 'diff', '--cached', '--name-only').stdout.split()
        if sorted(changed) != sorted(files):
            raise RuntimeError(f'the commit would touch {changed}, not only {sorted(files)}')
        _git(work, 'commit', '--quiet', '-m', message[0], '-m', message[1])
    finally:
        _git(repo, 'worktree', 'remove', '--force', str(work), check=False)
        shutil.rmtree(tmp, ignore_errors=True)
    out: dict[str, Any] = {'branch': branch, 'path': rel, 'pushed_to': None, 'pr': None,
                           'next': '', 'owner': owner_rel if owner_rel in files else None}
    if not push:
        out['next'] = (f'not pushed; push with: git -C {repo} push {remote} {branch}, '
                       f'then open a pull request into {upstream} {base}.' + note)
        return out
    has_gh = shutil.which('gh') is not None
    head_owner = ''
    push_remote = remote
    if has_gh:
        perm = run(['gh', 'repo', 'view', upstream, '--json', 'viewerPermission',
                    '-q', '.viewerPermission'])
        if perm.returncode == 0 and perm.stdout.strip() not in ('ADMIN', 'MAINTAIN', 'WRITE'):
            # No write access: the branch goes to the account's fork.
            fork = run(['gh', 'repo', 'fork', upstream, '--clone=false'])
            if fork.returncode and 'already exists' not in (fork.stderr or ''):
                out['next'] = f'could not fork {upstream}: {fork.stderr.strip()}'
                return out
            if _git(repo, 'remote', 'get-url', 'fork', check=False).returncode:
                login = run(['gh', 'api', 'user', '-q', '.login']).stdout.strip()
                _git(repo, 'remote', 'add', 'fork',
                     f'https://github.com/{login}/{upstream.split("/")[-1]}.git', check=False)
            push_remote = 'fork'
            head_owner = run(['gh', 'api', 'user', '-q', '.login']).stdout.strip()
    pushed = _git(repo, 'push', '--quiet', push_remote, f'{branch}:{branch}', check=False)
    if pushed.returncode:
        out['next'] = (f'the push to {push_remote} failed ({pushed.stderr.strip()}); '
                       f'the commit is on the local branch {branch}. Send '
                       f'{Path(rel).name} to the coordinator instead.')
        return out
    out['pushed_to'] = push_remote
    if not has_gh:
        out['next'] = (f'pushed {branch} to {push_remote}; open a pull request into '
                       f'{upstream} {base} on GitHub (gh is not installed).' + note)
        return out
    head = f'{head_owner}:{branch}' if head_owner else branch
    pr = run(['gh', 'pr', 'create', '--repo', upstream, '--base', base, '--head', head,
              '--title', title, '--body', body])
    if pr.returncode:
        out['next'] = f'pushed {branch}; gh pr create failed: {pr.stderr.strip()}'
    else:
        out['pr'] = pr.stdout.strip()
        out['next'] = note.strip()
    return out


def _runner(repo):
    return lambda cmd: subprocess.run(cmd, capture_output=True, text=True, cwd=str(repo or '.'))


def _version_for(assignment: dict[str, Any] | None, repo: Path, set_name: str) -> str:
    """The set version a file sent back names: from the assignment, else from the
    set in the clone (an assignment unpacked before the version was recorded)."""
    version = (assignment or {}).get('set_version')
    if version:
        return version
    return set_version(repo / 'review_sets' / set_name) or ''


def _stamp() -> str:
    return datetime.now(timezone.utc).strftime('%Y%m%d-%H%M%S')


def submit(folder, name: str, *, repo=None, set_name: str | None = None,
           remote: str = 'origin', base: str = 'main', upstream: str = UPSTREAM,
           push: bool = True, fetch: bool = True, runner=None,
           credit: str | None = None, github: str | None = None) -> dict[str, Any]:
    """Add ``<folder>/review_<name>.json`` as ``ratings/<set>/review_<name>.json`` in a
    commit of its own on a new branch from ``<remote>/<base>``, push it and open a
    pull request (to a fork when the account cannot push to *upstream*).

    The commit adds or changes only that one file, and a file of that name
    already on *base* is replaced only when it carries this reviewer's token.
    *credit* (a name to be listed under, or ``''`` for none) is stored in the
    file as ``credit``.  *github* is the account the name is bound to on its
    first pull request (default: the one ``gh`` is logged in with).  Returns
    what was done: ``branch``, ``path``,
    ``pushed_to``, ``pr`` and ``next`` (what to do by hand when something
    could not be done here).
    """
    run = runner or _runner(repo)
    folder = Path(folder).expanduser()
    found = find_assignment(folder)
    assignment = load_assignment(found) if found else None
    if assignment is not None:
        # An assigned name as the assignment spells it; any other name is a
        # helper who drew a block of their own.
        name = find_rater(assignment, name) or name.strip()
    source = rv.default_review_path(folder, name)
    if not source.is_file():
        raise FileNotFoundError(f'no review file {source}: rate in the Review tab first, '
                                'with the same folder and name')
    set_name = set_name or (assignment or {}).get('set')
    if not set_name:
        raise ValueError('which set is this? give --set (the folder name under review_sets/)')
    repo = _repo_top(repo)
    data = rv.load_review(source)
    rated = sum(1 for *_, stale in rv.iter_records(data) if not stale)
    if not rated:
        raise ValueError(f'{source} holds no rating yet')
    if credit is not None:
        # Kept in the reviewer's own file too, so a later submit keeps it.
        data['credit'] = _clean_credit(credit)
        here = rv.load_review(source)
        if here.get('credit') != data['credit']:
            here['credit'] = data['credit']
            rv.write_json_atomic(source, here)
    token = _token(folder, name)
    data['submitted'] = {'rater': name, 'token': token, 'time': _now(), 'rated': rated,
                         'set': set_name,
                         'set_version': _version_for(assignment, repo, set_name)}
    rel = f'{RATINGS_DIR}/{set_name}/{source.name}'
    stem = source.name[len('review_'):-len('.json')]
    out = _send(repo, rel, data, token=token, name=name,
                branch=f'ratings/{set_name}/{stem}-{_stamp()}',
                message=(f'Ratings of {name} for the review set {set_name}',
                         f'{rated} rated frames; adds only {rel}.'),
                title=f'Ratings of {name} for the review set {set_name}',
                body=(f'Adds {rel} ({rated} rated frames), written by the Review tab and sent '
                      'with "delfin review submit". Only this file changes (the first time also '
                      'the owner file of the name, with the GitHub account).'),
                token_of=lambda d: (d.get('submitted') or {}).get('token', ''),
                remote=remote, base=base, upstream=upstream, push=push, fetch=fetch, run=run,
                set_name=set_name, stem=stem, login=github)
    out['rated'] = rated
    return out


def _clean_credit(text: str) -> str:
    """A display name for the list of contributors: one line, at most 80 characters."""
    text = ' '.join(str(text or '').split())
    if len(text) > 80:
        raise ValueError('the credit name is longer than 80 characters')
    return text


def claim(folder, name: str, *, repo=None, set_name: str | None = None,
          remote: str = 'origin', base: str = 'main', upstream: str = UPSTREAM,
          push: bool = True, fetch: bool = True, runner=None,
          github: str | None = None) -> dict[str, Any]:
    """Reserve the block a helper draws for themselves, at once, for everyone.

    Draws (or reads back) the block of :func:`self_share`, keeps it in
    ``review_<name>.json`` as ``unpack`` does, and sends
    ``ratings/<set>/claims/<name>.json`` (set, set version, name, files, date)
    as a pull request of its own, the same way as :func:`submit`.  Others who
    draw after it is merged leave the block out for ``CLAIM_DAYS`` days from
    its date; claiming again renews the date and keeps the block.  A name
    listed in ``assignment.json`` is reserved there already.
    """
    run = runner or _runner(repo)
    folder = Path(folder).expanduser()
    found = find_assignment(folder)
    if found is None:
        raise ValueError(f'no {ASSIGNMENT_NAME} in {folder}: unpack the set first')
    assignment = load_assignment(found)
    rater = find_rater(assignment, name)
    if rater is not None:
        raise ValueError(f'{rater} has a block in {ASSIGNMENT_NAME} already; nothing to claim')
    own = self_share(folder, assignment, name)
    name = own['name']
    if not own['files']:
        raise ValueError('nothing is left to claim in this set: every structure is taken')
    path = rv.default_review_path(folder, name)
    data = rv.load_review(path)
    if data.get('assignment') != own:
        data['assignment'] = own
        rv.write_json_atomic(path, data)
    set_name = set_name or assignment.get('set')
    if not set_name:
        raise ValueError('which set is this? give --set (the folder name under review_sets/)')
    repo = _repo_top(repo)
    token = _token(folder, name)
    stem = path.name[len('review_'):-len('.json')]
    rel = f'{RATINGS_DIR}/{set_name}/{CLAIMS_DIR}/{stem}.json'
    body = {'version': SCHEMA_VERSION, 'set': set_name,
            'set_version': _version_for(assignment, repo, set_name), 'name': name,
            'files': own['files'], 'date': _now(), 'days': CLAIM_DAYS, 'token': token}
    out = _send(repo, rel, body, token=token, name=name,
                branch=f'ratings/{set_name}/claim-{stem}-{_stamp()}',
                message=(f'Claim of {name} on {len(own["files"])} structures of the review set '
                         f'{set_name}', f'Reserves the block for {CLAIM_DAYS} days; adds only '
                         f'{rel}.'),
                title=f'Claim of {name} on {len(own["files"])} structures of {set_name}',
                body=(f'Adds {rel}: reserves {len(own["files"])} structures for {name} for '
                      f'{CLAIM_DAYS} days ("delfin review claim"). Only this file changes (the '
                      'first time also the owner file of the name).'),
                token_of=lambda d: d.get('token', ''),
                remote=remote, base=base, upstream=upstream, push=push, fetch=fetch, run=run,
                set_name=set_name, stem=stem, login=github)
    out['files'] = len(own['files'])
    return out


# ---------------------------------------------------------------------------
# Command line: delfin review assign | merge | unpack | submit | claim | status | new-set
# ---------------------------------------------------------------------------

def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def _names(text: Sequence[str]) -> list[str]:
    return [n for part in text for n in re.split(r'[,\s]+', part) if n]


def _send_args(p) -> None:
    p.add_argument('--set', dest='set_name', default=None,
                   help='the set (default: from assignment.json)')
    p.add_argument('--repo', default=None, help='your clone of DELFIN (default: this one)')
    p.add_argument('--remote', default='origin')
    p.add_argument('--base', default='main')
    p.add_argument('--upstream', default=UPSTREAM)
    p.add_argument('--no-push', action='store_true', help='only commit on a local branch')
    p.add_argument('--github', default=None,
                   help='your GitHub login (default: the one gh is logged in with); the first '
                        'pull request of a name binds the name to it')


def add_parsers(sub) -> None:
    a = sub.add_parser('assign', help='Share the structures of a folder out to reviewers '
                                      '(writes assignment.json)')
    a.add_argument('folder', help='the folder with the .xyz files')
    a.add_argument('raters', nargs='*', help='names of the reviewers (commas or spaces); '
                                             'names already assigned are kept as they are')
    a.add_argument('--overlap', type=int, default=None,
                   help=f'structures every reviewer rates, for the agreement '
                        f'(new assignment only; default {DEFAULT_OVERLAP}, 0 allowed)')
    a.add_argument('--block', type=int, default=DEFAULT_BLOCK,
                   help=f'own structures per new reviewer (default {DEFAULT_BLOCK}); '
                        f'the rest stays for later')
    a.add_argument('--self-block', type=int, default=None,
                   help=f'the block a name not in the assignment draws for itself '
                        f'(default {DEFAULT_BLOCK}; kept when extending)')
    a.add_argument('--redundancy', type=float, default=None,
                   help=f'share of the set meant for two or three reviewers, spread over '
                        f'the set (default {DEFAULT_REDUNDANCY:g}; kept when extending)')
    a.add_argument('--seed', type=int, default=None, help='seed of a new assignment')
    a.add_argument('--test-fraction', type=float, default=None,
                   help='mark this fraction of the files as a held-out test part '
                        '(new assignment only)')
    a.add_argument('--ratings', action='append', default=[],
                   help='further folder(s) with review_<name>.json and claims/ (default: '
                        'ratings/<set> of the clone for review_sets/<set>)')
    a.add_argument('--groups', default=None,
                   help='JSON {file: group}: files of one group stay in one part of the '
                        'test split (e.g. a structure shown twice)')
    a.add_argument('--set', dest='set_name', default=None, help='name of the set')
    a.add_argument('--out', default=None,
                   help='where to write (default: assignment.json in the folder, or the '
                        'one found there or above); an existing file there is extended')
    a.add_argument('--regroup-split', action='store_true',
                   help='before any rating only: move groups of geometrically identical '
                        'structures (and --groups) that straddle the test split into one part')
    a.add_argument('--new', action='store_true',
                   help='ignore an existing assignment.json and draw a new one')
    m = sub.add_parser('merge', help='All review_<name>.json of a folder in one CSV, the '
                                     'agreement, the reliability and the consensus')
    m.add_argument('folder', help='the folder with the structures (and assignment.json)')
    m.add_argument('--ratings', action='append', default=[],
                   help='further folder(s) with review_<name>.json (default: ratings/<set> '
                        'of the clone for review_sets/<set>)')
    m.add_argument('--csv', default=None,
                   help='CSV to write (default: review_merged_<folder name>.csv here)')
    m.add_argument('--consensus', default=None,
                   help='consensus CSV (default: <csv stem>_consensus.csv)')
    m.add_argument('--raters-csv', default=None,
                   help='reliability of every reviewer per question (default: '
                        '<csv stem>_raters.csv)')
    m.add_argument('--split', choices=('all', 'calibration', 'test'), default='all',
                   help='only this part of the split, for every number and file written '
                        '(calibration: the test part never enters)')
    m.add_argument('--bootstrap', type=int, default=None,
                   help='bootstrap draws for the interval of alpha (default 1000)')
    m.add_argument('--repeats', default=None,
                   help='private file of the structures shown twice (JSON {file: group} or '
                        '[[a, b], ...]); default: found from the geometry of the set')
    m.add_argument('--controls', default=None,
                   help='private key of control structures (JSON or CSV, e.g. the crystal '
                        'structures); only printed')
    m.add_argument('--no-csv', action='store_true', help='only print, write no CSV')
    u = sub.add_parser('unpack', help='Write the structures of a published set as single '
                                      'files into a folder to open in the Review tab')
    u.add_argument('set_dir', help='e.g. review_sets/kalib_2026_10')
    u.add_argument('dest', help='an empty folder, e.g. ~/review_kalib_2026_10')
    u.add_argument('--name', default=None,
                   help='your name: a name not in assignment.json draws its own block')
    s = sub.add_parser('submit', help='Send your review_<name>.json back as a pull request')
    s.add_argument('folder', help='the folder you rated in')
    s.add_argument('--name', required=True, help='your name, as in the Review tab')
    s.add_argument('--credit', default=None,
                   help='the name to be listed under as a contributor ("" for none)')
    _send_args(s)
    c = sub.add_parser('claim', help='Reserve your block for 28 days, as a pull request '
                                     'of its own')
    c.add_argument('folder', help='the folder you unpacked the set into')
    c.add_argument('--name', required=True, help='your name, as in the Review tab')
    _send_args(c)
    t = sub.add_parser('status', help='How far a set is: how often each structure is '
                                      'rated, reservations, reviewers, alpha')
    t.add_argument('folder', help='e.g. review_sets/kalib_2026_10')
    t.add_argument('--ratings', action='append', default=[],
                   help='further folder(s) with review_<name>.json and claims/')
    t.add_argument('--bootstrap', type=int, default=None,
                   help='bootstrap draws for the interval of alpha (default 1000)')
    n = sub.add_parser('new-set', help='Make the public part of a new review set from a '
                                       'multi-frame XYZ file')
    n.add_argument('source', help='multi-frame XYZ, one frame per structure')
    n.add_argument('out_dir', help='e.g. review_sets/<name> (must be empty or new)')
    n.add_argument('--categories', required=True, help='kategorien.txt, one category per line')
    n.add_argument('--set', dest='set_name', default=None,
                   help='name of the set (default: the folder name)')
    n.add_argument('--seed', type=int, default=None)
    n.add_argument('--test-fraction', type=float, default=DEFAULT_TEST_FRACTION,
                   help=f'the held-out test part, drawn now (default {DEFAULT_TEST_FRACTION:g})')
    n.add_argument('--overlap', type=int, default=DEFAULT_OVERLAP)
    n.add_argument('--self-block', type=int, default=DEFAULT_BLOCK)
    n.add_argument('--redundancy', type=float, default=DEFAULT_REDUNDANCY)
    n.add_argument('--groups', default=None,
                   help='JSON {id or file: group}: one group stays in one part of the split')
    n.add_argument('--renumber', action='store_true',
                   help='give the structures the ids T0001 ... in file order')


def _print_sent(res: dict[str, Any], what: str) -> None:
    print(f'committed {res["path"]} ({what}) on branch {res["branch"]}')
    if res['pr']:
        print(f'pull request: {res["pr"]}')
    if res['next']:
        print(res['next'])


def run(args) -> int:
    try:
        if args.cmd == 'assign':
            return _run_assign(args)
        if args.cmd == 'merge':
            return _run_merge(args)
        if args.cmd == 'unpack':
            written, same, msg = unpack_for(args.set_dir, args.dest, args.name)
            print(f'{written} file(s) written, {same} unchanged in {args.dest}; '
                  'open this folder in the Review tab')
            if msg:
                print(msg)
            return 0
        if args.cmd in ('submit', 'claim'):
            kw = dict(repo=args.repo, set_name=args.set_name, remote=args.remote,
                      base=args.base, upstream=args.upstream, push=not args.no_push,
                      github=args.github)
            if args.cmd == 'submit':
                res = submit(args.folder, args.name, credit=args.credit, **kw)
                _print_sent(res, f'{res["rated"]} rated')
            else:
                res = claim(args.folder, args.name, **kw)
                _print_sent(res, f'{res["files"]} structures for {CLAIM_DAYS} days')
            return 0 if (args.no_push or res['pushed_to']) else 1
        if args.cmd == 'status':
            print(format_status(status(args.folder, args.ratings, n_boot=args.bootstrap)))
            return 0
        if args.cmd == 'new-set':
            groups = None
            if args.groups:
                raw = json.loads(Path(args.groups).read_text(encoding='utf-8'))
                groups = {(k if k.endswith('.xyz') else f'{k}.xyz'): v for k, v in raw.items()}
            _, msgs = new_set(args.source, args.out_dir, args.categories,
                              set_name=args.set_name, seed=args.seed,
                              test_fraction=args.test_fraction, overlap=args.overlap,
                              self_block=args.self_block, redundancy=args.redundancy,
                              groups=groups, renumber=args.renumber)
            print('\n'.join(msgs))
            return 0
    except (OSError, ValueError, RuntimeError, KeyError, subprocess.CalledProcessError) as exc:
        detail = getattr(exc, 'stderr', '') or ''
        print(f'delfin review {args.cmd}: {exc} {detail}'.rstrip())
        return 1
    return 2


def _run_assign(args) -> int:
    folder = Path(args.folder).expanduser()
    if args.new:
        found = None
    elif args.out and Path(args.out).is_file():
        found = Path(args.out)           # extend the assignment kept elsewhere
    else:
        found = find_assignment(folder)
    existing = load_assignment(found) if found else None
    data, msgs = assign(folder, _names(args.raters), overlap=args.overlap, block=args.block,
                        seed=args.seed, test_fraction=args.test_fraction,
                        existing=existing, set_name=args.set_name, ratings=args.ratings,
                        self_block=args.self_block, redundancy=args.redundancy,
                        regroup=args.regroup_split,
                        groups=(json.loads(Path(args.groups).read_text(encoding='utf-8'))
                                if args.groups else None))
    out = Path(args.out) if args.out else (found or folder / ASSIGNMENT_NAME)
    write_assignment(out, data)
    print('\n'.join(msgs))
    if data.get('split'):
        print(f'test part: {len(data["split"]["test"])} structures (never use their ratings '
              'for tuning)')
    print(f'written: {out}')
    return 0


def _run_merge(args) -> int:
    from delfin import review_stats as rs
    folder = Path(args.folder).expanduser()
    found = find_assignment(folder)
    assignment = load_assignment(found) if found else None
    reviews = collect_reviews(*ratings_dirs(folder, args.ratings))
    if not reviews:
        print(f'no review_<name>.json in {folder}')
        return 1
    files = structure_files(folder)
    items = split_items(assignment, args.split, files)
    print(f'{len(reviews)} review file(s): {", ".join(sorted(reviews))}'
          + ('' if items is None else f'; only the {args.split} part ({len(items)} structures)'))
    stem = folder.resolve().name
    print(format_agreement(agreement(reviews, items)))
    res = reliability(reviews, items, categories=rv.read_category_file(folder) or (),
                      n_boot=args.bootstrap)
    print()
    print(format_reliability(res))
    if args.repeats:
        pairs, source = load_repeats(args.repeats), 'pairs from the private file'
    else:
        pairs, source = rs.find_repeats(set_frames(folder)), 'pairs found in the set'
    if items is not None:
        pairs = [(a, b) for a, b in pairs if a in items and b in items]
    print()
    print(format_self_consistency(rs.self_consistency(verdicts_by_file(reviews, items), pairs),
                                  source))
    if args.controls:
        p_block = {f: row[VERDICT_QUESTION] for (f, i), row in res['consensus'].items() if i == 0}
        print()
        print(format_controls(rs.control_report(verdicts_by_file(reviews, items),
                                                load_controls(args.controls) & set(files),
                                                p_block)))
    if not args.no_csv:
        print()
        base = Path(args.csv) if args.csv else Path(f'review_merged_{stem}.csv')
        n = write_merged_csv(reviews, base, assignment, args.split)
        print(f'{n} rows written: {base}')
        out = (Path(args.consensus) if args.consensus
               else base.with_name(f'{base.stem}_consensus.csv'))
        n = write_consensus_csv(res, out, assignment)
        print(f'consensus of {n} frames written: {out}')
        out = (Path(args.raters_csv) if args.raters_csv
               else base.with_name(f'{base.stem}_raters.csv'))
        n = write_raters_csv(res, out)
        print(f'reliability, {n} rows (reviewer x question), written: {out}')
    return 0


__all__ = [
    'ASSIGNMENT_NAME', 'CLAIM_DAYS', 'DEFAULT_OVERLAP', 'DEFAULT_REDUNDANCY', 'agreement',
    'assign', 'claim', 'cohen_kappa', 'collect_claims', 'collect_reviews', 'coverage',
    'draw_block', 'find_assignment', 'find_rater', 'fleiss_kappa', 'format_agreement',
    'format_reliability', 'format_status', 'load_assignment', 'merged_rows', 'new_set',
    'rater_order', 'raters_wanted', 'ratings_dirs', 'reliability', 'reservations',
    'self_share', 'set_version', 'share_of', 'status', 'structure_files', 'submit', 'unpack',
    'unpack_for', 'write_assignment', 'write_consensus_csv', 'write_merged_csv',
]
