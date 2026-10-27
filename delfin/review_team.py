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
  first ``overlap`` become the shared ones, and every new reviewer takes the
  next ``block`` files that nobody has yet (by default the unassigned files
  split evenly).  Reviewers added later get only files still unassigned, so
  nothing given out before changes; files a reviewer has already rated (in
  ``review_<name>.json``) stay theirs.  An optional split marks a fixed
  fraction of the files, drawn once with its own seed, as *test*: a
  pre-registered hold-out whose ratings are never used to tune anything.
* :func:`rater_order` is the order the Review tab shows a reviewer: the
  shared files mixed among as many of their own, then the rest of their own
  in file order, so the shared ones are rated early but do not stand out.
* :func:`agreement` measures, on the shared structures, how well the
  reviewers agree: Cohen's kappa for two, Fleiss' kappa for more, for the
  verdict (pass/block) and for every category, with the share of frames on
  which all agree and the number of frames.  :func:`write_merged_csv` writes
  every rating of every reviewer in one CSV.
* :func:`unpack` writes the structures of a published set
  (``structures.xyz.gz``, one frame per structure) as single files, and
  :func:`submit` brings one reviewer's ``review_<name>.json`` back into the
  repository as a pull request that adds only that one file.

``delfin review assign|merge|unpack|submit`` are the command line to all of it.
"""

from __future__ import annotations

import gzip
import json
import random
import re
import secrets
import shutil
import subprocess
import tempfile
from collections.abc import Iterable, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from delfin import review as rv

ASSIGNMENT_NAME = 'assignment.json'
SET_ARCHIVE_NAME = 'structures.xyz.gz'
DEFAULT_OVERLAP = 50
SCHEMA_VERSION = 1
UPSTREAM = 'ComPlat/DELFIN'
RATINGS_DIR = 'ratings'
COMMANDS = ('assign', 'merge', 'unpack', 'submit')


# ---------------------------------------------------------------------------
# The assignment
# ---------------------------------------------------------------------------

def structure_files(folder) -> list[str]:
    """The ``*.xyz`` file names of *folder*, in name order (the order of the Review tab)."""
    return sorted(p.name for p in Path(folder).glob('*.xyz') if p.is_file())


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


def rated_files(folder, name: str) -> set[str]:
    """The files *name* has rated so far in ``<folder>/review_<name>.json``."""
    path = rv.default_review_path(Path(folder), name)
    try:
        data = rv.load_review(path)
    except (OSError, ValueError):
        return set()
    return {fname for fname, _, _, stale in rv.iter_records(data) if not stale}


def _shuffled(files: Sequence[str], seed) -> list[str]:
    order = sorted(files)
    random.Random(seed).shuffle(order)
    return order


def assign(folder, raters: Iterable[str], *, overlap: int | None = None,
           block: int | None = None, seed: int | None = None,
           test_fraction: float | None = None,
           existing: dict[str, Any] | None = None,
           set_name: str | None = None) -> tuple[dict[str, Any], list[str]]:
    """Share the structures of *folder* out to *raters*; returns (assignment, messages).

    With *existing*, the reviewers there and their files stay as they are
    and only the new names get files, from those nobody has yet.
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
            data['split'] = {'seed': split_seed, 'test_fraction': frac,
                             'test': sorted(_shuffled(files, split_seed)[:n_test],
                                            key=place.get)}
    if set_name:
        data['set'] = set_name
    raters_now: dict[str, list[str]] = data['raters']
    stems = {_file_stem(n).lower(): n for n in raters_now}
    new: list[str] = []
    for name in names:
        if find_rater(data, name):
            msgs.append(f'{name}: already in the assignment, kept as it is')
            continue
        stem = _file_stem(name).lower()
        if stem in stems:
            raise ValueError(f'"{name}" and "{stems[stem]}" would write the same review file')
        stems[stem] = name
        new.append(name)
    shared = set(data['shared'])
    taken = shared | {f for v in raters_now.values() for f in v}
    pool = [f for f in _shuffled(files, data['seed']) if f not in taken]
    if new:
        per = int(block) if block is not None else len(pool) // len(new)
        if per < 0:
            raise ValueError('block must be 0 or more')
        for name in new:
            own = [f for f in pool if f in rated_files(folder, name) and f not in taken]
            taken.update(own)
            fresh = [f for f in pool if f not in taken][:max(0, per - len(own))]
            taken.update(fresh)
            raters_now[name] = sorted(own + fresh, key=place.get)
            note = f' ({len(own)} already rated)' if own else ''
            msgs.append(f'{name}: {len(raters_now[name])} own + {len(shared)} shared{note}')
    data['updated'] = _now()
    left = len(files) - len(taken)
    msgs.append(f'{len(files)} structures: {len(shared)} shared, '
                f'{len(taken) - len(shared)} given out, {left} not yet assigned')
    return data, msgs


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
    own = list(assignment['raters'][rater])
    mine = set(own)
    shared = [f for f in assignment.get('shared', []) if f not in mine]
    head = shared + own[:len(shared)]
    random.Random(f'{assignment.get("seed", 0)}:{rater}').shuffle(head)
    return head + own[len(shared):]


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


MERGED_FIELDS = ('rater',) + rv.CSV_FIELDS + ('stale', 'shared', 'split')


def merged_rows(reviews: dict[str, dict[str, Any]],
                assignment: dict[str, Any] | None = None, split: str = 'all'):
    """Every row of every reviewer (:func:`delfin.review.csv_rows`), with the
    reviewer, whether the file is shared, and its part of the split."""
    shared = set((assignment or {}).get('shared', []))
    test = set(((assignment or {}).get('split') or {}).get('test', []))
    has_split = bool((assignment or {}).get('split'))
    for name in sorted(reviews):
        for row in rv.csv_rows(reviews[name]):
            part = ('test' if row['file'] in test else 'calibration') if has_split else ''
            if split != 'all' and part != split:
                continue
            yield {'rater': name, **row, 'shared': 'yes' if row['file'] in shared else '',
                   'split': part}


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


def fleiss_kappa(counts: Sequence[int], raters: int) -> float | None:
    """Fleiss' kappa on a yes/no question: *counts* = yes-votes per item out of *raters*."""
    N, m = len(counts), int(raters)
    if not N or m < 2:
        return None
    p_items = [(k * k + (m - k) * (m - k) - m) / (m * (m - 1)) for k in counts]
    p_bar = sum(p_items) / N
    p = sum(counts) / (N * m)
    pe = p * p + (1 - p) * (1 - p)
    return None if pe >= 1 else (p_bar - pe) / (1 - pe)


def _measure(votes: list[list[bool]]) -> dict[str, Any]:
    """Agreement of one yes/no question; *votes* = one list per item, one vote per rater."""
    n = len(votes)
    m = len(votes[0]) if votes else 0
    if m == 2:
        kappa = cohen_kappa([v[0] for v in votes], [v[1] for v in votes])
    else:
        kappa = fleiss_kappa([sum(v) for v in votes], m)
    return {'kappa': kappa, 'agreement': (sum(len(set(v)) == 1 for v in votes) / n) if n else None,
            'n': n, 'yes': sum(sum(v) for v in votes)}


def agreement(reviews: dict[str, dict[str, Any]],
              items: Iterable[str] | None = None) -> dict[str, Any] | None:
    """How well the reviewers agree on the frames of the files *items* (all
    files rated by two or more when None), counting only frames every
    reviewer who rated any of them has rated.  None with fewer than two."""
    rated = {name: current_ratings(data) for name, data in reviews.items()}
    if items is None:
        seen: dict[str, int] = {}
        for recs in rated.values():
            for f in {f for f, _ in recs}:
                seen[f] = seen.get(f, 0) + 1
        wanted = {f for f, k in seen.items() if k >= 2}
    else:
        wanted = set(items)
    names = sorted(n for n, recs in rated.items() if any(f in wanted for f, _ in recs))
    if len(names) < 2:
        return None
    keys = sorted(set.intersection(*({k for k in rated[n] if k[0] in wanted} for n in names)))
    cats: list[str] = []
    for n in names:
        for c in reviews[n].get('categories') or []:
            if c not in cats:
                cats.append(c)
    for n in names:
        for k in keys:
            for c in rated[n][k].get('categories') or []:
                if c not in cats:
                    cats.append(c)
    out = {'raters': names, 'method': 'Cohen' if len(names) == 2 else 'Fleiss',
           'frames': len(keys),
           'verdict': _measure([[rated[n][k]['verdict'] == 'block' for n in names] for k in keys]),
           'unsure': sum(rv.is_unsure(rated[n][k].get('note', '')) for n in names for k in keys),
           'categories': {}}
    for c in cats:
        out['categories'][c] = _measure(
            [[c in (rated[n][k].get('categories') or []) for n in names] for k in keys])
    return out


def format_agreement(result: dict[str, Any] | None) -> str:
    if not result:
        return 'Agreement: fewer than two reviewers rated the shared structures.'

    def fmt(m):
        kappa = 'n/a' if m['kappa'] is None else f'{m["kappa"]:.3f}'
        agree = 'n/a' if m['agreement'] is None else f'{100 * m["agreement"]:.1f} %'
        return f'{kappa:>7}  {agree:>8}  {m["n"]:5d}  {m["yes"]:5d}'

    width = max([len('ok/defect')] + [len(c) for c in result['categories']])
    lines = [f'Agreement of {len(result["raters"])} reviewers ({result["method"]} kappa) '
             f'on {result["frames"]} frames rated by all: {", ".join(result["raters"])}',
             f'{"question".ljust(width)}    kappa  all agree      n    yes',
             f'{"ok/defect".ljust(width)}  {fmt(result["verdict"])}']
    for c, m in result['categories'].items():
        lines.append(f'{c.ljust(width)}  {fmt(m)}')
    lines.append('"yes" counts the votes for block / for the category over all reviewers; '
                 f'{result["unsure"]} of the votes are marked unsure.')
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# A published set: unpacking it, and sending the ratings back
# ---------------------------------------------------------------------------

def unpack(set_dir, dest) -> tuple[int, int]:
    """Write the frames of ``<set_dir>/structures.xyz.gz`` as ``<dest>/<id>.xyz``
    (the id is the comment line), next to the set's ``kategorien.txt`` and
    ``assignment.json``.  Returns (written, unchanged); a file that exists
    with other content is an error, so a review in progress is never altered.
    """
    set_dir, dest = Path(set_dir).expanduser(), Path(dest).expanduser()
    archive = set_dir / SET_ARCHIVE_NAME
    text = gzip.open(archive, 'rt', encoding='utf-8').read()
    blocks: list[tuple[str, str]] = []
    lines = text.splitlines(keepends=True)
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        block = ''.join(lines[i:i + 2 + n])
        ident = lines[i + 1].strip()
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', ident):
            raise ValueError(f'{archive}: frame {len(blocks) + 1} has no plain id ({ident!r})')
        blocks.append((ident, block))
        i += 2 + n
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
    return written, same


def _git(repo, *args, check=True) -> subprocess.CompletedProcess:
    return subprocess.run(['git', '-C', str(repo), *args], capture_output=True, text=True,
                          check=check)


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


def submit(folder, name: str, *, repo=None, set_name: str | None = None,
           remote: str = 'origin', base: str = 'main', upstream: str = UPSTREAM,
           push: bool = True, fetch: bool = True, runner=None) -> dict[str, Any]:
    """Add ``<folder>/review_<name>.json`` as ``ratings/<set>/review_<name>.json`` in a
    commit of its own on a new branch from ``<remote>/<base>``, push it and open a
    pull request (to a fork when the account cannot push to *upstream*).

    The commit adds or changes only that one file, and a file of that name
    already on *base* is replaced only when it carries this reviewer's token.
    Returns what was done: ``branch``, ``path``, ``pushed_to``, ``pr`` and
    ``next`` (what to do by hand when something could not be done here).
    """
    run = runner or (lambda cmd: subprocess.run(cmd, capture_output=True, text=True,
                                                cwd=str(repo or '.')))
    folder = Path(folder).expanduser()
    found = find_assignment(folder)
    assignment = load_assignment(found) if found else None
    if assignment is not None:
        rater = find_rater(assignment, name)
        if rater is None:
            raise ValueError(f'"{name}" is not in {found}')
        name = rater
    source = rv.default_review_path(folder, name)
    if not source.is_file():
        raise FileNotFoundError(f'no review file {source}: rate in the Review tab first, '
                                'with the same folder and name')
    set_name = set_name or (assignment or {}).get('set')
    if not set_name:
        raise ValueError('which set is this? give --set (the folder name under review_sets/)')
    if repo is None:
        repo = Path(__file__).resolve().parent.parent
    top = _git(repo, 'rev-parse', '--show-toplevel', check=False)
    if top.returncode:
        raise ValueError(f'{repo} is not a git checkout of DELFIN; give --repo <your clone>')
    repo = Path(top.stdout.strip())
    data = rv.load_review(source)
    rated = sum(1 for *_, stale in rv.iter_records(data) if not stale)
    if not rated:
        raise ValueError(f'{source} holds no rating yet')
    token = _token(folder, name)
    data['submitted'] = {'rater': name, 'token': token, 'time': _now(), 'rated': rated}
    rel = f'{RATINGS_DIR}/{set_name}/{source.name}'
    if fetch:
        res = _git(repo, 'fetch', remote, base, check=False)
        if res.returncode:
            raise RuntimeError(f'git fetch {remote} {base} failed: {res.stderr.strip()}')
    start = f'{remote}/{base}'
    there = _git(repo, 'show', f'{start}:{rel}', check=False)
    if there.returncode == 0:
        try:
            theirs = json.loads(there.stdout).get('submitted') or {}
        except ValueError:
            theirs = {}
        if theirs.get('token') and theirs['token'] != token:
            raise ValueError(f'{rel} on {base} was sent by another reviewer named "{name}"; '
                             'ask for a name of your own in assignment.json')
    branch = f'ratings/{set_name}/{_file_stem(name)[len("review_"):-len(".json")]}' \
             f'-{datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")}'
    stem, k = branch, 1
    while _git(repo, 'rev-parse', '--verify', '--quiet', f'refs/heads/{branch}',
               check=False).returncode == 0:
        k += 1
        branch = f'{stem}-{k}'
    tmp = Path(tempfile.mkdtemp(prefix='delfin-submit-'))
    work = tmp / 'wt'
    try:
        _git(repo, 'worktree', 'add', '--quiet', '-b', branch, str(work), start)
        target = work / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        rv.write_json_atomic(target, data)
        _git(work, 'add', '--', rel)
        changed = _git(work, 'diff', '--cached', '--name-only').stdout.split()
        if changed != [rel]:
            raise RuntimeError(f'the commit would touch {changed}, not only {rel}')
        _git(work, 'commit', '--quiet', '-m', f'Ratings of {name} for the review set {set_name}',
             '-m', f'{rated} rated frames; adds only {rel}.')
    finally:
        _git(repo, 'worktree', 'remove', '--force', str(work), check=False)
        shutil.rmtree(tmp, ignore_errors=True)
    out: dict[str, Any] = {'branch': branch, 'path': rel, 'rated': rated,
                           'pushed_to': None, 'pr': None, 'next': ''}
    if not push:
        out['next'] = (f'not pushed; push with: git -C {repo} push {remote} {branch}, '
                       f'then open a pull request into {upstream} {base}')
        return out
    title = f'Ratings of {name} for the review set {set_name}'
    body = (f'Adds {rel} ({rated} rated frames), written by the Review tab and sent with '
            '"delfin review submit". Only this one file changes.')
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
                       f'{source} to the coordinator instead.')
        return out
    out['pushed_to'] = push_remote
    if not has_gh:
        out['next'] = (f'pushed {branch} to {push_remote}; open a pull request into '
                       f'{upstream} {base} on GitHub (gh is not installed)')
        return out
    head = f'{head_owner}:{branch}' if head_owner else branch
    pr = run(['gh', 'pr', 'create', '--repo', upstream, '--base', base, '--head', head,
              '--title', title, '--body', body])
    if pr.returncode:
        out['next'] = f'pushed {branch}; gh pr create failed: {pr.stderr.strip()}'
    else:
        out['pr'] = pr.stdout.strip()
    return out


# ---------------------------------------------------------------------------
# Command line: delfin review assign | merge | unpack | submit
# ---------------------------------------------------------------------------

def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def _names(text: Sequence[str]) -> list[str]:
    return [n for part in text for n in re.split(r'[,\s]+', part) if n]


def add_parsers(sub) -> None:
    a = sub.add_parser('assign', help='Share the structures of a folder out to reviewers '
                                      '(writes assignment.json)')
    a.add_argument('folder', help='the folder with the .xyz files')
    a.add_argument('raters', nargs='+', help='names of the reviewers (commas or spaces); '
                                             'names already assigned are kept as they are')
    a.add_argument('--overlap', type=int, default=None,
                   help=f'structures every reviewer rates, for the agreement '
                        f'(new assignment only; default {DEFAULT_OVERLAP}, 0 allowed)')
    a.add_argument('--block', type=int, default=None,
                   help='own structures per new reviewer (default: the unassigned ones '
                        'split evenly; the rest stays for reviewers added later)')
    a.add_argument('--seed', type=int, default=None, help='seed of a new assignment')
    a.add_argument('--test-fraction', type=float, default=None,
                   help='mark this fraction of the files as a held-out test part '
                        '(new assignment only)')
    a.add_argument('--set', dest='set_name', default=None, help='name of the set')
    a.add_argument('--out', default=None,
                   help='where to write (default: assignment.json in the folder, or the '
                        'one found there or above)')
    a.add_argument('--new', action='store_true',
                   help='ignore an existing assignment.json and draw a new one')
    m = sub.add_parser('merge', help='All review_<name>.json of a folder in one CSV, '
                                     'and the agreement on the shared structures')
    m.add_argument('folder', help='the folder with the structures (and assignment.json)')
    m.add_argument('--ratings', action='append', default=[],
                   help='further folder(s) with review_<name>.json (e.g. ratings/<set>)')
    m.add_argument('--csv', default=None, help='CSV to write (default: <folder>/review_merged.csv)')
    m.add_argument('--split', choices=('all', 'calibration', 'test'), default='all',
                   help='only the rows of this part of the split')
    m.add_argument('--no-csv', action='store_true', help='only print the agreement')
    u = sub.add_parser('unpack', help='Write the structures of a published set as single '
                                      'files into a folder to open in the Review tab')
    u.add_argument('set_dir', help='e.g. review_sets/kalib_2026_10')
    u.add_argument('dest', help='an empty folder, e.g. ~/review_kalib_2026_10')
    s = sub.add_parser('submit', help='Send your review_<name>.json back as a pull request')
    s.add_argument('folder', help='the folder you rated in')
    s.add_argument('--name', required=True, help='your name, as in the Review tab')
    s.add_argument('--set', dest='set_name', default=None,
                   help='the set (default: from assignment.json)')
    s.add_argument('--repo', default=None, help='your clone of DELFIN (default: this one)')
    s.add_argument('--remote', default='origin')
    s.add_argument('--base', default='main')
    s.add_argument('--upstream', default=UPSTREAM)
    s.add_argument('--no-push', action='store_true', help='only commit on a local branch')


def run(args) -> int:
    try:
        if args.cmd == 'assign':
            return _run_assign(args)
        if args.cmd == 'merge':
            return _run_merge(args)
        if args.cmd == 'unpack':
            written, same = unpack(args.set_dir, args.dest)
            print(f'{written} file(s) written, {same} unchanged in {args.dest}; '
                  'open this folder in the Review tab')
            return 0
        if args.cmd == 'submit':
            res = submit(args.folder, args.name, repo=args.repo, set_name=args.set_name,
                         remote=args.remote, base=args.base, upstream=args.upstream,
                         push=not args.no_push)
            print(f'committed {res["path"]} ({res["rated"]} rated) on branch {res["branch"]}')
            if res['pr']:
                print(f'pull request: {res["pr"]}')
            if res['next']:
                print(res['next'])
            return 0 if (args.no_push or res['pushed_to']) else 1
    except (OSError, ValueError, RuntimeError, KeyError, subprocess.CalledProcessError) as exc:
        detail = getattr(exc, 'stderr', '') or ''
        print(f'delfin review {args.cmd}: {exc} {detail}'.rstrip())
        return 1
    return 2


def _run_assign(args) -> int:
    folder = Path(args.folder).expanduser()
    found = None if args.new else find_assignment(folder)
    existing = load_assignment(found) if found else None
    data, msgs = assign(folder, _names(args.raters), overlap=args.overlap, block=args.block,
                        seed=args.seed, test_fraction=args.test_fraction,
                        existing=existing, set_name=args.set_name)
    out = Path(args.out) if args.out else (found or folder / ASSIGNMENT_NAME)
    write_assignment(out, data)
    print('\n'.join(msgs))
    if data.get('split'):
        print(f'test part: {len(data["split"]["test"])} structures (never use their ratings '
              'for tuning)')
    print(f'written: {out}')
    return 0


def _run_merge(args) -> int:
    folder = Path(args.folder).expanduser()
    found = find_assignment(folder)
    assignment = load_assignment(found) if found else None
    reviews = collect_reviews(folder, *args.ratings)
    if not reviews:
        print(f'no review_<name>.json in {folder}')
        return 1
    print(f'{len(reviews)} review file(s): {", ".join(sorted(reviews))}')
    if not args.no_csv:
        out = Path(args.csv) if args.csv else folder / 'review_merged.csv'
        n = write_merged_csv(reviews, out, assignment, args.split)
        print(f'{n} rows written: {out}')
    shared = assignment.get('shared') if assignment else None
    print(format_agreement(agreement(reviews, shared if shared else None)))
    return 0


__all__ = [
    'ASSIGNMENT_NAME', 'DEFAULT_OVERLAP', 'agreement', 'assign', 'cohen_kappa',
    'collect_reviews', 'find_assignment', 'find_rater', 'fleiss_kappa', 'format_agreement',
    'load_assignment', 'merged_rows', 'rater_order', 'structure_files', 'submit', 'unpack',
    'write_assignment', 'write_merged_csv',
]
