"""Checks a helper's pull request into ``ratings/`` before anyone merges it.

Helpers send their ratings (``delfin review submit``) and their reservations
(``delfin review claim``) as pull requests that add or change exactly one
file: ``ratings/<set>/review_<name>.json`` or
``ratings/<set>/claims/<name>.json`` -- the first time under a name also
``ratings/<set>/owners/<name>.json``, which binds the name to the GitHub
account of the pull request.  :func:`check_pull_request` reads the change
from git (base and head commits, nothing from the working tree) and says
whether it is such a pull request and whether it holds:

* exactly one rating or claim file is added or changed (and at most its
  owner file added), nothing deleted or renamed, no other path touched;
* the file is valid JSON of the right shape: the set exists under
  ``review_sets/<set>/`` on the base, the set version matches its
  ``assignment.json``, every structure named is in the set, every category is
  in its ``kategorien.txt``, every marked atom lies within its structure;
* the name inside the file is the name of the file;
* the token rule: a file that already exists on the base is replaced only by
  the same reviewer (the same token), so nobody overwrites someone else's
  file by mistake;
* the account rule: a name bound on the base comes only from its account
  (the author's login, ``--author``); an unbound name must be bound by this
  pull request, to the author's account; an owner file on the base never
  changes;
* a claim reserves no more than one block (``self_block``) and is dated
  today or earlier.

The result says ``first`` for the pull request that binds a name (a
maintainer reviews it) and ``automerge`` for a valid one from the account a
name is bound to.  A pull request that touches no rating, claim or owner file
is not a helper's and is left to the ordinary review (``helper`` is False,
nothing fails).  The check runs from the base commit
(``.github/workflows/ratings.yml``), so a pull request cannot change the
rules it is checked by.

``python -m delfin.review_check --repo <checkout> --base <sha> --head <sha> --author <login>``
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import sys
from collections.abc import Callable, Sequence
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from delfin import review as rv
from delfin import review_team as rt
from delfin.agent import safe_git

RATING_RE = re.compile(r'ratings/(?P<set>[A-Za-z0-9_.-]+)/review_(?P<name>[A-Za-z0-9_.-]+)\.json')
CLAIM_RE = re.compile(r'ratings/(?P<set>[A-Za-z0-9_.-]+)/claims/(?P<name>[A-Za-z0-9_.-]+)\.json')
OWNER_RE = re.compile(r'ratings/(?P<set>[A-Za-z0-9_.-]+)/owners/(?P<name>[A-Za-z0-9_.-]+)\.json')
# A GitHub login: letters, digits and single hyphens, at most 39 characters.
LOGIN_RE = re.compile(r'[A-Za-z0-9](?:[A-Za-z0-9]|-(?=[A-Za-z0-9])){0,38}')
MAX_BYTES = 30_000_000
MAX_CREDIT = 80


class Problem(Exception):
    pass


def _git_blob(repo, rev: str, path: str) -> bytes | None:
    res = safe_git.run(['show', f'{rev}:{path}'], repo, text=False)
    return res.stdout if res.returncode == 0 else None


def pr_changes(repo, base: str, head: str) -> list[tuple[str, str]]:
    """(status letter, path) of every change from the merge base of *base* and *head*."""
    res = safe_git.run(['diff', '--name-status', '--no-renames', f'{base}...{head}'], repo)
    if res.returncode:
        raise Problem(f'git diff {base}...{head} failed: {res.stderr.strip()}')
    out = []
    for line in res.stdout.splitlines():
        parts = line.split('\t')
        if len(parts) >= 2:
            out.append((parts[0][:1], parts[-1]))
    return out


def base_set(read: Callable[[str], bytes | None], set_name: str) -> dict[str, Any]:
    """The set as the base holds it: ``atoms`` per file, ``categories``, ``assignment``."""
    where = f'review_sets/{set_name}'
    raw = read(f'{where}/{rt.SET_ARCHIVE_NAME}')
    if raw is None:
        raise Problem(f'there is no review set {where}/ on the base')
    text = gzip.decompress(raw)
    atoms = {f'{ident}.xyz': int(block.split(None, 1)[0])
             for ident, block in rt.frames_of_text(text.decode('utf-8'), where)}
    cats_raw = read(f'{where}/{rv.CATEGORY_FILE_NAME}')
    cats = rv.clean_categories(
        line for line in (cats_raw or b'').decode('utf-8').splitlines()
        if not line.strip().startswith('#'))
    if not cats:
        raise Problem(f'{where}/{rv.CATEGORY_FILE_NAME} is missing on the base')
    assignment_raw = read(f'{where}/{rt.ASSIGNMENT_NAME}')
    assignment = json.loads(assignment_raw) if assignment_raw else {}
    return {'name': set_name, 'atoms': atoms, 'categories': cats, 'assignment': assignment,
            'version': assignment.get('set_version') or rt.set_version_of(text)}


def _stem(name: str) -> str:
    return rt._file_stem(name)[len('review_'):-len('.json')]


def _atoms_ok(values, n: int) -> bool:
    return isinstance(values, list) and all(
        isinstance(i, int) and not isinstance(i, bool) and 0 <= i < n for i in values)


def _check_record(rec, where: str, n_atoms: int, cats: set[str], *, rated: bool) -> list[str]:
    bad: list[str] = []
    if not isinstance(rec, dict):
        return [f'{where}: not an object']
    if rated and rec.get('verdict') not in rv.VERDICTS:
        bad.append(f'{where}: verdict {rec.get("verdict")!r} is not pass or block')
    marked = rec.get('categories', [])
    if not isinstance(marked, list) or not all(isinstance(c, str) for c in marked):
        bad.append(f'{where}: categories must be a list of names')
    else:
        unknown = sorted(set(marked) - cats)
        if unknown:
            bad.append(f'{where}: categories not in kategorien.txt: {", ".join(unknown)}')
    if 'atoms' in rec and not _atoms_ok(rec['atoms'], n_atoms):
        bad.append(f'{where}: atoms outside 0..{n_atoms - 1}')
    by_cat = rec.get('atoms_by_category', {})
    if not isinstance(by_cat, dict):
        bad.append(f'{where}: atoms_by_category must be an object')
    else:
        for c, v in by_cat.items():
            if c not in cats:
                bad.append(f'{where}: atoms marked for the unknown category {c!r}')
            if not _atoms_ok(v, n_atoms):
                bad.append(f'{where}: atoms of {c} outside 0..{n_atoms - 1}')
    if not isinstance(rec.get('note', ''), str):
        bad.append(f'{where}: the note must be text')
    for k, old in enumerate(rec.get('history') or []):
        bad += _check_record(old, f'{where} (earlier rating {k + 1})', n_atoms, cats, rated=False)
    return bad


def _per_frame(section, label: str, the_set: dict[str, Any], *, rated: bool) -> list[str]:
    bad: list[str] = []
    if not isinstance(section, dict):
        return [f'{label} must be an object']
    cats = set(the_set['categories'])
    for fname, frames in section.items():
        if fname not in the_set['atoms']:
            bad.append(f'{label}: {fname} is not in the set')
            continue
        if not isinstance(frames, dict):
            bad.append(f'{label}: {fname} must be an object')
            continue
        for key, rec in frames.items():
            if key != '0':
                bad.append(f'{label}: {fname} has no frame {key!r} (one frame per structure)')
                continue
            bad += _check_record(rec, fname, the_set['atoms'][fname], cats, rated=rated)
    return bad


def _credit_problems(data: dict[str, Any]) -> list[str]:
    credit = data.get('credit')
    if credit is None:
        return []
    if not isinstance(credit, str) or '\n' in credit or len(credit) > MAX_CREDIT:
        return [f'credit must be one line of at most {MAX_CREDIT} characters']
    return []


def _files_in_set(files, the_set, label: str) -> list[str]:
    if not isinstance(files, list) or not all(isinstance(f, str) for f in files):
        return [f'{label} must be a list of file names']
    missing = [f for f in files if f not in the_set['atoms']]
    if missing:
        return [f'{label}: {len(missing)} file(s) not in the set (first: {missing[0]})']
    if len(set(files)) != len(files):
        return [f'{label}: a file is listed twice']
    return []


def check_review(data, name: str, the_set: dict[str, Any], base: dict | None) -> list[str]:
    """What is wrong with ``review_<name>.json`` (empty when nothing)."""
    if not isinstance(data, dict):
        return ['the file is not a JSON object']
    bad: list[str] = []
    sent = data.get('submitted')
    if not isinstance(sent, dict):
        return ['no "submitted" block: send the file with "delfin review submit"']
    rater = sent.get('rater')
    if not isinstance(rater, str) or _stem(rater) != name:
        bad.append(f'the reviewer inside ({rater!r}) is not the name of the file ({name})')
    if not isinstance(sent.get('token'), str) or not sent.get('token'):
        bad.append('no token in "submitted"')
    if sent.get('set') != the_set['name']:
        bad.append(f'the file names the set {sent.get("set")!r}, it lies in {the_set["name"]}')
    if sent.get('set_version') != the_set['version']:
        bad.append(f'set version {sent.get("set_version")!r} is not that of the set '
                   f'({the_set["version"]}): unpack the newest set and submit again')
    if 'reviews' not in data:
        bad.append('no "reviews"')
    else:
        bad += _per_frame(data['reviews'], 'reviews', the_set, rated=True)
    if 'notes' in data:
        bad += _per_frame(data['notes'], 'notes', the_set, rated=False)
    # "files" describes every .xyz file of the reviewer's folder (also ones of
    # their own); only its shape is checked, the ratings above name the set.
    if not isinstance(data.get('files', {}), dict):
        bad.append('"files" must be an object')
    own = data.get('assignment')
    if own is not None:
        if not isinstance(own, dict):
            bad.append('"assignment" must be an object')
        else:
            bad += _files_in_set(own.get('files', []), the_set, 'assignment.files')
    bad += _credit_problems(data)
    if base is not None:
        theirs = ((base.get('submitted') or {}).get('token') or '') if isinstance(base, dict) else ''
        if theirs and theirs != sent.get('token'):
            bad.append('the file on the base was sent by another reviewer of this name '
                       '(other token): pick a name of your own')
    return bad


def check_claim(data, name: str, the_set: dict[str, Any], base: dict | None, *,
                now: datetime | None = None) -> list[str]:
    """What is wrong with ``claims/<name>.json`` (empty when nothing)."""
    if not isinstance(data, dict):
        return ['the file is not a JSON object']
    now = now or datetime.now(timezone.utc)
    bad: list[str] = []
    who = data.get('name')
    if not isinstance(who, str) or _stem(who) != name:
        bad.append(f'the reviewer inside ({who!r}) is not the name of the file ({name})')
    if not isinstance(data.get('token'), str) or not data.get('token'):
        bad.append('no token')
    if data.get('set') != the_set['name']:
        bad.append(f'the claim names the set {data.get("set")!r}, it lies in {the_set["name"]}')
    if data.get('set_version') != the_set['version']:
        bad.append(f'set version {data.get("set_version")!r} is not that of the set '
                   f'({the_set["version"]})')
    files = data.get('files')
    bad += _files_in_set(files, the_set, 'files')
    limit = int(the_set['assignment'].get('self_block', rt.DEFAULT_BLOCK))
    if isinstance(files, list) and not files:
        bad.append('the claim reserves nothing')
    if isinstance(files, list) and len(files) > limit:
        bad.append(f'the claim reserves {len(files)} structures, more than one block ({limit})')
    when = rt._parse_time(data.get('date'))
    if when is None:
        bad.append(f'date {data.get("date")!r} is not a date')
    elif when > now + timedelta(days=1):
        bad.append('the claim is dated in the future')
    if data.get('days', rt.CLAIM_DAYS) != rt.CLAIM_DAYS:
        bad.append(f'a claim runs {rt.CLAIM_DAYS} days, not {data.get("days")!r}')
    if base is not None:
        theirs = (base.get('token') or '') if isinstance(base, dict) else ''
        if theirs and theirs != data.get('token'):
            bad.append("the claim on the base is someone else's (other token)")
    return bad


def _load_json(raw: bytes | None):
    if raw is None:
        return None
    try:
        return json.loads(raw.decode('utf-8'))
    except (UnicodeDecodeError, ValueError):
        return None


def check_owner(data, name: str, set_name: str, author: str | None) -> list[str]:
    """What is wrong with a new ``owners/<name>.json`` (empty when nothing)."""
    if not isinstance(data, dict):
        return ['the owner file is not a JSON object']
    bad: list[str] = []
    who, login = data.get('name'), data.get('login')
    if not isinstance(who, str) or _stem(who) != name:
        bad.append(f'the owner file names {who!r}, not {name}')
    if data.get('set') != set_name:
        bad.append(f'the owner file names the set {data.get("set")!r}, it lies in {set_name}')
    if not isinstance(login, str) or not LOGIN_RE.fullmatch(login):
        bad.append(f'the owner file holds no GitHub login ({login!r})')
    elif author and login.casefold() != author.casefold():
        bad.append(f'the owner file binds the name to {login}, but the pull request comes '
                   f'from {author}')
    return bad


def check_pull_request(repo, base: str, head: str, *, author: str | None = None,
                       now: datetime | None = None) -> dict[str, Any]:
    """``{'helper', 'path', 'kind', 'problems', 'changes', 'owner', 'first',
    'known', 'automerge'}`` of a pull request.

    *author* is the GitHub login the pull request comes from.  The first pull
    request of a name adds ``ratings/<set>/owners/<name>.json`` with that
    login (``first``); later ones must come from the same account (``known``
    when it does).  ``automerge`` is true only for a valid pull request of a
    known account."""
    changes = pr_changes(repo, base, head)
    ours = [(s, p) for s, p in changes if RATING_RE.fullmatch(p) or CLAIM_RE.fullmatch(p)]
    owners = [(s, p) for s, p in changes if OWNER_RE.fullmatch(p)]
    out: dict[str, Any] = {'helper': bool(ours or owners), 'path': None, 'kind': None,
                           'problems': [], 'changes': changes, 'owner': None,
                           'first': False, 'known': False, 'automerge': False}
    if not out['helper']:
        return out
    other = [p for s, p in changes if (s, p) not in ours and (s, p) not in owners]
    if len(ours) != 1 or other or len(owners) > 1:
        out['problems'].append(
            "a helper's pull request changes exactly one rating or claim file (the first "
            f'time also its owner file); this one changes {len(changes)}: '
            + ', '.join(p for _, p in changes[:10]))
        return out
    status, path = ours[0]
    out['path'] = path
    if status not in ('A', 'M'):
        out['problems'].append(f'{path} is {"deleted" if status == "D" else "changed in kind"}; '
                               "a helper's pull request only adds or changes its file")
        return out
    match = RATING_RE.fullmatch(path) or CLAIM_RE.fullmatch(path)
    out['kind'] = 'claim' if CLAIM_RE.fullmatch(path) else 'review'
    owner_rel = f'ratings/{match["set"]}/owners/{match["name"]}.json'
    if owners and owners[0][1] != owner_rel:
        out['problems'].append(f'the owner file {owners[0][1]} is not that of {path}')
        return out
    raw = _git_blob(repo, head, path)
    if raw is None:
        out['problems'].append(f'{path} cannot be read at {head}')
        return out
    if len(raw) > MAX_BYTES:
        out['problems'].append(f'{path} is larger than {MAX_BYTES // 1_000_000} MB')
        return out
    try:
        data = json.loads(raw.decode('utf-8'))
    except (UnicodeDecodeError, ValueError) as exc:
        out['problems'].append(f'{path} is not valid JSON: {exc}')
        return out
    old = _load_json(_git_blob(repo, base, path))
    try:
        the_set = base_set(lambda p: _git_blob(repo, base, p), match['set'])
    except (Problem, ValueError, OSError) as exc:
        out['problems'].append(str(exc))
        return out
    if out['kind'] == 'claim':
        out['problems'] += check_claim(data, match['name'], the_set, old, now=now)
    else:
        out['problems'] += check_review(data, match['name'], the_set, old)
    # The GitHub account the name belongs to.
    bound = _load_json(_git_blob(repo, base, owner_rel))
    if bound is not None:
        login = str(bound.get('login') or '') if isinstance(bound, dict) else ''
        out['owner'] = login
        if owners:
            out['problems'].append(f'{owner_rel} is on the base already; a pull request may '
                                   'not change it')
        elif author and login.casefold() != author.casefold():
            out['problems'].append(f'the name {match["name"]} belongs to the GitHub account '
                                   f'{login}; this pull request comes from {author}')
        else:
            out['known'] = bool(author)
    elif owners:
        out['first'] = True
        if owners[0][0] != 'A':
            out['problems'].append(f'{owner_rel} is not added')
        else:
            new = _load_json(_git_blob(repo, head, owner_rel))
            out['problems'] += check_owner(new, match['name'], match['set'], author)
            out['owner'] = (new or {}).get('login') if isinstance(new, dict) else None
    else:
        out['first'] = True
        out['problems'].append(
            f'first pull request of the name {match["name"]}: it must add {owner_rel} with '
            'your GitHub login ("delfin review submit" does it, or give --github <login>)')
    out['automerge'] = out['known'] and not out['problems']
    return out


def format_check(res: dict[str, Any]) -> str:
    if not res['helper']:
        return ("Not a helper's pull request (no ratings/<set>/review_<name>.json, "
                'claims/<name>.json or owners/<name>.json changed): left to the ordinary review.')
    if not res['problems']:
        if res['known']:
            who = f'from {res["owner"]}, the account the name belongs to; may merge by itself'
        elif res['first']:
            who = (f'first contribution under this name (binds it to {res["owner"]}): '
                   'maintainer review')
        else:
            who = 'account not checked'
        return f"Helper's pull request ({res['kind']}): {res['path']} -- all checks pass; {who}."
    lines = [f"Helper's pull request {res['path'] or ''}: {len(res['problems'])} problem(s)"]
    lines += [f'- {p}' for p in res['problems'][:50]]
    if len(res['problems']) > 50:
        lines.append(f'- ... and {len(res["problems"]) - 50} more')
    return '\n'.join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog='python -m delfin.review_check',
                                 description="Check a helper's pull request into ratings/.")
    ap.add_argument('--repo', default='.', help='the checkout with both commits')
    ap.add_argument('--base', required=True, help='the base commit (e.g. the target branch)')
    ap.add_argument('--head', required=True, help='the pull request commit')
    ap.add_argument('--author', default=None,
                    help='the GitHub login the pull request comes from (the workflow passes '
                         'github.event.pull_request.user.login)')
    args = ap.parse_args(list(argv) if argv is not None else None)
    try:
        res = check_pull_request(Path(args.repo), args.base, args.head,
                                 author=(args.author or '').strip() or None)
    except Problem as exc:
        print(f'review_check: {exc}')
        return 2
    text = format_check(res)
    print(text)
    ok = res['helper'] and not res['problems']

    def flag(v):
        return 'true' if v else 'false'

    outputs = (f'helper={flag(res["helper"])}\nvalid={flag(ok)}\n'
               f'known={flag(res["known"])}\nfirst={flag(res["first"])}\n'
               f'automerge={flag(res["automerge"])}\n')
    for var, body in (('GITHUB_STEP_SUMMARY', text + '\n'), ('GITHUB_OUTPUT', outputs)):
        target = os.environ.get(var)
        if target:
            with open(target, 'a', encoding='utf-8') as handle:
                handle.write(body)
    return 1 if res['problems'] else 0


if __name__ == '__main__':
    sys.exit(main())
