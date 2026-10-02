"""Reviewing structures: rate every frame of multi-frame XYZ files by eye.

A review walks through the frames of one multi-frame XYZ file, or of every
``*.xyz`` file in a folder, and stores for each frame a verdict (``pass`` or
``block``), the categories of what is wrong, a free comment and a timestamp.
Everything that belongs together is kept together in one JSON file next to
the structures (``<folder>/review_<name>.json``):

* ``files`` -- one entry per structure file: its name, sha256, frame count,
  frame labels, a comment on the file as a whole and, when the folder (or
  its parent) holds an ``index.tsv`` with an ``id`` column, that row
  (SMILES and whatever else it lists), matched by the file stem.
* ``reviews`` -- ``{file: {frame_index: record}}``; a record names its file,
  frame, label and the sha256 it was given against, so a record stays
  readable on its own. A frame rated again keeps the earlier ratings in
  ``history``; nothing is overwritten away.
* ``order`` -- in blinded mode the shuffled order, so the anonymous number
  shown to the reviewer maps back to the file and frame.

The file is written after every rating (atomically, through a temporary file
and a rename), so an interrupted session loses nothing, and opening the same
folder again resumes where it stopped.  A record only counts as rated while
the file still has the sha256 it was rated against.

Findings from elsewhere can be shown next to a frame in the normal mode, read
from a JSON ``{"<file>": {"<frame_index>": ["text", ...]}}`` (by default
``findings.json`` in the folder).  The blinded mode never shows them.

``delfin review summary <review.json>`` prints the counts per verdict and per
category and writes the records as CSV.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import tempfile
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

DEFAULT_CATEGORIES: tuple[str, ...] = (
    'bond length',
    'angle/hybridisation',
    'contact/clash',
    'coordination polyhedron',
    'planarity/aromaticity',
    'torn/broken ligand',
    'wrong isomer',
    'other',
)

VERDICTS = ('pass', 'block')
SCHEMA_VERSION = 1
FINDINGS_NAME = 'findings.json'
INDEX_NAME = 'index.tsv'
CSV_FIELDS = (
    'file', 'frame', 'label', 'verdict', 'categories', 'note',
    'timestamp', 'sha256', 'blind_number', 'smiles', 'file_comment',
)


# ---------------------------------------------------------------------------
# Reading structures
# ---------------------------------------------------------------------------

@dataclass
class Frame:
    """One frame of a multi-frame XYZ file."""

    file: str            # file name, relative to the review folder
    index: int           # 0-based position of the frame in its file
    comment: str         # the comment line as written
    n_atoms: int
    coords: str          # the atom lines, without count and comment

    @property
    def key(self) -> tuple[str, int]:
        return (self.file, self.index)

    @property
    def label(self) -> str:
        """The frame label: the comment without a leading ``<id> frame<k>``."""
        parts = self.comment.split()
        if len(parts) >= 2 and parts[1].startswith('frame'):
            return ' '.join(parts[2:])
        return self.comment

    def xyz(self, comment: str | None = None) -> str:
        """The frame as a stand-alone XYZ block."""
        text = self.comment if comment is None else comment
        return f'{self.n_atoms}\n{text}\n{self.coords}\n'


def parse_frames(text: str, file: str = '') -> list[Frame]:
    """Split multi-frame XYZ text into frames ("N / comment / N atom lines")."""
    lines = text.splitlines()
    frames: list[Frame] = []
    i = 0
    while i < len(lines):
        head = lines[i].strip()
        if not head:
            i += 1
            continue
        try:
            n_atoms = int(head.split()[0])
        except ValueError:
            break
        body = lines[i + 2:i + 2 + n_atoms]
        if n_atoms <= 0 or len(body) < n_atoms:
            break
        comment = lines[i + 1].strip() if i + 1 < len(lines) else ''
        frames.append(Frame(file=file, index=len(frames), comment=comment,
                            n_atoms=n_atoms, coords='\n'.join(body)))
        i += 2 + n_atoms
    return frames


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read_index(folder: Path) -> dict[str, dict[str, str]]:
    """Rows of an ``index.tsv`` in *folder* or its parent, keyed by ``id``."""
    for candidate in (folder / INDEX_NAME, folder.parent / INDEX_NAME):
        if not candidate.is_file():
            continue
        try:
            with open(candidate, newline='', encoding='utf-8') as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
        except (OSError, csv.Error, UnicodeDecodeError):
            continue
        if rows and 'id' in rows[0]:
            return {row['id']: {k: v for k, v in row.items() if k and k != 'id'}
                    for row in rows if row.get('id')}
    return {}


def read_findings(path: Path | None) -> dict[str, dict[str, list[str]]]:
    """Read ``{"<file>": {"<frame>": ["text", ...]}}``; anything else is empty."""
    if path is None or not Path(path).is_file():
        return {}
    try:
        data = json.loads(Path(path).read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return {}
    if not isinstance(data, dict):
        return {}
    out: dict[str, dict[str, list[str]]] = {}
    for fname, frames in data.items():
        if not isinstance(frames, dict):
            continue
        out[str(fname)] = {
            str(k): [str(t) for t in (v if isinstance(v, list) else [v])]
            for k, v in frames.items()
        }
    return out


def blinded_order(keys: Sequence[tuple[str, int]], seed: int) -> list[tuple[str, int]]:
    """The keys shuffled by *seed*; the same keys and seed give the same order."""
    order = sorted(keys)
    random.Random(int(seed)).shuffle(order)
    return order


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def default_review_path(source: Path, name: str = '') -> Path:
    """``<folder>/review_<name>.json``; *name* defaults to the folder or file stem."""
    source = Path(source)
    folder = source if source.is_dir() else source.parent
    stem = (name or '').strip() or (source.name if source.is_dir() else source.stem)
    safe = ''.join(c if c.isalnum() or c in '-_.' else '_' for c in stem)
    return folder / f'review_{safe}.json'


# ---------------------------------------------------------------------------
# A review session
# ---------------------------------------------------------------------------

@dataclass
class ReviewSession:
    """The frames under review, the ratings so far and where the reader is."""

    source: Path
    review_path: Path
    frames: list[Frame]
    data: dict[str, Any]
    blind: bool = False
    seed: int = 0
    findings: dict[str, dict[str, list[str]]] = field(default_factory=dict)
    position: int = 0
    order: list[tuple[str, int]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    # -- opening -----------------------------------------------------------
    @classmethod
    def open(cls, source, *, name: str = '', review_path=None, blind: bool = False,
             seed: int = 0, findings_path=None,
             categories: Iterable[str] | None = None) -> ReviewSession:
        """Open a multi-frame XYZ file or a folder of them and resume its review."""
        source = Path(source).expanduser()
        if source.is_dir():
            folder = source
            paths = sorted(p for p in source.glob('*.xyz') if p.is_file())
        elif source.is_file():
            folder = source.parent
            paths = [source]
        else:
            raise FileNotFoundError(f'No such file or folder: {source}')
        review_path = Path(review_path) if review_path else default_review_path(source, name)
        data = load_review(review_path)
        index = read_index(folder)
        frames: list[Frame] = []
        warnings: list[str] = []
        files_meta = data.setdefault('files', {})
        for path in paths:
            rel = path.name
            try:
                text = path.read_text(encoding='utf-8', errors='replace')
            except OSError as exc:
                warnings.append(f'{rel}: not readable ({exc})')
                continue
            these = parse_frames(text, rel)
            if not these:
                warnings.append(f'{rel}: no XYZ frames')
                continue
            frames.extend(these)
            sha = file_sha256(path)
            meta = files_meta.setdefault(rel, {})
            if meta.get('sha256') and meta['sha256'] != sha:
                warnings.append(f'{rel}: changed since it was reviewed; '
                                'its earlier ratings are kept but no longer count')
            meta.update({'sha256': sha, 'n_frames': len(these),
                         'labels': [f.label for f in these]})
            meta.setdefault('comment', '')
            row = index.get(path.stem)
            if row:
                meta['index'] = row
        if categories is not None:
            data['categories'] = [c for c in (s.strip() for s in categories) if c]
        data.setdefault('categories', list(DEFAULT_CATEGORIES))
        data.setdefault('reviews', {})
        data['source'] = source.name
        data.setdefault('created', _now())
        if findings_path is None:
            findings_path = folder / FINDINGS_NAME
        session = cls(source=source, review_path=review_path, frames=frames,
                      data=data, findings=read_findings(findings_path),
                      warnings=warnings)
        session.set_mode(blind, seed)
        return session

    # -- order and mode -----------------------------------------------------
    def set_mode(self, blind: bool, seed: int = 0) -> None:
        """Switch between the file order and the blinded, shuffled order."""
        current = self.current.key if self.frames and self.order else None
        self.blind = bool(blind)
        self.seed = int(seed)
        keys = [f.key for f in self.frames]
        self.order = blinded_order(keys, self.seed) if self.blind else keys
        self._by_key = {f.key: f for f in self.frames}
        self._pos_of = {k: i for i, k in enumerate(self.order)}
        self.data['blind'] = self.blind
        if self.blind:
            self.data['seed'] = self.seed
            self.data['order'] = [[f, i] for f, i in self.order]
        if self.blind or current is None:
            self.position = 0
        else:
            self.position = self._pos_of.get(current, 0)

    def __len__(self) -> int:
        return len(self.order)

    @property
    def current(self) -> Frame:
        return self._by_key[self.order[self.position]]

    def frame_at(self, position: int) -> Frame:
        return self._by_key[self.order[position]]

    def display_name(self, position: int | None = None) -> str:
        """What the reader is told about a frame: anonymous when blinded."""
        pos = self.position if position is None else position
        if self.blind:
            return f'structure #{pos + 1}'
        frame = self.frame_at(pos)
        n = self.data['files'].get(frame.file, {}).get('n_frames', 0)
        return f'{frame.file}  frame {frame.index + 1}/{n}  {frame.label}'.rstrip()

    # -- navigation ----------------------------------------------------------
    def go(self, position: int) -> Frame:
        if self.order:
            self.position = max(0, min(int(position), len(self.order) - 1))
        return self.current

    def next(self) -> Frame:
        return self.go(self.position + 1)

    def prev(self) -> Frame:
        return self.go(self.position - 1)

    def next_unrated(self) -> Frame | None:
        """Move to the next frame without a rating, wrapping around; None when all are rated."""
        n = len(self.order)
        for step in range(1, n + 1):
            pos = (self.position + step) % n
            if self.record(*self.order[pos]) is None:
                return self.go(pos)
        return None

    def go_to_file(self, file: str) -> Frame | None:
        """Move to the first frame of *file* (only in the normal mode)."""
        for pos, (fname, _) in enumerate(self.order):
            if fname == file:
                return self.go(pos)
        return None

    # -- ratings -------------------------------------------------------------
    @property
    def categories(self) -> list[str]:
        return list(self.data.get('categories') or DEFAULT_CATEGORIES)

    def set_categories(self, categories: Iterable[str]) -> None:
        self.data['categories'] = [c for c in (s.strip() for s in categories) if c]
        self.save()

    def record(self, file: str, index: int) -> dict[str, Any] | None:
        """The current rating of a frame, or None (also when the file changed since)."""
        rec = self.data['reviews'].get(file, {}).get(str(index))
        if not rec or rec.get('verdict') not in VERDICTS:
            return None
        sha = self.data['files'].get(file, {}).get('sha256')
        if sha and rec.get('sha256') != sha:
            return None
        return rec

    def rate(self, verdict: str, categories: Iterable[str] = (), note: str = '',
             *, advance: bool = True) -> dict[str, Any]:
        """Rate the current frame, save at once and (by default) move on."""
        verdict = str(verdict).strip().lower()
        if verdict not in VERDICTS:
            raise ValueError(f'verdict must be one of {VERDICTS}, not {verdict!r}')
        frame = self.current
        meta = self.data['files'].get(frame.file, {})
        rec: dict[str, Any] = {
            'verdict': verdict,
            'categories': [c for c in categories if c],
            'note': str(note or '').strip(),
            'timestamp': _now(),
            'file': frame.file,
            'frame': frame.index,
            'label': frame.label,
            'sha256': meta.get('sha256', ''),
        }
        if self.blind:
            rec['blind_number'] = self.position + 1
            rec['seed'] = self.seed
        smiles = (meta.get('index') or {}).get('smiles')
        if smiles:
            rec['smiles'] = smiles
        shown = self.findings_for(frame)
        if shown and not self.blind:
            rec['findings_shown'] = shown
        per_file = self.data['reviews'].setdefault(frame.file, {})
        earlier = per_file.get(str(frame.index))
        if earlier:
            history = list(earlier.pop('history', []))
            history.append(earlier)
            rec['history'] = history
        per_file[str(frame.index)] = rec
        self.save()
        if advance:
            self.advance()
        return rec

    def advance(self) -> None:
        """After a rating: the next unrated frame ahead, else simply the next one."""
        n = len(self.order)
        for pos in range(self.position + 1, n):
            if self.record(*self.order[pos]) is None:
                self.go(pos)
                return
        if self.next_unrated() is None:
            self.next()

    def set_file_comment(self, file: str, comment: str) -> None:
        self.data['files'].setdefault(file, {})['comment'] = str(comment or '').strip()
        self.save()

    def set_session_comment(self, comment: str) -> None:
        self.data['comment'] = str(comment or '').strip()
        self.save()

    def findings_for(self, frame: Frame) -> list[str]:
        if self.blind:
            return []
        per_file = (self.findings.get(frame.file)
                    or self.findings.get(Path(frame.file).stem) or {})
        return list(per_file.get(str(frame.index), []))

    # -- progress, saving, export ---------------------------------------------
    def progress(self) -> dict[str, int]:
        counts = {'rated': 0, 'total': len(self.order), 'pass': 0, 'block': 0}
        for key in self.order:
            rec = self.record(*key)
            if rec:
                counts['rated'] += 1
                counts[rec['verdict']] += 1
        return counts

    def save(self) -> Path:
        self.data['version'] = SCHEMA_VERSION
        self.data['updated'] = _now()
        write_json_atomic(self.review_path, self.data)
        return self.review_path

    def export_csv(self, path=None) -> Path:
        path = Path(path) if path else self.review_path.with_suffix('.csv')
        write_csv(self.data, path)
        return path


# ---------------------------------------------------------------------------
# The review file
# ---------------------------------------------------------------------------

def load_review(path: Path) -> dict[str, Any]:
    """The stored review, or an empty one when there is none yet."""
    path = Path(path)
    if path.is_file():
        data = json.loads(path.read_text(encoding='utf-8'))
        if isinstance(data, dict):
            return data
        raise ValueError(f'{path} is not a review file')
    return {}


def write_json_atomic(path: Path, data: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f'.{path.name}.', suffix='.tmp', dir=path.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as handle:
            json.dump(data, handle, indent=1, ensure_ascii=False)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def iter_records(data: dict[str, Any]):
    """Every current rating, in file and frame order (stale ones included, flagged)."""
    files = data.get('files', {})
    for fname in sorted(data.get('reviews', {})):
        frames = data['reviews'][fname]
        for key in sorted(frames, key=lambda k: int(k) if str(k).isdigit() else 0):
            rec = frames[key]
            if rec.get('verdict') not in VERDICTS:
                continue
            sha = files.get(fname, {}).get('sha256')
            yield fname, int(key), rec, bool(sha and rec.get('sha256') != sha)


def write_csv(data: dict[str, Any], path: Path) -> None:
    files = data.get('files', {})
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS + ('stale',))
        writer.writeheader()
        for fname, idx, rec, stale in iter_records(data):
            meta = files.get(fname, {})
            writer.writerow({
                'file': fname,
                'frame': idx,
                'label': rec.get('label', ''),
                'verdict': rec.get('verdict', ''),
                'categories': ';'.join(rec.get('categories', [])),
                'note': rec.get('note', ''),
                'timestamp': rec.get('timestamp', ''),
                'sha256': rec.get('sha256', ''),
                'blind_number': rec.get('blind_number', ''),
                'smiles': rec.get('smiles') or (meta.get('index') or {}).get('smiles', ''),
                'file_comment': meta.get('comment', ''),
                'stale': 'yes' if stale else '',
            })


def summarize(data: dict[str, Any]) -> dict[str, Any]:
    """Counts per verdict and per category (overall and among blocked frames)."""
    out: dict[str, Any] = {'rated': 0, 'pass': 0, 'block': 0, 'stale': 0,
                           'total': sum(int(m.get('n_frames', 0))
                                        for m in data.get('files', {}).values()),
                           'categories': {}, 'block_categories': {}}
    for _, _, rec, stale in iter_records(data):
        if stale:
            out['stale'] += 1
            continue
        out['rated'] += 1
        out[rec['verdict']] += 1
        for cat in rec.get('categories', []):
            out['categories'][cat] = out['categories'].get(cat, 0) + 1
            if rec['verdict'] == 'block':
                out['block_categories'][cat] = out['block_categories'].get(cat, 0) + 1
    return out


def format_summary(summary: dict[str, Any]) -> str:
    lines = [
        f"rated {summary['rated']} / {summary['total']}   "
        f"pass {summary['pass']}   block {summary['block']}",
    ]
    if summary['stale']:
        lines.append(f"stale (file changed since rating): {summary['stale']}")
    if summary['categories']:
        width = max(len(c) for c in summary['categories'])
        lines.append('')
        lines.append(f"{'category'.ljust(width)}  all  block")
        for cat, n in sorted(summary['categories'].items(), key=lambda kv: (-kv[1], kv[0])):
            lines.append(f'{cat.ljust(width)}  {n:3d}  {summary["block_categories"].get(cat, 0):5d}')
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# CLI: delfin review summary <review.json>
# ---------------------------------------------------------------------------

def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog='delfin review',
                                     description='Evaluate a structure review file.')
    sub = parser.add_subparsers(dest='cmd')
    s = sub.add_parser('summary', help='Print counts per verdict and category and write CSV')
    s.add_argument('review', help='review_<name>.json written by the Review tab')
    s.add_argument('--csv', default=None,
                   help='CSV to write (default: next to the review file)')
    s.add_argument('--no-csv', action='store_true', help='Only print, write no CSV')
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.cmd != 'summary':
        parser.print_help()
        return 2
    path = Path(args.review)
    try:
        data = load_review(path)
    except (OSError, ValueError) as exc:
        print(f'Cannot read {path}: {exc}')
        return 1
    if not data:
        print(f'No review file at {path}')
        return 1
    print(format_summary(summarize(data)))
    if not args.no_csv:
        out = Path(args.csv) if args.csv else path.with_suffix('.csv')
        write_csv(data, out)
        print(f'\nCSV written: {out}')
    return 0


__all__ = [
    'DEFAULT_CATEGORIES', 'Frame', 'ReviewSession', 'blinded_order',
    'default_review_path', 'format_summary', 'load_review', 'main',
    'parse_frames', 'read_findings', 'read_index', 'summarize', 'write_csv',
]
