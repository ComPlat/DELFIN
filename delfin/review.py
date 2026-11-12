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
* ``notes`` -- ``{file: {frame_index: record}}`` for a comment typed on a
  frame that has no rating yet (comment, timestamp, file, frame, label,
  sha256; in blinded mode also the anonymous number).  It is not a rating
  (progress and summary ignore it), it is exported as a CSV row with an
  empty verdict, and a later rating of that frame takes it over.  Kept out
  of ``reviews`` so that readers which count every entry there as a
  rating, and the ``history`` of earlier ratings, stay as they were.

The file is written after every rating (atomically, through a temporary file
and a rename), so an interrupted session loses nothing, and opening the same
folder again resumes where it stopped.  A record only counts as rated while
the file still has the sha256 it was rated against.

Findings from elsewhere can be shown next to a frame in the normal mode, read
from a JSON ``{"<file>": {"<frame_index>": ["text", ...]}}`` (by default
``findings.json`` in the folder).  The blinded mode never shows them.

The categories of a new review can be given in a ``kategorien.txt`` next to
the structures (or in the folder above): one category per line, blank lines
and lines starting with ``#`` ignored (:func:`read_category_file`).  A
review stores its own list; categories added later are appended, so the
categories already in use keep their place (and their key in the tab).

A frame the reviewer is not sure about carries the word ``unsicher`` in its
comment (:func:`with_unsure`); readers of the review look for that word (or
``unsure``) in the comment, so the file format stays as it was.

``delfin review summary <review.json>`` prints the counts per verdict and per
category and writes the records as CSV.  ``delfin review assign | merge |
unpack | submit`` share a set out to a team of reviewers and bring their
ratings back together (:mod:`delfin.review_team`).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import re
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

# What a category code means, for the reviewer: a short label for the button and
# one sentence (the table of docs/REVIEW_TEAM_GUIDE.md).  The code stays what is
# stored; a category without an entry shows its code.
CATEGORY_INFO: dict[str, tuple[str, str]] = {
    'H': ('H atoms', 'a hydrogen is missing, stretched, at a wrong angle, points into '
          'the metal or sits on the metal'),
    'LAENGE_ML': ('M–L length', 'a metal-ligand or metal-metal distance visibly too '
                  'long or too short (the ligand is still bound)'),
    'DURCHDRINGUNG': ('Interpenetration', 'ligands or rings stuck through each other, '
                      'a bond threaded through a ring'),
    'CLASH': ('Clash', 'atoms too close without interpenetration (also H...H)'),
    'METALLKONTAKT': ('Metal contact', 'a heavy atom that is not a donor sits at the '
                      'metal (a hydrogen at the metal: H atoms)'),
    'TOPO': ('Topology', 'a bond missing or one too many, a ligand clearly torn off, a '
             'wrong donor atom, a fragment or counter-ion floating away or stuck '
             'inside another'),
    'HYB': ('Angle / hybrid.', 'a wrong angle, sp3 flat, sp2 pyramidal, a linear group '
            'bent (also one angle at a donor)'),
    'PLANAR': ('Planarity', 'an aromatic ring or sp2 framework not flat, an amide, '
               'ester or conjugated bond twisted'),
    'POLYEDER': ('Polyhedron', 'the coordination geometry at the metal is implausible'),
    'SONST': ('Other (comment)', 'something else: always write a comment'),
    'DONOR': ('Donor direction', 'a whole donor turned away from the metal, its lone '
              'pair not pointing at the metal, the metal out of the donor plane'),
    'LAENGE_KOV': ('Bond length', 'a covalent bond too long or too short'),
    'KOLLAPS': ('Collapse', 'a ring or ligand balled up, the molecule too compact'),
    'HAPTO': ('Hapticity', 'a pi ligand slipped, wrong hapticity'),
    'STEREO': ('Stereo', 'wrong configuration (R/S, E/Z, atropisomer)'),
    'KONFORMER': ('Conformer', 'an unrealistic conformation (eclipsed, a ring bent, a '
                  'boat without reason)'),
    'FORMEL': ('Formula', 'atoms missing or too many, an H on the wrong atom '
               '(tautomer, protonation)'),
    'SYMMETRIE': ('Symmetry', 'chemically equal parts visibly built unequal (equal '
                  'ligands at different distances, a Cp ring tilted)'),
}

# Where two categories are close (the guide's list), shown with the meanings.
CATEGORY_BOUNDARIES: tuple[str, ...] = (
    'stuck through each other = Interpenetration; only too close = Clash',
    'one angle at a donor = Angle / hybrid.; the whole donor turned away = Donor direction',
    'still bound but too long = M–L length; clearly gone = Topology',
    'a hydrogen at the metal = H atoms; a heavy atom at the metal = Metal contact',
)


def category_label(cat: str) -> str:
    """The button label of category *cat* (its code when it has no entry)."""
    return CATEGORY_INFO.get(cat, (cat, ''))[0]


def category_meaning(cat: str) -> str:
    """One sentence on what category *cat* means ('' when it has no entry)."""
    return CATEGORY_INFO.get(cat, (cat, ''))[1]


VERDICTS = ('pass', 'block')
SCHEMA_VERSION = 1
FINDINGS_NAME = 'findings.json'
INDEX_NAME = 'index.tsv'
CATEGORY_FILE_NAME = 'kategorien.txt'
UNSURE_WORD = 'unsicher'
_UNSURE_RE = re.compile(r'(?<!\w)(?:unsicher|unsure)(?!\w)', re.IGNORECASE)
CSV_FIELDS = (
    'file', 'frame', 'label', 'verdict', 'categories', 'note',
    'timestamp', 'sha256', 'blind_number', 'smiles', 'file_comment', 'atoms',
    'atoms_by_category',
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


def clean_categories(categories: Iterable[str]) -> list[str]:
    """Stripped, non-empty, each once, in the given order."""
    out: list[str] = []
    for c in categories:
        c = str(c).strip()
        if c and c not in out:
            out.append(c)
    return out


def read_category_file(source) -> list[str] | None:
    """The categories of a ``kategorien.txt`` beside *source* (or one folder up).

    *source* is the folder under review or a file in it.  One category per
    line; blank lines and lines starting with ``#`` are skipped.  None when
    there is no such file or it lists nothing.
    """
    source = Path(source).expanduser()
    folder = source if source.is_dir() else source.parent
    for candidate in (folder / CATEGORY_FILE_NAME, folder.parent / CATEGORY_FILE_NAME):
        if not candidate.is_file():
            continue
        try:
            lines = candidate.read_text(encoding='utf-8', errors='replace').splitlines()
        except OSError:
            continue
        cats = clean_categories(line for line in lines if not line.strip().startswith('#'))
        if cats:
            return cats
    return None


def is_unsure(note: str) -> bool:
    """True when the comment carries the word ``unsicher`` (or ``unsure``)."""
    return bool(_UNSURE_RE.search(note or ''))


def with_unsure(note: str, unsure: bool) -> str:
    """The comment with the word ``unsicher`` added (at the end) or removed.

    Removing takes ``unsure`` along as well, so that afterwards no reader of
    the comment counts the frame as unsure; the rest of the text is kept.
    """
    note = str(note or '')
    if unsure:
        if is_unsure(note):
            return note
        return f'{note.rstrip()} {UNSURE_WORD}'.strip()
    if not is_unsure(note):
        return note
    text = _UNSURE_RE.sub(' ', note)
    text = re.sub(r'[ \t]+', ' ', text)
    text = re.sub(r' *\n *', '\n', text)
    return text.strip(' \t\n,;')


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
    # The files this reviewer was given by a team assignment, in the order
    # they are shown (None: every file of the folder, in name order).
    assigned: list[str] | None = None

    # -- opening -----------------------------------------------------------
    @classmethod
    def open(cls, source, *, name: str = '', review_path=None, blind: bool = False,
             seed: int = 0, findings_path=None,
             categories: Iterable[str] | None = None,
             default_categories: Iterable[str] | None = None,
             files: Sequence[str] | None = None) -> ReviewSession:
        """Open a multi-frame XYZ file or a folder of them and resume its review.

        *categories* replaces the stored list; *default_categories* is used
        only by a review that has none stored yet, so a list set up before
        a new review is opened is kept and a resumed review keeps its own.
        *files* (a folder only) restricts the review to these file names,
        shown in this order: the share of one reviewer in a team
        (:mod:`delfin.review_team`).
        """
        source = Path(source).expanduser()
        missing: list[str] = []
        if source.is_dir() and files is not None:
            folder = source
            paths = []
            for fname in dict.fromkeys(str(f) for f in files):
                path = source / fname
                if path.is_file():
                    paths.append(path)
                else:
                    missing.append(fname)
        elif source.is_dir():
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
        if missing:
            warnings.append(f'{len(missing)} assigned file(s) not in the folder '
                            f'(first: {missing[0]})')
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
        if not data.get('categories') and default_categories is not None:
            data['categories'] = [c for c in (s.strip() for s in default_categories) if c]
        if not data.get('categories'):
            data['categories'] = list(DEFAULT_CATEGORIES)
        data.setdefault('reviews', {})
        data['source'] = source.name
        data.setdefault('created', _now())
        if findings_path is None:
            findings_path = folder / FINDINGS_NAME
        session = cls(source=source, review_path=review_path, frames=frames,
                      data=data, findings=read_findings(findings_path),
                      warnings=warnings,
                      assigned=([p.name for p in paths]
                                if files is not None and source.is_dir() else None))
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

    def has_ratings(self) -> bool:
        """True once any frame was rated (also against an older version of its file)."""
        return any(True for _ in iter_records(self.data))

    def append_categories(self, categories: Iterable[str]) -> list[str]:
        """Add the categories not yet listed at the end, save, and return them.

        The categories already listed keep their order, so a rating and the
        key of each category stay what they were.
        """
        have = self.categories
        new = [c for c in clean_categories(categories) if c not in have]
        if new:
            self.set_categories(have + new)
        return new

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
             *, atoms: Iterable[int] = (),
             atoms_by_category: dict[str, Iterable[int]] | None = None,
             advance: bool = True) -> dict[str, Any]:
        """Rate the current frame, save at once and (by default) move on.

        *atoms_by_category* are the atoms the reviewer marked as the place of
        each defect, *atoms* further ones without a category: 0-based indices
        into the frame's atom lines (hydrogens counted).  ``atoms`` in the
        record is the union of both.
        """
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
        by_cat = {str(c): sorted({int(i) for i in v})
                  for c, v in (atoms_by_category or {}).items() if v}
        marked = sorted({int(i) for i in atoms} | {i for v in by_cat.values() for i in v})
        if marked:
            rec['atoms'] = marked
        if by_cat:
            rec['atoms_by_category'] = by_cat
        if self.blind:
            rec['blind_number'] = self.position + 1
            rec['seed'] = self.seed
        smiles = (meta.get('index') or {}).get('smiles')
        if smiles:
            rec['smiles'] = smiles
        shown = self.findings_for(frame)
        if shown and not self.blind:
            rec['findings_shown'] = shown
        notes = self.data.get('notes') or {}
        if frame.file in notes:
            notes[frame.file].pop(str(frame.index), None)
            if not notes[frame.file]:
                notes.pop(frame.file)
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

    def pending_note(self, file: str, index: int) -> str:
        """The comment stored for an unrated frame ('' when none, or the file changed)."""
        rec = (self.data.get('notes') or {}).get(file, {}).get(str(index))
        if not rec:
            return ''
        sha = self.data['files'].get(file, {}).get('sha256')
        if sha and rec.get('sha256') != sha:
            return ''
        return str(rec.get('note', ''))

    def set_pending_note(self, note: str) -> None:
        """Store (or, when empty, drop) the comment on the current, unrated frame and save."""
        frame = self.current
        text = str(note or '').strip()
        notes = self.data.setdefault('notes', {})
        per_file = notes.setdefault(frame.file, {})
        if text:
            rec: dict[str, Any] = {
                'note': text, 'timestamp': _now(), 'file': frame.file,
                'frame': frame.index, 'label': frame.label,
                'sha256': self.data['files'].get(frame.file, {}).get('sha256', ''),
            }
            if self.blind:
                rec['blind_number'] = self.position + 1
                rec['seed'] = self.seed
            per_file[str(frame.index)] = rec
        else:
            per_file.pop(str(frame.index), None)
        if not per_file:
            notes.pop(frame.file, None)
        self.save()

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


def csv_rows(data: dict[str, Any]):
    """The rows :func:`write_csv` writes (``CSV_FIELDS`` + ``stale``), in file and frame order."""
    files = data.get('files', {})
    rows = list(iter_records(data))
    rated = {(f, i) for f, i, _, stale in rows if not stale}
    for fname, frames in (data.get('notes') or {}).items():
        sha = files.get(fname, {}).get('sha256')
        for key, rec in frames.items():
            if str(key).isdigit() and (fname, int(key)) not in rated:
                rows.append((fname, int(key), rec, bool(sha and rec.get('sha256') != sha)))
    rows.sort(key=lambda row: (row[0], row[1]))
    for fname, idx, rec, stale in rows:
        meta = files.get(fname, {})
        yield {
            'file': fname,
            'frame': idx,
            'label': rec.get('label', ''),
            'verdict': rec.get('verdict', ''),
            'categories': ';'.join(rec.get('categories') or []),
            'note': rec.get('note', ''),
            'timestamp': rec.get('timestamp', ''),
            'sha256': rec.get('sha256', ''),
            'blind_number': rec.get('blind_number', ''),
            'smiles': rec.get('smiles') or (meta.get('index') or {}).get('smiles', ''),
            'file_comment': meta.get('comment', ''),
            'atoms': ';'.join(str(i) for i in rec.get('atoms') or []),
            'atoms_by_category': ';'.join(
                f'{c}:' + ','.join(str(i) for i in v)
                for c, v in (rec.get('atoms_by_category') or {}).items()),
            'stale': 'yes' if stale else '',
        }


def write_csv(data: dict[str, Any], path: Path) -> None:
    """One row per rating, and one with an empty verdict per comment on an unrated frame."""
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS + ('stale',))
        writer.writeheader()
        for row in csv_rows(data):
            writer.writerow(row)


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
# CLI: delfin review summary <review.json>  (and assign | merge | unpack | submit)
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
    # Sharing a set out to a team of reviewers and bringing it back together.
    from delfin import review_team
    review_team.add_parsers(sub)
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.cmd in review_team.COMMANDS:
        return review_team.run(args)
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
    'CATEGORY_FILE_NAME', 'DEFAULT_CATEGORIES', 'Frame', 'ReviewSession', 'UNSURE_WORD',
    'blinded_order', 'clean_categories', 'csv_rows', 'default_review_path', 'format_summary',
    'is_unsure', 'load_review', 'main', 'parse_frames', 'read_category_file',
    'read_findings', 'read_index', 'summarize', 'with_unsure', 'write_csv',
]
