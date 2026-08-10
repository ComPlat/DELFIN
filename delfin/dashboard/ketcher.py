"""Ketcher in the dashboard: draw a structure, get a SMILES back.

The Submit tab takes a SMILES or a block of coordinates.  Typing a SMILES is
fine for ethanol and hopeless for anything a chemist would actually draw, so
the editor everyone already knows is offered beside the box: draw it, press
Convert, and the SMILES lands in the field.

Ketcher is a browser application, not a Python one, and a large one -- the
standalone build is a single 29.5 MB bundle because the Indigo chemistry engine
is compiled into it.  That is too much to carry in this repository, so it is
fetched when it is first wanted, the way xtb and g-xTB are: an offer, a
confirmation, and then it is there and works without a network.

Four of those bundles ship in the archive, one per demo page (index, popup,
duo, closable), and they are near enough identical.  Only the one the main page
loads is kept, which is the difference between 115 MB on disk and 32.

It has to be reachable by the browser, which means a URL rather than a path.
Voila serves everything below its root directory at ``/voila/files/``; that is
the same route the literature tab already uses for PDFs, and it goes through
whatever tunnel the dashboard itself came down.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Callable, Dict, Optional

__all__ = ['app_directory', 'app_url', 'install', 'installed_version',
           'latest_release', 'is_installed', 'smiles_from_molfile']

#: Where the releases come from.  The versioned asset rather than the
#: unversioned one beside it: two files of the same content and only one of
#: them says which build it is.
RELEASES = 'https://api.github.com/repos/epam/ketcher/releases/latest'
_ASSET = re.compile(r'^ketcher-standalone-(\d[\w.]*)\.zip$')

#: Written beside the app so an installed copy can say which build it is.
_STAMP = '.delfin-ketcher-version'


def _voila_root() -> Optional[Path]:
    """The directory Voila serves, or None when nothing is being served.

    Everything the browser is to load has to live under it -- an absolute path
    on the machine the kernel runs on means nothing to a browser that may be
    on the other side of an SSH tunnel.
    """
    raw = os.environ.get('DELFIN_VOILA_ROOT_DIR', '').strip()
    if not raw:
        return None
    try:
        return Path(raw).expanduser().resolve()
    except OSError:
        return None


def app_directory() -> Optional[Path]:
    """Where the drawing editor lives, once it has been fetched."""
    root = _voila_root()
    return None if root is None else root / '.delfin' / 'ketcher'


def installed_version() -> Optional[str]:
    """Which build is installed, or None if none is."""
    folder = app_directory()
    if folder is None:
        return None
    stamp = folder / _STAMP
    page = folder / 'index.html'
    if not (stamp.is_file() and page.is_file()):
        return None
    try:
        return stamp.read_text(encoding='utf-8').strip() or None
    except OSError:
        return None


def is_installed() -> bool:
    return installed_version() is not None


def app_url() -> Optional[str]:
    """The URL an iframe can load, or None if it is not installed.

    Same origin as the dashboard, so the page may reach into the editor and
    ask it for the structure; a cross-origin frame could only be talked to
    through messages, and Ketcher does not speak any.
    """
    folder = app_directory()
    root = _voila_root()
    if folder is None or root is None or not is_installed():
        return None
    from urllib.parse import quote

    rel = folder.relative_to(root).as_posix()
    # Cache-busted by version: an update that kept the same URL would be shown
    # from the browser cache, and the new build would appear to have changed
    # nothing.
    return (f'/voila/files/{quote(rel)}/index.html'
            f'?v={quote(str(installed_version()))}')


def latest_release(timeout: float = 20.0) -> Dict[str, Any]:
    """What the newest published build is, asked of the release feed.

    Returns ``{'ok', 'version', 'url', 'size', 'status'}``.  Nothing is pinned
    here: the newest is whatever is newest when asked, so keeping up to date is
    a button rather than an edit to this file.
    """
    request = urllib.request.Request(
        RELEASES, headers={'Accept': 'application/vnd.github+json',
                           'User-Agent': 'delfin-dashboard'})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as answer:
            payload = json.loads(answer.read().decode('utf-8'))
    except (urllib.error.URLError, ValueError, OSError) as exc:
        return {'ok': False, 'version': None, 'url': None, 'size': 0,
                'status': f'The release list could not be read: {exc}'}
    for asset in payload.get('assets') or []:
        found = _ASSET.match(str(asset.get('name') or ''))
        if found:
            return {'ok': True, 'version': found.group(1),
                    'url': asset.get('browser_download_url'),
                    'size': int(asset.get('size') or 0),
                    'status': f'Ketcher {found.group(1)} is the newest build.'}
    return {'ok': False, 'version': None, 'url': None, 'size': 0,
            'status': ('That release carries no standalone build; nothing to '
                       'install.')}


def _wanted(names: list, page: str) -> set:
    """Which files of the archive the main page actually needs.

    Four bundles ship, one per demo page, and each is 29.5 MB of the same
    chemistry engine compiled in.  Keeping the one the main page loads is the
    difference between 115 MB on disk and 32.
    """
    referenced = set(re.findall(r'(?:\./)?static/[\w./-]+', page))
    keep = set()
    for name in names:
        tail = name.split('standalone/', 1)[-1]
        # The folder entries themselves, and the resource forks a zip made on
        # a Mac carries beside every file.  Opening a folder entry as a file
        # is where this first went wrong.
        if not tail or tail.endswith('/') or '__MACOSX' in name:
            continue
        if Path(tail).name.startswith('._'):
            continue
        if not tail.startswith('static/'):
            # The other demo pages -- popup, duo, closable -- each load a
            # bundle of their own that is not being kept, so keeping the page
            # would leave a link into nothing.  Nobody opens them; index.html
            # is the editor.
            if tail.endswith('.html') and tail != 'index.html':
                continue
            keep.add(name)          # index.html, the icons, the manifest
            continue
        # every chunk, and the one main bundle the page names
        if tail.endswith('.chunk.js') or tail.endswith('.css'):
            keep.add(name)
        elif any(tail == ref.lstrip('./') for ref in referenced):
            keep.add(name)
    return keep


def install(
    on_line: Optional[Callable[[str], None]] = None,
    timeout: float = 900.0,
) -> Dict[str, Any]:
    """Fetch the newest build and put it where the browser can load it.

    Never called by itself: it is thirty-odd megabytes over the network, and
    on a machine that has no network it is a wait ending in nothing.  The
    caller asks first.
    """
    def say(text: str) -> None:
        if on_line is not None:
            try:
                on_line(text)
            except Exception:
                pass

    folder = app_directory()
    if folder is None:
        return {'ok': False, 'version': None,
                'status': ('Voila is not serving a directory, so a drawing '
                           'editor could not be reached by the browser even '
                           'once it was fetched.')}
    newest = latest_release()
    if not newest['ok']:
        return {'ok': False, 'version': None, 'status': newest['status']}

    started = time.perf_counter()
    say(f'fetching Ketcher {newest["version"]} '
        f'({newest["size"] / 1e6:.0f} MB)')
    staging = folder.parent / ('ketcher-download-%d' % os.getpid())
    archive = staging / 'ketcher.zip'
    try:
        staging.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(
            newest['url'], headers={'User-Agent': 'delfin-dashboard'})
        with urllib.request.urlopen(request, timeout=timeout) as answer:
            with open(archive, 'wb') as sink:
                shutil.copyfileobj(answer, sink)
        say('unpacking')
        with zipfile.ZipFile(archive) as bundle:
            names = bundle.namelist()
            page = ''
            for name in names:
                if name.endswith('standalone/index.html'):
                    page = bundle.read(name).decode('utf-8', 'replace')
                    break
            if not page:
                return {'ok': False, 'version': None,
                        'status': 'The archive has no standalone/index.html.'}
            keep = _wanted(names, page)
            unpacked = staging / 'app'
            for name in sorted(keep):
                tail = name.split('standalone/', 1)[-1]
                target = unpacked / tail
                if not tail or tail.endswith('/'):
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                with bundle.open(name) as source, open(target, 'wb') as sink:
                    shutil.copyfileobj(source, sink)
        (unpacked / _STAMP).write_text(newest['version'], encoding='utf-8')
        # In place at the end, not written over piece by piece: a half-updated
        # editor is worse than an old one.  shutil.move rather than replace --
        # renaming a directory onto a path refuses with EISDIR on some
        # filesystems, and this one is going into whatever the user launched
        # the dashboard in.
        shutil.rmtree(folder, ignore_errors=True)
        folder.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(unpacked), str(folder))
    except (OSError, zipfile.BadZipFile, urllib.error.URLError) as exc:
        return {'ok': False, 'version': None,
                'status': f'Ketcher could not be installed: {exc}'}
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    on_disk = sum(p.stat().st_size for p in folder.rglob('*') if p.is_file())
    return {
        'ok': True, 'version': newest['version'],
        'status': (f'Ketcher {newest["version"]} is ready '
                   f'({on_disk / 1e6:.0f} MB on disk, '
                   f'{time.perf_counter() - started:.0f} s).'),
    }


def smiles_from_molfile(molfile: str) -> Dict[str, Any]:
    """Turn what was drawn into a SMILES the rest of DELFIN will accept.

    The editor can write a SMILES itself, and it would be a different dialect:
    everything downstream here reads structures with RDKit, so the SMILES is
    made by RDKit from the drawing.  One that RDKit wrote is one RDKit will
    certainly read back.
    """
    # Not stripped: a molfile begins with three header lines and the first of
    # them is usually blank.  Taking the leading newline off shifts every line
    # up by one, the counts line is read as a title, and a perfectly good
    # drawing comes back unreadable.
    text = str(molfile or '')
    if not text.strip():
        return {'ok': False, 'smiles': '', 'status': 'Nothing was drawn yet.'}
    try:
        from rdkit import Chem
        from rdkit import RDLogger

        RDLogger.DisableLog('rdApp.*')
    except ImportError:
        return {'ok': False, 'smiles': '',
                'status': 'RDKit is not installed, so a drawing cannot be read.'}
    mol = Chem.MolFromMolBlock(text, sanitize=True, removeHs=False)
    if mol is None:
        # Drawn structures are often chemically loose -- an open valence, a
        # metal with a bond count no table agrees with.  Reading it without
        # sanitising says what is there rather than refusing it.
        mol = Chem.MolFromMolBlock(text, sanitize=False, removeHs=False)
        if mol is None:
            return {'ok': False, 'smiles': '',
                    'status': 'That drawing could not be read as a structure.'}
        try:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL
                             ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            pass
    if mol.GetNumAtoms() == 0:
        return {'ok': False, 'smiles': '', 'status': 'The drawing is empty.'}
    try:
        smiles = Chem.MolToSmiles(mol)
    except Exception as exc:
        return {'ok': False, 'smiles': '',
                'status': f'That drawing could not be written as SMILES: {exc}'}
    if not smiles:
        return {'ok': False, 'smiles': '',
                'status': 'That drawing came out as an empty SMILES.'}
    return {'ok': True, 'smiles': smiles,
            'status': f'{mol.GetNumAtoms()} atoms drawn: {smiles}'}
