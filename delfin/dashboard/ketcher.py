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
           'latest_release', 'is_installed', 'smiles_from_molfile',
           'reaction_smiles_from_rxnfile', 'smiles_from_drawing',
           'DRAWINGS_FOLDER', 'DRAWING_SUFFIXES', 'drawings_directory',
           'is_drawing', 'list_drawings', 'save_drawing', 'read_drawing',
           'delete_drawing', 'frame_html', 'focus_js', 'load_js']

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
    """Where the browser loads the editor from, which is under the served root.

    That root is whatever the dashboard was started in, so this path moves
    with the launch directory -- which is why a Ketcher that had been fetched
    was gone the next time somebody started the dashboard from somewhere else,
    and why one fetched into the cache directory disappeared when /tmp was
    swept.  It is where the editor has to be *served* from; it is not where it
    should be *kept*.  See :func:`stored_directory`.
    """
    root = _voila_root()
    return None if root is None else root / '.delfin' / 'ketcher'


def stored_directory() -> Path:
    """Where the editor is kept, which does not move.

    One place per user, beside everything else DELFIN installs for them.  The
    served copy is made from this one, and making it is a local file copy of
    thirty megabytes -- under a second -- against fetching the same thirty
    over the network again.
    """
    return Path.home() / '.delfin' / 'ketcher'


def _version_in(folder: Optional[Path]) -> Optional[str]:
    if folder is None:
        return None
    stamp, page = folder / _STAMP, folder / 'index.html'
    if not (stamp.is_file() and page.is_file()):
        return None
    try:
        return stamp.read_text(encoding='utf-8').strip() or None
    except OSError:
        return None


def place_from_store() -> Optional[str]:
    """Put the kept copy where the browser can load it, if it is not there.

    Returns the version now being served, or None when there is nothing to
    place.  Copied rather than linked: a server that refuses to follow a
    symlink out of the directory it serves is doing the right thing, and
    thirty megabytes of local copy is not worth arguing with it about.
    """
    served = app_directory()
    if served is None:
        return None
    there = _version_in(served)
    if there:
        # Already being served, and if it is not kept anywhere yet then this
        # is the copy that was fetched before there was a place to keep it.
        # Taken in rather than left to be lost with the next launch directory.
        if not _version_in(stored_directory()):
            try:
                stored_directory().parent.mkdir(parents=True, exist_ok=True)
                shutil.copytree(served, stored_directory())
            except OSError:
                pass
        return there
    kept = stored_directory()
    if not _version_in(kept):
        return None
    try:
        served.parent.mkdir(parents=True, exist_ok=True)
        shutil.rmtree(served, ignore_errors=True)
        shutil.copytree(kept, served)
    except OSError:
        return None
    return _version_in(served)


def installed_version() -> Optional[str]:
    """Which build is there, placing the kept copy first if it has to.

    Asking used to mean asking the served directory alone, so an editor that
    had been fetched once looked absent the moment the dashboard was started
    somewhere else -- and the answer to that was to fetch thirty megabytes
    again, into a place that would lose it just as easily.
    """
    # Always through place_from_store: it is the one that both serves what is
    # kept and takes in what is only served, and short-circuiting on the
    # served copy meant an existing install was never taken in at all.
    return place_from_store()


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
        # Kept where it will still be next time, and served from a copy.
        kept = stored_directory()
        shutil.rmtree(kept, ignore_errors=True)
        kept.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(unpacked), str(kept))
        shutil.rmtree(folder, ignore_errors=True)
        folder.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(kept, folder)
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


def _charge_separate_dative(mol: Any) -> int:
    """Turn dative bonds into the charge-separated form DELFIN works in.

    A coordination bond can be written two ways.  RDKit writes it as an arrow,
    ``[n]->[Cd]``; the structures DELFIN is built around write it as a plain
    bond with the charge shown on both ends -- the donor as ``[N+]`` and the
    metal carrying one minus for each bond it accepts, ``[Cd-3]`` for three.
    Both describe the same molecule and only the second one is understood by
    everything downstream here, so a drawing is brought into it.

    A drawn structure arrives as arrows either way: a dative bond drawn as one
    is read as dative, and a plain bond to a metal that would leave the donor
    hypervalent is read as dative too.  So this is the one conversion that
    turns what a chemist draws into what the rest of DELFIN reads.

    Returns how many bonds were converted.
    """
    from rdkit import Chem

    editable = Chem.RWMol(mol)
    moved = 0
    for bond in editable.GetBonds():
        if bond.GetBondType() != Chem.BondType.DATIVE:
            continue
        # The arrow points from the donor to what accepts it.
        donor = bond.GetBeginAtom()
        acceptor = bond.GetEndAtom()
        bond.SetBondType(Chem.BondType.SINGLE)
        donor.SetFormalCharge(donor.GetFormalCharge() + 1)
        acceptor.SetFormalCharge(acceptor.GetFormalCharge() - 1)
        # The hydrogens the drawing showed are the hydrogens meant.  They have
        # to be written down before the count is frozen, or an ammonia donor
        # comes back as a bare nitrogen: [NH3]->[Pt] became [N+][Pt-].
        donor.SetNumExplicitHs(donor.GetTotalNumHs())
        donor.SetNoImplicit(True)
        moved += 1
    if not moved:
        return 0
    mol.__init__(editable.GetMol())
    return moved


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
    dative = 0
    try:
        dative = _charge_separate_dative(mol)
        if dative:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL
                             ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    except Exception:
        pass
    try:
        smiles = Chem.MolToSmiles(mol)
    except Exception as exc:
        return {'ok': False, 'smiles': '',
                'status': f'That drawing could not be written as SMILES: {exc}'}
    if not smiles:
        return {'ok': False, 'smiles': '',
                'status': 'That drawing came out as an empty SMILES.'}
    said = f'{mol.GetNumAtoms()} atoms drawn: {smiles}'
    if dative:
        said += (f' ({dative} coordination bond(s) written with the charge on '
                 'both ends, which is the form the rest of DELFIN reads.)')
    return {'ok': True, 'smiles': smiles, 'dative': dative, 'status': said}


def reaction_smiles_from_rxnfile(rxnfile: str) -> Dict[str, Any]:
    """Turn a drawn reaction into a reaction SMILES, the way RDKit writes one.

    An arrow on the canvas makes the drawing a reaction, and a reaction is not
    a molfile: asked for one, Ketcher refuses with "The structure cannot be
    saved as *.MOL due to reaction".  So a drawing with an arrow is fetched as
    an RXN file instead, and this is what reads it.

    The result is the ordinary form, ``reactants>agents>products`` with the
    sides joined by dots -- ``CCO.CC(=O)O>>CCOC(C)=O`` -- which is what RDKit
    writes and what the reaction SMARTS elsewhere in DELFIN already read.

    The reaction is rebuilt rather than edited in place.  Each side is read,
    tidied and added to a new reaction, because the molecules a parsed reaction
    hands out are its own and replacing one of them under it is not something
    to rely on.
    """
    text = str(rxnfile or '')
    if not text.strip():
        return {'ok': False, 'smiles': '', 'status': 'Nothing was drawn yet.'}
    try:
        from rdkit import RDLogger
        from rdkit.Chem import rdChemReactions

        RDLogger.DisableLog('rdApp.*')
    except ImportError:
        return {'ok': False, 'smiles': '',
                'status': 'RDKit is not installed, so a drawing cannot be read.'}
    # sanitize is off by default here, which is the right default for a
    # drawing: a reaction between two things a valence table dislikes is still
    # a reaction, and refusing it would refuse the reason for drawing it.
    try:
        drawn = rdChemReactions.ReactionFromRxnBlock(text)
    except Exception:                                   # noqa: BLE001
        drawn = None
    if drawn is None:
        return {'ok': False, 'smiles': '',
                'status': 'That drawing could not be read as a reaction.'}

    made = rdChemReactions.ChemicalReaction()
    dative = atoms = 0
    for taken, add in ((drawn.GetReactants(), made.AddReactantTemplate),
                       (drawn.GetAgents(), made.AddAgentTemplate),
                       (drawn.GetProducts(), made.AddProductTemplate)):
        for mol in taken:
            fresh = _tidied(mol)
            dative += fresh[1]
            atoms += fresh[0].GetNumAtoms()
            add(fresh[0])
    if atoms == 0:
        return {'ok': False, 'smiles': '', 'status': 'The drawing is empty.'}
    if not made.GetNumProductTemplates():
        return {'ok': False, 'smiles': '',
                'status': ('The arrow has nothing after it yet, so there is no '
                           'reaction to write.')}
    try:
        smiles = rdChemReactions.ReactionToSmiles(made)
    except Exception as exc:                            # noqa: BLE001
        return {'ok': False, 'smiles': '',
                'status': f'That reaction could not be written as SMILES: {exc}'}
    if not smiles:
        return {'ok': False, 'smiles': '',
                'status': 'That reaction came out as an empty SMILES.'}
    said = (f'{made.GetNumReactantTemplates()} drawn into '
            f'{made.GetNumProductTemplates()}: {smiles}')
    if dative:
        said += (f' ({dative} coordination bond(s) written with the charge on '
                 'both ends, which is the form the rest of DELFIN reads.)')
    return {'ok': True, 'smiles': smiles, 'dative': dative, 'status': said}


def _tidied(mol: Any):
    """A copy of *mol* read as far as it can be, with its arrows charged.

    Returns ``(molecule, converted)``.  A full sanitisation first, because that
    is what gives the clean SMILES; the partial one that :func:`smiles_from_molfile`
    falls back to when it fails, because a drawn metal complex frequently does
    fail it and is still the structure the chemist means.
    """
    from rdkit import Chem

    fresh = Chem.RWMol(mol).GetMol()
    try:
        Chem.SanitizeMol(fresh)
    except Exception:                                   # noqa: BLE001
        try:
            Chem.SanitizeMol(fresh, Chem.SanitizeFlags.SANITIZE_ALL
                             ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:                               # noqa: BLE001
            pass
    moved = 0
    try:
        moved = _charge_separate_dative(fresh)
        if moved:
            Chem.SanitizeMol(fresh, Chem.SanitizeFlags.SANITIZE_ALL
                             ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    except Exception:                                   # noqa: BLE001
        pass
    return fresh, moved


def smiles_from_drawing(payload: str) -> Dict[str, Any]:
    """Whatever the editor handed back, read as the thing it actually is.

    One call site for both: the editor is asked for an RXN file when there is
    an arrow on the canvas and for a molfile when there is not, and an RXN file
    says so on its first line.
    """
    text = str(payload or '')
    if text.lstrip().startswith('$RXN'):
        outcome = reaction_smiles_from_rxnfile(text)
        outcome['reaction'] = True
        return outcome
    outcome = smiles_from_molfile(text)
    outcome['reaction'] = False
    return outcome


# ---------------------------------------------------------------------------
# Keeping a drawing, and opening it again
# ---------------------------------------------------------------------------

#: The folder drawings are kept in, inside the calculation directory the
#: Calculations tab already shows.  A drawing belongs with the work it was
#: drawn for, not in a store of its own that nothing else can see.
DRAWINGS_FOLDER = 'Ketcher'

#: What can be written and read back.  ``.ket`` first and by default: it is
#: Ketcher's own format and the only one of these that keeps an arrow, a text
#: label and the layout, so a drawing opened again is the drawing that was
#: saved.  The others are for handing a structure to something that is not
#: Ketcher.
DRAWING_SUFFIXES = ('.ket', '.mol', '.rxn', '.smi', '.cdxml')


def drawings_directory(calc_dir: Any) -> Path:
    """Where drawings are kept for this calculation directory."""
    return Path(calc_dir) / DRAWINGS_FOLDER


def is_drawing(path: Any) -> bool:
    """Whether the editor can open this file."""
    return Path(path).suffix.lower() in DRAWING_SUFFIXES


def list_drawings(calc_dir: Any) -> list:
    """Every drawing kept here, by name, or nothing when there is no folder."""
    try:
        found = [item for item in drawings_directory(calc_dir).iterdir()
                 if item.is_file() and is_drawing(item)]
    except OSError:
        return []
    return sorted(found, key=lambda item: item.name.lower())


def save_drawing(calc_dir: Any, name: str, text: str,
                 suffix: str = '.ket') -> Dict[str, Any]:
    """Keep what was drawn, under the name it was given.

    The name is made safe the same way a reaction graph's folder name is, so
    ``../../etc/passwd`` becomes ``etc_passwd`` and stays in the folder; the
    resolved path is checked against the folder afterwards as well, because a
    rule that is only applied is a rule that can be got round.
    """
    from .reaction_graph import safe_name

    wanted = str(suffix or '.ket').strip().lower()
    if wanted and not wanted.startswith('.'):
        wanted = '.' + wanted
    if wanted not in DRAWING_SUFFIXES:
        return {'ok': False, 'path': None,
                'status': f'{wanted or "That"} is not a format the editor writes.'}
    body = str(text or '')
    if not body.strip():
        return {'ok': False, 'path': None,
                'status': 'There is nothing drawn to save.'}
    raw = str(name or '').strip()
    # Typing "acetone.ket" into the name box means the drawing is called
    # acetone, not acetone.ket.ket.
    if raw.lower().endswith(wanted):
        raw = raw[:-len(wanted)].strip()
    if not raw:
        return {'ok': False, 'path': None,
                'status': 'Give the drawing a name first.'}
    folder = drawings_directory(calc_dir)
    target = folder / f'{safe_name(raw)}{wanted}'
    try:
        folder.mkdir(parents=True, exist_ok=True)
        if target.resolve().parent != folder.resolve():
            return {'ok': False, 'path': None,
                    'status': 'That name would leave the Ketcher folder.'}
        target.write_text(body, encoding='utf-8')
    except OSError as exc:
        return {'ok': False, 'path': None,
                'status': f'It could not be saved: {exc}'}
    return {'ok': True, 'path': target,
            'status': f'Saved as {target.name} in {DRAWINGS_FOLDER}.'}


def read_drawing(path: Any) -> Dict[str, Any]:
    """A kept drawing, as the text the editor takes back."""
    where = Path(path)
    if not where.is_file():
        return {'ok': False, 'text': '', 'name': where.name,
                'status': f'{where.name} is not there any more.'}
    if not is_drawing(where):
        return {'ok': False, 'text': '', 'name': where.name,
                'status': f'{where.name} is not something the editor opens.'}
    try:
        text = where.read_text(encoding='utf-8')
    except (OSError, UnicodeDecodeError) as exc:
        return {'ok': False, 'text': '', 'name': where.name,
                'status': f'{where.name} could not be read: {exc}'}
    if not text.strip():
        return {'ok': False, 'text': '', 'name': where.name,
                'status': f'{where.name} is empty.'}
    return {'ok': True, 'text': text, 'name': where.name,
            'status': f'{where.name} is in the editor.'}


def delete_drawing(path: Any) -> Dict[str, Any]:
    """Throw a kept drawing away."""
    where = Path(path)
    if not is_drawing(where):
        return {'ok': False,
                'status': f'{where.name} is not a drawing, so it is left alone.'}
    try:
        where.unlink()
    except FileNotFoundError:
        return {'ok': False, 'status': f'{where.name} was already gone.'}
    except OSError as exc:
        return {'ok': False, 'status': f'{where.name} could not be deleted: {exc}'}
    return {'ok': True, 'status': f'{where.name} is gone.'}


# ---------------------------------------------------------------------------
# The frame the editor lives in, and reaching into it
# ---------------------------------------------------------------------------

def frame_html(url: str, *, height: str = '560px') -> str:
    """The editor in a frame of its own, reachable by keyboard and clipboard.

    A frame rather than the page: Ketcher is a React application that owns its
    own globals, and the dashboard is another one.  Same origin, so the page
    may reach in and ask it for the drawing -- across origins there would be
    nothing to ask with, because Ketcher speaks no messages.

    ``tabindex`` so the frame can hold the focus at all, and ``allow`` so the
    clipboard permission is handed into it.  Same-origin frames are already
    inside the default allowlist for those two, so this changes nothing today;
    it is what keeps Paste working if the editor is ever served from somewhere
    else.
    """
    import html as _html

    return (
        "<iframe src='" + _html.escape(str(url), quote=True) + "' "
        "tabindex='0' allow='clipboard-read; clipboard-write' "
        "style='width:100%; height:" + _html.escape(str(height), quote=True) +
        "; border:1px solid #d0d0d0; "
        "border-radius:6px; background:#fff;' "
        "title='Ketcher'></iframe>"
    )


def focus_js(host_selector: str) -> str:
    """Give the frame the focus when the pointer is on it, and not before.

    Ketcher listens for keys, and for copy, cut and paste, on the document
    *inside* the frame.  Those handlers run only while that document holds the
    focus, so a dashboard that never hands it over swallows every hotkey and
    makes ``navigator.clipboard.read()`` throw "Document is not focused".  One
    cause with two names: the keyboard does nothing, and Paste does nothing.

    On the pointer rather than on load or on a timer.  A frame that takes the
    focus by itself a couple of seconds after it appears is what scrolls the
    page away under the user -- which is the whole reason the scroll hold in
    the structure editor exists -- and taking it while somebody is typing in a
    field of the dashboard would be a worse bug than the one being fixed.
    Moving the pointer onto the canvas is the moment the editor is meant.

    The timer here looks for the element, not for a moment to steal focus: a
    widget's value is set from the kernel and the DOM catches up afterwards, so
    the host may not be there yet when this runs.
    """
    return (
        "(function(){\n"
        "  var tries=0;\n"
        "  function reach(){\n"
        "    var host=document.querySelector(" + json.dumps(host_selector) + ");\n"
        "    var frame=host&&host.querySelector('iframe');\n"
        "    try{ if(frame&&frame.contentWindow) frame.contentWindow.focus(); }\n"
        "    catch(e){}\n"
        "  }\n"
        "  function bind(){\n"
        "    var host=document.querySelector(" + json.dumps(host_selector) + ");\n"
        "    if(!host){ if(++tries<40) window.setTimeout(bind,150); return; }\n"
        "    if(host.__delfinKetcherFocus) return;\n"
        "    host.__delfinKetcherFocus=true;\n"
        "    host.addEventListener('pointerenter',reach);\n"
        "    host.addEventListener('pointerdown',reach,true);\n"
        "  }\n"
        "  bind();\n"
        "})();"
    )


def load_js(host_selector: str, text: str) -> str:
    """Put a kept drawing back into the editor.

    ``setMolecule`` and not one call per format: Indigo reads what it is given,
    so the same line takes a .ket, a .mol, a .rxn or a SMILES.

    It waits for the editor rather than giving up on it.  The frame may have
    only just been made, and Ketcher is a 30 MB application: the global this
    reaches for does not exist for the first second or two after the frame is
    put on the page, and a file opened in that window used to vanish.
    """
    return (
        "(function(){\n"
        "  var tries=0;\n"
        "  function put(){\n"
        "    var host=document.querySelector(" + json.dumps(host_selector) + ");\n"
        "    var frame=host&&host.querySelector('iframe');\n"
        "    var api=null;\n"
        "    try{ api=frame&&frame.contentWindow&&frame.contentWindow.ketcher; }\n"
        "    catch(e){ api=null; }\n"
        "    if(!api){ if(++tries<60) window.setTimeout(put,200); return; }\n"
        "    try{ Promise.resolve(api.setMolecule(" + json.dumps(str(text or '')) +
        ")); }catch(e){}\n"
        "  }\n"
        "  put();\n"
        "})();"
    )
