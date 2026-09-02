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
           'is_drawing', 'list_drawings', 'list_in', 'save_drawing',
           'save_into', 'read_drawing',
           'delete_drawing', 'frame_html', 'focus_js', 'load_js', 'KET_MARK',
           'files_js', 'wire_js']

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


#: What separates the two things a drawn reaction is fetched as.
#:
#: An RXN file is what RDKit reads, and it holds one arrow -- Indigo writes a
#: canvas with three molecules and two arrows as "one into two", because the
#: format has nowhere to put the second one.  The arrows themselves survive
#: only in Ketcher's own KET, so both come back and are read together.
KET_MARK = '\n$DELFIN-KET$\n'


def _arrows(ket: str) -> list:
    """Where the arrows are, as ``{x0, x1, y}`` sorted left to right.

    Everything about how a scheme is read comes off these: what is before an
    arrow and what is after it, and what is written over or under one.
    """
    try:
        document = json.loads(str(ket or ''))
        nodes = document['root']['nodes']
    except (ValueError, KeyError, TypeError):
        return []
    found = []
    for node in nodes:
        if not isinstance(node, dict) or node.get('type') != 'arrow':
            continue
        where = (node.get('data') or {}).get('pos') or []
        points = [point for point in where if isinstance(point, dict)
                  and point.get('x') is not None]
        if not points:
            continue
        xs = [point['x'] for point in points]
        ys = [point.get('y', 0.0) for point in points]
        found.append({'x0': min(xs), 'x1': max(xs),
                      'y': sum(ys) / len(ys)})
    return sorted(found, key=lambda one: one['x0'])


def _extent(mol: Any) -> Optional[Dict[str, float]]:
    """Where a component sits on the canvas, or None if it says nothing."""
    try:
        if not mol.GetNumConformers() or not mol.GetNumAtoms():
            return None
        frame = mol.GetConformer()
        spots = [frame.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
        xs = [spot.x for spot in spots]
        ys = [spot.y for spot in spots]
    except Exception:                                   # noqa: BLE001
        return None
    if not xs:
        return None
    return {'x0': min(xs), 'x1': max(xs), 'y0': min(ys), 'y1': max(ys),
            'cx': (min(xs) + max(xs)) / 2.0}


def _as_drawn(drawn: Any, arrows: list) -> Dict[str, Any]:
    """The scheme read off the canvas, which is where it is actually written.

    Four places mean four things around an arrow, and only two of them survive
    a trip through an RXN file:

    * **before** it -- what goes in;
    * **after** it -- what comes out;
    * **over** it -- what is added, the reagent or solvent;
    * **under** it -- what is given off, which comes out with the products.

    Indigo keeps neither the second arrow nor the difference between over and
    under.  A three-step scheme comes back as "the first thing, into
    everything else", not even in the drawn order, and a molecule under the
    arrow is handed back as an agent exactly like one over it -- measured:
    benzene, cyclobutane over, cyclopropane under, cyclohexane after came back
    as ``C1C=CC=CC=1>C1CCC1.C1CC1>C1CCCCC1``.

    What it does keep is every component and its coordinates, and the KET
    keeps the arrows and theirs in the same frame.  So the reading is done
    from the geometry instead, and the result is the ordinary
    ``reactants>agents>products`` -- with a further ``>agents>products`` for
    every arrow after the first.

    Over and under are told apart strictly: the component has to sit inside
    the arrow's span and clear of its line altogether.  Reactants and products
    straddle that line -- measured at y -8.18..-6.17 against an arrow at
    -7.18 -- so anything that merely reaches across it is not a reagent.
    """
    from rdkit import Chem

    befores: Dict[int, list] = {}
    overs: Dict[int, list] = {}
    unders: Dict[int, list] = {}
    dative = 0
    placed = 0
    for mol in list(drawn.GetReactants()) + list(drawn.GetAgents()) \
            + list(drawn.GetProducts()):
        where = _extent(mol)
        if where is None:
            continue
        fresh, moved = _tidied(mol)
        dative += moved
        placed += 1
        on = None
        for index, arrow in enumerate(arrows):
            if not (arrow['x0'] <= where['cx'] <= arrow['x1']):
                continue
            if where['y0'] > arrow['y']:                # clear of it, above
                on = ('over', index)
            elif where['y1'] < arrow['y']:              # clear of it, below
                on = ('under', index)
            break
        if on is None:
            step = sum(1 for arrow in arrows
                       if (arrow['x0'] + arrow['x1']) / 2.0 < where['cx'])
            befores.setdefault(step, []).append(fresh)
        elif on[0] == 'over':
            overs.setdefault(on[1], []).append(fresh)
        else:
            # Given off by that step.  Kept apart from what is drawn between
            # the arrows, because those two are not the same thing once there
            # is a second arrow: the water a step gives off is a product of
            # that step and not something the next one is made from.
            unders.setdefault(on[1], []).append(fresh)
    if placed == 0:
        return {'ok': False, 'smiles': '', 'status': 'The drawing is empty.'}

    def written(*groups):
        parts = [one for group in groups for one in (group or [])]
        return '.'.join(sorted(Chem.MolToSmiles(one) for one in parts))

    # One reaction per step, one to a line.  A scheme is a sequence of
    # reactions and not one long one: written as a single chain,
    # ``A>reagent>B.HCl>>C``, the field holding B.HCl is at once the products
    # of the first step and the reactants of the second, so the hydrogen
    # chloride the first step gives off is read as something the second one is
    # made from.  A line each says what each step is, and every one of them is
    # a reaction SMILES anything can read on its own.
    try:
        lines = []
        for index in range(len(arrows)):
            lines.append('>'.join([
                written(befores.get(index)),
                written(overs.get(index)),
                written(befores.get(index + 1), unders.get(index)),
            ]))
    except Exception as exc:                            # noqa: BLE001
        return {'ok': False, 'smiles': '',
                'status': f'That reaction could not be written as SMILES: {exc}'}
    if not lines or not lines[-1].rsplit('>', 1)[-1]:
        return {'ok': False, 'smiles': '',
                'status': ('The last arrow has nothing after it yet, so there '
                           'is no reaction to write.')}
    smiles = '\n'.join(lines)
    steps = len(arrows)
    said = (f'{steps} step{"" if steps == 1 else "s"} drawn: '
            + ' / '.join(lines))
    if dative:
        said += (f' ({dative} coordination bond(s) written with the charge on '
                 'both ends, which is the form the rest of DELFIN reads.)')
    return {'ok': True, 'smiles': smiles, 'dative': dative, 'steps': steps,
            'status': said}


def reaction_smiles_from_rxnfile(rxnfile: str, ket: str = '') -> Dict[str, Any]:
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

    arrows = _arrows(ket)
    if arrows:
        return _as_drawn(drawn, arrows)

    # No KET to read the canvas from, so the RXN file's own split is all there
    # is: one arrow, and whatever Indigo decided was an agent.
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
    body, _, ket = text.partition(KET_MARK)
    if body.lstrip().startswith('$RXN'):
        outcome = reaction_smiles_from_rxnfile(body, ket)
        outcome['reaction'] = True
    else:
        outcome = smiles_from_molfile(body)
        outcome['reaction'] = False
    # The drawing itself travels with the answer, so it can be kept beside the
    # job it was drawn for rather than fetched again from a frame that may by
    # then hold something else.
    outcome['ket'] = ket
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


def list_in(folder: Any) -> list:
    """Every drawing in *folder*, by name, or nothing when there is no folder.

    A folder rather than the store, because a drawing opened from a job's own
    directory is saved back into that directory.  The store is one particular
    folder and not a different kind of thing.
    """
    try:
        found = [item for item in Path(folder).iterdir()
                 if item.is_file() and is_drawing(item)]
    except OSError:
        return []
    return sorted(found, key=lambda item: item.name.lower())


def list_drawings(calc_dir: Any) -> list:
    """Every drawing kept in the Ketcher folder of *calc_dir*."""
    return list_in(drawings_directory(calc_dir))


def save_drawing(calc_dir: Any, name: str, text: str,
                 suffix: str = '.ket') -> Dict[str, Any]:
    """Keep what was drawn in the Ketcher folder of *calc_dir*."""
    return save_into(drawings_directory(calc_dir), name, text, suffix)


def save_into(folder: Any, name: str, text: str,
              suffix: str = '.ket') -> Dict[str, Any]:
    """Keep what was drawn in *folder*, under the name it was given.

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
    where = Path(folder)
    target = where / f'{safe_name(raw)}{wanted}'
    try:
        where.mkdir(parents=True, exist_ok=True)
        if target.resolve().parent != where.resolve():
            return {'ok': False, 'path': None,
                    'status': f'That name would leave {where.name}.'}
        target.write_text(body, encoding='utf-8')
    except OSError as exc:
        return {'ok': False, 'path': None,
                'status': f'It could not be saved: {exc}'}
    return {'ok': True, 'path': target,
            'status': f'Saved as {target.name} in {where.name}.'}


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

    # A little narrower than the box it sits in, and its border counted
    # inside its own width.  At a plain 100% the frame stood 4 px past its
    # parent -- measured: parent right edge 1475, frame right edge 1479 -- and
    # the right-hand end of Ketcher's toolbar was clipped by it.
    return (
        "<iframe src='" + _html.escape(str(url), quote=True) + "' "
        "tabindex='0' allow='clipboard-read; clipboard-write' "
        "style='display:block; box-sizing:border-box; "
        "width:calc(100% - 8px); max-width:100%; height:" +
        _html.escape(str(height), quote=True) +
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


def files_js(host_selector: str, names: list) -> str:
    """Tell the frame which drawings are kept, so its Open list can show them."""
    return (
        "(function(){\n"
        "  var host=document.querySelector(" + json.dumps(host_selector) + ");\n"
        "  var frame=host&&host.querySelector('iframe');\n"
        "  var w=null; try{ w=frame&&frame.contentWindow; }catch(e){ w=null; }\n"
        "  if(!w) return;\n"
        "  w.__delfinKetcherFiles=" + json.dumps([str(n) for n in names]) + ";\n"
        "  if(w.__delfinKetcherRedraw) w.__delfinKetcherRedraw();\n"
        "})();"
    )


def wire_js(host_selector: str, sync_selector: str) -> str:
    """Wire Ketcher's own Save and Open buttons to the calculation directory.

    The editor already has the two buttons a chemist looks for, with a name
    field and a format list behind Save.  A second pair beside the frame is a
    second answer to a question the editor had already answered, so these are
    the ones that are made to work.

    **Save.**  Ketcher writes a file the way every browser application does:
    it builds a Blob, makes an object URL for it, hangs that on a detached
    ``<a download>`` and dispatches a click at it.  Nothing of that reaches the
    page's own listeners -- the anchor is never in the document -- so the two
    ends are taken instead: ``URL.createObjectURL`` remembers which Blob each
    URL stands for, and ``dispatchEvent`` recognises the click when it comes.
    The name comes off the anchor, which is the name typed into Ketcher's own
    dialog, and the download never happens.

    Only the formats that can be opened again are taken.  Save as SVG or PNG
    is a picture for somewhere else and still downloads the way it always did.

    **Open.**  Ketcher's dialog offers the clipboard and a file picker, and a
    file picker shows the machine the *browser* is on -- which, down an SSH
    tunnel, is not the machine the calculations are on.  So the button is
    answered with the drawings that are actually kept, and the picker is still
    one click away for a file that is genuinely local.
    """
    return (
        "(function(){\n"
        "  var tries=0;\n"
        "  function hand(kind,payload){\n"
        "    var box=document.querySelector(" + json.dumps(sync_selector) + ");\n"
        "    var input=box&&box.querySelector('textarea, input');\n"
        "    if(!input) return;\n"
        "    var proto=(input.tagName==='TEXTAREA')\n"
        "      ? window.HTMLTextAreaElement.prototype\n"
        "      : window.HTMLInputElement.prototype;\n"
        "    var setter=Object.getOwnPropertyDescriptor(proto,'value');\n"
        "    var line=(Date.now())+'\\n'+kind+'\\n'+payload;\n"
        "    if(setter&&setter.set) setter.set.call(input,line);\n"
        "    else input.value=line;\n"
        "    input.dispatchEvent(new Event('input',{bubbles:true}));\n"
        "    input.dispatchEvent(new Event('change',{bubbles:true}));\n"
        "  }\n"
        "  var KEEP=" + json.dumps(list(DRAWING_SUFFIXES)) + ";\n"
        "  function wire(){\n"
        "    var host=document.querySelector(" + json.dumps(host_selector) + ");\n"
        "    var frame=host&&host.querySelector('iframe');\n"
        "    var w=null,d=null;\n"
        "    try{ w=frame&&frame.contentWindow; d=w&&w.document; }catch(e){ w=null; }\n"
        "    if(!w||!d||!w.ketcher){ if(++tries<80) window.setTimeout(wire,250); return; }\n"
        "    if(w.__delfinKetcherWired) return;\n"
        "    w.__delfinKetcherWired=true;\n"
        "\n"
        "    /* Save: the two ends of a download nothing else can see. */\n"
        "    var blobs=new Map();\n"
        "    var makeURL=w.URL.createObjectURL.bind(w.URL);\n"
        "    w.URL.createObjectURL=function(thing){\n"
        "      var url=makeURL(thing);\n"
        "      try{ if(thing instanceof w.Blob) blobs.set(url,thing); }catch(e){}\n"
        "      return url;\n"
        "    };\n"
        "    var fire=w.EventTarget.prototype.dispatchEvent;\n"
        "    w.EventTarget.prototype.dispatchEvent=function(evt){\n"
        "      try{\n"
        "        if(evt&&evt.type==='click'&&this instanceof w.HTMLAnchorElement\n"
        "           &&this.hasAttribute('download')&&blobs.has(this.href)){\n"
        "          var name=this.getAttribute('download')||'drawing';\n"
        "          var dot=name.lastIndexOf('.');\n"
        "          var ext=dot<0?'':name.slice(dot).toLowerCase();\n"
        "          if(KEEP.indexOf(ext)>=0){\n"
        "            var blob=blobs.get(this.href); blobs.delete(this.href);\n"
        "            blob.text().then(function(text){ hand('save',name+'\\n'+text); },\n"
        "                             function(err){ hand('save-failed',''+err); });\n"
        "            /* What was just saved is what there is no reason to warn\n"
        "               about losing. */\n"
        "            try{ Promise.resolve(w.ketcher.getKet()).then(\n"
        "              function(k){ w.__delfinKetcherClean=k; }); }catch(e){}\n"
        "            return true;\n"
        "          }\n"
        "        }\n"
        "      }catch(e){}\n"
        "      return fire.apply(this,arguments);\n"
        "    };\n"
        "\n"
        "    /* Save: opened on Ket, which is the only one of these formats\n"
        "       that keeps an arrow, a text label and the layout, so a drawing\n"
        "       saved and opened again is the drawing that was saved.  The\n"
        "       control is a Material select: its own hidden input can be set\n"
        "       without React noticing, so the list is opened and the entry is\n"
        "       pressed, which is what a person would do. */\n"
        "    d.addEventListener('click',function(ev){\n"
        "      var press=ev.target&&ev.target.closest\n"
        "        ? ev.target.closest('[data-testid=\"save-file-button\"]') : null;\n"
        "      if(!press) return;\n"
        "      var waited=0;\n"
        "      function pick(){\n"
        "        var box=d.querySelector('[data-testid=\"save-dialog\"]');\n"
        "        if(!box){ if(++waited<40) window.setTimeout(pick,100); return; }\n"
        "        var label=d.querySelector('label[data-testid=\"file-format-list\"]');\n"
        "        if(!label) return;\n"
        "        var now=label.querySelector('input.MuiSelect-nativeInput');\n"
        "        if(now&&now.value==='ket') return;\n"
        "        var face=label.querySelector('[role=\"combobox\"]');\n"
        "        if(!face) return;\n"
        "        face.dispatchEvent(new w.MouseEvent('mousedown',{bubbles:true}));\n"
        "        var tries=0;\n"
        "        function press_it(){\n"
        "          var all=d.querySelectorAll('[role=\"option\"]');\n"
        "          for(var i=0;i<all.length;i++){\n"
        "            if(/^ket/i.test((all[i].textContent||'').trim())){\n"
        "              all[i].click(); return;\n"
        "            }\n"
        "          }\n"
        "          if(++tries<30) window.setTimeout(press_it,100);\n"
        "        }\n"
        "        press_it();\n"
        "      }\n"
        "      pick();\n"
        "    },true);\n"
        "\n"
        "    /* Open: answered with what is kept, not with the browser's disk. */\n"
        "    var sheet=d.createElement('style');\n"
        "    sheet.textContent=''\n"
        "      +'.delfin-open{position:fixed;inset:0;z-index:2147483000;'\n"
        "      +'background:rgba(0,0,0,.35);display:flex;align-items:center;'\n"
        "      +'justify-content:center;font-family:Arial,sans-serif;}'\n"
        "      +'.delfin-open-box{background:#fff;border-radius:8px;padding:18px;'\n"
        "      +'min-width:380px;max-width:520px;max-height:70vh;overflow:auto;'\n"
        "      +'box-shadow:0 8px 32px rgba(0,0,0,.3);}'\n"
        "      +'.delfin-open h3{margin:0 0 4px;font-size:15px;color:#333;}'\n"
        "      +'.delfin-open p{margin:0 0 12px;font-size:12px;color:#777;}'\n"
        "      +'.delfin-open button.row{display:block;width:100%;text-align:left;'\n"
        "      +'padding:8px 10px;margin:2px 0;border:1px solid #e0e0e0;'\n"
        "      +'border-radius:4px;background:#fafafa;cursor:pointer;font-size:13px;}'\n"
        "      +'.delfin-open button.row:hover{background:#e8f0fe;border-color:#1976d2;}'\n"
        "      +'.delfin-open .feet{margin-top:14px;display:flex;gap:8px;'\n"
        "      +'justify-content:flex-end;}'\n"
        "      +'.delfin-open .feet button{padding:6px 12px;border-radius:4px;'\n"
        "      +'border:1px solid #bbb;background:#fff;cursor:pointer;font-size:12px;}';\n"
        "    d.head.appendChild(sheet);\n"
        "\n"
        "    var shown=null;\n"
        "    function shut(){ if(shown&&shown.parentNode) shown.parentNode.removeChild(shown);\n"
        "                     shown=null; }\n"
        "    function draw(){\n"
        "      shut();\n"
        "      var kept=w.__delfinKetcherFiles||[];\n"
        "      shown=d.createElement('div');\n"
        "      shown.className='delfin-open';\n"
        "      var box=d.createElement('div'); box.className='delfin-open-box';\n"
        "      var head=d.createElement('h3'); head.textContent='Open a drawing';\n"
        "      var note=d.createElement('p');\n"
        "      note.textContent=kept.length\n"
        "        ? 'Kept beside the calculations, in the Ketcher folder.'\n"
        "        : 'Nothing is kept in the Ketcher folder yet. Save one first.';\n"
        "      box.appendChild(head); box.appendChild(note);\n"
        "      kept.forEach(function(name){\n"
        "        var row=d.createElement('button');\n"
        "        row.className='row'; row.textContent=name;\n"
        "        row.onclick=function(){ shut(); hand('open',name); };\n"
        "        box.appendChild(row);\n"
        "      });\n"
        "      var feet=d.createElement('div'); feet.className='feet';\n"
        "      var local=d.createElement('button');\n"
        "      local.textContent='From this computer...';\n"
        "      local.onclick=function(){\n"
        "        shut();\n"
        "        /* Ketcher's own dialog, let through for one press. */\n"
        "        w.__delfinKetcherPass=true;\n"
        "        /* The one on screen. Ketcher carries two, one per mode, and\n"
        "           the first in the document is the hidden macromolecule\n"
        "           editor's -- clicking that one does nothing at all. */\n"
        "        var all=d.querySelectorAll('[data-testid=\"open-file-button\"]');\n"
        "        var real=null;\n"
        "        for(var i=0;i<all.length;i++){\n"
        "          if(all[i].offsetParent!==null){ real=all[i]; break; }\n"
        "        }\n"
        "        if(!real) real=all[0];\n"
        "        if(real) real.click();\n"
        "      };\n"
        "      var stop=d.createElement('button');\n"
        "      stop.textContent='Cancel'; stop.onclick=shut;\n"
        "      feet.appendChild(local); feet.appendChild(stop);\n"
        "      box.appendChild(feet);\n"
        "      shown.appendChild(box);\n"
        "      shown.onclick=function(ev){ if(ev.target===shown) shut(); };\n"
        "      d.body.appendChild(shown);\n"
        "    }\n"
        "    w.__delfinKetcherRedraw=function(){ if(shown) draw(); };\n"
        "    d.addEventListener('click',function(ev){\n"
        "      var button=ev.target&&ev.target.closest\n"
        "        ? ev.target.closest('[data-testid=\"open-file-button\"]') : null;\n"
        "      if(!button) return;\n"
        "      if(w.__delfinKetcherPass){ w.__delfinKetcherPass=false; return; }\n"
        "      ev.preventDefault(); ev.stopPropagation();\n"
        "      hand('open-list','');\n"
        "      draw();\n"
        "    },true);\n"
        "  }\n"
        "  wire();\n"
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

    And it fits the drawing to the frame afterwards.  ``setMolecule`` does
    that itself, but against the size the frame has at the moment it runs --
    which, for a drawing opened into a tab that is not on screen yet, or into
    an editor that is still folded away, is not the size it will be looked at.
    That is how a structure ends up at 10% zoom somewhere off the side.  Fitted
    again a moment later, when the frame has the size it is going to keep.
    """
    return (
        "(function(){\n"
        "  var tries=0;\n"
        "  function settled(api, go){\n"
        "    /* What is on the canvas that nobody has kept.  The clean mark is\n"
        "       set when a drawing is opened and when one is saved, so a\n"
        "       difference from it is work that replacing this would lose.\n"
        "       Unknown means unknown, not dirty: an editor nobody has opened\n"
        "       anything into has nothing to warn about. */\n"
        "    var w=api.__delfinWindow||window;\n"
        "    var clean=w.__delfinKetcherClean;\n"
        "    if(clean===undefined||clean===null){ go(true); return; }\n"
        "    try{\n"
        "      Promise.resolve(api.getKet()).then(function(now){\n"
        "        if(now===clean){ go(true); return; }\n"
        "        go(w.confirm('This drawing has changes that were never '\n"
        "          + 'saved. Open the other one and lose them?'));\n"
        "      }, function(){ go(true); });\n"
        "    }catch(e){ go(true); }\n"
        "  }\n"
        "  function fit(api){\n"
        "    try{\n"
        "      var ed=api.editor;\n"
        "      ed.zoomAccordingContent(ed.struct());\n"
        "      ed.centerStruct();\n"
        "    }catch(e){}\n"
        "  }\n"
        "  function put(){\n"
        "    var host=document.querySelector(" + json.dumps(host_selector) + ");\n"
        "    var frame=host&&host.querySelector('iframe');\n"
        "    var api=null;\n"
        "    try{ api=frame&&frame.contentWindow&&frame.contentWindow.ketcher; }\n"
        "    catch(e){ api=null; }\n"
        "    if(!api){ if(++tries<60) window.setTimeout(put,200); return; }\n"
        "    api.__delfinWindow=frame.contentWindow;\n"
        "    settled(api, function(go){\n"
        "      if(!go) return;\n"
        "      try{\n"
        "        Promise.resolve(api.setMolecule(" + json.dumps(str(text or '')) +
        ")).then(function(){\n"
        "          fit(api);\n"
        "          /* Again once the frame is the size it will be read at. */\n"
        "          window.setTimeout(function(){ fit(api); }, 400);\n"
        "          window.setTimeout(function(){ fit(api); }, 1200);\n"
        "          /* What was just opened is the clean state from here on. */\n"
        "          try{ Promise.resolve(api.getKet()).then(function(k){\n"
        "            api.__delfinWindow.__delfinKetcherClean=k; }); }catch(e){}\n"
        "        });\n"
        "      }catch(e){}\n"
        "    });\n"
        "  }\n"
        "  put();\n"
        "})();"
    )
