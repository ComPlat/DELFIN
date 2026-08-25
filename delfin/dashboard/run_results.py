"""What a calculation folder says about itself, read without running anything.

The reaction graph starts calculations it will not see finish.  A DFT
optimisation of a forty-atom complex runs for two days; the dashboard is
restarted twice in that time, the person goes home, and what has to survive all
of it is a folder on disk and a line in a document pointing at it.  So the
question this module answers is asked of the folder and of nothing else: not of
a scheduler, which forgets; not of a process, which is gone; not of a handle in
memory, which never lasted the night.

Four answers, and the useful part is that there are four rather than two.

``nothing``   Nothing has been written yet.  The job may be queued, it may
              never have been submitted, and from here those look the same --
              which is why this is not called "failed".
``running``   An output exists and has not finished.  A number read out of it
              now is a number from the middle of an optimisation.
``done``      It finished, and what it found can be read.
``failed``    It stopped and said why.  This is a *result*: it is a fact about
              this structure at this level, and a graph that quietly showed a
              gap instead would invite the same two days to be spent again next
              week by somebody reading that gap.

Nothing here starts a process.  It opens files, and only files inside the
folder it was given.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

#: What ORCA writes when it is content, and when it is not.  Both are looked
#: for in the tail: an output is megabytes and the verdict is at the end of it,
#: so reading the whole file to learn one word is a cost paid on every poll of
#: every pending calculation.
_FINISHED = '****ORCA TERMINATED NORMALLY****'
_FINISHED_LOOSE = 'ORCA TERMINATED NORMALLY'
_BROKE = (
    'error termination',
    'aborting the run',
    'ORCA finished by error termination',
)

#: How much of the end of an output is enough to find the verdict in.  ORCA's
#: closing block is a few hundred bytes; a crash message and the line before it
#: are less.  Sixty-four kilobytes is generous by two orders of magnitude and
#: still nothing against a file that is often a hundred megabytes.
_TAIL_BYTES = 64 * 1024

#: How much of the beginning to keep when quoting a failure back to the user.
_WHY_CHARS = 300


def _tail(path: Path, size: int = _TAIL_BYTES) -> str:
    """The last *size* bytes of a file, as text, or '' when unreadable."""
    try:
        with open(path, 'rb') as handle:
            try:
                handle.seek(max(0, path.stat().st_size - size))
            except OSError:
                pass
            return handle.read().decode('utf-8', errors='replace')
    except OSError:
        return ''


def outputs_in(folder: Any) -> List[Path]:
    """Every ORCA output in *folder*, largest last.

    Largest rather than newest.  A job that restarts leaves a short output
    beside a long one and the long one is the run that got somewhere; a
    timestamp only says which was touched last, which after a copy or a
    filesystem with coarse times is not a fact about the chemistry at all.
    """
    root = Path(folder)
    try:
        found = [p for p in root.glob('*.out') if p.is_file()]
    except OSError:
        return []
    return sorted(found, key=lambda p: (p.stat().st_size, p.name))


def what_a_run_says(folder: Any) -> Dict[str, Any]:
    """Read the folder and say where the calculation stands.

    Returns ``{'state', 'why', 'output', 'energy', 'free_energy', 'zpe',
    'imaginary', 'frequency', 'xyz'}``.  The numbers are ``None`` unless the
    state is ``done``: a free energy taken off an output that is still being
    written is a number from the middle of an optimisation, and it would be
    indistinguishable from a finished one once it was in the document.
    """
    blank: Dict[str, Any] = {
        'state': 'nothing', 'why': '', 'output': None,
        'energy': None, 'free_energy': None, 'zpe': None,
        'imaginary': None, 'frequency': None, 'xyz': None,
    }
    root = Path(folder)
    if not root.is_dir():
        return dict(blank, why='there is no such folder')
    found = outputs_in(root)
    if not found:
        return dict(blank, why='nothing has been written yet')
    target = found[-1]
    tail = _tail(target)
    said = dict(blank, output=target.name)

    lowered = tail.lower()
    if _FINISHED in tail or _FINISHED_LOOSE in tail:
        said['state'] = 'done'
    elif any(one.lower() in lowered for one in _BROKE):
        said['state'] = 'failed'
        said['why'] = _why_it_broke(tail)
        return said
    else:
        said['state'] = 'running'
        said['why'] = 'it has not finished'
        return said

    said.update(_numbers_from(target))
    said['xyz'] = _geometry_beside(target)
    return said


def _why_it_broke(tail: str) -> str:
    """The line ORCA stopped on, quoted rather than summarised.

    Summarised, every failure reads as "the calculation failed", and the user
    has to open the output to learn whether it was an SCF that would converge
    with a different guess or a basis set that does not exist.  Both are two
    days either way; only one of them is worth starting again.
    """
    lines = [one.strip() for one in tail.splitlines() if one.strip()]
    for n, one in enumerate(lines):
        if any(mark.lower() in one.lower() for mark in _BROKE):
            around = lines[max(0, n - 2):n + 1]
            said = ' / '.join(around)
            return said[:_WHY_CHARS]
    return 'it stopped without saying why'


def _numbers_from(output: Path) -> Dict[str, Any]:
    """The energies and the modes, out of the readers DELFIN already has.

    Imported here rather than at the top: ``delfin.api`` pulls in a great deal
    that a dashboard tab has no use for until a calculation has actually
    landed, and a graph is opened far more often than a job finishes.
    """
    out: Dict[str, Any] = {}
    try:
        from delfin import energies as _energies
        out['energy'] = _energies.find_electronic_energy(str(output))
        out['free_energy'] = _energies.find_gibbs_energy(str(output))
        out['zpe'] = _energies.find_ZPE(str(output))
    except Exception:                                # noqa: BLE001
        pass
    try:
        from delfin.api import extract_imaginary_frequencies
        found = extract_imaginary_frequencies(str(output.parent))
        if not getattr(found, 'error', None):
            out['imaginary'] = int(getattr(found, 'n_imag', 0) or 0)
            worst = getattr(found, 'most_negative', None)
            out['frequency'] = float(worst) if worst is not None else None
    except Exception:                                # noqa: BLE001
        pass
    return out


def _geometry_beside(output: Path) -> Optional[str]:
    """The geometry this output produced, or None if it produced none.

    ``<stem>.xyz`` beside ``<stem>.out``, which is where ORCA writes the
    optimised structure and what every other part of DELFIN treats as the
    geometry handle.  Nothing else in the folder is considered, and that is
    deliberate: the graph writes ``from_graph.xyz`` there when it sets the job
    up, so a rule like "the newest .xyz" would hand back the structure that was
    sent in and call it the result -- a record saying an optimisation changed
    nothing, on every job that failed to write one.
    """
    candidate = output.with_suffix('.xyz')
    try:
        if not candidate.is_file():
            return None
        text = candidate.read_text(encoding='utf-8', errors='replace')
    except OSError:
        return None
    return text if text.strip() else None


__all__ = ['what_a_run_says', 'outputs_in']
