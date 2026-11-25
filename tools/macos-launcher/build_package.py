"""Rebuild the macOS launcher ZIP and its checksums (developer-side, stdlib only).

    python tools/macos-launcher/build_package.py

Output: SHA256SUMS.txt here and ../delfin-macos-launcher.zip, with one folder
``macos-launcher/`` inside, as Finder's Archive Utility extracts it.

The package is this folder's own files plus files that live elsewhere in the
repository and are taken from there at build time, never kept as copies here:

* ``remote_bootstrap.sh`` and ``remote_launcher.py``: the server half of the
  Windows launcher. The Mac uses it byte for byte, so both apps share the
  server's ``~/.cache/delfin-launcher/<version>/`` copy.
* ``fetch_micromamba.sh``: DELFIN's pinned, digest-checked micromamba fetch.
* ``DELFIN_logo.png``: the logo (the installer makes the app icon from it).
* ``VERSION``: the DELFIN release that ``setup_local.sh`` installs.

Unsigned and not notarized: Install.command explains Gatekeeper and offers
to remove the download mark (START-HERE.txt has the manual way).
"""
from __future__ import annotations

import hashlib
import re
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

OWN = ['DELFIN.sh', 'Install.command', 'Uninstall.command', 'setup_local.sh',
       'README.md', 'START-HERE.txt']
SHARED = {
    'remote_bootstrap.sh': ROOT / 'tools' / 'windows-launcher' / 'remote_bootstrap.sh',
    'remote_launcher.py': ROOT / 'tools' / 'windows-launcher' / 'remote_launcher.py',
    'fetch_micromamba.sh': ROOT / 'delfin' / 'installers' / 'fetch_micromamba.sh',
    'DELFIN_logo.png': ROOT / 'tools' / 'windows-launcher' / 'DELFIN_logo.png',
}
EXECUTABLE = {'DELFIN.sh', 'Install.command', 'Uninstall.command', 'setup_local.sh',
              'fetch_micromamba.sh', 'remote_bootstrap.sh'}
#: The order of the checksum file (DELFIN.sh's PACKAGE_FILES lists the same set).
NAMES = ['DELFIN.sh', 'Install.command', 'Uninstall.command', 'setup_local.sh',
         'fetch_micromamba.sh', 'remote_bootstrap.sh', 'remote_launcher.py',
         'README.md', 'START-HERE.txt', 'DELFIN_logo.png', 'VERSION']
#: Fixed entry time: the same sources give the same ZIP on every machine.
STAMP = (2026, 1, 1, 0, 0, 0)
DESTINATION = ROOT / 'tools' / 'delfin-macos-launcher.zip'


def version() -> str:
    text = (ROOT / 'pyproject.toml').read_text(encoding='utf-8')
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if not match:
        raise SystemExit('No version in pyproject.toml')
    return match.group(1)


def contents() -> dict[str, bytes]:
    """Every file of the package, by its name in the package."""
    files = {name: (HERE / name).read_bytes() for name in OWN}
    files.update({name: path.read_bytes() for name, path in SHARED.items()})
    files['VERSION'] = (version() + '\n').encode('ascii')
    assert sorted(files) == sorted(NAMES), sorted(files)
    for name, data in files.items():
        # A checkout with CRLF line ends (Windows, core.autocrlf) gives
        # scripts that bash on the Mac cannot run: build from an LF copy.
        if not name.endswith('.png') and b'\r' in data:
            raise SystemExit(f'{name} has CR line ends; build from an LF checkout (git archive).')
    return files


def checksums(files: dict[str, bytes]) -> str:
    return ''.join(f'{hashlib.sha256(files[name]).hexdigest()}  {name}\n' for name in NAMES)


def main() -> None:
    files = contents()
    sums = checksums(files)
    (HERE / 'SHA256SUMS.txt').write_bytes(sums.encode('ascii'))
    with zipfile.ZipFile(DESTINATION, 'w', zipfile.ZIP_DEFLATED) as archive:
        folder = zipfile.ZipInfo('macos-launcher/', date_time=STAMP)
        folder.create_system = 3
        folder.external_attr = (0o40755 << 16) | 0x10
        archive.writestr(folder, b'')
        for name in NAMES + ['SHA256SUMS.txt']:
            data = sums.encode('ascii') if name == 'SHA256SUMS.txt' else files[name]
            item = zipfile.ZipInfo('macos-launcher/' + name, date_time=STAMP)
            item.compress_type = zipfile.ZIP_DEFLATED
            # Unix modes, which Archive Utility and unzip apply: the
            # .command files and scripts stay executable after extraction.
            item.create_system = 3
            item.external_attr = (0o100755 if name in EXECUTABLE else 0o100644) << 16
            archive.writestr(item, data)
    with zipfile.ZipFile(DESTINATION) as archive:
        if archive.testzip() is not None:
            raise RuntimeError('Corrupt ZIP')
    print(f'{DESTINATION.name}: {len(NAMES) + 1} files; CRC checked')


if __name__ == '__main__':
    main()
