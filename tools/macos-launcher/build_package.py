"""Rebuild the macOS launcher ZIP, its checksums and the one-line install (stdlib only).

    python tools/macos-launcher/build_package.py      (from an LF checkout)

Output:
  ../delfin-macos-launcher.zip   one folder ``macos-launcher/`` inside, as
                                 Finder's Archive Utility extracts it:
        Install.command, Uninstall.command, README.md, START-HERE.txt,
        SHA256SUMS.txt, and DELFIN.app -- a finished app bundle, never
        assembled on the Mac, so it can be signed and notarized as it is
        (sign_and_notarize.sh).
  SHA256SUMS.txt here            the same list, for review.
  README.md of the repository    the one-line install between its
                                 ``macos-install-line`` markers: the release
                                 URL of this version and the ZIP's SHA-256.
                                 It cannot live inside the ZIP (it names the
                                 ZIP's own digest).

Taken from elsewhere in the repository at build time, never kept as copies:
``remote_bootstrap.sh`` and ``remote_launcher.py`` (the Windows launcher's
server half, byte for byte, so both apps share the server's
``~/.cache/delfin-launcher/<version>/``), ``fetch_micromamba.sh`` (DELFIN's
pinned micromamba fetch), and the version from pyproject.toml.
"""
from __future__ import annotations

import hashlib
import re
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DESTINATION = ROOT / 'tools' / 'delfin-macos-launcher.zip'
README = ROOT / 'README.md'
REPOSITORY = 'ComPlat/DELFIN'
RESOURCES = 'DELFIN.app/Contents/Resources/'

#: Package path -> source (a path, or None for what the build writes).
SOURCES = {
    'Install.command': HERE / 'Install.command',
    'Uninstall.command': HERE / 'Uninstall.command',
    'README.md': HERE / 'README.md',
    'START-HERE.txt': HERE / 'START-HERE.txt',
    'DELFIN.app/Contents/Info.plist': HERE / 'app' / 'Info.plist',
    'DELFIN.app/Contents/MacOS/DELFIN': HERE / 'app' / 'DELFIN',
    RESOURCES + 'DELFIN.sh': HERE / 'DELFIN.sh',
    RESOURCES + 'setup_local.sh': HERE / 'setup_local.sh',
    RESOURCES + 'fetch_micromamba.sh': ROOT / 'delfin' / 'installers' / 'fetch_micromamba.sh',
    RESOURCES + 'remote_bootstrap.sh': ROOT / 'tools' / 'windows-launcher' / 'remote_bootstrap.sh',
    RESOURCES + 'remote_launcher.py': ROOT / 'tools' / 'windows-launcher' / 'remote_launcher.py',
    RESOURCES + 'DELFIN.icns': HERE / 'app' / 'DELFIN.icns',
    RESOURCES + 'VERSION': None,
}
NAMES = list(SOURCES)
EXECUTABLE = {'Install.command', 'Uninstall.command', 'DELFIN.app/Contents/MacOS/DELFIN',
              RESOURCES + 'DELFIN.sh', RESOURCES + 'setup_local.sh',
              RESOURCES + 'fetch_micromamba.sh', RESOURCES + 'remote_bootstrap.sh'}
BINARY = {RESOURCES + 'DELFIN.icns'}
#: Fixed entry time: the same sources give the same ZIP on every machine.
STAMP = (2026, 1, 1, 0, 0, 0)
LINE_START = '<!-- macos-install-line: written by tools/macos-launcher/build_package.py -->'
LINE_END = '<!-- /macos-install-line -->'


def version() -> str:
    text = (ROOT / 'pyproject.toml').read_text(encoding='utf-8')
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if not match:
        raise SystemExit('No version in pyproject.toml')
    return match.group(1)


def contents() -> dict[str, bytes]:
    """Every file of the package, by its path in the package."""
    files = {}
    for name, source in SOURCES.items():
        files[name] = source.read_bytes() if source else (version() + '\n').encode('ascii')
    files['DELFIN.app/Contents/Info.plist'] = files['DELFIN.app/Contents/Info.plist'].replace(
        b'@VERSION@', version().encode('ascii'))
    for name, data in files.items():
        # A checkout with CRLF line ends (Windows, core.autocrlf) gives
        # scripts bash on the Mac cannot run: build from an LF copy.
        if name not in BINARY and b'\r' in data:
            raise SystemExit(f'{name} has CR line ends; build from an LF checkout (git archive).')
    return files


def checksums(files: dict[str, bytes]) -> str:
    return ''.join(f'{hashlib.sha256(files[name]).hexdigest()}  {name}\n' for name in NAMES)


def folders() -> list[str]:
    seen = []
    for name in NAMES:
        parts = name.split('/')[:-1]
        for depth in range(1, len(parts) + 1):
            folder = '/'.join(parts[:depth]) + '/'
            if folder not in seen:
                seen.append(folder)
    return seen


def install_line(digest: str) -> str:
    """The one-line install: a fixed release asset, checked against its
    SHA-256 before anything of it runs; never piped into a shell."""
    url = f'https://github.com/{REPOSITORY}/releases/download/v{version()}/delfin-macos-launcher.zip'
    return ('d=$(mktemp -d) && curl -fL --proto =https -o "$d/delfin.zip" ' + url
            + f' && echo "{digest}  $d/delfin.zip" | shasum -a 256 -c'
            + ' && ditto -x -k "$d/delfin.zip" "$d" && bash "$d/macos-launcher/Install.command"')


def readme_block(digest: str) -> str:
    return (f'{LINE_START}\n```bash\n{install_line(digest)}\n```\n{LINE_END}')


def main() -> None:
    files = contents()
    sums = checksums(files)
    (HERE / 'SHA256SUMS.txt').write_bytes(sums.encode('ascii'))
    with zipfile.ZipFile(DESTINATION, 'w', zipfile.ZIP_DEFLATED) as archive:
        for folder in ['macos-launcher/'] + ['macos-launcher/' + f for f in folders()]:
            entry = zipfile.ZipInfo(folder, date_time=STAMP)
            entry.create_system = 3
            entry.external_attr = (0o40755 << 16) | 0x10
            archive.writestr(entry, b'')
        for name in NAMES + ['SHA256SUMS.txt']:
            data = sums.encode('ascii') if name == 'SHA256SUMS.txt' else files[name]
            item = zipfile.ZipInfo('macos-launcher/' + name, date_time=STAMP)
            item.compress_type = zipfile.ZIP_DEFLATED
            # Unix modes, which Archive Utility, ditto and unzip apply.
            item.create_system = 3
            item.external_attr = (0o100755 if name in EXECUTABLE else 0o100644) << 16
            archive.writestr(item, data)
    with zipfile.ZipFile(DESTINATION) as archive:
        if archive.testzip() is not None:
            raise RuntimeError('Corrupt ZIP')
    digest = hashlib.sha256(DESTINATION.read_bytes()).hexdigest()
    text = README.read_text(encoding='utf-8')
    start, end = text.index(LINE_START), text.index(LINE_END) + len(LINE_END)
    README.write_text(text[:start] + readme_block(digest) + text[end:], encoding='utf-8', newline='')
    print(f'{DESTINATION.name}: {len(NAMES) + 1} files, sha256 {digest}; CRC checked; README line updated')


if __name__ == '__main__':
    main()
