"""Rebuild the portable Windows ZIP and checksums (developer-side, stdlib only).

Input: the files in NAMES, in this folder (an LF `git archive` copy).
Output: SHA256SUMS.txt here and ../delfin-windows-launcher.zip.

Optional Authenticode signing of the PowerShell scripts, before the
checksums are computed (a signature changes the file). Unset, the build is
unsigned and works the same. Never put a key or password in the repository.

  Linux (osslsigncode >= 2.6, which signs PowerShell scripts):
    DELFIN_SIGN_PFX            path to the PKCS#12 code-signing certificate
    DELFIN_SIGN_PASSWORD_FILE  file holding its password (optional)
  Windows (certificate in the user's store, Set-AuthenticodeSignature):
    DELFIN_SIGN_THUMBPRINT     thumbprint under Cert:\\CurrentUser\\My
  Both:
    DELFIN_SIGN_TIMESTAMP_URL  RFC 3161 timestamp server (required when
                               signing: a timestamped signature stays valid
                               after the certificate expires)

.cmd files cannot carry an Authenticode signature (Windows has no signature
format for batch files); they are two short plain-text wrappers.
"""
from pathlib import Path
import hashlib
import os
import shutil
import subprocess
import sys
import time
import zipfile

NAMES = ['Install.cmd', 'Uninstall.cmd', 'DELFIN.ps1', 'Install.ps1', 'Uninstall.ps1', 'Test-Launcher.ps1',
         'remote_launcher.py', 'remote_bootstrap.sh', 'README.md', 'DELFIN_logo.png', 'DELFIN.ico']
SIGNED = [name for name in NAMES if name.endswith('.ps1')]

#: Fixed text: the file list and settings reach it through the environment,
#: never through the command line.
_WINDOWS_SIGN = (
    "$ErrorActionPreference = 'Stop'; "
    "$certificate = Get-Item -LiteralPath ('Cert:\\CurrentUser\\My\\' + $env:DELFIN_SIGN_THUMBPRINT); "
    "foreach ($file in $env:DELFIN_SIGN_FILES.Split('|')) { "
    "$result = Set-AuthenticodeSignature -LiteralPath $file -Certificate $certificate "
    "-TimestampServer $env:DELFIN_SIGN_TIMESTAMP_URL -HashAlgorithm SHA256; "
    "if ($result.Status -ne 'Valid') { throw ($file + ': ' + $result.StatusMessage) } }"
)


def sign(folder):
    """Sign SIGNED in place when a certificate is configured; return the mode."""
    pfx = os.environ.get('DELFIN_SIGN_PFX')
    thumbprint = os.environ.get('DELFIN_SIGN_THUMBPRINT')
    if not pfx and not thumbprint:
        return 'unsigned'
    timestamp = os.environ.get('DELFIN_SIGN_TIMESTAMP_URL')
    if not timestamp:
        raise SystemExit('DELFIN_SIGN_TIMESTAMP_URL is required for signing.')
    files = [folder / name for name in SIGNED]
    if pfx:
        tool = shutil.which('osslsigncode')
        if not tool:
            raise SystemExit('osslsigncode (>= 2.6) is required for DELFIN_SIGN_PFX.')
        for path in files:
            signed = path.with_name(path.name + '.signed')
            command = [tool, 'sign', '-pkcs12', pfx, '-h', 'sha256', '-n', 'DELFIN Launcher',
                       '-ts', timestamp, '-in', str(path), '-out', str(signed)]
            if os.environ.get('DELFIN_SIGN_PASSWORD_FILE'):
                command[4:4] = ['-readpass', os.environ['DELFIN_SIGN_PASSWORD_FILE']]
            subprocess.run(command, check=True)
            signed.replace(path)
        return 'signed (osslsigncode)'
    if sys.platform != 'win32':
        raise SystemExit('DELFIN_SIGN_THUMBPRINT signs on Windows only; use DELFIN_SIGN_PFX here.')
    powershell = Path(os.environ['SystemRoot']) / 'System32/WindowsPowerShell/v1.0/powershell.exe'
    env = dict(os.environ, DELFIN_SIGN_FILES='|'.join(str(path) for path in files))
    subprocess.run([str(powershell), '-NoProfile', '-NonInteractive', '-Command', _WINDOWS_SIGN],
                   env=env, check=True)
    return 'signed (Set-AuthenticodeSignature)'


def main():
    folder = Path(__file__).resolve().parent
    mode = sign(folder)
    checksum = folder / 'SHA256SUMS.txt'
    checksum.write_text(''.join(
        f'{hashlib.sha256((folder / name).read_bytes()).hexdigest()}  {name}\n'
        for name in NAMES), encoding='ascii')
    destination = folder.parent / 'delfin-windows-launcher.zip'
    # Ordinary entries, as Windows' own "Send to > Compressed folder" writes
    # them: a folder entry, DEFLATE, MS-DOS attributes (archive bit), and each
    # file's own modification time (from `git archive`: the commit's time,
    # written in UTC so every builder gets the same bytes).
    with zipfile.ZipFile(destination, 'w', zipfile.ZIP_DEFLATED) as archive:
        newest = max((folder / name).stat().st_mtime for name in NAMES)
        entry = zipfile.ZipInfo('windows-launcher/', date_time=time.gmtime(max(newest, 315532800))[:6])
        entry.create_system = 0
        entry.external_attr = 0x10
        archive.writestr(entry, b'')
        for name in NAMES + [checksum.name]:
            path = folder / name
            stamp = time.gmtime(max(path.stat().st_mtime, 315532800))[:6]
            item = zipfile.ZipInfo('windows-launcher/' + name, date_time=stamp)
            item.compress_type = zipfile.ZIP_DEFLATED
            item.create_system = 0
            item.external_attr = 0x20
            archive.writestr(item, path.read_bytes())
    with zipfile.ZipFile(destination) as archive:
        if archive.testzip() is not None:
            raise RuntimeError('Corrupt ZIP')
    print(f'{destination.name}: {len(NAMES) + 1} files, {mode}; CRC checked')


if __name__ == '__main__':
    main()
