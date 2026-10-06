"""Rebuild the portable Windows ZIP and checksums (developer-side, stdlib only)."""
from pathlib import Path
import hashlib
import zipfile


def main():
    folder = Path(__file__).resolve().parent
    names = ['Install.cmd', 'Uninstall.cmd', 'DELFIN.ps1', 'Install.ps1', 'Uninstall.ps1', 'Test-Launcher.ps1', 'remote_launcher.py', 'remote_bootstrap.sh',
             'README.md', 'DELFIN_logo.png', 'DELFIN.ico']
    checksum = folder / 'SHA256SUMS.txt'
    checksum.write_text(''.join(
        f'{hashlib.sha256((folder / name).read_bytes()).hexdigest()}  {name}\n'
        for name in names), encoding='ascii')
    destination = folder.parent / 'delfin-windows-launcher.zip'
    with zipfile.ZipFile(destination, 'w', zipfile.ZIP_DEFLATED) as archive:
        for name in names + [checksum.name]:
            item = zipfile.ZipInfo('windows-launcher/' + name, date_time=(1980, 1, 1, 0, 0, 0))
            item.compress_type = zipfile.ZIP_DEFLATED
            item.external_attr = 0o100644 << 16
            archive.writestr(item, (folder / name).read_bytes())
    with zipfile.ZipFile(destination) as archive:
        if archive.testzip() is not None:
            raise RuntimeError('Corrupt ZIP')
    print(f'{destination.name}: {len(names) + 1} files; CRC checked')


if __name__ == '__main__':
    main()
