"""The release distributions carry no third-party program binaries.

Input: the paths of the built distributions (``dist/*``): one sdist
(``.tar.gz``) and one wheel (``.whl``), as ``python -m build`` writes them.

Output: exit 0 and one line per distribution with its member count, or
exit 1 and one line per problem:
- not exactly one sdist and one wheel;
- a member under ``qm_tools/bin/`` or ``qm_tools/share/xtb4stda/`` (the
  xtb4stda/stda binaries and parameters that delfin_complat 1.1.12 shipped
  without their LGPL notice; compliance audit item B2);
- a member that is an executable image whatever its name: ELF
  (``7f 45 4c 46``), Mach-O (thin or fat, either byte order), or a PE file
  (``MZ`` with ``PE\\0\\0`` at the offset the DOS header names).

Semantics: the check reads bytes, not file names, so a renamed binary is
still found; DELFIN is a ``py3-none-any`` wheel, and any native image in it
is someone else's program redistributed. The installers download those at
run time, with checksums (delfin/qm_tools/install_qm_tools.sh).

Rejected: ``unzip -l dist/*.whl | grep -q qm_tools/bin && exit 1`` (the
audit's sketch): it reads names only, the wheel only, and a listing that
fails (no wheel, a corrupt one) makes grep find nothing, which reads as a
pass.
"""

from __future__ import annotations

import sys
import tarfile
import zipfile
from pathlib import Path

_FORBIDDEN_DIRS = ("qm_tools/bin/", "qm_tools/share/xtb4stda/")
_MACHO = {bytes.fromhex(h) for h in ("feedface", "feedfacf", "cefaedfe", "cffaedfe", "cafebabe")}


def native_image(head: bytes) -> str | None:
    """The kind of executable image ``head`` (the first 4096 bytes) starts, or None."""
    if head[:4] == b"\x7fELF":
        return "ELF"
    if head[:4] in _MACHO:
        return "Mach-O"
    if head[:2] == b"MZ" and len(head) >= 0x40:
        offset = int.from_bytes(head[0x3C:0x40], "little")
        if head[offset:offset + 4] == b"PE\0\0":
            return "PE"
    return None


def _members(path: Path):
    """(name, first 4096 bytes) for every regular file in a distribution."""
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            for info in archive.infolist():
                if not info.is_dir():
                    with archive.open(info) as handle:
                        yield info.filename, handle.read(4096)
    else:
        with tarfile.open(path, "r:gz") as archive:
            for info in archive.getmembers():
                if info.isfile():
                    handle = archive.extractfile(info)
                    yield info.name, handle.read(4096) if handle else b""


def check(paths: list[Path]) -> tuple[list[str], list[str]]:
    problems, summary = [], []
    wheels = [p for p in paths if p.name.endswith(".whl")]
    sdists = [p for p in paths if p.name.endswith(".tar.gz")]
    if len(wheels) != 1 or len(sdists) != 1:
        problems.append(f"expected one wheel and one sdist, got {[p.name for p in paths]}")
    for path in wheels + sdists:
        count = 0
        for name, head in _members(path):
            count += 1
            if any(d in name for d in _FORBIDDEN_DIRS):
                problems.append(f"{path.name}: {name} is a third-party tool file")
            kind = native_image(head)
            if kind:
                problems.append(f"{path.name}: {name} is an executable image ({kind})")
        if count == 0:
            problems.append(f"{path.name}: no files")
        summary.append(f"{path.name}: {count} files, no program binaries"
                       if not any(p.startswith(path.name) for p in problems) else
                       f"{path.name}: {count} files")
    return problems, summary


def main(argv: list[str] | None = None) -> int:
    paths = [Path(a) for a in (sys.argv[1:] if argv is None else argv)]
    problems, summary = check(paths)
    for line in summary:
        print(line)
    for line in problems:
        print(f"ERROR: {line}", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
