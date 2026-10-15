#!/usr/bin/env python3
"""Check a built wheel or sdist for what it must and must not contain.

Input: one or more paths to ``*.whl`` and ``*.tar.gz`` files, as written by
``python -m build``. Output: one line per finding on stdout and exit status
1 when there is any, 0 when every archive passes, 2 on a usage error.

What it refuses:

- any member under an installer's output directory
  (``delfin/<x>_tools/{bin,share,downloads,build,third_party,logs,...}``):
  ``install_qm_tools.sh`` and its siblings download third-party programs
  there, and a release built on a machine where they ran used to ship them
  (the 1.1.12 wheel carried xtb4stda and stda_v1.6.1, LGPL-3.0, without
  their licence text or source);
- any member below ``delfin/`` whose file name is one of those programs or
  their parameter files (xtb4stda, stda, xtb, crest, ...), wherever it sits;
- native code of any kind: an ELF, Mach-O or PE header, or a
  ``.so``/``.dylib``/``.dll``/``.exe``/``.pyd`` member. The wheel is tagged
  ``py3-none-any``; native code in it is wrong whatever its licence.

What it requires:

- the licence files in :data:`REQUIRED_LICENCE_FILES`, in the wheel's
  ``*.dist-info/licenses/`` and in the sdist's root directory;
- ``License-Expression: LGPL-3.0-or-later`` and one ``License-File`` line per
  required file in the wheel's METADATA and the sdist's PKG-INFO (PEP 639).

Rejected alternatives: ``unzip -l | grep qm_tools/bin`` (sees one path and
no binary placed anywhere else); ``check-wheel-contents`` (a further
dependency, and it knows nothing about this project's installer
directories or its licence files).
"""

from __future__ import annotations

import argparse
import re
import sys
import tarfile
import zipfile
from dataclasses import dataclass
from email.parser import HeaderParser
from pathlib import Path, PurePosixPath

LICENSE_EXPRESSION = "LGPL-3.0-or-later"

#: Files that carry the project's licence and the notices of what it bundles.
REQUIRED_LICENCE_FILES: tuple[str, ...] = ("LICENSE", "COPYING", "COPYING.LESSER")

#: Directories the installers write into (see delfin/qm_tools/install_qm_tools.sh
#: and the csp/mlp/analysis/ai installers).
_INSTALLER_OUTPUT = re.compile(
    r"(^|/)delfin/(qm|csp|mlp|analysis|ai)_tools/"
    r"(bin|share|downloads|build|third_party|logs|\.mamba_env|\.build)(/|$)"
)

#: File names of third-party programs and data the installers fetch.
_THIRD_PARTY_NAMES = frozenset({
    "xtb4stda", "stda", "stda_v1.6.1", "std2", "xtb", "xtb-gxtb", "crest",
    "dftb+", "mopac", "orca", "gnrs",
    ".xtb4stdarc", ".param_stda1.xtb", ".param_stda2.xtb",
})

_NATIVE_SUFFIXES = (".so", ".dylib", ".dll", ".exe", ".pyd")

_MACHO_MAGIC = (
    b"\xfe\xed\xfa\xce", b"\xfe\xed\xfa\xcf",
    b"\xce\xfa\xed\xfe", b"\xcf\xfa\xed\xfe",
    b"\xca\xfe\xba\xbe",
)


@dataclass(frozen=True)
class Member:
    name: str      # path inside the archive, POSIX separators
    head: bytes    # first bytes of a regular file, b"" for anything else


def native_kind(head: bytes) -> str:
    """The executable format *head* starts with, or ``""``."""
    if head.startswith(b"\x7fELF"):
        return "ELF"
    if head[:4] in _MACHO_MAGIC:
        return "Mach-O"
    if head.startswith(b"MZ") and len(head) >= 64:
        offset = int.from_bytes(head[60:64], "little")
        if head[offset:offset + 4] == b"PE\x00\x00":
            return "PE"
    return ""


def content_findings(members: list[Member]) -> list[str]:
    """Members that must not be in a DELFIN distribution."""
    out = []
    for m in members:
        path = PurePosixPath(m.name)
        if _INSTALLER_OUTPUT.search(m.name):
            out.append(f"{m.name}: inside an installer's output directory")
            continue
        if path.parts[:1] == ("delfin",) and path.name in _THIRD_PARTY_NAMES:
            out.append(f"{m.name}: named like a third-party program the "
                       "installers download")
            continue
        kind = native_kind(m.head)
        if kind:
            out.append(f"{m.name}: native executable code ({kind})")
        elif path.suffix.lower() in _NATIVE_SUFFIXES:
            out.append(f"{m.name}: native library suffix {path.suffix}")
    return out


def metadata_findings(text: str, where: str) -> list[str]:
    """PEP 639 fields of a METADATA or PKG-INFO document."""
    headers = HeaderParser().parsestr(text)
    out = []
    expression = (headers.get("License-Expression") or "").strip()
    if expression != LICENSE_EXPRESSION:
        out.append(f"{where}: License-Expression is {expression!r}, "
                   f"expected {LICENSE_EXPRESSION!r}")
    listed = {v.strip() for v in headers.get_all("License-File") or []}
    for name in REQUIRED_LICENCE_FILES:
        if name not in listed:
            out.append(f"{where}: no License-File entry for {name}")
    return out


def _wheel(path: Path) -> tuple[list[Member], list[str]]:
    members, findings = [], []
    with zipfile.ZipFile(path) as zf:
        infos = zf.infolist()
        for info in infos:
            head = b""
            if not info.is_dir():
                with zf.open(info) as fh:
                    head = fh.read(4096)
            members.append(Member(info.filename, head))
        dist_info = sorted({n.split("/", 1)[0] for n in zf.namelist()
                            if n.split("/", 1)[0].endswith(".dist-info")})
        if len(dist_info) != 1:
            return members, [f"{path.name}: expected one .dist-info, "
                             f"found {dist_info}"]
        names = set(zf.namelist())
        for lic in REQUIRED_LICENCE_FILES:
            if f"{dist_info[0]}/licenses/{lic}" not in names:
                findings.append(f"{path.name}: {dist_info[0]}/licenses/{lic} "
                                "is missing")
        meta = f"{dist_info[0]}/METADATA"
        if meta in names:
            findings += metadata_findings(
                zf.read(meta).decode("utf-8"), f"{path.name}:{meta}")
        else:
            findings.append(f"{path.name}: {meta} is missing")
    return members, findings


def _sdist(path: Path) -> tuple[list[Member], list[str]]:
    members, findings = [], []
    with tarfile.open(path, "r:*") as tf:
        infos = tf.getmembers()
        for info in infos:
            head = b""
            if info.isfile():
                fh = tf.extractfile(info)
                head = fh.read(4096) if fh else b""
            members.append(Member(info.name, head))
        roots = sorted({i.name.split("/", 1)[0] for i in infos})
        if len(roots) != 1:
            return members, [f"{path.name}: expected one root directory, "
                             f"found {roots}"]
        names = {i.name for i in infos}
        for lic in REQUIRED_LICENCE_FILES:
            if f"{roots[0]}/{lic}" not in names:
                findings.append(f"{path.name}: {roots[0]}/{lic} is missing")
        pkg_info = f"{roots[0]}/PKG-INFO"
        if pkg_info in names:
            fh = tf.extractfile(pkg_info)
            findings += metadata_findings(
                fh.read().decode("utf-8") if fh else "", f"{path.name}:PKG-INFO")
        else:
            findings.append(f"{path.name}: PKG-INFO is missing")
    # The sdist's own root prefix is not part of the package path.
    stripped = [Member(m.name.split("/", 1)[1] if "/" in m.name else "", m.head)
                for m in members]
    return stripped, findings


def check(path: Path) -> list[str]:
    """Every finding for one archive; an empty list means it passes."""
    if path.name.endswith(".whl"):
        members, findings = _wheel(path)
    elif path.name.endswith((".tar.gz", ".tgz")):
        members, findings = _sdist(path)
    else:
        return [f"{path}: neither a wheel nor an sdist"]
    return findings + [f"{path.name}: {f}" for f in content_findings(members)]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("dists", nargs="+", type=Path)
    args = parser.parse_args(argv)
    total = 0
    for dist in args.dists:
        if not dist.is_file():
            print(f"{dist}: not a file", file=sys.stderr)
            return 2
        findings = check(dist)
        total += len(findings)
        for line in findings:
            print(line)
        if not findings:
            print(f"{dist.name}: ok")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
