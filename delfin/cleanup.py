from pathlib import Path
from typing import Optional, Sequence

DEFAULT_PATTERNS: Sequence[str] = [
    "*.cpcm",
    "*.cpcm_corr",
    "*.densitiesinfo",
    "*.tmp",
    "*.tmp*",
    "*.bas*",
    "*_D0*",
]

def cleanup(folder: str = ".",
            recursive: bool = False,
            dry_run: bool = False,
            patterns: Optional[Sequence[str]] = None) -> int:
    root = Path(folder).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"{root} ist kein Ordner.")

    pats = list(patterns or DEFAULT_PATTERNS)

    def iter_files():
        for pat in pats:
            yield from (root.rglob(pat) if recursive else root.glob(pat))

    to_delete = sorted({p for p in iter_files() if p.is_file()})
    count = 0
    for f in to_delete:
        if dry_run:
            # don't spend anything
            pass
        else:
            try:
                f.unlink()
                count += 1
            except Exception:
                # Ignore errors or log them if necessary
                pass
    return count

def cleanup_all(folder: str = ".",
                dry_run: bool = False,
                patterns: Optional[Sequence[str]] = None) -> int:    
    return cleanup(folder=folder, recursive=True, dry_run=dry_run, patterns=patterns)
