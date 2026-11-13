"""Whether a bypass shell command would destroy what nothing can restore.

Input: a shell command, the directory it runs in and the session's
permissions. Output: ``refusal(cmd, cwd, perms)`` -- None, or the text that
refuses the command and names the way that keeps the file: move it into the
workspace's trash (``.delfin-trash/<UTC date>/``), which DELFIN empties after
``TRASH_DAYS`` days.

Semantics. Bypass asks nothing, so a command that deletes or overwrites a
file is undone only by what kept a copy. A file counts as restorable when:

* it is empty, or small enough for the shell change journal to hold its
  pre-image (``shell_changes.MAX_PRE_IMAGE_BYTES``; the journal's total
  budget is counted over the command's targets);
* git holds it: tracked, and identical to HEAD;
* its current bytes are the ones this session wrote (the change journal
  says so): destroying it loses this session's own output, nothing that
  was there before;
* it lies in the session's scratch, or it is a symbolic link being removed
  (the link goes, the file it names stays).

Anything else that a command names for deletion or overwrite -- ``rm``,
``unlink``, ``shred``, ``truncate``, ``find ... -delete`` / ``-exec rm``,
``dd of=``, a ``>`` redirection onto an existing file, a ``cp``/``mv``/
``install``/``ln -f`` destination that exists -- refuses the command.
A directory target is walked (bounded); a tree too large to walk counts as
not restorable. Unattended runs keep going: the refusal is the answer, with
the command that does the job without the loss.

Not seen, said here rather than implied: a deletion inside an interpreter
or a script (``python -c``, ``./clean.sh``), a variable or a command
substitution as the target, and an overwrite by a program that opens its
output itself. The overlay workspace is the answer for those (every write
of the command lands in a layer DELFIN applies afterwards); this module is
the deterministic floor until then.

Rejected: asking the person -- bypass is the mode nobody answers, and the
user's rule is that safety must not stop unattended work. Rejected:
hard-linking every file before each command (cheap pre-images for any
size) -- a link raises st_nlink, which DELFIN's own checks and other tools
read as a planted file, and it does not protect an in-place overwrite.
Rejected: a trash outside the workspace -- the cage lets the agent write
only inside its roots, so the move would itself be refused.
"""

from __future__ import annotations

import glob
import os
import shlex
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional

#: The workspace's trash directory, relative to a root.
TRASH = ".delfin-trash"
#: Days a trashed file is kept before DELFIN removes it.
TRASH_DAYS = 7
#: A directory target is walked up to this many entries / seconds; beyond,
#: it counts as not restorable.
WALK_MAX_ENTRIES = 5_000
WALK_MAX_SECONDS = 2.0

_WRAPPERS = {"nice", "nohup", "time", "command", "env", "ionice", "stdbuf",
             "sudo", "doas", "exec"}
_TIMEOUT_LIKE = {"timeout"}
#: Wrapper options that take the next word as their value (nice -n 5,
#: ionice -c 3, stdbuf -o L, timeout -s KILL, env -u NAME, sudo -u user).
_WRAPPER_VALUED = {"-n", "-c", "-p", "-s", "-k", "-o", "-e", "-i", "-u",
                   "--signal", "--kill-after", "--adjustment", "--class",
                   "--classdata", "--unset", "--user"}
_DELETERS = {"rm", "unlink", "shred"}
_COPIERS = {"cp", "mv", "install", "ln"}


def _journal_limits() -> tuple[int, int]:
    from . import shell_changes as _sc
    return _sc.MAX_PRE_IMAGE_BYTES, _sc.MAX_HELD_BYTES


def _strip_prefix(argv: list[str]) -> list[str]:
    """The command proper: leading VAR=value words and wrappers dropped."""
    i = 0
    while i < len(argv):
        word = argv[i]
        if "=" in word and not word.startswith("-") and word.split("=", 1)[0].isidentifier():
            i += 1
            continue
        base = os.path.basename(word)
        if base in _WRAPPERS or base in _TIMEOUT_LIKE:
            i += 1
            while i < len(argv) and argv[i].startswith("-"):
                i += 2 if argv[i] in _WRAPPER_VALUED else 1
            if base in _TIMEOUT_LIKE:
                i += 1                 # the duration
            continue
        break
    return argv[i:]


def _operands(args: list[str], takes_value: Iterable[str] = ()) -> list[str]:
    """Non-option words; after ``--`` every word is an operand."""
    out, i, takes = [], 0, set(takes_value)
    while i < len(args):
        word = args[i]
        if word == "--":
            out.extend(args[i + 1:])
            break
        if word.startswith("-") and word != "-":
            if word in takes:
                i += 1
            i += 1
            continue
        out.append(word)
        i += 1
    return out


def _expand(word: str, cwd: Path) -> list[Path]:
    """A word as the shell would name files: ~ and globs; a word holding a
    variable or a substitution names nothing this reading can know."""
    if "$" in word or "`" in word:
        return []
    word = os.path.expanduser(word)
    path = word if os.path.isabs(word) else str(cwd / word)
    if any(ch in word for ch in "*?["):
        return [Path(p) for p in sorted(glob.glob(path))]
    return [Path(path)]


#: find predicates whose meaning this reading follows exactly: with only
#: these, the files a ``-delete`` removes are the ones whose name matches.
_FIND_PLAIN = {"-type": 1, "-maxdepth": 1, "-mindepth": 1, "-name": 1,
               "-iname": 1, "-delete": 0, "-print": 0, "-exec": None}


def _find_names(preds: list[str]) -> Optional[list[tuple[str, bool]]]:
    """``[(pattern, ignore case)]`` when the predicates are ``-name``/
    ``-iname`` (ANDed with -type f / depth limits) and ``-delete`` or
    ``-exec rm``; None for anything else (then the whole tree counts)."""
    names: list[tuple[str, bool]] = []
    i = 0
    while i < len(preds):
        word = preds[i]
        if word not in _FIND_PLAIN:
            return None
        if word == "-exec":
            end = next((j for j in range(i + 1, len(preds)) if preds[j] in (";", "+")),
                       None)
            if end is None:
                return None
            i = end + 1
            continue
        if word == "-type" and (i + 1 >= len(preds) or preds[i + 1] != "f"):
            return None
        if word in ("-name", "-iname"):
            if i + 1 >= len(preds):
                return None
            names.append((preds[i + 1], word == "-iname"))
        i += 1 + _FIND_PLAIN[word]
    return names or None


class _Matching(type(Path())):
    """A find start path, with the name patterns its -delete obeys."""

    def __new__(cls, path, names):
        self = super().__new__(cls, path)
        self.names = names
        return self

    def __init__(self, path, names):
        try:
            super().__init__(path)
        except TypeError:
            pass
        self.names = names


def _redirect_targets(segment: str) -> list[str]:
    """Words after ``>``, ``>|``, ``&>`` and ``N>`` (not ``>>``, not ``>&N``)."""
    try:
        lex = shlex.shlex(segment, posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        toks = list(lex)
    except ValueError:
        return []
    out = []
    for i, tok in enumerate(toks):
        if tok in (">", ">|", "&>", "1>", "2>") and i + 1 < len(toks):
            nxt = toks[i + 1]
            if nxt.startswith("&") or nxt == "/dev/null":
                continue
            out.append(nxt)
    return out


def destroyed(cmd: str, cwd: Path) -> list[tuple[str, Path, bool]]:
    """``(verb, path, recursive)`` for each existing path the command names
    for deletion or overwrite. Never raises."""
    from .api_client import _split_shell_segments
    found: list[tuple[str, Path, bool]] = []
    try:
        segments = _split_shell_segments(cmd)
    except Exception:
        return found
    for seg in segments:
        for word in _redirect_targets(seg):
            for p in _expand(word, cwd):
                if p.is_file():
                    found.append(("overwrite", p, False))
        try:
            argv = shlex.split(seg, comments=False)
        except ValueError:
            continue
        # The redirections are not operands.
        clean, skip = [], False
        for word in argv:
            if skip:
                skip = False
                continue
            if word in (">", ">>", ">|", "&>", "1>", "2>", "<", "2>>"):
                skip = True
                continue
            if word.startswith((">", "<", "2>", "1>", "&>")):
                continue
            clean.append(word)
        argv = _strip_prefix(clean)
        if not argv:
            continue
        verb, args = os.path.basename(argv[0]), argv[1:]
        if verb in _DELETERS:
            recursive = verb == "rm" and any(
                a in ("-r", "-R", "--recursive") or (a.startswith("-") and not a.startswith("--")
                                                     and ("r" in a or "R" in a))
                for a in args)
            for word in _operands(args):
                for p in _expand(word, cwd):
                    if p.exists() or p.is_symlink():
                        found.append((verb, p, recursive))
        elif verb == "truncate":
            for word in _operands(args, ("-s", "--size", "-r", "--reference")):
                for p in _expand(word, cwd):
                    if p.is_file():
                        found.append(("truncate", p, False))
        elif verb == "dd":
            for a in args:
                if a.startswith("of="):
                    for p in _expand(a[3:], cwd):
                        if p.is_file():
                            found.append(("overwrite", p, False))
        elif verb == "find":
            if "-delete" in args or any(
                    a == "-exec" and i + 1 < len(args)
                    and os.path.basename(args[i + 1]) in _DELETERS
                    for i, a in enumerate(args)):
                starts = []
                for a in args:
                    if a.startswith("-") or a in ("(", "!", ")"):
                        break
                    starts.append(a)
                names = _find_names(args[len(starts):])
                for word in starts or ["."]:
                    for p in _expand(word, cwd):
                        if not p.exists():
                            continue
                        if names is None:
                            found.append(("delete", p, True))
                        else:
                            found.append(("delete", _Matching(p, names), True))
        elif verb in _COPIERS:
            if verb == "ln" and not any(a in ("-f", "--force") or (
                    a.startswith("-") and not a.startswith("--") and "f" in a) for a in args):
                continue
            target_dir = None
            for i, a in enumerate(args):
                if a in ("-t", "--target-directory") and i + 1 < len(args):
                    target_dir = args[i + 1]
                elif a.startswith("--target-directory="):
                    target_dir = a.split("=", 1)[1]
            ops = _operands(args, ("-t", "--target-directory", "-S", "--suffix",
                                   "-m", "--mode", "-o", "--owner", "-g", "--group"))
            if target_dir is not None:
                sources, dest = ops, target_dir
            elif len(ops) >= 2:
                sources, dest = ops[:-1], ops[-1]
            else:
                continue
            for dp in _expand(dest, cwd):
                if dp.is_dir():
                    for src in sources:
                        cand = dp / os.path.basename(src.rstrip("/"))
                        if cand.is_file():
                            found.append(("overwrite", cand, False))
                elif dp.is_file():
                    found.append(("overwrite", dp, False))
    return found


#: git runs outside the cage here: nothing a repository configures may run
#: (fsmonitor, hooks), and no filter or diff driver is involved -- the
#: blob id of HEAD's copy is compared with the file's own, hashed raw.
_GIT = ["git", "-c", "core.fsmonitor=false", "-c", "core.hooksPath=/dev/null"]


def _git_holds(path: Path) -> bool:
    """Tracked by git and byte-identical to HEAD's copy (git gives it back)."""
    try:
        top = subprocess.run([*_GIT, "-C", str(path.parent), "rev-parse", "--show-toplevel"],
                             capture_output=True, text=True, timeout=10)
        if top.returncode != 0:
            return False
        root = top.stdout.strip()
        rel = os.path.relpath(path.resolve(), root)
        head = subprocess.run([*_GIT, "-C", root, "rev-parse", "--verify", "-q", f"HEAD:{rel}"],
                              capture_output=True, text=True, timeout=10)
        if head.returncode != 0:
            return False
        mine = subprocess.run([*_GIT, "-C", root, "hash-object", "--no-filters", "--", rel],
                              capture_output=True, text=True, timeout=60)
        return mine.returncode == 0 and mine.stdout.strip() == head.stdout.strip()
    except Exception:
        return False


def _own_output(path: Path, perms) -> bool:
    try:
        from .api_client import _DocToolExecutor
        return _DocToolExecutor._own_written_file(path, perms) is not None
    except Exception:
        return False


def _scratch(path: Path, perms) -> bool:
    try:
        from .api_client import _is_ephemeral_sink
        return bool(_is_ephemeral_sink(path, getattr(perms, "workspace", None)))
    except Exception:
        return False


def _files_under(root: Path) -> Optional[list[Path]]:
    """Regular files under ``root``; None when the walk ran past its bound."""
    out: list[Path] = []
    seen = 0
    deadline = time.monotonic() + WALK_MAX_SECONDS
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        seen += len(dirnames) + len(filenames)
        if seen > WALK_MAX_ENTRIES or time.monotonic() > deadline:
            return None
        for name in filenames:
            p = Path(dirpath) / name
            if not p.is_symlink():
                out.append(p)
    return out


def unrestorable(found: list[tuple[str, Path, bool]], perms) -> list[tuple[str, Path, str]]:
    """``(verb, path, size text)`` for each target nothing could restore."""
    per_file, budget = _journal_limits()
    held = 0
    bad: list[tuple[str, Path, str]] = []
    for verb, path, recursive in found:
        if path.is_symlink() and verb != "overwrite":
            continue
        if path.is_dir():
            if not recursive:
                continue                       # rm without -r refuses a directory itself
            files = _files_under(path)
            if files is None:
                bad.append((verb, path, f"more than {WALK_MAX_ENTRIES} entries"))
                continue
            names = getattr(path, "names", None)
            if names:
                import fnmatch
                files = [f for f in files if any(
                    fnmatch.fnmatch(f.name.lower(), pat.lower()) if fold
                    else fnmatch.fnmatchcase(f.name, pat) for pat, fold in names)]
        else:
            files = [path]
        for f in files:
            try:
                size = f.stat().st_size
            except OSError:
                continue
            in_git_dir = ".git" in f.parts
            if in_git_dir and f.name.endswith(".lock"):
                continue
            if size == 0 or _scratch(f, perms):
                continue
            if not in_git_dir and TRASH not in f.parts and size <= per_file \
                    and held + size <= budget:
                held += size
                continue
            if TRASH not in f.parts and (_git_holds(f) or _own_output(f, perms)):
                continue
            bad.append((verb, f, _size_text(size)))
    return bad


def _size_text(size: int) -> str:
    for unit in ("B", "KiB", "MiB", "GiB"):
        if size < 1024 or unit == "GiB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} GiB"


def _root_of(path: Path, perms) -> Path:
    try:
        roots = [Path(r).resolve() for r in perms.all_workspace_roots()]
    except Exception:
        roots = [Path(str(getattr(perms, "workspace", ".") or ".")).resolve()]
    rp = path.resolve()
    for r in sorted(roots, key=lambda r: -len(str(r))):
        if rp == r or str(rp).startswith(str(r).rstrip("/") + "/"):
            return r
    return roots[0] if roots else rp.parent


def prepare_trash(root: Path) -> Path:
    """``root/.delfin-trash/<UTC date>``, created, kept out of git, and the
    trash's expired days removed. Run by DELFIN, outside the cage."""
    trash = root / TRASH
    day = trash / datetime.now(timezone.utc).strftime("%Y-%m-%d")
    day.mkdir(parents=True, exist_ok=True)
    try:
        info = subprocess.run([*_GIT, "-C", str(root), "rev-parse", "--git-path", "info/exclude"],
                              capture_output=True, text=True, timeout=10)
        if info.returncode == 0:
            exclude = Path(info.stdout.strip())
            if not exclude.is_absolute():
                exclude = root / exclude
            text = exclude.read_text(encoding="utf-8") if exclude.exists() else ""
            if f"/{TRASH}/" not in text.splitlines():
                exclude.parent.mkdir(parents=True, exist_ok=True)
                with exclude.open("a", encoding="utf-8") as fh:
                    fh.write(("" if text.endswith("\n") or not text else "\n") + f"/{TRASH}/\n")
    except Exception:
        pass
    sweep(trash)
    return day


def sweep(trash: Path, days: int = TRASH_DAYS) -> None:
    """Remove the trash's day folders older than ``days``. Never raises."""
    try:
        cutoff = datetime.now(timezone.utc).timestamp() - days * 86400
        for child in trash.iterdir():
            try:
                stamp = datetime.strptime(child.name, "%Y-%m-%d").replace(
                    tzinfo=timezone.utc).timestamp()
            except ValueError:
                continue
            if stamp + 86400 < cutoff and child.is_dir() and not child.is_symlink():
                shutil.rmtree(child, ignore_errors=True)
    except Exception:
        pass


def refusal(cmd: str, cwd, perms) -> Optional[str]:
    """None, or why this bypass command is refused and what keeps the file."""
    try:
        found = destroyed(cmd, Path(cwd))
        if not found:
            return None
        bad = unrestorable(found, perms)
    except Exception:
        return None
    if not bad:
        return None
    verb, path, size = bad[0]
    root = _root_of(path, perms)
    try:
        day = prepare_trash(root)
    except Exception:
        day = root / TRASH / datetime.now(timezone.utc).strftime("%Y-%m-%d")
    more = f" (and {len(bad) - 1} more)" if len(bad) > 1 else ""
    try:
        from .api_client import _record_security_event
        _record_security_event("unrestorable_change", "bash",
                               f"{verb} {path} ({size}){more}"[:200], blocked=True)
    except Exception:
        pass
    return (f"refused in bypassPermissions: this command would {verb} "
            f"'{path}' ({size}){more}, and nothing could restore it afterwards "
            "(too large for the undo journal, not identical to a commit in "
            "git, not written by this session). Keep it instead: move it "
            f"aside with  mv {shlex.quote(str(path))} {shlex.quote(str(day))}/  "
            f"-- the trash keeps it {TRASH_DAYS} days, then DELFIN removes "
            "it -- and then run the command. To write a new version, give "
            "it a new name. Do not remove the file another way.")


__all__ = ["refusal", "destroyed", "unrestorable", "prepare_trash", "sweep",
           "TRASH", "TRASH_DAYS"]
