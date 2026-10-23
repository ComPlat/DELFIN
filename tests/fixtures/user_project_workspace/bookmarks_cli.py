"""Entry point of the toy bookmark keeper.

Prints how many bookmarks the store holds; ``--tags`` adds the tag list. The
script inserts its own directory into ``sys.path`` so it runs from any
working directory, which keeps a hand-written launcher script trivial.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from bookmark_store import collect_tags, load_bookmarks  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    """Print the bookmark count and, on ``--tags``, the tags. Returns 0."""
    args = list(sys.argv[1:] if argv is None else argv)
    bookmarks = load_bookmarks()
    print(f"bookmarks: {len(bookmarks)}")
    if "--tags" in args:
        print("tags: " + ", ".join(collect_tags(bookmarks)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
