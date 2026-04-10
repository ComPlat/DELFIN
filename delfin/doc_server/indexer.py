"""Document indexer — parses PDFs, Markdown, and text into a searchable JSON index."""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# PDF extraction
# ---------------------------------------------------------------------------

def _extract_pdf_text(path: Path) -> list[dict[str, Any]]:
    """Extract text from a PDF file, one entry per page.

    Returns a list of ``{"page": int, "text": str}`` dicts.
    """
    try:
        from pypdf import PdfReader  # type: ignore
    except ImportError:
        print(f"  [skip] pypdf not installed — cannot index {path.name}", file=sys.stderr)
        return []

    reader = PdfReader(str(path))
    pages: list[dict[str, Any]] = []
    for i, page in enumerate(reader.pages):
        text = (page.extract_text() or "").strip()
        if text:
            pages.append({"page": i + 1, "text": text})
    return pages


def _chunk_pdf_into_sections(pages: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Split extracted PDF pages into sections based on numbered headings.

    Detects patterns like ``1 Introduction``, ``9.3.2 DFT Grids``,
    ``A.1 Appendix``, etc.
    """
    # Heading pattern: start of line, optional number(s), title text
    heading_re = re.compile(
        r"^(\d+(?:\.\d+)*(?:\.\d+)?|[A-Z](?:\.\d+)*)\s+"  # number like 9.3.2 or A.1
        r"([A-Z][A-Za-z0-9\s\-\/:,()&]+)$",                # title in mixed case
        re.MULTILINE,
    )

    full_text = "\n\n".join(
        f"--- PAGE {p['page']} ---\n{p['text']}" for p in pages
    )

    # Find all heading positions
    matches = list(heading_re.finditer(full_text))

    if not matches:
        # No headings found — treat entire document as one section
        clean_text = re.sub(r"--- PAGE \d+ ---\n", "", full_text)
        return [{
            "section_id": "full_document",
            "title": "Full Document",
            "level": 0,
            "text": clean_text[:12000],
        }]

    sections: list[dict[str, Any]] = []
    for i, match in enumerate(matches):
        number = match.group(1)
        title = match.group(2).strip()
        level = number.count(".") + 1

        start = match.start()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(full_text)
        section_text = full_text[start:end].strip()

        # Remove page markers
        section_text = re.sub(r"--- PAGE \d+ ---\n?", "", section_text)

        # Generate a clean section ID
        section_id = re.sub(r"[^a-z0-9]+", "_", f"ch{number}_{title}".lower()).strip("_")
        section_id = section_id[:80]

        # Cap section size
        if len(section_text) > 12000:
            section_text = section_text[:12000] + "\n[... truncated]"

        sections.append({
            "section_id": section_id,
            "title": f"{number} {title}",
            "level": level,
            "text": section_text,
        })

    return sections


# ---------------------------------------------------------------------------
# Markdown / text extraction
# ---------------------------------------------------------------------------

def _extract_markdown_sections(path: Path) -> list[dict[str, Any]]:
    """Split a Markdown file into sections based on ``#`` headings."""
    text = path.read_text(encoding="utf-8", errors="replace")
    heading_re = re.compile(r"^(#{1,4})\s+(.+)$", re.MULTILINE)

    matches = list(heading_re.finditer(text))
    if not matches:
        return [{
            "section_id": "full_document",
            "title": path.stem,
            "level": 0,
            "text": text[:12000],
        }]

    sections: list[dict[str, Any]] = []
    for i, match in enumerate(matches):
        level = len(match.group(1))
        title = match.group(2).strip()
        start = match.start()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        section_text = text[start:end].strip()

        section_id = re.sub(r"[^a-z0-9]+", "_", title.lower()).strip("_")
        section_id = section_id[:80]

        if len(section_text) > 12000:
            section_text = section_text[:12000] + "\n[... truncated]"

        sections.append({
            "section_id": section_id,
            "title": title,
            "level": level,
            "text": section_text,
        })

    return sections


def _extract_text_sections(path: Path) -> list[dict[str, Any]]:
    """Read a plain text or RST file as a single section."""
    text = path.read_text(encoding="utf-8", errors="replace")
    if len(text) > 12000:
        text = text[:12000] + "\n[... truncated]"
    return [{
        "section_id": "full_document",
        "title": path.stem,
        "level": 0,
        "text": text,
    }]


# ---------------------------------------------------------------------------
# Index builder
# ---------------------------------------------------------------------------

_SUPPORTED_EXTENSIONS = {
    ".pdf": "pdf",
    ".md": "markdown",
    ".txt": "text",
    ".rst": "text",
}


def _discover_documents(literature_dir: Path) -> list[dict[str, str]]:
    """Discover indexable documents in the literature directory."""
    docs: list[dict[str, str]] = []
    for ext, doc_type in _SUPPORTED_EXTENSIONS.items():
        for path in sorted(literature_dir.rglob(f"*{ext}")):
            if path.name.startswith("."):
                continue
            # Generate doc_id from relative path
            rel = path.relative_to(literature_dir)
            doc_id = re.sub(r"[^a-z0-9]+", "_", str(rel.with_suffix("")).lower()).strip("_")
            docs.append({
                "path": str(path),
                "doc_id": doc_id,
                "title": path.stem.replace("_", " ").replace("-", " "),
                "type": doc_type,
            })
    return docs


def _discover_repo_docs(literature_dir: Path) -> list[dict[str, str]]:
    """Auto-discover DELFIN's own docs/ folder relative to the literature dir."""
    docs_dir = literature_dir.parent / "docs"
    if not docs_dir.is_dir():
        return []
    extra: list[dict[str, str]] = []
    for path in sorted(docs_dir.glob("*.md")):
        if path.name.startswith("."):
            continue
        doc_id = f"delfin_{path.stem}"
        extra.append({
            "path": str(path),
            "doc_id": doc_id,
            "title": f"DELFIN: {path.stem.replace('_', ' ')}",
            "type": "markdown",
        })
    return extra


def build_index(literature_dir: Path, extra_paths: list[dict[str, str]] | None = None) -> dict:
    """Build the search index from documents in the literature directory.

    Automatically includes DELFIN's own ``docs/`` folder as well.

    Parameters
    ----------
    literature_dir : Path
        Root of the literature folder (e.g., ``DELFIN/literature/``).
    extra_paths : list, optional
        Additional documents outside the literature folder.
        Each entry: ``{"path": ..., "doc_id": ..., "title": ...}``.

    Returns
    -------
    dict
        The complete index structure.
    """
    doc_specs = _discover_documents(literature_dir)

    # Auto-include DELFIN's own docs/ folder
    doc_specs.extend(_discover_repo_docs(literature_dir))

    if extra_paths:
        for spec in extra_paths:
            p = Path(spec["path"]).expanduser()
            ext = p.suffix.lower()
            doc_type = _SUPPORTED_EXTENSIONS.get(ext, "text")
            doc_specs.append({
                "path": str(p),
                "doc_id": spec.get("doc_id", p.stem),
                "title": spec.get("title", p.stem),
                "type": doc_type,
            })

    documents: dict[str, Any] = {}
    for spec in doc_specs:
        path = Path(spec["path"])
        if not path.exists():
            print(f"  [skip] not found: {path}", file=sys.stderr)
            continue

        doc_type = spec.get("type", _SUPPORTED_EXTENSIONS.get(path.suffix.lower(), "text"))
        print(f"  indexing: {path.name} ({doc_type})", file=sys.stderr)

        if doc_type == "pdf":
            pages = _extract_pdf_text(path)
            sections_list = _chunk_pdf_into_sections(pages)
        elif doc_type == "markdown":
            sections_list = _extract_markdown_sections(path)
        else:
            sections_list = _extract_text_sections(path)

        # Deduplicate section IDs
        seen_ids: set[str] = set()
        for s in sections_list:
            base_id = s["section_id"]
            sid = base_id
            counter = 2
            while sid in seen_ids:
                sid = f"{base_id}_{counter}"
                counter += 1
            s["section_id"] = sid
            seen_ids.add(sid)

        sections = {s["section_id"]: s for s in sections_list}
        total_chars = sum(len(s["text"]) for s in sections.values())

        documents[spec["doc_id"]] = {
            "title": spec["title"],
            "source_path": str(path),
            "source_type": doc_type,
            "section_count": len(sections),
            "total_chars": total_chars,
            "sections": sections,
        }

    return {
        "version": 1,
        "built_at": datetime.now(timezone.utc).isoformat(),
        "document_count": len(documents),
        "documents": documents,
    }


def get_default_index_path() -> Path:
    """Return the default index file path: ``~/.delfin/doc_index.json``."""
    return Path.home() / ".delfin" / "doc_index.json"


def get_default_literature_dir() -> Path | None:
    """Try to find the literature/ directory relative to the DELFIN repo."""
    # Walk up from this file to find the repo root
    current = Path(__file__).resolve()
    for parent in current.parents:
        candidate = parent / "literature"
        if candidate.is_dir():
            return candidate
    return None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    """CLI entry point: ``delfin-docs-index``."""
    parser = argparse.ArgumentParser(
        prog="delfin-docs-index",
        description="Build the DELFIN documentation search index from the literature/ folder.",
    )
    parser.add_argument(
        "--path",
        type=str,
        default=None,
        help="Path to the literature directory (default: auto-detect from repo).",
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None,
        help="Output path for the index JSON (default: ~/.delfin/doc_index.json).",
    )
    parser.add_argument(
        "--extra",
        type=str,
        nargs="*",
        default=[],
        help="Additional file paths to index (e.g., ~/my_papers/review.pdf).",
    )
    args = parser.parse_args(argv)

    # Resolve literature directory
    if args.path:
        lit_dir = Path(args.path).expanduser().resolve()
    else:
        lit_dir = get_default_literature_dir()
        if lit_dir is None:
            print("Error: Could not find literature/ directory. Use --path to specify.", file=sys.stderr)
            sys.exit(1)

    if not lit_dir.is_dir():
        print(f"Error: {lit_dir} is not a directory.", file=sys.stderr)
        sys.exit(1)

    # Resolve output path
    output_path = Path(args.output).expanduser() if args.output else get_default_index_path()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Extra paths
    extra = [{"path": p, "doc_id": Path(p).stem, "title": Path(p).stem} for p in (args.extra or [])]

    print(f"Indexing documents from: {lit_dir}", file=sys.stderr)
    index = build_index(lit_dir, extra_paths=extra or None)

    output_path.write_text(json.dumps(index, indent=2, ensure_ascii=False), encoding="utf-8")
    doc_count = index["document_count"]
    total_sections = sum(d["section_count"] for d in index["documents"].values())
    print(f"Index written to: {output_path}", file=sys.stderr)
    print(f"  {doc_count} document(s), {total_sections} section(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
