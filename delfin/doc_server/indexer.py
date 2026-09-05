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

def _extract_pdf_text(
    path: Path, quiet: bool = False, report: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Extract text from a PDF file, one entry per page.

    Returns a list of ``{"page": int, "text": str}`` dicts. When nothing
    could be extracted, *report* (if given) is filled with ``reason`` — the
    caller needs to tell "no pypdf" from "scanned or DRM-protected PDF"
    apart, because both produce the same empty list and neither is a
    document.
    """
    try:
        from pypdf import PdfReader  # type: ignore
    except ImportError:
        if not quiet:
            print(f"  [skip] pypdf not installed — cannot index {path.name}", file=sys.stderr)
        if report is not None:
            report["reason"] = (
                "pypdf is not installed, so no text could be read from this "
                "PDF (pip install pypdf)"
            )
        return []

    try:
        reader = PdfReader(str(path))
        page_count = len(reader.pages)
    except Exception as exc:            # unreadable / encrypted / corrupt
        if report is not None:
            report["reason"] = f"PDF could not be opened: {type(exc).__name__}: {exc}"
        return []

    pages: list[dict[str, Any]] = []
    for i, page in enumerate(reader.pages):
        try:
            text = (page.extract_text() or "").strip()
        except Exception:
            text = ""
        if text:
            pages.append({"page": i + 1, "text": text})
    if not pages and report is not None:
        report["reason"] = (
            f"no text on any of the {page_count} page(s) — a scanned or "
            "image-only PDF needs OCR before it can be searched"
        )
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

    if not pages:
        # No pages with text is NOT a one-section document with an empty
        # body. It used to be: the chunker emitted ``full_document`` with
        # ``text=""``, the document reported ``section_count: 1``, the "is
        # the manual indexed" check passed on the title alone, and every
        # search against it returned nothing, forever. The caller records
        # this as an indexing failure instead.
        return []

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


def build_index(
    literature_dir: Path,
    extra_paths: list[dict[str, str]] | None = None,
    quiet: bool = False,
) -> dict:
    """Build the search index from documents in the literature directory.

    Automatically includes DELFIN's own ``docs/`` folder as well.

    Parameters
    ----------
    literature_dir : Path
        Root of the literature folder (e.g., ``DELFIN/literature/``).
    extra_paths : list, optional
        Additional documents outside the literature folder.
        Each entry: ``{"path": ..., "doc_id": ..., "title": ...}``.
    quiet : bool
        Suppress progress output (used when called from the dashboard).

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
    failed: list[dict[str, Any]] = []

    def _fail(spec: dict, path: Path, reason: str) -> None:
        failed.append({
            "doc_id": spec.get("doc_id", ""),
            "title": spec.get("title", ""),
            "source_path": str(path),
            "reason": reason,
        })
        if not quiet:
            print(f"  [fail] {path.name}: {reason}", file=sys.stderr)

    for spec in doc_specs:
        path = Path(spec["path"])
        if not path.exists():
            if not quiet:
                print(f"  [skip] not found: {path}", file=sys.stderr)
            _fail(spec, path, "file not found")
            continue

        doc_type = spec.get("type", _SUPPORTED_EXTENSIONS.get(path.suffix.lower(), "text"))
        if not quiet:
            print(f"  indexing: {path.name} ({doc_type})", file=sys.stderr)

        report: dict[str, Any] = {}
        if doc_type == "pdf":
            pages = _extract_pdf_text(path, quiet=quiet, report=report)
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

        # Zero extracted characters is an indexing FAILURE, not a document
        # with nothing in it. Recorded WITH the reason, and kept out of
        # ``documents`` so nothing downstream can report it as present:
        # a present-and-empty entry answered "yes, indexed" to the manual
        # check and returned nothing to every search, permanently.
        if total_chars == 0:
            _fail(spec, path, report.get("reason")
                  or "no text could be extracted from this file")
            continue

        documents[spec["doc_id"]] = {
            "title": spec["title"],
            "source_path": str(path),
            "source_type": doc_type,
            "section_count": len(sections),
            "total_chars": total_chars,
            "source_mtime": _mtime(path),
            "sections": sections,
        }

    return {
        "version": 1,
        "built_at": datetime.now(timezone.utc).isoformat(),
        "document_count": len(documents),
        "documents": documents,
        # Named, with the reason. A file that produced nothing has to be
        # visible somewhere or the only symptom is a search that never
        # matches.
        "failed_documents": failed,
    }


def _mtime(path: Path) -> float:
    """Source mtime, or 0.0 when it cannot be read."""
    try:
        return float(path.stat().st_mtime)
    except OSError:
        return 0.0


def stale_documents(index: dict) -> list[dict[str, Any]]:
    """Indexed documents whose source changed or vanished since the build.

    ``built_at`` was written by both indexers and compared to nothing, so
    an answer from a section deleted from the source last month was
    byte-identical to a fresh hit. Every consumer that hands index content
    to a reader stamps this alongside it.
    """
    out: list[dict[str, Any]] = []
    for doc_id, doc in (index or {}).get("documents", {}).items():
        if not isinstance(doc, dict):
            continue
        src = doc.get("source_path") or ""
        if not src:
            continue
        path = Path(src)
        if not path.exists():
            out.append({"doc_id": doc_id, "source_path": src,
                        "reason": "source file no longer exists"})
            continue
        recorded = doc.get("source_mtime")
        if not recorded:
            # Indexed before mtimes were recorded — unknowable, not stale.
            continue
        current = _mtime(path)
        if current > float(recorded) + 1.0:
            out.append({"doc_id": doc_id, "source_path": src,
                        "reason": "source file changed after the index was built"})
    return out


def index_provenance(index: dict) -> dict[str, Any]:
    """What the caller needs to judge an answer from this index.

    Carried on every search and listing payload, INCLUDING an empty one:
    "nothing matched" and "there is nothing to match against" are the same
    sentence otherwise.
    """
    index = index or {}
    docs = index.get("documents", {}) or {}
    failures = [f for f in (index.get("failed_documents", []) or [])
                if isinstance(f, dict)]
    return {
        "built_at": index.get("built_at", ""),
        "document_count": len(docs),
        "section_count": sum(
            len((d or {}).get("sections", {}) or {}) for d in docs.values()
            if isinstance(d, dict)),
        "failed_documents": [
            {"doc_id": f.get("doc_id", ""), "reason": f.get("reason", "")}
            for f in failures
        ],
        "stale_documents": stale_documents(index),
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
    for entry in index.get("failed_documents", []):
        print(f"  NOT indexed: {entry['source_path']} — {entry['reason']}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
