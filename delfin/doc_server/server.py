"""MCP server for the DELFIN documentation search.

Exposes four tools via the Model Context Protocol (stdio transport):

- ``search_docs``   — TF-IDF search across all indexed documents
- ``read_section``   — read the full text of a specific section
- ``list_docs``      — list all indexed documents
- ``list_sections``  — list sections / table of contents of a document
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from .indexer import get_default_index_path
from .search import DocSearchEngine


def _load_index(index_path: str | Path) -> dict:
    """Load the JSON index from disk."""
    p = Path(index_path)
    if not p.exists():
        print(f"Error: Index file not found: {p}", file=sys.stderr)
        print("Run 'delfin-docs-index' to build the index.", file=sys.stderr)
        sys.exit(1)
    return json.loads(p.read_text(encoding="utf-8"))


def run_server(argv: list[str] | None = None) -> None:
    """Start the MCP server (stdio transport)."""
    parser = argparse.ArgumentParser(prog="delfin-doc-server")
    parser.add_argument(
        "--index",
        type=str,
        default=None,
        help="Path to the doc index JSON (default: ~/.delfin/doc_index.json).",
    )
    args = parser.parse_args(argv)

    index_path = Path(args.index) if args.index else get_default_index_path()
    index = _load_index(index_path)
    engine = DocSearchEngine(index)

    # --- MCP server setup via FastMCP ---
    from mcp.server.fastmcp import FastMCP

    mcp = FastMCP(
        "delfin-docs",
        instructions=(
            "DELFIN documentation server. Search and read ORCA manuals, "
            "xTB/CREST documentation, methodology docs, and scientific papers "
            "that have been indexed from the literature/ folder."
        ),
    )

    @mcp.tool()
    def search_docs(
        query: str,
        doc_filter: str = "",
        max_results: int = 10,
    ) -> str:
        """Search indexed documentation for sections matching a query.

        Use this to find relevant information about ORCA keywords, DFT methods,
        basis sets, xTB options, CREST workflows, and other computational
        chemistry topics.

        Args:
            query: Free-text search query (e.g., "RIJCOSX approximation",
                   "broken symmetry DFT", "def2-TZVP basis set")
            doc_filter: Optional doc_id to restrict search to one document
            max_results: Maximum number of results (default: 10)

        Returns:
            JSON array of matching sections with doc_id, section_id, title,
            score, and a text snippet.
        """
        results = engine.search(query, doc_filter=doc_filter, max_results=max_results)
        return json.dumps(results, indent=2, ensure_ascii=False)

    @mcp.tool()
    def read_section(doc_id: str, section_id: str) -> str:
        """Read the full text of a specific section from an indexed document.

        Use this after search_docs to read a section in detail.

        Args:
            doc_id: Document identifier (from search_docs results)
            section_id: Section identifier (from search_docs results)

        Returns:
            The full section text, or an error message if not found.
        """
        doc = index.get("documents", {}).get(doc_id)
        if not doc:
            available = list(index.get("documents", {}).keys())
            return f"Document '{doc_id}' not found. Available: {available}"

        section = doc.get("sections", {}).get(section_id)
        if not section:
            available = list(doc.get("sections", {}).keys())[:20]
            return (
                f"Section '{section_id}' not found in '{doc_id}'. "
                f"First sections: {available}"
            )

        return (
            f"# {section.get('title', section_id)}\n"
            f"Source: {doc.get('title', doc_id)}\n\n"
            f"{section.get('text', '')}"
        )

    @mcp.tool()
    def list_docs() -> str:
        """List all indexed documents with metadata.

        Returns:
            JSON array of documents with doc_id, title, source_path,
            section_count, and total_chars.
        """
        docs = []
        for doc_id, doc in index.get("documents", {}).items():
            docs.append({
                "doc_id": doc_id,
                "title": doc.get("title", ""),
                "source_path": doc.get("source_path", ""),
                "source_type": doc.get("source_type", ""),
                "section_count": doc.get("section_count", 0),
                "total_chars": doc.get("total_chars", 0),
            })
        return json.dumps(docs, indent=2, ensure_ascii=False)

    @mcp.tool()
    def list_sections(doc_id: str) -> str:
        """List all sections (table of contents) of a specific document.

        Args:
            doc_id: Document identifier (from list_docs)

        Returns:
            JSON array of sections with section_id, title, level, and char_count.
        """
        doc = index.get("documents", {}).get(doc_id)
        if not doc:
            available = list(index.get("documents", {}).keys())
            return f"Document '{doc_id}' not found. Available: {available}"

        sections = []
        for section_id, section in doc.get("sections", {}).items():
            sections.append({
                "section_id": section_id,
                "title": section.get("title", ""),
                "level": section.get("level", 0),
                "char_count": len(section.get("text", "")),
            })
        return json.dumps(sections, indent=2, ensure_ascii=False)

    # -- Calculation search tools ---
    # Build calc index lazily on first call
    _calc_engine_cache: dict[str, Any] = {}

    def _get_calc_engine():
        if "engine" not in _calc_engine_cache:
            from .calc_indexer import build_calc_index
            from .calc_search import CalcSearchEngine
            calc_dir = Path.home() / "calc"
            archive_dir = Path.home() / "archive"
            idx = build_calc_index(
                calc_dir=calc_dir if calc_dir.is_dir() else None,
                archive_dir=archive_dir if archive_dir.is_dir() else None,
                quiet=True,
            )
            _calc_engine_cache["engine"] = CalcSearchEngine(idx)
        return _calc_engine_cache["engine"]

    @mcp.tool()
    def search_calcs(
        query: str = "",
        source: str = "",
        functional: str = "",
        basis_set: str = "",
        solvent: str = "",
        module: str = "",
        max_results: int = 20,
    ) -> str:
        """Search DELFIN calculations across calc/, archive/, and remote_archive/.

        Find calculations by keyword or structured filters. Searches method,
        basis set, solvent, molecule name, DELFIN modules, and more.

        Args:
            query: Free-text keyword query (e.g. 'PBE0 def2-TZVP', 'TDDFT toluene')
            source: Filter by source: 'calc', 'archive', or 'remote_archive'
            functional: Filter by DFT functional (e.g. 'PBE0', 'B3LYP')
            basis_set: Filter by basis set (e.g. 'def2-TZVP')
            solvent: Filter by solvent (e.g. 'toluene', 'DMF')
            module: Filter by DELFIN module (e.g. 'ESD', 'GUPPY', 'IMAG')
            max_results: Maximum results to return (default: 20)

        Returns:
            JSON array of matching calculations with metadata.
        """
        eng = _get_calc_engine()
        results = eng.search(
            query=query, source=source, functional=functional,
            basis_set=basis_set, solvent=solvent, module=module,
            max_results=max_results,
        )
        return json.dumps(results, indent=2, ensure_ascii=False)

    @mcp.tool()
    def get_calc_info(calc_id: str) -> str:
        """Get detailed information about a specific DELFIN calculation.

        Returns functional, basis set, solvent, charge, SMILES, energies,
        modules, output files, completion status, and more.

        Args:
            calc_id: Calculation name or ID (e.g. 'Emitter8_CAMB3LYP_ma-def2-TZVP')

        Returns:
            JSON object with all available metadata, or error if not found.
        """
        eng = _get_calc_engine()
        info = eng.get_calc_info(calc_id)
        if info is None:
            return json.dumps({"error": f"Calculation '{calc_id}' not found."})
        return json.dumps(info, indent=2, ensure_ascii=False)

    @mcp.tool()
    def calc_summary() -> str:
        """Get a summary of all indexed DELFIN calculations.

        Returns total count, breakdown by source, most-used functionals,
        basis sets, solvents, and DELFIN modules.
        """
        eng = _get_calc_engine()
        return json.dumps(eng.summary(), indent=2, ensure_ascii=False)

    # Run the server on stdio
    mcp.run(transport="stdio")
