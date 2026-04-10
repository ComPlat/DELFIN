"""DELFIN documentation server — MCP-based search over local literature.

Provides a Model Context Protocol (MCP) server that indexes PDF, Markdown,
and text documents from the ``literature/`` folder and makes them searchable
for the DELFIN agent.

Quick start::

    delfin-docs-index          # build the search index
    python -m delfin.doc_server   # run the MCP server (normally spawned by the agent)
"""

from __future__ import annotations
