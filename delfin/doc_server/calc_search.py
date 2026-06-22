"""Search engine over the DELFIN calculation index.

Provides keyword search and structured filtering across all indexed
calculations (calc/, archive/, remote_archive/).

**Read-only**: never modifies any calculation data.
"""

from __future__ import annotations

import re
from typing import Any


class CalcSearchEngine:
    """Search indexed DELFIN calculations by keyword or structured query.

    The engine works on the index built by :mod:`calc_indexer`.  No
    external dependencies — pure Python string matching.
    """

    def __init__(self, index: dict) -> None:
        self._index = index
        self._records: list[dict[str, Any]] = index.get("records", [])

    @property
    def record_count(self) -> int:
        return len(self._records)

    # ------------------------------------------------------------------
    # Main search
    # ------------------------------------------------------------------

    def search(
        self,
        query: str = "",
        *,
        source: str = "",
        functional: str = "",
        basis_set: str = "",
        solvent: str = "",
        module: str = "",
        has_data: str = "",
        completed: bool | None = None,
        max_results: int = 20,
    ) -> list[dict[str, Any]]:
        """Search calculations by keyword and/or structured filters.

        Parameters
        ----------
        query : str
            Free-text keyword query (matched against all fields).
        source : str
            Filter by source: 'calc', 'archive', or 'remote_archive'.
        functional : str
            Filter by DFT functional (e.g. 'PBE0', 'B3LYP').
        basis_set : str
            Filter by basis set (e.g. 'def2-TZVP').
        solvent : str
            Filter by solvent name (e.g. 'toluene', 'DMF').
        module : str
            Filter by DELFIN module (e.g. 'ESD', 'GUPPY', 'IMAG').
        has_data : str
            Filter by result section (e.g. 'excited_states', 'esd_results').
        completed : bool or None
            Filter by completion status (True/False/None for any).
        max_results : int
            Maximum results to return.

        Returns
        -------
        list of dict
            Matching calculation records with a ``score`` field.
        """
        results: list[tuple[float, dict]] = []

        # Normalise query tokens
        query_tokens = query.lower().split() if query.strip() else []

        for rec in self._records:
            # --- Structured filters (hard filters) ---
            if source and rec.get("source", "") != source:
                continue
            if functional and not _ci_match(
                rec.get("functional", ""), functional
            ):
                continue
            if basis_set and not _ci_match(
                rec.get("basis_set", ""), basis_set
            ):
                continue
            if solvent and not _ci_match(rec.get("solvent", ""), solvent):
                continue
            if module:
                mods = rec.get("modules", [])
                if not any(_ci_match(m, module) for m in mods):
                    continue
            if has_data:
                hd = rec.get("has_data", [])
                if not any(_ci_match(h, has_data) for h in hd):
                    continue
            if completed is not None and rec.get("completed") != completed:
                continue

            # --- Keyword scoring ---
            if query_tokens:
                search_text = rec.get("_search_text", "")
                score = _score_match(query_tokens, search_text)
                if score < 0.01:
                    continue
            else:
                score = 1.0  # No query → all pass with equal score

            results.append((score, rec))

        # Sort by score descending
        results.sort(key=lambda x: x[0], reverse=True)

        # Build output
        output: list[dict[str, Any]] = []
        for score, rec in results[:max_results]:
            entry = _format_result(rec)
            entry["score"] = round(score, 4)
            output.append(entry)

        return output

    def get_calc_info(self, calc_id: str) -> dict[str, Any] | None:
        """Get detailed info for a specific calculation by name/ID.

        Searches by exact calc_id first, then by substring match.
        """
        # Exact match
        for rec in self._records:
            if rec.get("calc_id", "") == calc_id:
                return _format_result_detailed(rec)

        # Substring match
        q = calc_id.lower()
        for rec in self._records:
            if q in rec.get("calc_id", "").lower():
                return _format_result_detailed(rec)

        return None

    def list_functionals(self) -> list[tuple[str, int]]:
        """List all functionals with their occurrence count."""
        return _count_field(self._records, "functional")

    def list_basis_sets(self) -> list[tuple[str, int]]:
        """List all basis sets with their occurrence count."""
        return _count_field(self._records, "basis_set")

    def list_solvents(self) -> list[tuple[str, int]]:
        """List all solvents with their occurrence count."""
        return _count_field(self._records, "solvent")

    def list_modules(self) -> list[tuple[str, int]]:
        """List all DELFIN modules with their occurrence count."""
        counts: dict[str, int] = {}
        for rec in self._records:
            for m in rec.get("modules", []):
                if m:
                    counts[m] = counts.get(m, 0) + 1
        return sorted(counts.items(), key=lambda x: x[1], reverse=True)

    def summary(self) -> dict[str, Any]:
        """Return a summary of the index contents."""
        sources: dict[str, int] = {}
        for rec in self._records:
            s = rec.get("source", "unknown")
            sources[s] = sources.get(s, 0) + 1
        return {
            "total_calculations": len(self._records),
            "by_source": sources,
            "functionals": self.list_functionals()[:10],
            "basis_sets": self.list_basis_sets()[:10],
            "solvents": self.list_solvents()[:10],
            "modules": self.list_modules(),
        }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ci_match(value: str, pattern: str) -> bool:
    """Case-insensitive substring match."""
    return pattern.lower() in value.lower()


def _score_match(tokens: list[str], text: str) -> float:
    """Score how well query tokens match the search text."""
    if not tokens or not text:
        return 0.0
    matched = sum(1 for t in tokens if t in text)
    return matched / len(tokens)


def _count_field(records: list[dict], field: str) -> list[tuple[str, int]]:
    counts: dict[str, int] = {}
    for rec in records:
        val = rec.get(field, "")
        if val and isinstance(val, str) and val.strip():
            counts[val.strip()] = counts.get(val.strip(), 0) + 1
    return sorted(counts.items(), key=lambda x: x[1], reverse=True)


def _format_result(rec: dict[str, Any]) -> dict[str, Any]:
    """Format a record for search result output (compact)."""
    return {
        "calc_id": rec.get("calc_id", ""),
        "source": rec.get("source", ""),
        "rel_path": rec.get("rel_path", ""),
        "functional": rec.get("functional", ""),
        "basis_set": rec.get("basis_set", ""),
        "solvent": rec.get("solvent", ""),
        "charge": rec.get("charge", ""),
        "modules": rec.get("modules", []),
        "completed": rec.get("completed"),
        "smiles": rec.get("smiles", "")[:80],
    }


def _format_result_detailed(rec: dict[str, Any]) -> dict[str, Any]:
    """Format a record with all available details."""
    result = {
        "calc_id": rec.get("calc_id", ""),
        "path": rec.get("path", ""),
        "source": rec.get("source", ""),
        "index_source": rec.get("index_source", ""),
        "functional": rec.get("functional", ""),
        "basis_set": rec.get("basis_set", ""),
        "aux_basis": rec.get("aux_basis", ""),
        "ri_method": rec.get("ri_method", ""),
        "dispersion": rec.get("dispersion", ""),
        "solvation": rec.get("solvation", ""),
        "solvent": rec.get("solvent", ""),
        "charge": rec.get("charge", ""),
        "smiles": rec.get("smiles", ""),
        "name": rec.get("name", ""),
        "xtb_method": rec.get("xtb_method", ""),
        "modules": rec.get("modules", []),
        "has_data": rec.get("has_data", []),
        "out_files": rec.get("out_files", []),
        "completed": rec.get("completed"),
        "exit_code": rec.get("exit_code"),
    }
    # Energies
    for key in ("gibbs_energy", "electronic_energy", "zpe", "run_time_seconds"):
        val = rec.get(key)
        if val is not None:
            result[key] = val
    # ESD
    if rec.get("esd_modus"):
        result["esd_modus"] = rec["esd_modus"]
    # ORCA keywords
    if rec.get("orca_keywords"):
        result["orca_keywords"] = rec["orca_keywords"]
    return result
