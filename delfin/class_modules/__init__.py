"""Class-specialist modules — chemistry-domain-specific SMILES→XYZ converters.

Each specialist handles one (coord_class, metal_block) region of the
SMILES space.  Specialists are cherry-picked / forensik-ported from
historical champion commits per empirical evidence in
quality_framework/results/champion_routing.md.

Public API:
    get_specialist(features: ClassFeatures) -> Specialist
    list_registered() -> List[SpecialistID]

Per nature_project/15_HYBRID_PATH_FINAL.md Phase 2.
"""
from __future__ import annotations

from typing import Dict, Optional, Protocol, Tuple

from delfin.classify import ClassFeatures


class Specialist(Protocol):
    """Protocol every class-specialist must implement."""

    id: str
    """Unique identifier (e.g. 'multi_sigma_123a130')."""

    source_commit: str
    """Champion git commit this specialist was forensik-ported from."""

    def convert(
        self,
        smiles: str,
        **kwargs,
    ) -> Tuple[Optional[str], Optional[str]]:
        """Convert SMILES → XYZ.

        Returns (xyz, error) tuple matching smiles_converter.smiles_to_xyz API.
        """
        ...


# Registry: SpecialistID → Specialist instance
_REGISTRY: Dict[str, Specialist] = {}


def register(specialist: Specialist) -> None:
    """Register a specialist in the global registry."""
    if specialist.id in _REGISTRY:
        raise ValueError(f"Specialist {specialist.id} already registered")
    _REGISTRY[specialist.id] = specialist


def get_specialist(specialist_id: str) -> Optional[Specialist]:
    """Return registered specialist by ID, or None if not found."""
    return _REGISTRY.get(specialist_id)


def list_registered() -> list[str]:
    """Return list of all registered specialist IDs."""
    return sorted(_REGISTRY.keys())
