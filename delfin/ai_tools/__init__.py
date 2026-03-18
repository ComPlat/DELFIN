"""AI/ML tools for DELFIN.

Provides unified availability checks and lazy access to:

Foundation Models:
  - MoLFormer    — SMILES-based molecular embeddings & property prediction
  - Uni-Mol      — 3D-aware molecular representation & property prediction
  - ChemBERTa    — SMILES-based QSAR/ADMET classification

Generative Models:
  - REINVENT4    — goal-directed molecular design (RL + Transformer)
  - SyntheMol    — synthesisable molecule generation
  - GeoMol       — DL-based 3D conformer generation
  - torsional-diffusion — diffusion-based conformer generation
  - MatterGen    — crystal structure generation (diffusion)
  - CDVAE        — crystal diffusion variational autoencoder

Retrosynthesis & Reactions:
  - AiZynthFinder — retrosynthetic route planning
  - RXNMapper     — atom-mapping in chemical reactions
  - LocalRetro    — template-based retrosynthesis

Screening & ADMET:
  - DeepChem     — broad ML platform for life sciences
  - ADMETlab     — ADMET property prediction

Metal Complex ML:
  - molSimplify  — ML-based transition metal complex design
  - architector  — automated metal complex structure generation

Visualization:
  - plotly       — interactive plots for dashboards
"""

from __future__ import annotations

import importlib
import shutil
from typing import Optional


def _has_spec(module: str) -> bool:
    try:
        return importlib.util.find_spec(module) is not None
    except (ModuleNotFoundError, ValueError):
        return False


def _pkg_version(dist_name: str) -> Optional[str]:
    try:
        from importlib.metadata import version
        return version(dist_name)
    except Exception:
        return None


# ── Foundation Models ─────────────────────────────────────────────────

def molformer_available() -> bool:
    return _has_spec("transformers") and _has_spec("torch")

def get_molformer_version() -> Optional[str]:
    if not molformer_available():
        return None
    return _pkg_version("transformers") or "installed"

def unimol_available() -> bool:
    return _has_spec("unimol_tools")

def get_unimol_version() -> Optional[str]:
    if not unimol_available():
        return None
    return _pkg_version("unimol_tools") or "installed"

def chemberta_available() -> bool:
    return _has_spec("transformers") and _has_spec("torch")

def get_chemberta_version() -> Optional[str]:
    if not chemberta_available():
        return None
    return _pkg_version("transformers") or "installed"


# ── Generative Models ────────────────────────────────────────────────

def reinvent_available() -> bool:
    return _has_spec("reinvent")

def get_reinvent_version() -> Optional[str]:
    if not reinvent_available():
        return None
    return _pkg_version("reinvent") or _pkg_version("REINVENT4") or "installed"

def synthemol_available() -> bool:
    return _has_spec("synthemol")

def get_synthemol_version() -> Optional[str]:
    if not synthemol_available():
        return None
    return _pkg_version("synthemol") or "installed"

def geomol_available() -> bool:
    return _has_spec("geomol")

def get_geomol_version() -> Optional[str]:
    if not geomol_available():
        return None
    return _pkg_version("geomol") or "installed"

def torsional_diffusion_available() -> bool:
    return _has_spec("torsional_diffusion") or _has_spec("diffusion_hopping")

def get_torsional_diffusion_version() -> Optional[str]:
    if not torsional_diffusion_available():
        return None
    return _pkg_version("torsional-diffusion") or "installed"

def mattergen_available() -> bool:
    return _has_spec("mattergen")

def get_mattergen_version() -> Optional[str]:
    if not mattergen_available():
        return None
    return _pkg_version("mattergen") or "installed"

def cdvae_available() -> bool:
    return _has_spec("cdvae")

def get_cdvae_version() -> Optional[str]:
    if not cdvae_available():
        return None
    return _pkg_version("cdvae") or "installed"


# ── Retrosynthesis & Reactions ───────────────────────────────────────

def aizynthfinder_available() -> bool:
    return _has_spec("aizynthfinder")

def get_aizynthfinder_version() -> Optional[str]:
    if not aizynthfinder_available():
        return None
    return _pkg_version("aizynthfinder") or "installed"

def rxnmapper_available() -> bool:
    return _has_spec("rxnmapper")

def get_rxnmapper_version() -> Optional[str]:
    if not rxnmapper_available():
        return None
    return _pkg_version("rxnmapper") or "installed"

def localretro_available() -> bool:
    return _has_spec("localretro")

def get_localretro_version() -> Optional[str]:
    if not localretro_available():
        return None
    return _pkg_version("localretro") or "installed"


# ── Screening & ADMET ────────────────────────────────────────────────

def deepchem_available() -> bool:
    return _has_spec("deepchem")

def get_deepchem_version() -> Optional[str]:
    if not deepchem_available():
        return None
    return _pkg_version("deepchem") or "installed"

def admetlab_available() -> bool:
    return _has_spec("admetlab3") or _has_spec("admet")

def get_admetlab_version() -> Optional[str]:
    if not admetlab_available():
        return None
    return _pkg_version("admetlab3") or _pkg_version("admet") or "installed"


# ── Metal Complex ML ─────────────────────────────────────────────────

def molsimplify_available() -> bool:
    return _has_spec("molSimplify")

def get_molsimplify_version() -> Optional[str]:
    if not molsimplify_available():
        return None
    return _pkg_version("molSimplify") or "installed"

def architector_available() -> bool:
    return _has_spec("architector")

def get_architector_version() -> Optional[str]:
    if not architector_available():
        return None
    return _pkg_version("architector") or "installed"


# ── Visualization ────────────────────────────────────────────────────

def plotly_available() -> bool:
    return _has_spec("plotly")

def get_plotly_version() -> Optional[str]:
    if not plotly_available():
        return None
    return _pkg_version("plotly") or "installed"


# ── Summary ──────────────────────────────────────────────────────────

_TOOL_REGISTRY = [
    # (label, category, avail_fn, ver_fn, description, install_hint)
    ("MoLFormer", "Foundation Models", molformer_available, get_molformer_version,
     "SMILES-based molecular embeddings & property prediction",
     "pip install transformers torch"),
    ("Uni-Mol", "Foundation Models", unimol_available, get_unimol_version,
     "3D-aware molecular representation & property prediction",
     "pip install unimol_tools"),
    ("ChemBERTa", "Foundation Models", chemberta_available, get_chemberta_version,
     "SMILES-based QSAR/ADMET classification",
     "pip install transformers torch"),
    ("REINVENT4", "Generative", reinvent_available, get_reinvent_version,
     "Goal-directed molecular design (RL + Transformer)",
     "pip install reinvent"),
    ("SyntheMol", "Generative", synthemol_available, get_synthemol_version,
     "Synthesisable molecule generation",
     "pip install synthemol"),
    ("GeoMol", "Conformers", geomol_available, get_geomol_version,
     "DL-based 3D conformer generation",
     "pip install geomol"),
    ("torsional-diffusion", "Conformers", torsional_diffusion_available, get_torsional_diffusion_version,
     "Diffusion-based torsional conformer generation",
     "pip install torsional-diffusion"),
    ("MatterGen", "Crystal Generation", mattergen_available, get_mattergen_version,
     "Crystal structure generation via diffusion model",
     "pip install mattergen"),
    ("CDVAE", "Crystal Generation", cdvae_available, get_cdvae_version,
     "Crystal diffusion variational autoencoder",
     "pip install cdvae"),
    ("AiZynthFinder", "Retrosynthesis", aizynthfinder_available, get_aizynthfinder_version,
     "Retrosynthetic route planning (AstraZeneca)",
     "pip install aizynthfinder"),
    ("RXNMapper", "Retrosynthesis", rxnmapper_available, get_rxnmapper_version,
     "Atom-mapping in chemical reactions",
     "pip install rxnmapper"),
    ("LocalRetro", "Retrosynthesis", localretro_available, get_localretro_version,
     "Template-based retrosynthesis prediction",
     "pip install localretro"),
    ("DeepChem", "Screening", deepchem_available, get_deepchem_version,
     "Broad ML platform for life sciences",
     "pip install deepchem"),
    ("ADMETlab", "Screening", admetlab_available, get_admetlab_version,
     "ADMET property prediction",
     "pip install admetlab3"),
    ("molSimplify", "Metal Complex ML", molsimplify_available, get_molsimplify_version,
     "ML-based transition metal complex design & property prediction",
     "pip install molSimplify"),
    ("architector", "Metal Complex ML", architector_available, get_architector_version,
     "Automated metal complex structure generation",
     "pip install architector"),
    ("plotly", "Visualization", plotly_available, get_plotly_version,
     "Interactive plots for dashboards",
     "pip install plotly"),
]


def available_tools() -> list[str]:
    """Return names of all available AI tools."""
    return [label for label, _, avail_fn, *_ in _TOOL_REGISTRY if avail_fn()]


def collect_ai_summary() -> dict:
    """Return a dict summarising all AI tool status for dashboard/diagnostics."""
    tools_info = []
    for label, category, avail_fn, ver_fn, description, install_hint in _TOOL_REGISTRY:
        ok = avail_fn()
        tools_info.append({
            "name": label,
            "category": category,
            "installed": ok,
            "version": ver_fn() or "" if ok else "",
            "description": description,
            "install_hint": install_hint,
        })
    return {"tools": tools_info}
