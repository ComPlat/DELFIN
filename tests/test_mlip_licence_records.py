"""Licence records for third-party ML potentials (LICENSE_GUARD.md family 4).

Guarantees:
1. Every registry backend carries a code-licence AND a weight-licence record,
   each with a source URL, and each marked verified/unverified.
2. The tracked tree contains no weight/checkpoint file and no token-like
   string (HF token, api key patterns).

Red-first: written before the registry gained the licence fields.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from delfin.mlp_tools import registry  # noqa: E402

WEIGHT_SUFFIXES = (".pt", ".pth", ".ckpt", ".safetensors", ".npz")

TOKEN_RE_STRINGS = (
    "hf_[A-Za-z0-9]{30,}",       # HuggingFace user access token literal
    "github_pat_[A-Za-z0-9_]{20,}",
)


def test_every_backend_has_code_and_weight_licence():
    missing = []
    for name in registry.backend_names():
        spec = registry.get_backend(name)
        for field, what in (("code_licence", "code"), ("weight_licence", "weight")):
            record = getattr(spec, field)
            if not isinstance(record, dict):
                missing.append(f"{name}: {field} is not a dict")
                continue
            for key in ("name", "url", "verified"):
                if key not in record:
                    missing.append(f"{name}: {field}.{key} missing")
            if record.get("name") == "unknown":
                missing.append(f"{name}: {field} name unknown")
    assert not missing, "incomplete licence records:\n" + "\n".join(missing)


def test_licence_records_carry_source_urls():
    bad = []
    for name in registry.backend_names():
        for field in ("code_licence", "weight_licence"):
            record = getattr(registry.get_backend(name), field)
            url = record.get("url", "")
            if not url.startswith("http"):
                bad.append(f"{name}/{field}: {url!r}")
    assert not bad, "non-URL licence sources:\n" + "\n".join(bad)


def test_weight_licence_flags_noncommercial():
    """Non-commercial weight licences must be flagged so a filter can refuse them."""
    omol = registry.get_backend("mace_omol")
    assert omol.weight_licence["non_commercial"] is True
    uma = registry.get_backend("fairchem_uma")
    assert uma.weight_licence["non_commercial"] is False
    assert uma.weight_licence["gated"] is True


def test_no_weight_or_checkpoint_files_tracked():
    out = subprocess.run(
        ["git", "ls-files"], cwd=REPO, capture_output=True, text=True, check=True
    )
    bad = [
        line for line in out.stdout.splitlines()
        if line.lower().endswith(WEIGHT_SUFFIXES)
    ]
    assert not bad, "weight/checkpoint files tracked in the tree:\n" + "\n".join(bad)


def test_no_token_like_strings_tracked():
    import re

    out = subprocess.run(
        ["git", "ls-files"], cwd=REPO, capture_output=True, text=True, check=True
    )
    patterns = [re.compile(p) for p in TOKEN_RE_STRINGS]
    bad = []
    for rel in out.stdout.splitlines():
        path = REPO / rel
        if not path.is_file() or path.stat().st_size > 2_000_000:
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for pat in patterns:
            if pat.search(text):
                bad.append(f"{rel}: matches {pat.pattern}")
    assert not bad, "token-like strings tracked:\n" + "\n".join(bad)
