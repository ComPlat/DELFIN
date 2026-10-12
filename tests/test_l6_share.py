"""Package L6: an approved skill can be shared, and arrives as a proposal.

L6 exports an ACCEPTED skill (one that lives in ``~/.delfin/skills/<name>/``)
to a local shared folder as one JSON bundle and imports one back as a
PROPOSAL (source "shared") that is never active without a human accept.

Design rules (task + QS ruling nacht-s29):
* no PR path, no network;
* the export scrub REUSES ``skill_safety.check`` (the command + secret scan
  DELFIN already has) on top of a *targeted* absolute-path / user@host /
  job-id scan (those are not covered by the generic scan);
* an import MUST pass non-empty evidence so the safety gate inside propose()
  runs and unsafe text is STORED as status "blocked" (never active). The
  default evidence is a NON-fabricated provenance record -- kind "recipe",
  ref = the bundle's sha256 -- not a fake green test/calc/job run;
* the accepted skill home and the proposal store are resolved per call from
  ``Path.home()`` so a redirected HOME redirects them and no test ever
  touches the real ``~/.delfin``.
"""
import json
from pathlib import Path

import pytest

from delfin.agent import skill_proposals as sp
import delfin.agent.skill_share as ss
from delfin.agent.skills import discover_skills


@pytest.fixture
def home(tmp_path, monkeypatch):
    """Redirect HOME so the live skills tree and proposal store are ours."""
    monkeypatch.setenv("HOME", str(tmp_path))
    return tmp_path


def _live_dir(home):
    return Path(home) / ".delfin" / "skills"


def _install_accepted(home, name="shared-skill", text="# A skill\nbody",
                      version="v1"):
    """Simulate a human-accepted skill: SKILL.md + proposal.json in the live
    skills tree exactly as ``accept`` leaves them."""
    d = _live_dir(home) / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "SKILL.md").write_text(text, encoding="utf-8")
    meta = {
        "name": name, "text": text, "evidence": [], "source": "local",
        "status": "accepted", "findings": [], "created": "2026-10-01",
        "base_version": version,
    }
    (d / "proposal.json").write_text(json.dumps(meta), encoding="utf-8")
    return d


def _bundle_path(tmp_path, name="shared-skill"):
    return Path(tmp_path) / f"{name}.skill.json"


def _sha256(payload):
    import hashlib
    return hashlib.sha256(json.dumps(payload, sort_keys=True,
                                     ensure_ascii=False).encode()).hexdigest()


def _make_valid_bundle(dirpath, *, name="shared-skill", text="# A skill\nbody",
                       version="v1", author="kb"):
    """Build a valid bundle the way export_bundle does."""
    payload = {"name": name, "text": text, "version": version,
               "author": author}
    bundle = {"sha256": _sha256(payload), **payload}
    target = dirpath / f"{name}.skill.json"
    target.write_text(json.dumps(bundle, ensure_ascii=False), encoding="utf-8")
    return target


# --------------------------------------------------------------------------
# export_bundle
# --------------------------------------------------------------------------

def test_export_writes_one_json_bundle_with_text_version_author_sha(
        home, tmp_path):
    _install_accepted(home, text="# A skill\nbody", version="v3")
    out = ss.export_bundle("shared-skill", str(tmp_path), author="kb")
    bundle = _bundle_path(tmp_path)
    assert bundle.is_file()
    data = json.loads(bundle.read_text(encoding="utf-8"))
    assert data["name"] == "shared-skill"
    assert data["text"] == "# A skill\nbody"
    assert data["version"] == "v3"
    assert data["author"] == "kb"
    assert isinstance(data["sha256"], str) and len(data["sha256"]) == 64
    # returns exactly the text that left this machine
    assert out == "# A skill\nbody"


def test_export_of_a_non_accepted_skill_is_refused(home, tmp_path):
    # a name with no live SKILL.md is not an accepted skill
    with pytest.raises(ValueError):
        ss.export_bundle("ghost", str(tmp_path))
    assert not _bundle_path(tmp_path).exists()


def test_export_carries_no_author_when_omitted(home, tmp_path):
    _install_accepted(home, text="# A skill\nbody")
    ss.export_bundle("shared-skill", str(tmp_path))
    data = json.loads(_bundle_path(tmp_path).read_text(encoding="utf-8"))
    assert data["author"] == ""


def test_export_detects_a_version_even_without_proposal_json(home, tmp_path):
    d = _live_dir(home) / "plain"
    d.mkdir(parents=True, exist_ok=True)
    (d / "SKILL.md").write_text("# plain\nbody", encoding="utf-8")
    ss.export_bundle("plain", str(tmp_path))
    data = json.loads(_bundle_path(tmp_path, "plain").read_text(encoding="utf-8"))
    assert data["version"] == ""


# --------------------------------------------------------------------------
# export scrub -- blocks, names the line, writes nothing
# --------------------------------------------------------------------------

def test_export_blocked_on_absolute_home_path_names_the_line(home, tmp_path):
    _install_accepted(home, text="# S\nuses /home/alice/refs\nmore\n")
    with pytest.raises(Exception) as exc:
        ss.export_bundle("shared-skill", str(tmp_path))
    msg = str(exc.value)
    assert "/home/alice/refs" in msg and "line 2" in msg
    assert not _bundle_path(tmp_path).exists()


def test_export_blocked_on_a_job_id(home, tmp_path):
    _install_accepted(home, text="# S\nslurm job 9876543 finished\n")
    with pytest.raises(Exception) as exc:
        ss.export_bundle("shared-skill", str(tmp_path))
    assert not _bundle_path(tmp_path).exists()


def test_export_routes_through_the_reused_skill_safety_scan(
        home, tmp_path, monkeypatch):
    """The scrub reuses DELFIN's own skill_safety.check (the command+secret
    scan) rather than a second scan: a finding it returns BLOCKS the export.
    Proved by watching the call count and the blocking effect."""
    from delfin.agent import skill_safety as real_ss
    calls = []

    def fake_check(text):
        calls.append(text)
        return ["safety finding: blocked by the reused scan"]

    monkeypatch.setattr(real_ss, "check", fake_check)
    _install_accepted(home, text="# A skill\nbody")
    with pytest.raises(Exception) as exc:
        ss.export_bundle("shared-skill", str(tmp_path))
    assert calls, "export must consult the reused skill_safety scan"
    assert not _bundle_path(tmp_path).exists()


# --------------------------------------------------------------------------
# import_bundle -- bounds, no proposal of its own
# --------------------------------------------------------------------------

def test_import_clean_bundle_is_verified_but_creates_no_proposal(
        home, tmp_path):
    """QS corrected ruling: import_bundle does NOT create a proposal. A clean,
    hash-verified bundle returns status "verified" and nothing is stored in
    the proposal store."""
    bundle = _make_valid_bundle(dirpath=tmp_path)
    info = ss.import_bundle(bundle)
    assert info["status"] == "verified"
    assert info["name"] == "shared-skill"
    assert info["text"] == "# A skill\nbody"
    assert sp.get_proposal("shared-skill") is None
    assert "shared-skill" not in [s.name for s in discover_skills()]


def test_import_unsafe_text_is_stored_blocked_naming_the_line(home, tmp_path):
    """Unsafe text (git push) is recorded by import_bundle itself as status
    "blocked", naming the line -- the package's "stored blocked" outcome.
    No proposal is created and nothing is active."""
    bundle = _make_valid_bundle(
        dirpath=tmp_path, text="# S\n```bash\ngit push origin main\n```\n")
    info = ss.import_bundle(bundle)
    assert info["status"] == "blocked"
    assert "git push" in info["reason"].lower()
    assert sp.get_proposal("shared-skill") is None
    assert "shared-skill" not in [s.name for s in discover_skills()]


def test_import_clean_text_with_an_absolute_path_is_stored_blocked(
        home, tmp_path):
    bundle = _make_valid_bundle(dirpath=tmp_path,
                                text="# S\nuses /home/alice/refs\n")
    info = ss.import_bundle(bundle)
    assert info["status"] == "blocked"
    assert "/home/alice/refs" in info["reason"] and "line 2" in info["reason"]


def test_import_never_activates_without_accept(home, tmp_path):
    """A verified import stays inactive; only a human attaching REAL evidence
    via skill_proposals.propose(source="shared") makes a proposal, and that
    proposal is pending (never active) until accept."""
    bundle = _make_valid_bundle(dirpath=tmp_path)
    verified = ss.import_bundle(bundle)
    assert verified["status"] == "verified"
    assert "shared-skill" not in [s.name for s in discover_skills()]
    # the human path: attach real evidence and propose with source "shared"
    prop = sp.propose(
        verified["name"], verified["text"],
        evidence=[sp.Evidence(kind="test", ref="tests/test_l6_share.py::x",
                              detail="importer's own green run")],
        source="shared", base_version=verified["version"])
    assert prop.source == "shared" and prop.status == "pending"
    assert "shared-skill" not in [s.name for s in discover_skills()]


def test_import_detects_a_tampered_hash(home, tmp_path):
    bundle = _make_valid_bundle(dirpath=tmp_path)
    bundle.write_text(bundle.read_text(encoding="utf-8").replace(
        "A skill", "A setup"), encoding="utf-8")
    with pytest.raises(Exception):
        ss.import_bundle(bundle)
    assert sp.get_proposal("shared-skill") is None


def test_import_refuses_a_symlink(home, tmp_path):
    real = _make_valid_bundle(dirpath=tmp_path, name="real-target")
    link = tmp_path / "linked.skill.json"
    try:
        link.symlink_to(real)
    except OSError:
        pytest.skip("symlinks unavailable on this filesystem")
    with pytest.raises(Exception):
        ss.import_bundle(link)


def test_import_refuses_an_oversize_bundle(home, tmp_path, monkeypatch):
    bundle = _make_valid_bundle(dirpath=tmp_path, text="x" * 4000)
    # an explicit tiny cap proves the size guard is honoured
    monkeypatch.setattr(ss, "_MAX_BUNDLE_BYTES", 512)
    with pytest.raises(Exception):
        ss.import_bundle(bundle)
