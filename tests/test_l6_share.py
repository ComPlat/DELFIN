"""skill_share.export_bundle / import_bundle through the module API.

* export writes one bundle for an accepted skill and refuses any other;
* the scrub blocks, naming field and line, an absolute or ``~/`` path, a
  user@host and a job id next to a scheduler word, in text, author and
  version, while method strings (``B3LYP/def2-TZVP``, ``kcal/mol``) pass;
* a symlink planted at the temp name is not followed;
* import refuses a symlink, a changed hash and an oversize bundle, and a
  clean bundle is stored as a pending proposal that accept() refuses.

HOME is redirected so the skills tree and proposal store are in tmp.
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


def test_export_chemistry_method_strings_and_plain_numbers_pass(home, tmp_path):
    """A path is a token that BEGINS with / or ~/, never a "/"
    inside a word, and a job id is only next to a scheduler word. So a
    chemistry skill mentioning B3LYP/def2-TZVP, kcal/mol, and/or and a bare
    6-digit number must EXPORT, while a real /home path, a ~/ path, and a
    scheduler job id must BLOCK."""
    _install_accepted(home, text=(
        "# DFT\n"
        "PBE0 with the def2-TZVP basis (B3LYP/def2-SVP reference)\n"
        "energy 12.4 kcal/mol, and/or correction; mass 987654 u\n"))
    out = ss.export_bundle("shared-skill", str(tmp_path))
    assert out.startswith("# DFT")
    assert _bundle_path(tmp_path).exists()

    _install_accepted(home, text="# X\nworkdir is /home/alice/x\n",
                      name="blocked-path")
    with pytest.raises(Exception) as exc:
        ss.export_bundle("blocked-path", str(tmp_path))
    assert "line 2" in str(exc.value)
    assert not _bundle_path(tmp_path, "blocked-path").exists()

    _install_accepted(home, text="# Y\nscratch under ~/calc\n",
                      name="blocked-tilde")
    with pytest.raises(Exception) as exc:
        ss.export_bundle("blocked-tilde", str(tmp_path))
    assert "~/calc" in str(exc.value)
    assert not _bundle_path(tmp_path, "blocked-tilde").exists()


def test_export_blocks_a_job_id_only_next_to_a_scheduler_word(home, tmp_path):
    _install_accepted(home, text="# S\nslurm job 9876543 finished\n")
    with pytest.raises(Exception) as exc:
        ss.export_bundle("shared-skill", str(tmp_path))
    assert "9876543" in str(exc.value)
    assert not _bundle_path(tmp_path).exists()
    # a bare 6-digit number (not a scheduler job id) still exports
    _install_accepted(home, text="# S\nvalue 123456 normal\n")
    assert ss.export_bundle("shared-skill", str(tmp_path)) == "# S\nvalue 123456 normal\n"


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
    with pytest.raises(Exception):
        ss.export_bundle("shared-skill", str(tmp_path))
    assert calls, "export must consult the reused skill_safety scan"
    assert not _bundle_path(tmp_path).exists()


def test_export_refuses_a_planted_symlink_at_the_temp_name(home, tmp_path):
    """The export temp is created with O_EXCL|O_NOFOLLOW (0600),
    so a symlink planted at the temp name is NOT followed -- the write would
    otherwise land outside (os.replace'd into an attacker's path). A planted
    symlink at the temp name must be refused, never followed."""
    _install_accepted(home, text="# S\nA clean body\n", name="sym-skill")
    tmp = tmp_path / "sym-skill.skill.json.tmp"
    tmp.symlink_to(tmp_path / "outside.txt")
    with pytest.raises(Exception):
        ss.export_bundle("sym-skill", str(tmp_path))
    # refusal: nothing was written THROUGH the symlink (outside.txt stays
    # absent) and no bundle landed; the planted symlink itself is left
    # untouched (not followed, not deleted).
    assert not tmp_path.joinpath("outside.txt").exists()
    assert not tmp_path.joinpath("sym-skill.skill.json").exists()
    assert tmp_path.joinpath("sym-skill.skill.json.tmp").is_symlink()


def test_import_refuses_a_symlink_via_the_open_descriptor(home, tmp_path):
    """The bundle read decides on the OPEN descriptor
    (O_RDONLY|O_NOFOLLOW|O_NONBLOCK + fstat), not on a check-then-use of the
    path -- a symlink is refused even if the path stats as a regular file."""
    target = tmp_path / "real_bundle.txt"
    target.write_text("not a bundle", encoding="utf-8")
    link = tmp_path / "bundle.skill.json"
    link.symlink_to(target)
    with pytest.raises(Exception) as exc:
        ss.import_bundle(str(link))
    assert "symlink" in str(exc.value).lower() or "not a regular" in str(exc.value)


def test_export_blocked_on_a_userhost_author_label(home, tmp_path):
    """The author label is also machine identity and must
    be blocked, naming the author field -- not just the skill text."""
    _install_accepted(home, text="# A skill\nbody")
    with pytest.raises(Exception) as exc:
        ss.export_bundle("shared-skill", str(tmp_path), author="alice@srv")
    msg = str(exc.value)
    assert "author" in msg and "alice@srv" in msg
    assert not _bundle_path(tmp_path).exists()


def test_export_blocked_on_an_absolute_path_author(home, tmp_path):
    _install_accepted(home, text="# A skill\nbody")
    with pytest.raises(Exception) as exc:
        ss.export_bundle("shared-skill", str(tmp_path),
                         author="/home/x/me")
    assert "author" in str(exc.value)
    assert not _bundle_path(tmp_path).exists()


def test_export_blocked_on_a_machine_ref_in_version(home, tmp_path):
    """A base_version that names a machine path is blocked on export -- the
    version field leaves too and is scrubbed."""
    _install_accepted(home, text="# A skill\nbody", version="/home/x/base")
    with pytest.raises(Exception) as exc:
        ss.export_bundle("shared-skill", str(tmp_path))
    assert "version" in str(exc.value)
    assert not _bundle_path(tmp_path).exists()


# --------------------------------------------------------------------------
# import_bundle -- bounds, no proposal of its own
# --------------------------------------------------------------------------

def test_import_clean_bundle_is_stored_pending_never_active(home, tmp_path):
    """A clean, hash-verified bundle is stored as a pending proposal of
    source "shared"; it is not discovered, and accept() refuses it."""
    bundle = _make_valid_bundle(dirpath=tmp_path)
    info = ss.import_bundle(bundle)
    assert info["status"] == "pending"
    assert info["name"] == "shared-skill"
    assert info["text"] == "# A skill\nbody"
    prop = sp.get_proposal("shared-skill")
    assert prop is not None and prop.status == "pending"
    assert prop.source == "shared"
    assert "shared-skill" not in [s.name for s in discover_skills()]
    with pytest.raises(ValueError, match="no verifying evidence"):
        sp.accept("shared-skill", by="tester")


def test_import_unsafe_text_is_blocked_naming_the_line(home, tmp_path):
    """Unsafe text (git push) is recorded by import_bundle itself as status
    "blocked", naming the line.
    No proposal is created and nothing is active."""
    bundle = _make_valid_bundle(
        dirpath=tmp_path, text="# S\n```bash\ngit push origin main\n```\n")
    info = ss.import_bundle(bundle)
    assert info["status"] == "blocked"
    assert "git push" in info["reason"].lower()
    assert sp.get_proposal("shared-skill") is None
    assert "shared-skill" not in [s.name for s in discover_skills()]


def test_import_clean_text_with_an_absolute_path_is_blocked(
        home, tmp_path):
    bundle = _make_valid_bundle(dirpath=tmp_path,
                                text="# S\nuses /home/alice/refs\n")
    info = ss.import_bundle(bundle)
    assert info["status"] == "blocked"
    assert "/home/alice/refs" in info["reason"] and "line 2" in info["reason"]


def test_import_re_scrubs_author_and_version_before_proposal(home, tmp_path):
    """Import side: the author and version are re-scrubbed
    before a bundle is proposed -- a machine ref in either is "blocked",
    not silently accepted into a proposal."""
    bundle = _make_valid_bundle(dirpath=tmp_path, author="alice@srv")
    info = ss.import_bundle(bundle)
    assert info["status"] == "blocked"
    assert "author" in info["reason"]
    # a clean author still stores a pending proposal
    clean = _make_valid_bundle(dirpath=tmp_path, name="shared-skill-clean",
                               author="kb", version="v1")
    assert ss.import_bundle(clean)["status"] == "pending"


def test_import_detects_a_tampered_hash(home, tmp_path):
    bundle = _make_valid_bundle(dirpath=tmp_path)
    bundle.write_text(bundle.read_text(encoding="utf-8").replace(
        "A skill", "A setup"), encoding="utf-8")
    with pytest.raises(Exception):
        ss.import_bundle(bundle)
    assert sp.get_proposal("shared-skill") is None


def test_import_refuses_an_oversize_bundle(home, tmp_path, monkeypatch):
    bundle = _make_valid_bundle(dirpath=tmp_path, text="x" * 4000)
    # an explicit tiny cap proves the size guard is honoured
    monkeypatch.setattr(ss, "_MAX_BUNDLE_BYTES", 512)
    with pytest.raises(Exception):
        ss.import_bundle(bundle)
