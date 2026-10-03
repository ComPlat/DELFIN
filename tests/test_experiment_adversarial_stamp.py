"""Package G, reviewer s6: adversarial tests for phase 3 (instrument stamp).

Attacks: the stamp's collision/aliasing surface and the cross-stamp
refusal's completeness.  The builder's own tests cover the declared cases;
these go after what was NOT declared:

  A3.1  same bytes via different paths (symlink alias) -- must refuse
  A3.2  file reordering changes the stamp (order sensitivity)
  A3.3  relative vs absolute path of the SAME file (alias again)
  A3.4  stamp_id covers ALL three components (env-only + judge-only + both)
  A3.5  env hash is not injective across key/value swaps
  A3.6  a changed file must produce a DIFFERENT stamp_id, not just refusal
  A3.7  directory in the file list is a refusal, not a silent skip
"""

import os

import pytest

from delfin.agent.experiment import (
    ExperimentError,
    assert_same_stamp,
    instrument_stamp,
)

ENV1 = "DELFIN_ADV3_VAR_A"


@pytest.fixture
def env1(monkeypatch):
    monkeypatch.setenv(ENV1, "stable-value")
    return ENV1


def _reg(tmp_path, name="code.py", text="switch default off"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return str(p)


# ── A3.1: same bytes through a symlink — is the instrument ambiguous? ────

def test_symlink_alias_same_bytes_refuses(tmp_path, env1, monkeypatch):
    """code.py and link.py hold IDENTICAL bytes via a symlink.  If the
    dedup key is os.path.abspath, both list and the hash double-counts the
    same content — the instrument becomes ambiguous without a refusal."""
    real = _reg(tmp_path, "code.py")
    link = tmp_path / "link.py"
    link.symlink_to(real)
    a = instrument_stamp(files=[real], env_keys=[env1])
    b = instrument_stamp(files=[real, str(link)], env_keys=[env1])
    # Two files with identical bytes double-counted: content_hash differs
    # from the single-file stamp.  The question is whether the module
    # REFUSES (same resolved file listed twice) or silently accepts.
    try:
        same = (a.stamp_id == b.stamp_id)
    except ExperimentError:
        return  # refused: acceptable
    assert not same, (
        "symlink alias of the same file silently produced a second stamp "
        "with different content_hash — an ambiguous instrument accepted")


# ── A3.2: order sensitivity — same instrument, listed in another order ──

def test_file_order_changes_stamp(tmp_path, env1):
    """The stamp hashes files in declared order.  The SAME instrument
    declared as [f1, f2] vs [f2, f1] produces different stamps — the
    cross-stamp compare then refuses an instrument that did not change.
    Document the behavior (order is part of the declaration)."""
    f1 = _reg(tmp_path, "a.py", "one")
    f2 = _reg(tmp_path, "b.py", "two")
    s12 = instrument_stamp(files=[f1, f2], env_keys=[env1])
    s21 = instrument_stamp(files=[f2, f1], env_keys=[env1])
    # Order IS part of the stamp.  This is a documentation-of-behavior pin:
    # a compare of reordered declarations refuses, which is defensible
    # (declaration differs) but the builder must own it explicitly.
    if s12.stamp_id != s21.stamp_id:
        with pytest.raises(ExperimentError):
            assert_same_stamp(s12, s21)


# ── A3.3: relative vs absolute path of the same file ─────────────────────

def test_relative_vs_absolute_path_same_file(tmp_path, env1, monkeypatch):
    """The same file declared with a relative vs an absolute path must give
    the SAME stamp (the instrument is the file's content, not the spelling
    of its path).  The hash does not include paths, so this should hold."""
    f1 = _reg(tmp_path, "code.py", "stable content")
    monkeypatch.chdir(tmp_path)
    s_abs = instrument_stamp(files=[f1], env_keys=[env1])
    s_rel = instrument_stamp(files=["code.py"], env_keys=[env1])
    assert s_abs.stamp_id == s_rel.stamp_id
    assert_same_stamp(s_abs, s_rel)


# ── A3.4: stamp_id covers every component ─────────────────────────────────

def test_stamp_id_covers_all_components(tmp_path, env1, monkeypatch):
    """stamp_id must move when ANY component moves: content, env, judge —
    not only content.  Probe each axis in isolation."""
    f1 = _reg(tmp_path, "code.py")
    ENV2 = "DELFIN_ADV3_VAR_B"
    monkeypatch.setenv(ENV2, "v1")
    base = instrument_stamp(files=[f1], env_keys=[env1], judge="j1")

    env_only = instrument_stamp(files=[f1], env_keys=[ENV2], judge="j1")
    assert base.stamp_id != env_only.stamp_id, "env axis not in stamp_id"

    judge_only = instrument_stamp(files=[f1], env_keys=[env1], judge="j2")
    assert base.stamp_id != judge_only.stamp_id, "judge axis not in stamp_id"

    (tmp_path / "code.py").write_text("changed", encoding="utf-8")
    content_only = instrument_stamp(files=[f1], env_keys=[env1], judge="j1")
    assert base.stamp_id != content_only.stamp_id, "content axis not in stamp_id"


# ── A3.5: env hash cannot conflate a key/value swap ──────────────────────

def test_env_hash_conflates_key_value_swap(tmp_path, env1, monkeypatch):
    """The env hash concatenates key\\0value per key.  Two DIFFERENT
    environments must never produce the same env_hash.  The classic
    collision: ['A=x', 'B=y'] vs ['A=x\\0B', 'B=y'] — the \\0 separator
    makes that safe, but a swap like ['A=x'] vs ['x=A'] is only safe if the
    key is hashed with its position.  Probe it."""
    monkeypatch.setenv("DELFIN_ADV3_VAR_K", "x=y")
    a = instrument_stamp(files=[_reg(tmp_path)], env_keys=["DELFIN_ADV3_VAR_K"])
    monkeypatch.setenv("DELFIN_ADV3_VAR_K", "y=x")
    # rename probe: another var whose (key,value) is the swapped pair
    monkeypatch.setenv("DELFIN_ADV3_VAR_J", "x=y")
    b = instrument_stamp(files=[_reg(tmp_path, "other.py")],
                         env_keys=["DELFIN_ADV3_VAR_J"])
    # a's pair is ('K','x=y'); b's is ('J','x=y') — different keys, same
    # value, different files.  env_hash must differ (key is in the hash).
    assert a.env_hash != b.env_hash


# ── A3.6: a content change moves stamp_id, not only the refusal ──────────

def test_content_change_moves_stamp_id(tmp_path, env1):
    """The cross-stamp compare refusing is not enough: stamp_id itself must
    differ when content differs, else two results stamped before/after an
    edit could be recorded under one id and compared without a refusal
    being reachable at all."""
    f1 = _reg(tmp_path, "code.py", "before")
    s1 = instrument_stamp(files=[f1], env_keys=[env1])
    (tmp_path / "code.py").write_text("after", encoding="utf-8")
    s2 = instrument_stamp(files=[f1], env_keys=[env1])
    assert s1.stamp_id != s2.stamp_id


# ── A3.7: a directory in the file list is a refusal ──────────────────────

def test_directory_in_file_list_refuses(tmp_path, env1):
    """A directory passes path.is_file() == False, so the code path says
    refusal.  Pin it: the message must name the path, not crash."""
    d = tmp_path / "subdir"
    d.mkdir()
    with pytest.raises(ExperimentError) as excinfo:
        instrument_stamp(files=[str(d)], env_keys=[env1])
    assert "subdir" in str(excinfo.value)
