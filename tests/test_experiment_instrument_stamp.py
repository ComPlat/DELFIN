"""Package G, Phase 3: the instrument stamp and cross-stamp compare refusal.

Every measurement carries a content stamp -- a hash of the code/judge/tool
files and the relevant environment, taken once.  Two results measured under
different instruments must NOT be compared: that is how an agent silently
moves the goalposts between the arms of an experiment.  ``assert_same_stamp``
is the refusal that makes a cross-stamp compare a hard error.

The stamp covers EXACTLY what the caller declares: the listed files'
contents, the listed environment variables' values, and the judge version.
An environment variable that was never declared is out of scope BY
CONSTRUCTION -- it cannot change the stamp.  A declared environment value or
a file's content that changes absolutely does.
"""

import pytest

from delfin.agent.experiment import (
    ExperimentError,
    assert_same_stamp,
    instrument_stamp,
)

# Unique env keys so the tests never collide with anything real.
ENV1 = "DELFIN_EXP_P3_VAR_A"
ENV2 = "DELFIN_EXP_P3_VAR_B"


def _write(path, text: str):
    path.write_text(text, encoding="utf-8")
    return str(path)


def _stamp(tmp_path, monkeypatch, *, env=("x",), text="switch default off", judge="j1"):
    f1 = _write(tmp_path / "code_a.py", text)
    f2 = _write(tmp_path / "judge_a.py", "def judge(): return 1")
    monkeypatch.setenv(ENV1, env[0])
    if len(env) > 1:
        monkeypatch.setenv(ENV2, env[1])
    return instrument_stamp(files=[f1, f2], env_keys=[ENV1] + ([ENV2] if len(env) > 1 else []), judge=judge)


# ── A stamp is taken once and is immutable ───────────────────────────────


def test_stamp_is_immutable_after_taken(tmp_path, monkeypatch):
    st = _stamp(tmp_path, monkeypatch)
    with pytest.raises(Exception):
        st.content_hash = "different"  # frozen: assignment must fail


def test_same_inputs_same_stamp(tmp_path, monkeypatch):
    a = _stamp(tmp_path, monkeypatch)
    b = _stamp(tmp_path, monkeypatch)
    # Same declared instrument -> same instrument ID (taken_at is when the
    # stamp was RECORDED, not part of instrument identity; asserting the raw
    # dataclass would compare wall-clock values, which is machine-timing).
    assert a.stamp_id == b.stamp_id
    assert_same_stamp(a, b)


# ── A file-content change changes the stamp and the compare refuses ──────


def test_file_content_change_refuses_compare(tmp_path, monkeypatch):
    a = _stamp(tmp_path, monkeypatch, text="switch default off")
    b = _stamp(tmp_path, monkeypatch, text="switch default ON")  # same paths, new content
    assert a.stamp_id != b.stamp_id
    with pytest.raises(ExperimentError):
        assert_same_stamp(a, b)


# ── A declared environment change changes the stamp and refuses ──────────


def test_declared_env_change_refuses_compare(tmp_path, monkeypatch):
    a = _stamp(tmp_path, monkeypatch, env=("off",))
    b = _stamp(tmp_path, monkeypatch, env=("on",))
    assert a.stamp_id != b.stamp_id
    with pytest.raises(ExperimentError):
        assert_same_stamp(a, b)


# ── An UNDECLARED env change does not change the stamp ───────────────────


def test_undeclared_env_change_is_out_of_scope(tmp_path, monkeypatch):
    monkeypatch.setenv(ENV2, "some-unrelated-value")
    a = _stamp(tmp_path, monkeypatch, env=("x",))
    monkeypatch.setenv(ENV2, "now-changed")
    b = _stamp(tmp_path, monkeypatch, env=("x",))
    # ENV2 was never declared; the stamp must ignore it by construction.
    assert a.stamp_id == b.stamp_id
    assert_same_stamp(a, b)


def test_declaring_a_missing_env_var_refuses(tmp_path, monkeypatch):
    monkeypatch.delenv(ENV1, raising=False)
    f1 = _write(tmp_path / "code_a.py", "switch default off")
    with pytest.raises(ExperimentError):
        instrument_stamp(files=[f1], env_keys=[ENV1])


# ── Judge version is part of the instrument ──────────────────────────────


def test_different_judge_version_refuses_compare(tmp_path, monkeypatch):
    a = _stamp(tmp_path, monkeypatch, judge="judge-v1")
    b = _stamp(tmp_path, monkeypatch, judge="judge-v2")
    assert a.stamp_id != b.stamp_id
    with pytest.raises(ExperimentError):
        assert_same_stamp(a, b)


# ── Refusals on an unusable instrument ───────────────────────────────────


def test_stamp_refuses_no_files(tmp_path, monkeypatch):
    with pytest.raises(ExperimentError):
        instrument_stamp(files=[], env_keys=[ENV1])


def test_stamp_refuses_missing_file(tmp_path, monkeypatch):
    monkeypatch.setenv(ENV1, "x")
    missing = str(tmp_path / "does_not_exist.py")
    with pytest.raises(ExperimentError):
        instrument_stamp(files=[missing], env_keys=[ENV1])


def test_stamp_refuses_no_env_keys(tmp_path, monkeypatch):
    f1 = _write(tmp_path / "code_a.py", "x")
    with pytest.raises(ExperimentError):
        instrument_stamp(files=[f1], env_keys=[])


def test_stamp_message_is_clear(tmp_path, monkeypatch):
    a = _stamp(tmp_path, monkeypatch, judge="judge-v1")
    b = _stamp(tmp_path, monkeypatch, judge="judge-v2")
    try:
        assert_same_stamp(a, b)
        assert False, "assert_same_stamp must raise"
    except ExperimentError as e:
        msg = str(e)
        assert "judge" in msg or "stamp" in msg or "instrument" in msg


@pytest.mark.parametrize("key", ["KIT_TOOLBOX_API_KEY", "GITHUB_TOKEN", "db_password", "MY_SECRET", "AWS_CREDENTIALS"])
def test_a_secret_never_enters_the_stamp(tmp_path, monkeypatch, key):
    """A hash of a short key or password can be guessed offline, and a stamp
    travels with every result: secret-named variables are refused outright."""
    f = tmp_path / "judge.py"
    f.write_text("x = 1\n")
    monkeypatch.setenv(key, "hunter2")
    with pytest.raises(ExperimentError, match="never part of an instrument"):
        instrument_stamp(files=[str(f)], env_keys=[key])
