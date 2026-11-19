"""JC2 reviewer adversarial tests, part 3: sharp edges of the params fix.

Covers the cases I promised the builder plus one finding from reading
the fix itself:
- a Path-valued USER param must survive the write path (L7 json-safety
  must not be regressed by the spec.kwargs switch);
- a SMILES-string "geometry" (first step of a chain) hashes its CONTENT,
  never a failed file read silently;
- a pending LEGACY step is re-hashed at its resume submission (stale
  hash over polluted params does not survive);
- FINDING: a user who legitimately passes ``moread``/``hess_file``
  explicitly to Pipeline.add gets their param dropped on load —
  _drop_internal_params cannot distinguish rule-injected from explicit.
"""

import json
from pathlib import Path


from delfin.tools import chain_state
from delfin.tools.pipeline import Pipeline, StepResult, StepStatus


def _fake_run_step(step_name, geometry=None, work_dir=None, **kwargs):
    out = Path(work_dir) / f"{step_name}.out"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("fake run output\n")
    return StepResult(step_name=step_name, status=StepStatus.SUCCESS,
                      artifacts={"gbw": out.with_suffix(".gbw")})


def test_path_valued_user_param_survives_the_write_path(
    monkeypatch, tmp_path
):
    """L7: a Path-valued USER kwarg must be json_safe'd, not crash the
    save (the JC1 bug this code fixed once already)."""
    monkeypatch.setattr("delfin.tools._runner.run_step", _fake_run_step)
    wd = tmp_path / "c"
    pipe = Pipeline(name="p")
    pipe.add("xtb_opt", charge=0, some_file=Path(tmp_path) / "ref.xyz")
    pipe.run(work_dir=wd, resume_chain=True)
    state = chain_state.ChainState.load(wd / "chain_state.json")
    step = state.step("xtb_opt")
    assert isinstance(step.params.get("some_file"), str), (
        "a Path-valued user param must be persisted as a string"
    )


def test_smiles_first_step_hashes_string_content_not_failed_file_read(
    monkeypatch, tmp_path
):
    """A chain starting at smiles_to_xyz has a STRING 'geometry'. The
    hash must differ between two different SMILES (content matters) and
    must NOT log a failed-file-read warning per submission (a string is
    not a path — content is the string itself)."""
    monkeypatch.setattr("delfin.tools._runner.run_step", _fake_run_step)

    def run_with(smiles, wd):
        pipe = Pipeline(name="s")
        pipe.add("smiles_to_xyz", smiles=smiles)
        pipe.run(work_dir=wd, resume_chain=True, geometry=smiles)
        return chain_state.ChainState.load(wd / "chain_state.json")

    state_a = run_with("CCO", tmp_path / "a")
    state_b = run_with("CCC", tmp_path / "b")
    ha = state_a.step("smiles_to_xyz").inputs_hash
    hb = state_b.step("smiles_to_xyz").inputs_hash
    assert ha != hb, (
        "two different SMILES inputs must hash differently — the hash "
        "covers content, not the mere fact that no file existed"
    )


class _NullCtx:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_pending_legacy_step_rehashed_at_resume_submission(
    monkeypatch, tmp_path
):
    """A JC1-era state whose pending step carries a hash computed over
    polluted params: after resume, the step's submission rewrites the
    hash over clean user params (stale hash not trusted)."""
    wd = tmp_path / "legacy"
    wd.mkdir()
    legacy = {
        "name": "legacy-chain",
        "work_dir": str(wd),
        "steps": [
            {"name": "xtb_opt", "tool": "xtb_opt", "state": "pending",
             "params": {"charge": 0, "_prev_artifacts": {"gbw": "/old/x.gbw"}},
             "inputs_hash": "stale-legacy-hash"},
        ],
    }
    (wd / "chain_state.json").write_text(json.dumps(legacy), encoding="utf-8")

    monkeypatch.setattr("delfin.tools._runner.run_step", _fake_run_step)
    pipe = Pipeline(name="legacy-chain")
    pipe.add("xtb_opt", charge=0)
    pipe.run(work_dir=wd, resume_chain=True)
    state = chain_state.ChainState.load(wd / "chain_state.json")
    step = state.step("xtb_opt")
    assert step.inputs_hash != "stale-legacy-hash", (
        "the pending step's stale legacy hash must be replaced at its "
        "resume submission"
    )
    assert "_prev_artifacts" not in step.params
    # And the new hash equals a fresh chain's hash for the same inputs:
    fresh = tmp_path / "fresh"
    pipe2 = Pipeline(name="fresh")
    pipe2.add("xtb_opt", charge=0)
    pipe2.run(work_dir=fresh, resume_chain=True)
    fresh_state = chain_state.ChainState.load(fresh / "chain_state.json")
    assert step.inputs_hash == fresh_state.step("xtb_opt").inputs_hash, (
        "the re-hashed legacy step must equal a fresh chain's hash for "
        "the same user inputs"
    )


def test_explicit_moread_is_dropped_on_load_finding(monkeypatch, tmp_path):
    """FINDING (report to builder, red on the fix): a user who
    legitimately passes moread=... explicitly to Pipeline.add gets it
    silently dropped on load — _drop_internal_params cannot distinguish
    a rule-injected from an explicit kwarg."""
    monkeypatch.setattr("delfin.tools._runner.run_step", _fake_run_step)
    wd = tmp_path / "explicit"
    pipe = Pipeline(name="p")
    pipe.add("orca_sp", charge=0, method="B3LYP",
             moread=str(tmp_path / "my.gbw"))
    pipe.run(work_dir=wd, resume_chain=True)
    state = chain_state.ChainState.load(wd / "chain_state.json")
    step = state.step("orca_sp")
    assert step.params.get("moread") == str(tmp_path / "my.gbw"), (
        "an EXPLICIT user moread must survive the round trip; dropping it "
        "unconditionally by key silently changes the user's calculation"
    )


def test_legacy_state_without_explicit_set_still_loads(
    monkeypatch, tmp_path
):
    """QS ruling on JC2-R1, requirement 2: a legacy JC1 chain_state.json
    carries NO explicit-param set (the field does not exist yet). Load
    must treat the absent set as 'unknown': keep everything except the
    known internal keys (with the ONE log line the JC2 target already
    demands), and still resume. If load refused or dropped more user
    params on such a file, every JC1-written chain would lose data on
    first resume."""
    monkeypatch.setattr("delfin.tools._runner.run_step", _fake_run_step)
    wd = tmp_path / "legacy-unknown"
    wd.mkdir()
    legacy = {
        "name": "jc1-era",
        "work_dir": str(wd),
        "steps": [
            {
                "name": "orca_sp",
                "tool": "orca_sp",
                "state": "pending",
                # A user param whose fate depends on the explicit set:
                # with the set ABSENT (unknown provenance) it must stay.
                "params": {"method": "B3LYP", "user_scalar": 42},
            }
        ],
    }
    (wd / "chain_state.json").write_text(
        json.dumps(legacy), encoding="utf-8"
    )
    pipe = Pipeline(name="jc1-era")
    pipe.add("orca_sp", charge=0)
    pipe.run(work_dir=wd, resume_chain=True)
    state = chain_state.ChainState.load(wd / "chain_state.json")
    step = state.step("orca_sp")
    assert step.params.get("user_scalar") == 42, (
        "a legacy file without the explicit-param set keeps unknown-"
        "provenance user params (drop only the known internal keys)"
    )
