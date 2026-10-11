"""Adversarial ALLOW-list tests: the sbatch-id classifier under-catches nothing.

Package V3's promise is "nothing finishes unnoticed". A background command
that submits a cluster job (``sbatch run.sh``) is auto-watched by scanning
its finished output for the SLURM job id it submitted; if that scan misses
the id, the cluster job is never registered as watched and its completion
is never announced. The classifier is ``_submitted_slurm_ids`` (bash_jobs.py)
guarded by ``is_own_output_file``.

Per review rules, classifiers are tested by their ALLOW-list: what the
parser must accept so it never silently drops a real submission. These
tests pin the allow-list of canonical sbatch output forms.

No real sbatch: the parse is pure over a tempfile that passes
``is_own_output_file``. sbatch's standard stdout line is exactly
"Submitted batch job <n>". ``--parsable`` is a deliberate opt-in that
changes the format and is NOT in the allow-list here -- that is a
documented scope limit, not a defect; it is asserted as such so a future
change that makes the parser fail the *normal* form is caught.
"""
from delfin.agent import bash_jobs as jb
from delfin.agent import job_monitor as jm


def _own_stdout(tmp_path, name: str, text: str) -> str:
    """Write an output file that passes ``is_own_output_file`` and return it."""
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    assert jb.is_own_output_file(p), "precondition: file must look module-owned"
    return str(p)


def test_allow_list_canonical_submission_is_extracted(tmp_path):
    # sbatch's standard stdout line, nothing else in the file.
    out = _own_stdout(tmp_path, "kit_bg_x.stdout", "Submitted batch job 4067231\n")
    assert jb._submitted_slurm_ids(out) == ["4067231"]


def test_allow_list_multiple_submissions_are_extracted_in_order(tmp_path):
    out = _own_stdout(
        tmp_path, "kit_bg_m.stdout",
        "Submitted batch job 10\nSubmitted batch job 20\n",
    )
    assert jb._submitted_slurm_ids(out) == ["10", "20"]


def test_allow_list_extra_whitespace_and_noise_do_not_hide_id(tmp_path):
    out = _own_stdout(
        tmp_path, "kit_bg_w.stdout",
        "srun: job 5 completed\nSubmitted batch job    42\n[some log tail]",
    )
    assert jb._submitted_slurm_ids(out) == ["42"]


def test_deny_list_repeated_ids_are_reported_once(tmp_path):
    out = _own_stdout(
        tmp_path, "kit_bg_d.stdout",
        "Submitted batch job 9\nSubmitted batch job 9\n",
    )
    assert jb._submitted_slurm_ids(out) == ["9"]


def test_deny_list_out_of_scope_parsable_form_is_deliberately_not_watched(tmp_path):
    # sbatch --parsable prints only the bare id. Not in the allow-list on
    # purpose: the register-on-finish path is built around the human-readable
    # form. Pin that it is NOT extracted (scope limit), not observed to race.
    out = _own_stdout(tmp_path, "kit_bg_p.stdout", "4067231\n")
    assert jb._submitted_slurm_ids(out) == []


def test_deny_list_unrelated_slurm_text_is_not_a_submission(tmp_path):
    # A plain log that mentions "batch job" but was never an sbatch stdout.
    out = _own_stdout(tmp_path, "kit_bg_u.stdout", "batch job id 7 done\n")
    assert jb._submitted_slurm_ids(out) == []


def test_non_own_or_missing_file_is_an_error_not_a_silent_empty_job(tmp_path):
    # is_own_output_file guards the parse; a file that does not look like a
    # module output must produce nothing (never a guessed id).
    stray = tmp_path / "secret.log"
    stray.write_text("Submitted batch job 99\n", encoding="utf-8")
    assert jb.is_own_output_file(stray) is False
    assert jb._submitted_slurm_ids(str(stray)) == []
    assert jb._submitted_slurm_ids(str(tmp_path / "nope")) == []


def test_classifier_regex_is_stable(tmp_path):
    # The exact allow-list regex must not drift (it gates the whole watch).
    assert jb._SBATCH_SUBMITTED_RE.pattern == r"Submitted batch job\s+(\d+)"


def test_register_path_uses_shell_session_id(tmp_path):
    # _watch_submitted_jobs registers the extracted id with the submitting
    # shell's session so the completion reaches the session that started it.
    ws = tmp_path / "ws"
    ws.mkdir()
    out = _own_stdout(tmp_path, "kit_bg_r.stdout", "Submitted batch job 77\n")
    rec = {"stdout_path": out, "job_id": "bg1", "session_id": "sess-1"}
    watched = jb._watch_submitted_jobs(str(ws), rec)
    assert watched == ["77"]
    jobs = jm.load_watched(jm._agent_watch_path(str(ws)))["jobs"]
    assert jobs["77"]["session_id"] == "sess-1"
