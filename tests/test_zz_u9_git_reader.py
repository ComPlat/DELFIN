"""U9 (operator build): the helper's git reader takes refs, never options."""

from delfin.agent import pr_helper


def test_a_ref_that_looks_like_an_option_is_not_handed_to_git(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    target = tmp_path / "written.txt"
    assert pr_helper._default_git_log(f"--output={target}", "HEAD") == []
    assert pr_helper._default_git_log("main", f"--output={target}") == []
    assert not target.exists()


def test_the_reader_runs_in_the_workspace_it_is_given(tmp_path, monkeypatch):
    import subprocess
    from delfin.agent import safe_git
    seen = {}

    def fake_run(args, cwd, **kw):
        seen["argv"], seen["cwd"] = list(args), str(cwd)
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    monkeypatch.setattr(safe_git, "run", fake_run)
    pr_helper._default_git_log("origin/main", "HEAD", cwd=str(tmp_path))
    assert seen["cwd"] == str(tmp_path)
    assert "--end-of-options" in seen["argv"]
    assert seen["argv"].index("--end-of-options") < seen["argv"].index("origin/main..HEAD")
