"""The system-prompt head is byte-stable within a session.

The workspace decisions (which tool groups, which DELFIN blocks, which
stack advice) are taken once. A later turn whose text names chemistry
and documents, after the agent wrote an ``.xyz``, an ``.xlsx`` and a
``Cargo.toml`` into the workspace, must leave the head identical, or the
prefix cache is lost for the whole history behind it. Measured on GLM:
7-12 s warm against 199-266 s cold.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent import engine as E

# A model whose lazy modules follow the task text, so the text of turn
# two could move them; GLM keeps every module and would hide the effect.
MODEL = "kit.deepseek-v4-flash"


def _ts(base: Path) -> Path:
    base.mkdir(parents=True)
    (base / "package.json").write_text('{"name": "app"}')
    (base / "tsconfig.json").write_text("{}")
    return base


def _sections(eng: E.AgentEngine, task: str) -> dict[str, str]:
    eng.messages = [{"role": "user", "content": task}]
    captured: dict[str, str] = {}
    orig = eng.loader.compose_sections

    def _spy(**kw):
        secs = orig(**kw)
        captured.update({s.name: s.content for s in secs})
        return secs

    eng.loader.compose_sections = _spy
    try:
        eng._build_current_system_prompt(task_text=task)
    finally:
        eng.loader.compose_sections = orig
    return captured


def test_the_head_does_not_move_within_a_session(tmp_path):
    ws = _ts(tmp_path / "app")
    eng = E.AgentEngine(repo_dir=ws, backend="api", provider="kit",
                        api_key="stub-key", model=MODEL, mode="solo",
                        permission_mode="default")
    eng._load_mode("solo")
    first = _sections(eng, "fix the failing unit test in src/")
    (ws / "mol.xyz").write_text("1\n\nH 0 0 0\n")
    (ws / "sheet.xlsx").write_bytes(b"PK")
    (ws / "Cargo.toml").write_text('[package]\nname = "x"\n')
    second = _sections(eng, "the energy of the molecule complex looks "
                            "wrong; also check the column in the "
                            "spreadsheet")
    for name in ("role_prompt", "project_context"):
        assert first.get(name) == second.get(name), name


def test_the_decision_waits_for_the_first_turn(tmp_path):
    """Taken at the first turn, not at construction: the benchmark seeds
    a chemistry fixture between building the engine and sending the
    task (benchmark_runner: engine_factory, then run_setup)."""
    ws = tmp_path / "chem_ws"
    ws.mkdir()
    eng = E.AgentEngine(repo_dir=ws, backend="api", provider="kit",
                        api_key="stub-key", model=MODEL, mode="solo",
                        permission_mode="default")
    eng._load_mode("solo")
    (ws / "chem" / "benzene").mkdir(parents=True)
    (ws / "chem" / "benzene" / "start.xyz").write_text("1\n\nC 0 0 0\n")
    sections = _sections(eng, "optimise the geometry")
    assert eng.kit_permissions.toolsets.offered("chemistry")
    assert "Not loaded here: chemistry" not in sections["project_context"]
