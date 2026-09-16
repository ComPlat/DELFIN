"""Two writers saving memories at once both keep their entry.

Agreed with the user on 2026-09-16: parallel sessions share one memory
store (the home directory is shared by every session and every login
node), and saving was a read-merge-write cycle with no lock -- the second
index rewrite could drop the first one's line. And the index is read again
on every turn, so a fact one session saves reaches another mid-conversation.
"""
import inspect
import threading
import time

from delfin.agent import memory_store as MS


def test_the_save_holds_the_cross_process_lock():
    src = inspect.getsource(MS.save_typed_memory)
    assert "cross_process_lock(_MEMORY_WRITE_LOCK)" in src
    assert "_save_typed_memory_unlocked(text, **kwargs)" in src


def test_two_threads_saving_at_once_both_land(tmp_path, monkeypatch):
    repo = tmp_path / "repo"; repo.mkdir()
    slow = MS._update_memory_index

    def slow_index(*a, **k):
        time.sleep(0.05)          # widen the read-modify-write window
        return slow(*a, **k)
    monkeypatch.setattr(MS, "_update_memory_index", slow_index)
    errors = []

    facts = [
        "ORCA scratch files go to the node-local disk under /scratch",
        "The NMR training set holds eight hundred samples of thirteen shifts each",
        "Jerome works on carbon dioxide coordination at metal centres",
        "Continuous integration runs only on a push to the main branch",
        "The GLM endpoint answers slowly before noon on busy days",
        "Uploads from the dashboard land in agent_workspace/uploads",
    ]

    def save(i):
        try:
            MS.save_typed_memory(facts[i], repo_root=repo, memory_type="project",
                                 title=f"fact-{i}")
        except Exception as exc:
            errors.append(exc)
    threads = [threading.Thread(target=save, args=(i,)) for i in range(6)]
    for t in threads: t.start()
    for t in threads: t.join()
    assert not errors
    index = (MS._memory_dir_for_scope(repo, "project") / "MEMORY.md").read_text()
    for i in range(6):
        assert f"fact-{i}" in index, f"entry {i} lost to a concurrent save"


def test_the_prompt_reads_the_index_on_every_turn():
    """The memory index travels in the system prompt, and the prompt is
    built per turn -- so a fact saved by another session is seen on the
    next turn, not the next start."""
    from delfin.agent.engine import AgentEngine
    src = inspect.getsource(AgentEngine.stream_response)
    assert "self._build_current_system_prompt(" in src
