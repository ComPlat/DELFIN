"""Package G, Phase 6: the agent-experiments documentation is present and
grounded in the real API.

The doc is the deliverable that tells an agent how to run an experiment, so
it must exist and be non-trivial, and every public API name it cites must
resolve to a real callable in ``delfin.agent.experiment`` -- a doc that names
a function that doesn't exist is a doc that teaches the agent to fool itself.
"""

import re
from pathlib import Path

import delfin.agent.experiment as exp_mod

_DOC = Path(__file__).resolve().parents[1] / "docs" / "agent_experiments.md"


def test_doc_exists_and_is_non_trivial():
    assert _DOC.is_file(), "docs/agent_experiments.md must exist"
    text = _DOC.read_text(encoding="utf-8")
    assert len(text.split()) > 200, "the doc must be more than a stub"


def test_doc_covers_every_principle():
    text = _DOC.read_text(encoding="utf-8")
    for phrase in [
        "One switch",
        "reach before measuring",
        "Pre-registration",
        "content stamp",
        "Absolute counts",
        "null run",
        "real regression",
        "human approves",
    ]:
        assert phrase.lower() in text.lower(), f"doc must cover: {phrase}"


def test_doc_cites_only_real_api_names():
    text = _DOC.read_text(encoding="utf-8")
    # The doc may cite the experiment API and the one cross-module
    # statistics point it names (compare_runs in delfin.agent.benchmark).
    import delfin.agent.benchmark as bm_mod
    # Only API-STYLE code mentions must resolve to a callable: snake_case
    # function names or CamelCase class names.  Plain English words in
    # backquotes (like a status value) and dotted path-like references
    # (benchmark_fixtures._make_deterministic) are prose, not API.
    code_mentions = set(re.findall(r"`([A-Za-z_][A-Za-z0-9_.]*)`", text))
    for name in sorted(code_mentions):
        if name.startswith("delfin/agent") or " " in name:
            continue
        if "." in name:                # dotted module member reference
            continue
        if "_" not in name and not (name and name[0].isupper()):
            continue                   # prose word (status, heading), not API
        assert hasattr(exp_mod, name) or hasattr(bm_mod, name), (
            f"doc cites {name!r} which is not a public member of the "
            "experiment module or delfin.agent.benchmark")
