"""A side's worker must be able to authenticate at all.

``bench ab`` runs each side in a subprocess with a fresh ``HOME`` -- a process
boundary, so the attempt cannot read anything a previous attempt left behind.
Two things follow from that boundary, and only one of them was handled:

* the fresh ``HOME`` hides the credential FILE (``~/.delfin/credentials.json``
  lives under the real home), which is the point; and
* it must not hide the credential itself, or the verb can never measure.

``credentials.load_credential`` resolves env var > keyring > file. With the
file hidden by the redirected home, the env var is the only route that
survives the boundary -- and the child's environment was built from scratch
without one. So every attempt died in engine init and the comparison reported
"0/2 vs 0/2" for a run that called nothing.

WHAT IS DELIBERATELY NOT PASSED. The whole parent environment, or every key in
``credentials._WELL_KNOWN_KEYS``. A benchmark measures an agent that can run
shell commands; handing it a GitHub token would hand it the repository. Only
the key of the provider under test crosses the boundary, and only when the
parent actually holds it.
"""

from __future__ import annotations

import pytest


def _ab():
    from delfin.agent import bench_ab

    return bench_ab


def test_the_workers_environment_carries_the_providers_key(tmp_path):
    """The defect: no key reached the child, so no attempt could measure."""
    ab = _ab()
    parent = {"PATH": "/usr/bin:/bin", "ANTHROPIC_API_KEY": "sk-test-value"}
    env = ab._child_env(tmp_path, tmp_path / "home", "", "", parent=parent)

    assert env.get("ANTHROPIC_API_KEY") == "sk-test-value", (
        "the side's worker got no credential, so every attempt fails in "
        "engine init and a comparison reports a rate for a run that never "
        "happened")


@pytest.mark.parametrize("provider, name", [
    ("kit", "KIT_TOOLBOX_API_KEY"),
    ("blablador", "BLABLADOR_API_KEY"),
    ("openai", "OPENAI_API_KEY"),
    ("claude", "ANTHROPIC_API_KEY"),
])
def test_each_provider_gets_its_own_key(tmp_path, provider, name):
    """The key handed over follows ``--provider``, not a single default."""
    ab = _ab()
    parent = {"PATH": "/bin", name: "value-for-" + provider}
    env = ab._child_env(tmp_path, tmp_path / "home", "", "",
                        parent=parent, provider=provider)
    assert env.get(name) == "value-for-" + provider


def test_no_other_credential_crosses_the_boundary(tmp_path):
    """The control: a benchmark agent must not be handed the whole keyring.

    ``an agent that can run shell commands`` is what a benchmark measures. A
    GitHub token in its environment is a token it can spend.
    """
    ab = _ab()
    parent = {
        "PATH": "/bin",
        "ANTHROPIC_API_KEY": "the-one-named",
        "GITHUB_TOKEN": "must-not-cross",
        "GH_TOKEN": "must-not-cross",
        "BRAVE_SEARCH_API_KEY": "must-not-cross",
    }
    env = ab._child_env(tmp_path, tmp_path / "home", "", "", parent=parent)
    for name in ("GITHUB_TOKEN", "GH_TOKEN", "BRAVE_SEARCH_API_KEY"):
        assert name not in env, f"{name} was handed to the benchmark worker"


def test_a_missing_key_is_not_invented(tmp_path):
    """No key in the parent, no key in the child -- and no empty value either."""
    ab = _ab()
    env = ab._child_env(tmp_path, tmp_path / "home", "", "", parent={"PATH": "/bin"})
    assert "ANTHROPIC_API_KEY" not in env


def test_the_redirected_home_still_hides_the_credential_file(tmp_path):
    """The property the fresh home exists for survives the fix."""
    ab = _ab()
    env = ab._child_env(tmp_path, tmp_path / "home", "", "",
                        parent={"PATH": "/bin"})
    assert env["HOME"] != str(tmp_path)
    assert env["DELFIN_SCRATCH_STATE"].startswith(env["HOME"])


def test_every_key_name_the_worker_may_receive_is_a_real_credential():
    """A typo in the map would silently pass nothing.

    ``credentials._WELL_KNOWN_KEYS`` is the list of names the product reads as
    credentials, and a test guards it against drifting from the call sites. A
    name in ``bench_ab`` that is not in it is a name nothing will look up.
    """
    from delfin.agent.credentials import _WELL_KNOWN_KEYS

    ab = _ab()
    unknown = sorted(set(ab.PROVIDER_KEY_NAMES.values()) - set(_WELL_KNOWN_KEYS))
    assert not unknown, (
        "bench ab would hand the worker a variable the credential loader "
        f"never reads: {unknown}")
