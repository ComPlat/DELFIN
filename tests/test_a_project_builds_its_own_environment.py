"""A project builds and uses its own environment in the workspace.

Session 0db62cdf (2026-10-19, a browser game): three.js could not be fetched
(the CDN unknown to the egress proxy, refused under Bypass, curl exit 56) and
`npm install` was asked about every time, so the agent rebuilt a 3D renderer
by hand. A venv with pip worked; npm, yarn, pnpm and uv did not.

Pinned: the installs into the project are free and do not leave the
workspace; installs elsewhere (global, a prefix, a path, a URL) are judged
as before; the package indexes are known hosts; each package manager's
cache lands in the cage's temp dir; npm's 407 is named with its cause.
"""

from __future__ import annotations

import tempfile

import pytest

from delfin.agent import api_client as A
from delfin.agent import egress_proxy as EP
from delfin.agent import caged_test_runners as ctr
from delfin.agent import landlock_exec as LL


@pytest.fixture
def perms():
    return A.KitToolPermissions(workspace=tempfile.mkdtemp(),
                                mode="bypassPermissions")


def _stays(perms, cmd):
    return not A._leaves_workspace(perms, {"tool": "bash", "command": cmd},
                                   True)


def _admitted(cmd, ws):
    return bool(ctr.admitted_runners(cmd, str(ws), str(ws)))


@pytest.fixture
def project(tmp_path):
    (tmp_path / "package.json").write_text("{}")
    return tmp_path


INTO_THE_PROJECT = [
    "npm install three", "npm i three@0.160.1 --save", "npm ci", "npm install",
    "npm init -y", "yarn add three", "yarn", "pnpm add -D vite",
    "uv venv", "uv venv .venv --python 3.11", "uv pip install requests numpy",
    "uv pip install -r requirements.txt", "uv add httpx", "uv sync",
]


@pytest.mark.parametrize("cmd", INTO_THE_PROJECT)
def test_an_install_into_the_project_runs_in_the_cage(perms, project, cmd):
    # Admitted where the cage holds (caged_test_runners), not leaving the
    # workspace under Bypass, and NOT on the plain auto-allow list: without
    # a cage the user is asked, as for every stack runner.
    assert _admitted(cmd, project), cmd
    assert _stays(perms, cmd), cmd
    assert not perms.matches_bash_auto_allow(cmd), cmd


@pytest.mark.parametrize("cmd", [
    "python3 -m venv .venv", ".venv/bin/pip install requests",
    ".venv/bin/python -m pip install numpy",
])
def test_the_workspace_venv_is_free_as_before(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd) and _stays(perms, cmd), cmd


@pytest.mark.parametrize("cmd", [
    "npm install -g typescript", "npm i --global x",
    "npm i --location=global x", "yarn global add x",
    "npm install --prefix /usr/local x", "npm i ../other",
    "npm install https://example.org/x.tgz", "npm i git+ssh://host/x",
    "uv pip install --system requests", "uv pip install --target /tmp x",
    "npm install $(cat x)", "pip install requests", "npx serve",
    "npm run build", "cargo install ripgrep",
])
def test_an_install_elsewhere_is_judged_as_before(perms, project, cmd):
    assert not _admitted(cmd, project), cmd
    assert not _stays(perms, cmd), cmd


@pytest.mark.parametrize("host", [
    "registry.npmjs.org", "cdn.jsdelivr.net", "unpkg.com", "static.crates.io",
    "proxy.golang.org", "pypi.org",
])
def test_the_package_indexes_are_known_hosts(host):
    assert EP.domain_allowed(host, EP.DEFAULT_ALLOWED_DOMAINS)


def test_every_cage_gives_each_cache_a_writable_home(monkeypatch, tmp_path):
    names = dict(LL.PACKAGE_CACHE_DIRS)
    assert {"npm_config_cache", "UV_CACHE_DIR"} <= set(names)
    assert A._package_cache_dirs() == tuple(LL.PACKAGE_CACHE_DIRS)
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    monkeypatch.setattr(A, "_BWRAP_RUNTIME_FAILURE", None)
    p = A.KitToolPermissions(workspace=tmp_path, mode="default")
    argv = A._bash_isolation_argv("true", tmp_path, p, caged_unasked=True)
    env = {argv[i + 1]: argv[i + 2] for i, a in enumerate(argv)
           if a == "--setenv"}
    assert env.get("npm_config_cache") == "/tmp/npm-cache"
    assert env.get("UV_CACHE_DIR") == "/tmp/uv-cache"


def test_npms_407_is_named_with_its_cause():
    out = ("npm ERR! code E407\nnpm ERR! 407 Proxy Authentication Required - "
           "GET https://registry.npmjs.org/is-number")
    note = A._proxy_auth_note(1, out)
    assert "npm >= 10" in note and "Do not retry" in note
    assert A._proxy_auth_note(0, out) == ""
    assert A._proxy_auth_note(1, "npm ERR! code E404") == ""
