"""container_exec: each root becomes exactly the one mount it names.

Apptainer's ``--bind`` takes a comma-separated list of
``src[:dst[:opts]]`` specs (``ParseBindPath``); Docker's and Podman's
``-v`` take ``src:dst:opts`` and split on ':'. A path is written into
those specs unescaped, so a ',' in a path splits an Apptainer spec into
several binds (the first without ``ro``), and a ':' shifts the fields
of every spec. Such a path is refused, as is an image reference that
the runtime would read as an option.

Input: roots, masks, cwd and image for ``build_argv``; a session whose
root holds such a character for ``session_argv``. Output: ``ValueError``
(``UnsafeMountPath``) from ``build_argv``; a refusal argv naming the path
and the fix from ``session_argv``.
"""
from __future__ import annotations

import subprocess

import pytest

from delfin.agent import api_client as A
from delfin.agent import container_exec as CE
from delfin.agent.api_client import KitToolPermissions


@pytest.fixture(autouse=True)
def _fresh_choice():
    CE.reset()
    yield
    CE.reset()


def _apptainer_binds(value: str) -> list[tuple[str, str, str]]:
    """The binds Apptainer reads from one ``--bind`` value: comma-separated
    specs, each ``src[:dst[:opts]]``; a bare ``ro``/``rw`` joins the
    previous spec's options (simplified ParseBindPath)."""
    binds: list[list[str]] = []
    for part in value.split(","):
        if part in ("ro", "rw") and binds:
            binds[-1][2] = (binds[-1][2] + "," + part).strip(",")
            continue
        fields = part.split(":", 2)
        src = fields[0]
        dst = fields[1] if len(fields) > 1 else src
        opts = fields[2] if len(fields) > 2 else ""
        binds.append([src, dst, opts])
    return [tuple(b) for b in binds]


def _build(runtime, ws, *, read=(), image="img:1"):
    return CE.build_argv(
        runtime, f"/usr/bin/{runtime}", image, ["/bin/bash", "-c", "true"],
        cwd=str(ws), write_roots=[str(ws)], read_roots=list(read),
        network="none", uid=1000, gid=1000)


def _values(argv, flag):
    return [argv[i + 1] for i, a in enumerate(argv) if a == flag]


def test_a_plain_read_only_root_is_one_read_only_bind(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    ro = tmp_path / "archive"
    ro.mkdir()
    argv = _build("apptainer", ws, read=[ro])
    binds = [b for v in _values(argv, "--bind") for b in _apptainer_binds(v)]
    assert binds == [(str(ws), str(ws), "rw"), (str(ro), str(ro), "ro")]


def test_a_comma_in_a_read_only_root_is_refused_for_apptainer(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    ro = tmp_path / f"data,{tmp_path}"
    ro.mkdir(parents=True)
    # Before: the one --bind value read as three binds, the first of them
    # (tmp_path/data) without "ro".
    with pytest.raises(CE.UnsafeMountPath) as exc:
        _build("apptainer", ws, read=[ro])
    assert str(ro) in str(exc.value)


@pytest.mark.parametrize("runtime", ["apptainer", "podman", "docker"])
def test_a_colon_in_a_root_is_refused(runtime, tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    ro = tmp_path / "a:b"
    ro.mkdir()
    with pytest.raises(CE.UnsafeMountPath):
        _build(runtime, ws, read=[ro])


@pytest.mark.parametrize("runtime", ["apptainer", "podman", "docker"])
def test_a_control_character_in_a_root_is_refused(runtime, tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    ro = tmp_path / "a\nb"
    ro.mkdir()
    with pytest.raises(CE.UnsafeMountPath):
        _build(runtime, ws, read=[ro])


def test_a_backslash_in_a_root_is_refused_for_apptainer(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    ro = tmp_path / "a\\"
    ro.mkdir()
    with pytest.raises(CE.UnsafeMountPath):
        _build("apptainer", ws, read=[ro])


@pytest.mark.parametrize("runtime", ["podman", "docker"])
def test_a_comma_is_kept_for_an_oci_runtime(runtime, tmp_path):
    """``-v`` splits on ':' only; a comma in a path works there today."""
    ws = tmp_path / "ws,1"
    ws.mkdir()
    argv = _build(runtime, ws)
    assert f"{ws}:{ws}:rw" in _values(argv, "-v")


@pytest.mark.parametrize("runtime", ["apptainer", "podman", "docker"])
def test_an_image_that_reads_as_an_option_is_refused(runtime, tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    with pytest.raises(ValueError):
        _build(runtime, ws, image="--privileged")


def test_the_session_refuses_and_names_the_path_and_the_fix(tmp_path, monkeypatch):
    monkeypatch.setattr(CE, "choose", lambda: CE.Choice(
        "apptainer", "/usr/bin/apptainer", "img.sif", ""))
    events = []
    monkeypatch.setattr(A, "_record_security_event",
                        lambda kind, tool, detail, blocked=True:
                        events.append((kind, detail, blocked)))
    ws = tmp_path / "proj,x"
    ws.mkdir()
    perms = KitToolPermissions(workspace=str(ws))
    argv = A._bash_isolation_argv("echo ran", ws, perms, mode="container")
    assert "apptainer" not in argv[0]
    done = subprocess.run(argv, capture_output=True, text=True, timeout=20)
    assert done.returncode == 126
    assert "ran" not in done.stdout
    assert str(ws.resolve()) in done.stderr
    assert "workspace" in done.stderr
    assert any(blocked for _k, _d, blocked in events), events
