"""Every way out we know of, tried for real, and held.

"No escapes" cannot be proved: nobody can show the absence of a route
nobody has thought of. What can be done is to write the routes down, run
each one through the real executor, and assert the file outside the
workspace never appears. The list is the promise, and it grows whenever
somebody thinks of another one.

It is driven through ``execute``, not through the gate helper. Nine
helper cases passed here once while the real gate ran six of the seven
things it was meant to refuse, and the difference was the executor.

Each attempt below is a way of naming a path the GATE cannot read from
the command text — which is the point, because the gate reads text and
the text is not the act:

    a shell variable      echo x > "$HOME"/f
    command substitution  echo x > "$(echo /home/...)"/f
    a symlink             ln -s /home/... link; echo x > link/f
    leaving the folder    cd /home/...; echo x > f
    a nested shell        bash -c 'echo x > /home/.../f'
    an interpreter        python3 -c "open('/home/.../f','w')"
    a decoded path        echo x > "$(echo L2hvbWU... | base64 -d)"
    a copy                cp inside.txt "$HOME"/f
    a heredoc             cat > "$HOME"/f <<EOF

Measured on 539a0f09, in an ordinary attended session, the shell-variable
one wrote into the home directory and came back exit_code 0.

Two outcomes count as held, and they are different things. The gate may
refuse the command before it runs -- there is no exit code then. Or the
command runs and the sandbox denies the write. Either is containment;
what must never happen is the file appearing.
"""

from __future__ import annotations

import json
import os
import pathlib
import tempfile

import pytest

from delfin.agent import api_client as A


ESCAPE_NAME = "delfin-escape-probe.txt"


@pytest.fixture()
def outside(tmp_path):
    """A directory the session was never told about."""
    room = pathlib.Path(tempfile.mkdtemp(prefix="not-the-workspace-"))
    yield room
    for child in room.rglob("*"):
        try:
            child.unlink()
        except OSError:
            pass


@pytest.fixture()
def session(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "inside.txt").write_text("payload\n", encoding="utf-8")
    return ws


def _attempts(outside: pathlib.Path):
    target = outside / ESCAPE_NAME
    import base64
    encoded = base64.b64encode(str(outside).encode()).decode()
    return {
        "a shell variable":
            f'OUT="{outside}"; echo escaped > "$OUT"/{ESCAPE_NAME}',
        "command substitution":
            f'echo escaped > "$(echo {outside})"/{ESCAPE_NAME}',
        "a symlink":
            f'ln -sfn "{outside}" link && echo escaped > link/{ESCAPE_NAME}',
        "leaving the folder":
            f'cd "{outside}" && echo escaped > {ESCAPE_NAME}',
        "a nested shell":
            f"""bash -c 'echo escaped > "{outside}"/{ESCAPE_NAME}'""",
        "an interpreter":
            f"""python3 -c "open('{target}','w').write('escaped')" """,
        "a decoded path":
            f'echo escaped > "$(echo {encoded} | base64 -d)"/{ESCAPE_NAME}',
        "a copy":
            f'cp inside.txt "{outside}"/{ESCAPE_NAME}',
        "a heredoc":
            f'cat > "{outside}"/{ESCAPE_NAME} <<EOF\nescaped\nEOF',
        "an append":
            f'printf escaped >> "{outside}"/{ESCAPE_NAME}',
    }


@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
def test_no_attempt_reaches_outside_the_workspace(session, outside, mode):
    perms = A.KitToolPermissions(workspace=session, mode=mode)
    escaped = []
    for label, command in _attempts(outside).items():
        target = outside / ESCAPE_NAME
        if target.exists():
            target.unlink()
        try:
            A._doc_executor.execute("bash", {"command": command}, perms)
        except Exception as exc:                        # noqa: BLE001
            escaped.append(f"{label}: raised {type(exc).__name__}: {exc}")
            continue
        if target.exists():
            escaped.append(f"{label}: wrote {target}")
            target.unlink()
    assert not escaped, (
        f"a command in {mode} reached outside its workspace:\n  "
        + "\n  ".join(escaped))


def test_the_workspace_itself_still_works(session):
    """A cage that stops the work is not a cage, it is a wall."""
    perms = A.KitToolPermissions(workspace=session, mode="acceptEdits")
    out = json.loads(A._doc_executor.execute(
        "bash", {"command": "mkdir -p sub && echo hi > sub/a.txt && cat sub/a.txt"},
        perms))
    assert out.get("exit_code") == 0, out
    assert (session / "sub" / "a.txt").read_text().strip() == "hi"


@pytest.mark.parametrize("form", ["{d}/f", '"{d}/f"', '"{d}"/f', "'{d}'/f"])
def test_a_granted_directory_is_reachable_however_it_is_quoted(
        session, outside, form):
    """Containment must not make --add-dir a lie, and quoting is not a
    different intention.

    Found by the corpus above: the write-target reader stripped quotes
    only at the ENDS of a word, so `"dir"/file` came out as `dir"/file` --
    a path that exists nowhere, and a directory the user had granted was
    refused as "outside what you may modify". `dir/file` and `"dir/file"`
    passed; `"dir"/file` and `'dir'/file` did not. That punishes exactly
    the careful quoting a model writes around a variable.
    """
    perms = A.KitToolPermissions(workspace=session, mode="acceptEdits",
                                 extra_workspace_dirs=[str(outside)])
    name = "ok-" + str(abs(hash(form)))[:6] + ".txt"
    target = form.format(d=outside) .replace("/f", "/" + name)
    out = json.loads(A._doc_executor.execute(
        "bash", {"command": f"echo granted > {target}"}, perms))
    assert out.get("exit_code") == 0, out
    assert (outside / name).read_text().strip() == "granted"


@pytest.mark.parametrize("form", ['"{d}"/f', "'{d}'/f"])
def test_quoting_does_not_open_a_directory_nobody_granted(
        session, tmp_path, form):
    """The other half. Reading the word the way the shell will must make
    the gate see the REAL path, not accept more paths."""
    ungranted = tmp_path / "ungranted"
    ungranted.mkdir()
    perms = A.KitToolPermissions(workspace=session, mode="acceptEdits")
    target = form.format(d=ungranted).replace("/f", "/x.txt")
    A._doc_executor.execute("bash", {"command": f"echo a > {target}"}, perms)
    assert not (ungranted / "x.txt").exists()


def test_reading_a_credential_is_refused(session):
    """The deny list is about paths, and it holds under a variable too."""
    perms = A.KitToolPermissions(workspace=session, mode="acceptEdits")
    out = str(A._doc_executor.execute(
        "bash", {"command": 'P="$HOME/.ssh/id_rsa"; cat "$P"'}, perms))
    assert "BEGIN" not in out, out[:200]
