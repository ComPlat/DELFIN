"""orca-external-tools (OET) as an ExtOpt provider: one model server per job.

orca-external-tools (FACCTs, github.com/faccts/orca-external-tools, GPL-3.0)
ships ``oet_server <method>`` -- a long-lived HTTP server whose worker
processes each keep one model resident -- and ``oet_client``, the ProgExt
ORCA starts for every energy/gradient, which only forwards the request.
DELFIN starts the server once per job, before ORCA, and stops it afterwards;
every GOAT worker's client talks to that one server.  This removes the two
costs the orca-mlips route showed on an 81-atom chromophore
(docs/MLIP_GOAT_GPU.md): a Python start-up with torch on each call, and one
model server per GOAT worker (40 copies of the model, GPU out of memory).

OET is a separately installed external program (its own environment, built
by ``python -m delfin.installer --install orca-external-tools``); DELFIN
starts its commands and reads its distribution metadata, it does not import
or ship any of its code.

The client address is not written into the ORCA input: ProgExt names
``delfin-oet-client`` in the OET environment, a two-line shell wrapper that
reads the server address from ``.delfin_oet_server`` in the ORCA working
directory (written by :func:`oet_server`).  The input therefore stays the
same from run to run, which smart recalc relies on.
"""

from __future__ import annotations

import importlib.metadata
import os
import signal
import socket
import subprocess
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Mapping, Optional

OET_DIST = "oet"
OET_INSTALL = "python -m delfin.installer --install orca-external-tools"
#: Backends DELFIN may hand to OET, and the distribution that must be present
#: in the OET environment for each (the OET extras: aimnet2 -> aimnet,
#: mace -> mace-torch, uma -> fairchem-core).  OET has no ORB calculator.
OET_BACKEND_DISTS = {"aimnet2": "aimnet", "mace": "mace-torch", "uma": "fairchem-core"}
#: The file in the ORCA working directory that names the server address.
ADDRESS_FILE = ".delfin_oet_server"
WRAPPER_NAME = "delfin-oet-client"

_WRAPPER = """#!/bin/sh
# Written by DELFIN (delfin/mlp_tools/oet.py): ORCA's ProgExt for the
# orca-external-tools server of this job. The address is read from the ORCA
# working directory, so the ORCA input does not change from run to run.
addr=127.0.0.1:8888
[ -r {address_file} ] && addr=$(cat {address_file})
exec "{client}" -b "$addr" "$@"
"""


def _env_python(environ: Optional[Mapping[str, str]] = None) -> Optional[str]:
    from delfin.common import external_builders

    python = external_builders.tool_python("oet", environ)
    if python == sys.executable or not os.access(python, os.X_OK):
        return None
    return python


def _sites(python: str) -> list[str]:
    prefix = Path(python).resolve().parent.parent
    return [str(s) for s in sorted(prefix.glob("lib/python*/site-packages"))]


def detect_oet(environ: Optional[Mapping[str, str]] = None) -> dict:
    """The installed OET: version, client/server paths, the backends it can run.

    Distribution metadata only, read from its own environment; nothing of it
    is imported here.  ``backends`` are those whose model package is present
    (:data:`OET_BACKEND_DISTS`).
    """
    out = {"available": False, "version": "", "python": "", "client": "", "server": "",
           "backends": (), "install_hint": OET_INSTALL}
    python = _env_python(environ)
    if not python:
        return out
    sites = _sites(python)
    if not sites:
        return out
    found = list(importlib.metadata.distributions(name=OET_DIST, path=sites))
    if not found:
        return out
    bin_dir = Path(python).parent
    client, server = bin_dir / "oet_client", bin_dir / "oet_server"
    if not (client.is_file() and server.is_file()):
        return out
    backends = tuple(
        backend for backend, dist in OET_BACKEND_DISTS.items()
        if any(True for _ in importlib.metadata.distributions(name=dist, path=sites)))
    out.update(available=True, version=found[0].version or "", python=python,
               client=str(client), server=str(server), backends=backends)
    return out


#: The slim client (``mlip_fast_client=yes``): the same request oet_client
#: sends -- POST {"arguments": [...], "directory": <cwd>} to /calculate
#: (orca-external-tools src/oet/server_client/client.py) -- made with curl,
#: so no Python interpreter starts for each energy/gradient. Its answer is
#: checked the way oet_client checks it: "status" must not be "Error".
FAST_WRAPPER_NAME = "delfin-oet-client-curl"
_FAST_WRAPPER = r"""#!/bin/sh
# Written by DELFIN (delfin/mlp_tools/oet.py): the orca-external-tools request
# of oet_client, sent with curl instead of a Python start-up per call.
addr=127.0.0.1:8888
[ -r {address_file} ] && addr=$(cat {address_file})
esc() {{ printf '%s' "$1" | sed -e 's/\\/\\\\/g' -e 's/"/\\"/g'; }}
args=""
for a in "$@"; do args="$args${{args:+,}}\"$(esc "$a")\""; done
resp=$(curl -sS --max-time 3600 -H 'Content-Type: application/json' \
  -d "{{\"arguments\":[$args],\"directory\":\"$(esc "$(pwd)")\"}}" \
  "http://$addr/calculate") || {{ echo "oet request failed: $resp"; exit 1; }}
case "$resp" in
  *'"status":"Error"'*|*'"status": "Error"'*) echo "$resp"; exit 1;;
  *'"status"'*) exit 0;;
  *) echo "unexpected answer from the oet server: $resp"; exit 1;;
esac
"""


def _write_if_different(target: Path, text: str) -> str:
    try:
        if target.read_text(encoding="utf-8") == text and os.access(target, os.X_OK):
            return str(target)
    except OSError:
        pass
    tmp = target.with_name(f".{target.name}.{os.getpid()}")
    tmp.write_text(text, encoding="utf-8")
    tmp.chmod(0o755)
    os.replace(tmp, target)
    return str(target)


def client_wrapper(info: Optional[dict] = None, *, fast: bool = False) -> str:
    """The absolute ProgExt for OET in its environment's ``bin``, written if needed.

    ``delfin-oet-client`` runs the documented ``oet_client``;
    ``fast=True`` gives ``delfin-oet-client-curl``, the same request made with
    curl.  Rewritten only when the content differs (atomic rename), so
    concurrent jobs never see half a file.
    """
    info = info or detect_oet()
    if not info.get("available"):
        raise RuntimeError(f"orca-external-tools is not installed; {OET_INSTALL}")
    bin_dir = Path(info["client"]).parent
    if fast:
        return _write_if_different(bin_dir / FAST_WRAPPER_NAME,
                                   _FAST_WRAPPER.format(address_file=ADDRESS_FILE))
    return _write_if_different(bin_dir / WRAPPER_NAME,
                               _WRAPPER.format(address_file=ADDRESS_FILE, client=info["client"]))


def ext_params(backend: str, *, model: Optional[str] = None, device: str = "cuda") -> str:
    """``Ext_Params`` for an OET request: the model and device, per request.

    The server is started with the method only; the calculator's options
    (``-m``, ``-d``) travel with each request, as ``oet_client`` forwards them.
    """
    if backend not in OET_BACKEND_DISTS:
        raise ValueError(f"backend {backend!r} is not an orca-external-tools method "
                         f"({', '.join(OET_BACKEND_DISTS)})")
    parts = []
    if model:
        parts.append(f"-m {model}")
    parts.append(f"-d {device}")
    return " ".join(parts)


def _signal_group(proc, sig) -> None:
    try:
        os.killpg(proc.pid, sig)
    except (ProcessLookupError, PermissionError):
        try:
            proc.send_signal(sig)
        except ProcessLookupError:
            pass


def free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return int(s.getsockname()[1])


def _healthy(address: str, timeout: float = 2.0) -> bool:
    import urllib.request

    try:
        with urllib.request.urlopen(f"http://{address}/healthz", timeout=timeout) as resp:
            return resp.status == 200
    except Exception:
        return False


@contextmanager
def oet_server(backend: str, work_dir: Path | str, *, processes: int,
               info: Optional[dict] = None, startup_timeout: float = 600.0,
               log_name: str = "oet_server.log") -> Iterator[str]:
    """Run one ``oet_server <backend>`` for the duration of the block.

    ``processes`` is the server's ``-n``: worker processes, each holding one
    model, and the core budget a request's NCores is counted against.  The
    address is written to ``<work_dir>/.delfin_oet_server`` for the client
    wrapper; the server is stopped (SIGTERM, then SIGKILL) when the block ends,
    whatever ends it.  Yields the address.
    """
    info = info or detect_oet()
    if not info.get("available"):
        raise RuntimeError(f"orca-external-tools is not installed; {OET_INSTALL}")
    work_dir = Path(work_dir)
    port = free_port()
    address = f"127.0.0.1:{port}"
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME")}
    env["PYTHONNOUSERSITE"] = "1"
    log = open(work_dir / log_name, "a", encoding="utf-8")
    proc = subprocess.Popen(
        [info["server"], backend, "-b", address, "-n", str(int(processes))],
        cwd=str(work_dir), env=env, stdout=log, stderr=subprocess.STDOUT,
        stdin=subprocess.DEVNULL, start_new_session=True,
    )
    try:
        deadline = time.monotonic() + startup_timeout
        while not _healthy(address):
            if proc.poll() is not None:
                raise RuntimeError(f"oet_server exited with code {proc.returncode} before it "
                                   f"answered; see {work_dir / log_name}")
            if time.monotonic() > deadline:
                raise RuntimeError(f"oet_server did not answer within {startup_timeout:.0f} s")
            time.sleep(0.5)
        (work_dir / ADDRESS_FILE).write_text(address + "\n", encoding="utf-8")
        yield address
    finally:
        if proc.poll() is None:
            # The server and its model processes form one process group.
            _signal_group(proc, signal.SIGTERM)
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                _signal_group(proc, signal.SIGKILL)
                proc.wait(timeout=10)
        log.close()
        try:
            (work_dir / ADDRESS_FILE).unlink()
        except OSError:
            pass
