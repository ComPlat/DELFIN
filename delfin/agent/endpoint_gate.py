"""A host-wide, per-endpoint limit on concurrent model requests.

Problem measured in a 60-agent run (2026-10): many DELFIN processes on one
host called the same model endpoint (the KIT toolbox) independently.
Together they exceeded its capacity: requests waited 600 s without a byte
("No byte from the endpoint in 600 s"), some ended in errors, and every
process retried on its own, which made the congestion worse. Nothing in
DELFIN coordinated requests across sessions or processes.

The gate is a counting semaphore over FILE LOCKS, not a threads lock:

    <runtime_dir>/endpoint_gates/<key>/
        slot-0 .. slot-N     flock'd for the lifetime of the holder;
                             released automatically when the process dies
        state.json           the adaptive limit (written atomically)
        state.lock           flock serialising limit updates

The runtime directory is the service's own (delfin.agent.service.rundir:
/tmp/delfin-<uid>, 0700, refused if it is not safely the user's own), so
every DELFIN process of the user shares it and no other user can touch it.
A slot's flock is held by an open file descriptor; the kernel releases it
when the process dies -- crashed or killed, there is nothing to leak and
no stale-holder sweep to run.

The limit is adaptive (AIMD) and SHARED, per (endpoint, model): every
process reads the current limit before it waits and updates it after a
result.

    start   agent.endpoint_max_concurrency (default 16; 0 disables the gate)
            is the TOTAL cap of one endpoint: how many concurrent requests
            may hit it at all, across every model behind it.
    halve   an HTTP 429, a 503, or a request that saw no first byte within
            the stall timeout (the "no byte from the endpoint" deadline).
            Each (endpoint, model) halves ONLY its OWN limit, so a slow
            model (GLM: first byte often minutes) stops throttling the fast
            ones (DeepSeek, Mistral) behind the same endpoint.
    floor   the per-model halving floor is 2, not 1 (a model idling at 1
            queue is unusable: one slow request blocks everything).
    +1      (a) after _QUICK_SUCCESS_RUN consecutive requests that succeeded
            quickly (first byte within _QUICK_SUCCESS_S), OR (b) after
            _QUIET_RECOVER_S without a congestion signal -- time recovery,
            the safety net for load where no quick success ever occurs.
            Never above the configured maximum.

When per-model limits are on, every request ALSO takes one slot in the
endpoint's own gate (the total per-endpoint cap), so N models each at
their own high limit cannot together exceed the endpoint's total cap.

The gate never makes a request fail that would otherwise have succeeded:
every flock error, a missing runtime directory, a corrupt state file, any
unexpected exception disables the gate for that call (best effort, the
request goes out ungated).
"""

from __future__ import annotations

import errno
import json
import os
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Optional

#: Default for ``agent.endpoint_max_concurrency``. Zero disables the gate.
DEFAULT_MAX_CONCURRENCY = 16

#: Consecutive quick successes that raise the limit by one again.
_QUICK_SUCCESS_RUN = 8

#: A first byte within this time counts as a quick success -- the endpoint
#: was not saturated while the request ran. Roughly half of the measured
#: KIT median first-token time (18.5 s, 2026-10-19).
_QUICK_SUCCESS_S = 10.0

#: After this long with no congestion signal the limit climbs by one again
#: (time-based additive increase). This is the recovery that works under
#: load, where no quick success ever occurs to trigger the run-based climb:
#: without it a limit halved while the endpoint stayed saturated never rose
#: again (measured 2026-10: KIT AIMD sat at 1 for hours).
_QUIET_RECOVER_S = 120.0

#: The halving floor of ONE (endpoint, model) limit. 2, not 1: a limit of 1
#: makes a whole model nearly unusable (one slow request queues everything
#: behind it) and the model never produces the quick success that climbs
#: back. 1 request at a time is also below a usable queue for any real
#: retry ladder.
_HALVING_FLOOR = 2

#: How long one wait for a slot may park on the state lock before it gives
#: up and the request goes out ungated. A healthy holder holds it for
#: milliseconds; anything longer means the file is stuck (NFS hiccup) and
#: throttling the request would be worse than skipping the gate.
_STATE_LOCK_TIMEOUT_S = 5.0

#: Maximum wait for a free slot before raising :class:`SlotWaitTimeout`.
#: Generous on purpose: the point of the gate is that waiting here replaces
#: waiting inside the overloaded endpoint (which answered nothing for
#: 600 s in the measurement above). Still bounded, so a wedged endpoint
#: cannot hang a turn forever -- the request raises and the client's own
#: transient-retry ladder takes over.
_SLOT_WAIT_TIMEOUT_S = 900.0



def display_url(url: str) -> str:
    """``url`` without userinfo (``user:secret@``), for messages a user reads."""
    try:
        from urllib.parse import urlsplit, urlunsplit
        p = urlsplit(str(url or ""))
        host = p.hostname or ""
        if p.port:
            host = f"{host}:{p.port}"
        return urlunsplit((p.scheme, host, p.path, "", ""))
    except Exception:
        return ""

class SlotWaitTimeout(RuntimeError):
    """No slot became free within :data:`_SLOT_WAIT_TIMEOUT_S`."""


# ---------------------------------------------------------------------------
# The state: where it lives, what it says
# ---------------------------------------------------------------------------

def gate_root() -> Optional[Path]:
    """The directory the gates live in, or None when there is none.

    The service runtime directory is shared by every DELFIN process of the
    user on this host -- that sharing is the whole point. Its ``ensure``
    refuses a directory that is not safely ours (rundir.UnsafeRuntimeDir);
    a refusal disables the gate rather than failing the request.
    """
    try:
        from .service.rundir import ensure_runtime_dir
        return ensure_runtime_dir() / "endpoint_gates"
    except Exception:
        return None


def _endpoint_key(base_url: str, model: str = "",
                  per_model: bool = False) -> str:
    """The key one limit is shared under: the endpoint's base URL, lower
    cased, with a scheme marker so ``http://a:1`` never collides with
    ``https://a:1``. With ``per_model`` the model joins it -- for
    endpoints that limit per model, not per endpoint."""
    # Only scheme, host, port and path: user info (a password in the URL),
    # query and fragment never become part of a file name.
    from urllib.parse import urlsplit
    raw = (base_url or "").strip()
    try:
        parts = urlsplit(raw)
        if parts.scheme and parts.hostname:
            host = parts.hostname + (f":{parts.port}" if parts.port else "")
            raw = f"{parts.scheme}://{host}{parts.path}"
        else:
            raw = raw.split("?", 1)[0].split("#", 1)[0].rsplit("@", 1)[-1]
    except ValueError:
        raw = raw.split("?", 1)[0].split("#", 1)[0].rsplit("@", 1)[-1]
    u = raw.lower().rstrip("/")
    if not u:
        return ""
    key = u.replace("://", "_")
    safe = "".join(c if c.isalnum() or c in "-._" else "_" for c in key)
    if per_model:
        safe += "@" + "".join(
            c if c.isalnum() or c in "-._" else "_" for c in (model or ""))
    # Long keys overflow filename limits (255); a checksum suffix keeps
    # them unique. A checksum, not a cryptographic hash: it protects
    # nothing, it only tells two long keys apart.
    if len(safe) > 160:
        import zlib
        digest = f"{zlib.crc32(safe.encode('utf-8')) & 0xffffffff:08x}"
        safe = safe[:140] + "-" + digest
    return safe


def _config_limit() -> int:
    """``agent.endpoint_max_concurrency`` from the user's settings.
    Zero or missing -> the default; a non-number never raises."""
    limit = DEFAULT_MAX_CONCURRENCY
    try:
        from delfin.user_settings import load_settings
        raw = ((load_settings() or {}).get("agent") or {}).get(
            "endpoint_max_concurrency")
        if isinstance(raw, (int, float)) and not isinstance(raw, bool):
            limit = int(raw)
    except Exception:
        pass
    try:
        env = os.environ.get("DELFIN_ENDPOINT_MAX_CONCURRENCY", "")
        if env.strip():
            limit = int(float(env))
    except (TypeError, ValueError):
        pass
    return limit


def _read_state(dirpath: Path) -> dict:
    try:
        return json.loads((dirpath / "state.json").read_text())
    except (OSError, ValueError):
        return {}


def _write_state(dirpath: Path, state: dict) -> None:
    """Atomically replace the state file: a reader sees the old limit or
    the new one, never a torn file."""
    try:
        from . import state_paths
        state_paths.write_text_atomic(
            dirpath / "state.json",
            json.dumps(state, indent=2))
    except Exception:
        pass


class _FileLock:
    """A blocking flock with a deadline, for the state file."""

    def __init__(self, path: Path) -> None:
        self._path = path

    def __enter__(self) -> "_FileLock":
        import fcntl
        self._fd = os.open(str(self._path), os.O_RDWR | os.O_CREAT, 0o600)
        deadline = time.monotonic() + _STATE_LOCK_TIMEOUT_S
        while True:
            try:
                fcntl.flock(self._fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                return self
            except OSError as exc:
                if exc.errno not in (errno.EWOULDBLOCK, errno.EAGAIN,
                                     errno.EACCES):
                    os.close(self._fd)
                    raise
                if time.monotonic() >= deadline:
                    os.close(self._fd)
                    raise TimeoutError(
                        f"state lock busy: {self._path}") from None
                time.sleep(0.02)

    def __exit__(self, *exc) -> None:
        import fcntl
        try:
            fcntl.flock(self._fd, fcntl.LOCK_UN)
        finally:
            os.close(self._fd)


# ---------------------------------------------------------------------------
# The gate itself
# ---------------------------------------------------------------------------

class EndpointGate:
    """One endpoint's shared semaphore, resolved fresh per call.

    Constructing one is cheap and side-effect free; :meth:`acquire` does
    the work, so a limit the user changed applies without rebuilding the
    client.
    """

    def __init__(self, base_url: str, model: str = "") -> None:
        self.base_url = base_url or ""
        self.model = model or ""
        # Per-model is the DEFAULT: a slow model's congestion must not
        # throttle the fast models behind the same endpoint. The env var
        # now only DISABLES it (flip falsy -> legacy one-limit-per-endpoint);
        # an unset/empty value leaves the default ON.
        raw_pm = (os.environ.get("DELFIN_ENDPOINT_PER_MODEL_LIMIT", "")
                  .strip().lower())
        self._per_model = raw_pm not in ("0", "false", "no", "off")
        # The per-(endpoint, model) adaptive state: its own limit, its own
        # slots. ``_endpoint_key`` names the endpoint-wide TOTAL cap, the
        # slots every model behind this endpoint must also take.
        self._model_key = _endpoint_key(self.base_url, self.model,
                                        per_model=True)
        self._endpoint_key = _endpoint_key(self.base_url, "", per_model=False)
        if self._per_model:
            self._key = self._model_key
            self._cap_key = self._endpoint_key
        else:
            # Legacy mode: one limit per endpoint; there is no separate
            # endpoint cap (the single gate already spans the endpoint).
            self._key = self._endpoint_key
            self._cap_key = ""
        #: Emitted as a notice while waiting for a slot (the wait/notice
        #: mechanism the stream loop already shows the user).
        self.wait_notice = ""

    @property
    def enabled(self) -> bool:
        """Whether the gate is active: a keyed endpoint, a runtime
        directory, and a positive limit (0 disables it)."""
        if not self._key or _config_limit() <= 0:
            return False
        return gate_root() is not None

    # -- the shared limit ---------------------------------------------------

    def _dir(self) -> Optional[Path]:
        root = gate_root()
        if root is None or not self._key:
            return None
        try:
            from . import state_paths
            return state_paths.ensure_dir(root / self._key)
        except Exception:
            return None

    def _cap_dir(self) -> Optional[Path]:
        """The directory of the endpoint-wide TOTAL cap gate (empty when
        per-model is off and the single gate already spans the endpoint)."""
        if not self._cap_key:
            return None
        root = gate_root()
        if root is None:
            return None
        try:
            from . import state_paths
            return state_paths.ensure_dir(root / self._cap_key)
        except Exception:
            return None

    def current_limit(self, dirpath: Optional[Path] = None) -> int:
        """The limit in force right now: the configured maximum, or the
        lower adaptive value a peer halved after congestion."""
        config = _config_limit()
        if config <= 0:
            return 0
        if dirpath is None:
            dirpath = self._dir()
        if dirpath is None:
            return config
        try:
            state = _read_state(dirpath)
            adaptive = int(state.get("limit", config))
        except Exception:
            return config
        # A peer's halving survives only while it stays under the
        # configured maximum: raising the setting takes effect at once.
        return max(1, min(config, adaptive)) if adaptive > 0 else config

    def _note_result(self, dirpath: Path, *, halve: bool,
                     quick: bool) -> int:
        """Fold one result into the shared adaptive limit, returning the
        limit now in force. Atomic under the state lock; any failure is
        swallowed -- the limit is an optimisation, not a correctness
        requirement.

        Congestion (``halve``) halves the limit down to
        :data:`_HALVING_FLOOR` and records when it happened. Every
        non-congestion result that carries no signal -- quick or slow
        success -- may add one back: immediately once
        :data:`_QUICK_SUCCESS_RUN` quick successes have accumulated, or by
        one once :data:`_QUIET_RECOVER_S` have passed without a congestion
        signal (the recovery that works under load).
        """
        config = _config_limit()
        if config <= 0:
            return 0
        try:
            with _FileLock(dirpath / "state.lock"):
                state = _read_state(dirpath)
                try:
                    limit = int(state.get("limit", config))
                except (TypeError, ValueError):
                    limit = config
                limit = max(1, min(config, limit))
                now = time.monotonic()
                try:
                    last_cong = float(
                        state.get("last_congestion_ts", now))
                except (TypeError, ValueError):
                    last_cong = now
                if halve:
                    limit = max(_HALVING_FLOOR, limit // 2)
                    state["limit"] = limit
                    # A congestion signal resets BOTH climb paths: the quick
                    # run and the quiet-period baseline start over, or two
                    # successes after one halving would walk the limit back
                    # up while the endpoint is still saturated.
                    state["quick_run"] = 0
                    state["last_congestion_ts"] = now
                    _write_state(dirpath, state)
                    return limit
                # A non-congestion result: the endpoint answered, so it has
                # spare capacity and may recover one slot.
                incremented = False
                if quick:
                    run = int(state.get("quick_run", 0) or 0) + 1
                    if run >= _QUICK_SUCCESS_RUN and limit < config:
                        limit += 1
                        state["quick_run"] = 0
                        incremented = True
                    else:
                        state["quick_run"] = run
                else:
                    # A slow (but successful) first byte restarts the quick
                    # run: the endpoint was not fast. It still does not
                    # block the quiet-period recovery -- no congestion
                    # occurred.
                    state["quick_run"] = 0
                if (not incremented and now - last_cong >= _QUIET_RECOVER_S
                        and limit < config):
                    limit += 1
                    incremented = True
                if incremented:
                    # One recovery per quiet period: re-anchor the baseline
                    # so the next +1 needs another full quiet spell.
                    state["last_congestion_ts"] = now
                state["limit"] = limit
                _write_state(dirpath, state)
                return limit
        except Exception:
            return self.current_limit(dirpath)

    # -- acquire / release ---------------------------------------------------

    def acquire(self, *, timeout_s: float = _SLOT_WAIT_TIMEOUT_S
                ) -> Optional[int]:
        """Take one slot, or None when the gate is disabled.

        Polls the slot files with non-blocking flocks: slot ``i`` is taken
        by flocking ``slot-<i>``; a file whose flock is held by another
        process is in use. A holder that died had its descriptor closed by
        the kernel, so its slot is free -- no stale-holder sweep.  Never
        raises except :class:`SlotWaitTimeout` after ``timeout_s`` with no
        slot free anywhere (the caller turns that into the same
        transient-retry path a stalled request already takes).

        With per-model limits on, the request ALSO takes one slot in the
        endpoint's own total-cap gate ``_cap_dir``, so N models each at
        their own high limit cannot together exceed the endpoint's cap.
        Both slots are acquired under one shared deadline, in a fixed
        order (model first, then endpoint); releasing happens in the
        opposite order in :meth:`release`. On an endpoint-cap timeout the
        already-taken model slot is released before re-raising, so no
        request ever holds a model slot ungated.
        """
        if not self.enabled:
            return None
        dirpath = self._dir()
        if dirpath is None:
            return None
        limit = self.current_limit(dirpath)
        deadline = time.monotonic() + max(0.0, float(timeout_s))
        if (getattr(self, "_held_fd", None) is not None
                or getattr(self, "_cap_fd", None) is not None):
            # Re-entrant acquire would overwrite the outer slot's descriptor
            # and leak its flock on release. One set of slots per gate
            # instance; a second request in the same process gets its own
            # EndpointGate.
            raise RuntimeError(
                "EndpointGate already holds a slot; build a fresh "
                "EndpointGate for the second request")
        fd, index = self._take_slot(dirpath, limit, deadline)
        if fd is None:
            return None                       # gate broken, go ungated
        self._held_fd = fd
        if self._cap_key:
            cap_dir = self._cap_dir()
            if cap_dir is None:
                self._close(fd); self._held_fd = None
                return None                   # cap gate broken, ungated
            try:
                cap_fd, _ = self._take_slot(cap_dir, _config_limit(),
                                            deadline)
            except SlotWaitTimeout:
                # The endpoint is saturated by OTHER models: never leave a
                # model slot held ungated. Free it and let the caller treat
                # the wait as what it was -- waiting on the endpoint.
                self._close(fd); self._held_fd = None
                raise
            if cap_fd is None:
                self._close(fd); self._held_fd = None
                return None                   # cap gate broken, ungated
            self._cap_fd = cap_fd
        return index

    def _take_slot(self, dirpath: Path, n_slots: int,
                   deadline: float):
        """Flock one free slot file in ``dirpath`` within ``deadline``.

        Returns ``(fd, index)`` on success, ``(None, None)`` when the gate
        is broken (an open or flock error -- the request runs ungated),
        and raises :class:`SlotWaitTimeout` when ``deadline`` passes with
        every slot held. While it waits it reflects the in-use count into
        :attr:`wait_notice`, which the stream loop surfaces to the user.
        """
        announced = False
        while True:
            in_use = 0
            for i in range(n_slots):
                slot = dirpath / f"slot-{i}"
                try:
                    fd = os.open(str(slot), os.O_RDWR | os.O_CREAT, 0o600)
                except OSError:
                    return None, None         # gate broken, go ungated
                try:
                    import fcntl
                    fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except OSError as exc:
                    os.close(fd)
                    if exc.errno in (errno.EWOULDBLOCK, errno.EAGAIN,
                                     errno.EACCES):
                        in_use += 1
                        continue
                    return None, None         # unexpected: ungated
                # The fd now holds the flock for this process's lifetime;
                # release() undoes it. A process that dies without
                # releasing has the kernel close it: the slot frees itself.
                if in_use:
                    self.wait_notice = (
                        f"waiting for a slot on "
                        f"{self.base_url or 'endpoint'}: {in_use} in use")
                return fd, i
            if not announced:
                announced = True
                self.wait_notice = (
                    f"waiting for a slot on {self.base_url or 'endpoint'}:"
                    f" {in_use} in use")
            if time.monotonic() >= deadline:
                raise SlotWaitTimeout(
                    f"no slot on {display_url(self.base_url)} within the "
                    f"allotted time ({n_slots} slots, all in use)")
            time.sleep(0.05)

    def _close(self, fd) -> None:
        """Unflock and close one slot descriptor. Safe on any fd."""
        try:
            import fcntl
            fcntl.flock(fd, fcntl.LOCK_UN)
        except OSError:
            pass
        try:
            os.close(fd)
        except OSError:
            pass

    def release(self) -> None:
        """Free the slots this gate holds. Held in ``finally`` by the
        caller; safe to call when nothing was acquired. Releases the
        endpoint-cap slot first, then the per-model slot (the reverse of
        the acquire order)."""
        for attr in ("_cap_fd", "_held_fd"):
            fd = getattr(self, attr, None)
            if fd is None:
                continue
            setattr(self, attr, None)
            self._close(fd)

    @contextmanager
    def slot(self, *, timeout_s: float = _SLOT_WAIT_TIMEOUT_S,
             ) -> Iterator[Optional[int]]:
        """Acquire and always release: ::

            with gate.slot() as slot_i:
                ... the request ...

        ``slot_i`` is None when the gate is disabled (the request runs
        ungated). Releasing happens on any exit path, including an
        exception inside the block.
        """
        index = self.acquire(timeout_s=timeout_s)
        try:
            yield index
        finally:
            self.release()


# ---------------------------------------------------------------------------
# The result classification the stream loop feeds back
# ---------------------------------------------------------------------------

def halve_for_exception(exc: BaseException) -> bool:
    """Whether this exception is a congestion signal (halve the shared
    limit): an HTTP 429 (rate limit), a 503, or a stall. Kept next to
    api_client's ``_is_transient_api_error``, class-name + status based so
    the openai exception hierarchy is not imported here."""
    name = type(exc).__name__.lower()
    if "ratelimit" in name or "overloaded" in name:
        return True
    status = getattr(exc, "status_code", None)
    if not isinstance(status, int):
        status = getattr(exc, "status", None)
    if isinstance(status, int) and status in (429, 503):
        return True
    msg = str(exc).lower()
    if " 429" in msg or "rate limit" in msg or " 503" in msg \
            or "service unavailable" in msg:
        return True
    # The stream loop raises TimeoutError("no byte from the endpoint in
    # Ns") when the deadline passed without a first byte: the shared
    # symptom of an endpoint beyond its capacity.
    if isinstance(exc, TimeoutError) and "no byte from the endpoint" in msg:
        return True
    return False


def report_result(gate: Optional["EndpointGate"], *,
                  base_url: str = "", model: str = "",
                  exception: Optional[BaseException] = None,
                  first_byte_s: Optional[float] = None) -> None:
    """Fold one request's outcome into the shared adaptive limit.

    Called by the stream loop after each round: ``exception`` set means
    the round failed (halve when it was congestion); otherwise a first
    byte within :data:`_QUICK_SUCCESS_S` counts toward the climb.  Never
    raises.
    """
    if not isinstance(gate, EndpointGate) or not gate.enabled:
        return
    dirpath = gate._dir()
    if dirpath is None:
        return
    halve = exception is not None and halve_for_exception(exception)
    quick = (exception is None
             and first_byte_s is not None
             and first_byte_s <= _QUICK_SUCCESS_S)
    if exception is not None and not halve:
        return                                   # a deterministic failure: no signal
    try:
        gate._note_result(dirpath, halve=halve, quick=quick)
    except Exception:
        pass


def for_client(base_url: str, model: str) -> Optional[EndpointGate]:
    """The gate a client's request round uses, or None when disabled.
    Constructed per round so a changed setting applies without
    rebuilding the client."""
    try:
        gate = EndpointGate(base_url, model)
        return gate if gate.enabled else None
    except Exception:
        return None
