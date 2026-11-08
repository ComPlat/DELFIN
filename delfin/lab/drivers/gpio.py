"""A Raspberry Pi 5 GPIO relay/input driver through libgpiod.

The Pi 5's 40-pin header sits on the RP1 chip. RPi.GPIO does not work there
(kernel PR #6144 moved the lines from gpiochip4 to gpiochip0), so this
driver uses libgpiod v2 (the ``gpiod`` package on PyPI, ``request_lines`` /
``LineSettings``). It finds the chip by its label ``pinctrl-rp1`` and never
by number, requests only the lines the manifest declares, and releases them
when the request goes out of scope.

Safety rule inherited from :mod:`delfin.lab.drivers.base`: this driver never
decides whether an action is allowed -- the kernel does. It reports what the
device really did (the line read back) and raises :class:`DriverError` on
any transport failure so the kernel fails closed (drives the safe state and
latches).

There is no device-side timer on a GPIO line (a line keeps its value while
the chip is powered), so :attr:`GpioDriver.capabilities.device_watchdog` is
False; the kernel's dead-man plus a hardware watchdog relay cover that (design
doc, section 6). ``arm_watchdog`` is inherited as "unsupported" and returns
False without touching the chip.

Manifest ``address`` (flat mapping; the manifest validator rejects nested
values). ``chip_label`` is optional and defaults to ``pinctrl-rp1``. Each
declared write/read channel maps to a line offset under its own name, with an
optional ``<channel>_active_low`` boolean:

.. code-block:: yaml

    address:
      chip_label: pinctrl-rp1
      heater: 17            # relay output, active-low (energised by a LOW)
      heater_active_low: true
      overtemp: 22          # input, the over-temperature switch

A two-state ``[on, off]`` channel maps ``on``/``off`` to logical ACTIVE /
INACTIVE; ``active_low`` is handed to ``LineSettings`` so libgpiod does the
physical inversion. ``de-energise-to-trip`` means the safe state is ``off``
(the de-energised line), which is what ``safe_state`` writes. No secret ever
appears here: GPIO needs no credentials, and a manifest key that looks like a
secret is already refused by the manifest validator.
"""

from __future__ import annotations

import threading
import time
from glob import iglob
from typing import Any, Optional

from .base import (Driver, DriverCapabilities, DriverError, Health,  # noqa: F401
                   Reading)

DEFAULT_CHIP_LABEL = "pinctrl-rp1"
DEFAULT_TIMEOUT_S = 2.0

#: The ``[lab]`` pyproject extra a user must install for this driver.
EXTRA = "lab-gpio"
_IMPORT_HINT = (f"the GPIO driver needs the ``gpiod`` library; install it "
                f"with ``pip install delfin[{EXTRA}]`` (libgpiod v2).")


def _logical_value(state: Any) -> bool:
    """Map a manifest state value to a logical line level.

    ``on`` / truthy-first is ACTIVE (True); ``off`` / the fallback is INACTIVE
    (False). The manifest's channels are two-state relays; anything else is a
    driver error because there is no way to map it to one bit.
    """
    if isinstance(state, bool):
        return state
    s = str(state).strip().lower()
    if s in ("on", "high", "true", "1", "active"):
        return True
    if s in ("off", "low", "false", "0", "inactive"):
        return False
    raise DriverError(f"a GPIO channel is two-state (on/off); cannot map "
                      f"{state!r} to a line level")


class _ChipRequest:
    """One held line request per device, released on close."""

    __slots__ = ("request", "offsets", "lock")

    def __init__(self, request, offsets: dict):
        # offsets: line offset -> logical value name ("power") for readback
        self.request = request
        self.offsets = offsets
        self.lock = threading.Lock()

    def close(self) -> None:
        with self.lock:
            if self.request is not None:
                try:
                    self.request.close()
                finally:
                    self.request = None


class GpioDriver(Driver):
    name = "gpiod"
    capabilities = DriverCapabilities(device_watchdog=False, readback=True,
                                      observable_commands=False,
                                      device_lock=False)

    def __init__(self, *, chip_label: str = DEFAULT_CHIP_LABEL,
                 device_dir: str = "/dev", timeout_s: float = DEFAULT_TIMEOUT_S,
                 clock: Any = time.time):
        self.chip_label = chip_label
        self.device_dir = device_dir
        self.timeout_s = float(timeout_s)
        self.clock = clock
        self._gpiod: Optional[Any] = None
        self._opens: dict = {}        # device.id -> _ChipRequest
        self._init_lock = threading.Lock()
        self._io_lock = threading.Lock()  # serialise chip access per driver

    # -- private --------------------------------------------------------------

    def _gpiod_module(self):
        """Import libgpiod lazily; a clear message when the extra is absent."""
        if self._gpiod is None:
            with self._init_lock:
                if self._gpiod is None:
                    try:
                        import gpiod
                        from gpiod import line as _gpiod_line
                    except ImportError as exc:
                        raise DriverError(
                            f"{EXTRA}: {exc}; {_IMPORT_HINT}") from exc
                    self._gpiod = gpiod
                    self._gpiod_line = _gpiod_line
        return self._gpiod

    def _find_chip_path(self) -> str:
        """The path of the chip whose label is ``chip_label``; never a number.

        Scan ``/dev/gpiochip*`` and open each, matching by label. Raise
        ``DriverError`` (fail closed) when none matches.
        """
        gpiod = self._gpiod_module()
        pattern = self.device_dir.rstrip("/") + "/gpiochip*"
        for path in sorted(iglob(pattern)):
            chip = None
            try:
                chip = gpiod.Chip(path)
                if chip.get_label() == self.chip_label:
                    return path
            except Exception:
                # one unreachable chip must not stop the search; the error is
                # surfaced at request time if it was the one we needed
                continue
            finally:
                if chip is not None:
                    try:
                        chip.close()
                    except Exception:
                        pass
        raise DriverError(f"no {self.chip_label!r} gpio chip under "
                          f"{self.device_dir!r} (is the agent user in the "
                          f"gpio group?)")

    def _line_offsets(self, device) -> dict:
        """address: {channel: offset, <channel>_active_low: bool, ...}."""
        addr = device.address or {}
        offsets: dict = {}
        active_low: dict = {}
        for ch in list(device.writes) + list(device.reads):
            offset = addr.get(ch)
            if offset is None:
                raise DriverError(f"address has no line offset for channel "
                                  f"{ch!r} of {device.id}")
            try:
                offsets[ch] = int(offset)
            except (TypeError, ValueError) as exc:
                raise DriverError(
                    f"line offset for {ch!r} of {device.id} must be an "
                    f"integer, got {offset!r}") from exc
            al = addr.get(f"{ch}_active_low", False)
            if not isinstance(al, bool):
                raise DriverError(f"active_low for {ch!r} must be a bool")
            active_low[ch] = al
        return offsets, active_low

    def _request(self, device):
        """Open (once, held for the process) the line request for ``device``."""
        req = self._opens.get(device.id)
        if req is not None:
            return req
        gpiod = self._gpiod_module()
        dline = self._gpiod_line
        offsets, active_low = self._line_offsets(device)
        with self._init_lock:
            req = self._opens.get(device.id)
            if req is not None:
                return req
            path = self._find_chip_path()
            config = {}
            for ch, off in offsets.items():
                direction = (dline.Direction.INPUT if ch in device.reads
                             else dline.Direction.OUTPUT)
                # de-energise-to-trip: a fresh output line starts "off"
                settings = gpiod.LineSettings(
                    direction=direction,
                    active_low=bool(active_low[ch]),
                    output_value=dline.Value.INACTIVE,
                )
                config[off] = settings
            try:
                request = gpiod.request_lines(path, consumer="delfin-lab",
                                              config=config)
            except Exception as exc:
                raise DriverError(f"{device.id}: cannot request its lines on "
                                  f"{path!r}: {exc}") from exc
            self._opens[device.id] = _ChipRequest(request, offsets)
        return self._opens[device.id]

    def _closed(self, fn, *args):
        """Run a gpiod call with the driver's timeout; a hang becomes an error.

        GPIO ioctls are local and normally immediate, but a kernel/driver
        stall must not block the kernel forever: run the call in a daemon
        thread and fail closed on expiry (every call is bounded).
        """
        result: list = []
        error: list = []

        def run():
            try:
                result.append(fn(*args))
            except Exception as exc:                  # transport error: fail closed
                error.append(exc)

        t = threading.Thread(target=run, daemon=True)
        t.start()
        t.join(self.timeout_s)
        if t.is_alive():
            raise DriverError("gpio operation did not finish within "
                              f"{self.timeout_s}s")
        if error:
            raise DriverError(f"gpio transport error: {error[0]}")
        return result[0]

    # -- the interface --------------------------------------------------------

    def read(self, device) -> dict:
        with self._io_lock:
            creq = self._request(device)
            out = {}
            for ch, off in creq.offsets.items():
                val = self._closed(creq.request.get_value, off)
                on = bool(val == self._gpiod_line.Value.ACTIVE)
                spec = device.channel(ch)
                out[ch] = Reading("on" if on else "off",
                                  spec.unit if spec else None,
                                  float(self.clock()))
            return out

    def write(self, device, channel: str, value: Any) -> str:
        with self._io_lock:
            creq = self._request(device)
            offsets = creq.offsets
            if channel not in offsets:
                raise DriverError(f"{device.id} has no gpio line for "
                                  f"channel {channel!r}")
            logical = _logical_value(value)
            level = (self._gpiod_line.Value.ACTIVE if logical
                     else self._gpiod_line.Value.INACTIVE)
            self._closed(creq.request.set_value, offsets[channel], level)
            return f"gpio-{device.id}-{offsets[channel]}"

    def arm_watchdog(self, device, seconds: float) -> bool:
        # A GPIO line keeps its value while the chip is powered; there is no
        # device-side timer. The kernel's dead-man and a hardware watchdog
        # relay cover it. Report honestly and touch nothing.
        return False

    def disarm_watchdog(self, device) -> None:
        return None

    def health(self, device) -> Health:
        try:
            with self._io_lock:
                creq = self._request(device)
                off = next(iter(creq.offsets.values()))
                self._closed(creq.request.get_value, off)
            return Health(True, f"{self.chip_label} readback ok")
        except Exception as exc:
            return Health(False, str(exc))

    # -- release --------------------------------------------------------------

    def close(self) -> None:
        """Release every held line request (safe to call more than once)."""
        with self._init_lock:
            opens, self._opens = self._opens, {}
        for creq in opens.values():
            creq.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass
