"""A Raspberry Pi 5 GPIO relay/input driver through libgpiod.

The Pi 5's 40-pin header sits on the RP1 chip. RPi.GPIO does not work there
(kernel PR #6144 moved the lines from gpiochip4 to gpiochip0), so this
driver uses libgpiod v2 (the ``gpiod`` package on PyPI, ``request_lines`` /
``LineSettings``). It finds the chip by its label ``pinctrl-rp1``
(``Chip.get_info().label``) and never by number, and requests only the lines
the manifest declares.

One request per operation, released at once. Measured on a Pi 5: a released
output line keeps its direction and value, a new request with direction
``AS_IS`` reads it unchanged, and while one process holds a line every other
request for it fails with EBUSY. So a read requests outputs ``AS_IS`` (a read
never switches anything) and inputs as ``INPUT``; a write requests the line
as ``OUTPUT`` with the commanded value. Holding the request for the life of
the driver instead would lock the e-stop and dead-man daemons (other
processes) out of the very line they must switch off, and a fresh output
request with a default value would switch a running load off on every read.

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

import errno
import threading
import time
from glob import iglob
from typing import Any, Optional

from .base import (Driver, DriverCapabilities, DriverError, Health,  # noqa: F401
                   Reading)

DEFAULT_CHIP_LABEL = "pinctrl-rp1"
DEFAULT_TIMEOUT_S = 2.0
_BUSY_TRIES = 4
_BUSY_WAIT_S = 0.05

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
        self._chip_path: Optional[str] = None
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
        if self._chip_path is not None:
            return self._chip_path
        gpiod = self._gpiod_module()
        pattern = self.device_dir.rstrip("/") + "/gpiochip*"
        for path in sorted(iglob(pattern)):
            try:
                with gpiod.Chip(path) as chip:
                    label = chip.get_info().label
            except Exception:
                # one unreachable chip must not stop the search; the error is
                # surfaced at request time if it was the one we needed
                continue
            if label == self.chip_label:
                self._chip_path = path
                return path
        raise DriverError(f"no {self.chip_label!r} gpio chip under "
                          f"{self.device_dir!r} (is the agent user in the "
                          f"gpio group?)")

    def _line_offsets(self, device) -> tuple:
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

    def _directions(self, offsets: dict) -> dict:
        """offset -> the line's current direction (no request needed)."""
        gpiod = self._gpiod_module()
        path = self._find_chip_path()

        def op():
            with gpiod.Chip(path) as chip:
                return {off: chip.get_line_info(off).direction
                        for off in offsets.values()}
        return self._closed(op)

    def _with_lines(self, device, config: dict, fn):
        """Request ``config`` (offset -> LineSettings), run ``fn(request)``,
        release. EBUSY (another process holds a line for its own short
        operation) is retried a few times before it fails closed."""
        gpiod = self._gpiod_module()
        path = self._find_chip_path()

        def op():
            last = None
            for attempt in range(_BUSY_TRIES):
                try:
                    request = gpiod.request_lines(path, consumer="delfin-lab",
                                                  config=config)
                except OSError as exc:
                    if getattr(exc, "errno", None) == errno.EBUSY:
                        last = exc
                        if attempt + 1 < _BUSY_TRIES:
                            time.sleep(_BUSY_WAIT_S * (attempt + 1))
                        continue
                    raise
                try:
                    return fn(request)
                finally:
                    request.release()
            raise last
        try:
            return self._closed(op)
        except DriverError as exc:
            raise DriverError(f"{device.id}: {exc}") from None

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
        gpiod = self._gpiod_module()
        dline = self._gpiod_line
        offsets, active_low = self._line_offsets(device)
        directions = self._directions(offsets)
        config = {}
        for ch, off in offsets.items():
            al = bool(active_low[ch])
            if ch in device.reads and ch not in device.writes:
                config[off] = gpiod.LineSettings(
                    direction=dline.Direction.INPUT, active_low=al)
            elif directions.get(off) == dline.Direction.OUTPUT:
                # an output is read as it stands: a read never drives it
                config[off] = gpiod.LineSettings(
                    direction=dline.Direction.AS_IS, active_low=al)
            else:
                # A write line nobody has driven since boot is an input whose
                # level says nothing about the relay. It is taken over in the
                # device's safe state (de-energise-to-trip), then read.
                safe = (dline.Value.ACTIVE
                        if _logical_value(device.safe_state.get(ch, "off"))
                        else dline.Value.INACTIVE)
                config[off] = gpiod.LineSettings(
                    direction=dline.Direction.OUTPUT, active_low=al,
                    output_value=safe)
        with self._io_lock:
            values = self._with_lines(
                device, config,
                lambda req: {ch: req.get_value(off)
                             for ch, off in offsets.items()})
        out = {}
        for ch, val in values.items():
            on = bool(val == dline.Value.ACTIVE)
            spec = device.channel(ch)
            out[ch] = Reading("on" if on else "off",
                              spec.unit if spec else None,
                              float(self.clock()))
        return out

    def write(self, device, channel: str, value: Any) -> str:
        gpiod = self._gpiod_module()
        dline = self._gpiod_line
        offsets, active_low = self._line_offsets(device)
        if channel not in offsets or channel not in device.writes:
            raise DriverError(f"{device.id} has no writable gpio line for "
                              f"channel {channel!r}")
        level = (dline.Value.ACTIVE if _logical_value(value)
                 else dline.Value.INACTIVE)
        off = offsets[channel]
        config = {off: gpiod.LineSettings(
            direction=dline.Direction.OUTPUT,
            active_low=bool(active_low[channel]), output_value=level)}
        with self._io_lock:
            self._with_lines(device, config, lambda req: None)
        return f"gpio-{device.id}-{off}"

    def arm_watchdog(self, device, seconds: float) -> bool:
        # A GPIO line keeps its value while the chip is powered; there is no
        # device-side timer. The kernel's dead-man and a hardware watchdog
        # relay cover it. Report honestly and touch nothing.
        return False

    def disarm_watchdog(self, device) -> None:
        return None

    def health(self, device) -> Health:
        try:
            self.read(device)
            return Health(True, f"{self.chip_label} readback ok")
        except Exception as exc:
            return Health(False, str(exc))

    def close(self) -> None:
        """Nothing is held between operations; kept for callers."""
        return None
