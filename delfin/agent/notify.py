"""Local desktop notifications + optional remote-trigger webhooks.

Two pieces:

  - ``send_notification(title, body, urgency=...)``
        Fires a desktop notification via ``notify-send`` on Linux,
        ``terminal-notifier`` on macOS, or PowerShell on Windows.
        No network, no third-party deps. Failures silent.

  - ``send_remote_trigger(url, payload)``
        POSTs a small JSON blob to a user-configured webhook
        (so the user's own infra can react to agent events).
        Hard restrictions mirror web_tools._check_url:

          * https only
          * deny-list for localhost/RFC1918/cloud-metadata IPs

        ``url`` MUST be configured in ~/.delfin/settings.json under
        ``"remoteTrigger": {"url": "...", "secret": "..."}``. This
        prevents the agent from blasting arbitrary URLs.

Both calls are fire-and-forget — they never block the agent loop
for more than ~3 seconds and never raise.
"""

from __future__ import annotations

import ipaddress
import json
import os
import platform
import shutil
import socket
import subprocess
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any


_NOTIFY_TIMEOUT_S = 3
_REMOTE_TIMEOUT_S = 5
_BLOCKED_HOST_PATTERNS = (
    "localhost", "127.", "0.0.0.0",
    "169.254.",     # link-local / cloud-metadata
)


def _platform() -> str:
    return platform.system().lower()


def _send_linux(title: str, body: str, urgency: str) -> bool:
    if not shutil.which("notify-send"):
        return False
    try:
        subprocess.run(
            ["notify-send", "-u", urgency, title, body],
            timeout=_NOTIFY_TIMEOUT_S, check=False,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        return True
    except (FileNotFoundError, subprocess.SubprocessError):
        return False


def _send_macos(title: str, body: str, urgency: str) -> bool:
    if shutil.which("terminal-notifier"):
        try:
            subprocess.run(
                ["terminal-notifier", "-title", title, "-message", body],
                timeout=_NOTIFY_TIMEOUT_S, check=False,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
            return True
        except (FileNotFoundError, subprocess.SubprocessError):
            pass
    # Fall back to AppleScript
    if shutil.which("osascript"):
        script = (
            f'display notification {json.dumps(body)} '
            f'with title {json.dumps(title)}'
        )
        try:
            subprocess.run(
                ["osascript", "-e", script],
                timeout=_NOTIFY_TIMEOUT_S, check=False,
            )
            return True
        except (FileNotFoundError, subprocess.SubprocessError):
            return False
    return False


def _send_windows(title: str, body: str, urgency: str) -> bool:
    pw = shutil.which("powershell") or shutil.which("powershell.exe")
    if not pw:
        return False
    script = (
        '[Windows.UI.Notifications.ToastNotificationManager,'
        'Windows.UI.Notifications,ContentType=WindowsRuntime] | Out-Null;'
        f'$xml = [xml]"<toast><visual><binding template=\\"ToastText02\\">'
        f'<text id=\\"1\\">{title}</text>'
        f'<text id=\\"2\\">{body}</text>'
        '</binding></visual></toast>";'
    )
    try:
        subprocess.run(
            [pw, "-NoProfile", "-Command", script],
            timeout=_NOTIFY_TIMEOUT_S, check=False,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        return True
    except (FileNotFoundError, subprocess.SubprocessError):
        return False


def send_notification(
    title: str, body: str = "", *, urgency: str = "normal",
) -> bool:
    """Best-effort desktop notification. Returns True on dispatch."""
    title = str(title)[:120] or "delfin agent"
    body = str(body)[:400]
    if urgency not in ("low", "normal", "critical"):
        urgency = "normal"
    plat = _platform()
    if plat == "linux":
        return _send_linux(title, body, urgency)
    if plat == "darwin":
        return _send_macos(title, body, urgency)
    if plat == "windows":
        return _send_windows(title, body, urgency)
    return False


def _is_blocked_host(host: str) -> bool:
    h = host.lower()
    for pat in _BLOCKED_HOST_PATTERNS:
        if h.startswith(pat) or h == pat.rstrip("."):
            return True
    if h.endswith(".internal") or h.endswith(".local"):
        return True
    try:
        addr = ipaddress.ip_address(socket.gethostbyname(h))
        if addr.is_private or addr.is_loopback or addr.is_link_local:
            return True
    except (socket.gaierror, ValueError):
        pass
    return False


def _load_trigger_config(workspace: Path | None) -> dict:
    paths = [Path.home() / ".delfin" / "settings.json"]
    if workspace is not None:
        paths.append(workspace / ".delfin" / "settings.json")
    cfg: dict[str, Any] = {}
    for p in paths:
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if isinstance(data, dict) and isinstance(data.get("remoteTrigger"), dict):
            cfg.update(data["remoteTrigger"])
    return cfg


@dataclass
class TriggerResult:
    sent: bool
    status_code: int = 0
    error: str = ""


def send_remote_trigger(
    payload: dict[str, Any],
    *,
    workspace: Path | None = None,
    override_url: str = "",
    override_secret: str = "",
) -> TriggerResult:
    """POST a JSON payload to the configured remoteTrigger URL.

    - URL must come from settings.json (or ``override_url`` for tests).
    - Scheme must be https.
    - Host must not resolve to a private / loopback / link-local IP.
    """
    cfg = _load_trigger_config(workspace)
    url = override_url or str(cfg.get("url", ""))
    secret = override_secret or str(cfg.get("secret", ""))
    if not url:
        return TriggerResult(sent=False, error="remoteTrigger.url not configured")
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https":
        return TriggerResult(sent=False, error="https required")
    host = parsed.hostname or ""
    if not host or _is_blocked_host(host):
        return TriggerResult(sent=False, error=f"host blocked: {host}")
    body = json.dumps(payload, default=str).encode("utf-8")
    headers = {"Content-Type": "application/json"}
    if secret:
        headers["X-Delfin-Secret"] = secret
    req = urllib.request.Request(url, data=body, headers=headers, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=_REMOTE_TIMEOUT_S) as resp:
            return TriggerResult(sent=True, status_code=resp.status)
    except urllib.error.HTTPError as exc:
        return TriggerResult(sent=False, status_code=exc.code, error=str(exc))
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        return TriggerResult(sent=False, error=str(exc))


__all__ = [
    "TriggerResult",
    "send_notification",
    "send_remote_trigger",
]
