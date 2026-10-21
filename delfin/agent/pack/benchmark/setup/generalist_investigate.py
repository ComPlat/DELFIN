#!/usr/bin/env python3
"""Setup for generalist_investigate: a question about the code.

``ratelimit/`` -- a client configured with ``max_retries = 3`` sends a
request 5 times when the server keeps answering 503. Two causes, both
needed for the full count:

* retry.py: ``for attempt in range(policy.max_retries + 1)`` -- one
  first try plus 3 retries = 4 sends;
* client.py: when ``send_with_retry`` gives up, ``fetch`` makes one more
  ``self._send_once(request)`` = the 5th.

Distractors: config.py and settings.toml both say 3 (true and
irrelevant); backoff.py computes delays only.
The acceptance key is located by code text, not by line number, so an
edit to this fixture cannot silently move it.
Precondition: driving the client against an always-503 transport counts
5 sends.
argv[1] is the workspace.
"""
from __future__ import annotations

import subprocess
import sys

from generalist_repo import (GITIGNORE_PY, commit, expect, main, mark_base,
                             new_repo, write)

CONFIG = '''"""Client configuration."""

from dataclasses import dataclass


@dataclass
class RetryPolicy:
    """How often and how patiently to retry."""

    max_retries: int = 3
    base_delay: float = 0.0
    retry_on: tuple = (429, 502, 503, 504)


def load_policy(settings: dict) -> RetryPolicy:
    """A RetryPolicy from the ``[retry]`` table of settings.toml."""
    table = settings.get("retry", {})
    return RetryPolicy(
        max_retries=int(table.get("max_retries", 3)),
        base_delay=float(table.get("base_delay", 0.0)),
    )
'''

BACKOFF = '''"""Delays between attempts."""


def delay_for(attempt: int, base: float) -> float:
    """Exponential backoff: base, 2*base, 4*base, ..."""
    return base * (2 ** attempt)
'''

RETRY = '''"""Retrying a send."""

import time

from ratelimit.backoff import delay_for


class GaveUp(Exception):
    """Every attempt failed with a retryable status."""


def send_with_retry(send, request, policy):
    """Call ``send(request)`` until it returns a non-retryable status."""
    last = None
    for attempt in range(policy.max_retries + 1):
        last = send(request)
        if last.status not in policy.retry_on:
            return last
        time.sleep(delay_for(attempt, policy.base_delay))
    raise GaveUp(last)
'''

CLIENT = '''"""The HTTP-ish client."""

from ratelimit.config import RetryPolicy
from ratelimit.retry import GaveUp, send_with_retry


class Response:
    def __init__(self, status, body=""):
        self.status = status
        self.body = body


class Client:
    """Sends requests through a transport callable."""

    def __init__(self, transport, policy=None):
        self._transport = transport
        self.policy = policy or RetryPolicy()

    def _send_once(self, request):
        return self._transport(request)

    def fetch(self, request):
        """The response to ``request``, retried per the policy."""
        try:
            return send_with_retry(self._send_once, request, self.policy)
        except GaveUp:
            return self._send_once(request)
'''

FILES = {
    ".gitignore": GITIGNORE_PY,
    "README.md": ("# ratelimit\n\nA client that retries on 429/5xx.\n"
                  "Configured in settings.toml (`[retry] max_retries`).\n"),
    "settings.toml": "[retry]\nmax_retries = 3\nbase_delay = 0.0\n",
    "pyproject.toml": ("[project]\nname = \"ratelimit\"\nversion = \"0.9.0\"\n\n"
                       "[tool.pytest.ini_options]\npythonpath = [\".\"]\n"
                       "testpaths = [\"tests\"]\n"),
    "ratelimit/__init__.py": "\"\"\"A retrying client.\"\"\"\n",
    "ratelimit/config.py": CONFIG,
    "ratelimit/backoff.py": BACKOFF,
    "ratelimit/retry.py": RETRY,
    "ratelimit/client.py": CLIENT,
    "tests/test_client.py": (
        "from ratelimit.client import Client, Response\n\n\n"
        "def test_a_success_is_returned_at_once():\n"
        "    calls = []\n"
        "    c = Client(lambda r: calls.append(r) or Response(200, 'ok'))\n"
        "    assert c.fetch('GET /').body == 'ok'\n"
        "    assert len(calls) == 1\n"),
}

PROBE = (
    "import sys; sys.path.insert(0, '.')\n"
    "from ratelimit.client import Client, Response\n"
    "n = []\n"
    "Client(lambda r: n.append(r) or Response(503)).fetch('GET /')\n"
    "print(len(n))\n")


def build(workspace: str):
    repo = new_repo(workspace, "ratelimit")
    write(repo, FILES)
    mark_base(repo, commit(repo, "ratelimit 0.9.0"))
    proc = subprocess.run([sys.executable, "-c", PROBE], cwd=str(repo),
                          capture_output=True, text=True, timeout=60,
                          env={"PYTHONDONTWRITEBYTECODE": "1",
                               "PATH": "/usr/bin:/bin"})
    expect("an always-503 server sees 5 sends",
           proc.stdout.strip() == "5", proc.stdout + proc.stderr)
    return repo


if __name__ == "__main__":
    main(build)
