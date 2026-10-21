"""``set -o pipefail`` names a shell option, not a file to write.

The write-target parser read the ``-o`` of ``set`` like ``sort -o``, so the
documented test recipe ``set -o pipefail; gate ...`` reported a write to a
file called ``pipefail``. Harmless inside the workspace, it became a refusal
once a repository narrowed its write scope.
"""
from __future__ import annotations

import pytest

from delfin.agent.api_client import _bash_write_targets


@pytest.mark.parametrize("cmd", [
    "set -o pipefail; gate tests/test_x.py -q | tail -30",
    "set -o pipefail",
    "set -o errexit -o nounset; ls",
    "shopt -o -s nounset",
])
def test_a_shell_option_is_not_a_target(cmd):
    assert _bash_write_targets(cmd) == []


def test_sort_dash_o_is_still_a_target():
    assert _bash_write_targets("sort -o out.txt in.txt") == ["out.txt"]


def test_a_redirect_after_set_is_still_a_target():
    assert _bash_write_targets("set -o pipefail; ls > list.txt") == ["list.txt"]
