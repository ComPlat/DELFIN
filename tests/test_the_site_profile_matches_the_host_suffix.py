"""The site profile is chosen from the host name's domain suffix.

``detect_site_profile`` tested ``'scc.kit.edu' in fqdn``, which also
matches a host whose name merely contains the string
(``scc.kit.edu.example.org``, ``notscc.kit.edu``) -- CodeQL
py/incomplete-url-substring-sanitization. The domain must be the suffix.
"""

from __future__ import annotations

import pytest

from delfin import slurm_submit


@pytest.mark.parametrize("fqdn, expected", [
    ("uc3n990.localdomain.hpc.scc.kit.edu", "bwunicluster3"),
    ("login1.scc.kit.edu", "bwunicluster3"),
    ("scc.kit.edu", "bwunicluster3"),
    ("scc.kit.edu.example.org", ""),
    ("notscc.kit.edu", ""),
    ("host.example.org", ""),
])
def test_the_domain_must_be_the_suffix(monkeypatch, fqdn, expected):
    monkeypatch.setattr(slurm_submit.socket, "getfqdn", lambda: fqdn)
    monkeypatch.setenv("HOSTNAME", "workstation")
    assert slurm_submit.detect_site_profile() == expected


def test_a_uc3_host_name_still_counts(monkeypatch):
    monkeypatch.setattr(slurm_submit.socket, "getfqdn", lambda: "x.example.org")
    monkeypatch.setenv("HOSTNAME", "uc3n990")
    assert slurm_submit.detect_site_profile() == "bwunicluster3"
