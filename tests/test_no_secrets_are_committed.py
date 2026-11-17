"""No tracked file carries a generic secret.

Completes the guard family: personal data, account names, CCDC/CSD
licence material, MLIP tokens and home paths in commit messages are
covered elsewhere; this file covers credentials a scanner can recognise
by shape — private keys, provider API keys, bearer credentials,
credential assignments in config-like files, key-material file names,
high-entropy token literals, and private-network addresses in the
developer-facing directories (scripts/, docs/, tools/, README*).

The rules live in ``scripts/secret_guard.py`` (single source of truth,
shared with the pre-commit hook). This test:

1. proves each rule fires on a synthetic sample planted in a temp tree
   and stays silent on the legitimate stand-ins (RFC 5737 addresses,
   placeholder assignments, test data);
2. runs the whole tracked tree — clean today, and kept clean by the
   pre-commit entry ``python scripts/secret_guard.py --staged``.

Samples are assembled from parts so this file stays clean of the shapes
it guards. Everything here is synthetic.

Rejected: a directory-level allow-list (it grows quietly); scanning test
files with the entropy rule (test fixtures legitimately carry long token
shapes); a length-1 JWT rule (any "eyJ..." fragment would flood docs).
Measured cost: one pass over the tracked tree, a few seconds.
"""
from __future__ import annotations

import importlib.util
import pathlib
import subprocess
import sys

import pytest

_ROOT = pathlib.Path(__file__).resolve().parents[1]

_spec = importlib.util.spec_from_file_location(
    "secret_guard", _ROOT / "scripts" / "secret_guard.py")
secret_guard = importlib.util.module_from_spec(_spec)
sys.modules["secret_guard"] = secret_guard
_spec.loader.exec_module(secret_guard)


def _plant(tree: pathlib.Path, relpath: str, text: str) -> None:
    path = tree / relpath
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


# -- 1. every rule refuses its planted synthetic sample ---------------------
# The dangerous fragments are concatenated at run time so this file never
# contains a matching literal itself.

_PEM = ("-----BEGIN " + "RSA PRIVATE" + " KEY-----\n"
        + "MIICXAIBAAKBgQC" + "9xPlantedNotARealKey" + "\n"
        + "-----END " + "RSA PRIVATE" + " KEY-----")
_OPENSSH = ("-----BEGIN OPEN" + "SSH PRIVATE" + " KEY-----\n"
            + "b3BlbnNzaC1rZXktdjEAAAAA" + "planted" + "\n")
_PPK = "PuTTY-User-Key-File-" + "2\nEncryption: none\nComment: test"
_JWT = ("eyJ" + "hbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9"
        + "." + "eyJ" + "zdWIiOiJ0ZXN0In0" + "." + "SflKxwRJSMeKKF2QT4fwpMeJf36")

@pytest.mark.parametrize("rule,text,relpath", [
    ("private_key", _PEM, "docs/example.pem"),
    ("private_key", _OPENSSH, "notes.txt"),
    ("private_key", _PPK, "tools/putty.ppk"),
    ("jwt", f"Authorization: {_JWT}", "docs/api.md"),
    ("api_key_shape", "OPENAI_API_KEY=" + "sk-proj-" + "a1B2c3D4e5F6g7H8i9J0kKlM",
     "README.md"),
    ("api_key_shape", "token: " + "ghp_" + "A" * 36, "docs/ci.md"),
    ("api_key_shape", "xoxp-" + "123-456-abcdef", "docs/slack.md"),
    ("api_key_shape", "key " + "AKIA" + "IOSFODNN7EXAMPLE", "docs/aws.md"),
    ("api_key_shape", "hf_" + "A" * 34, "docs/mlp.md"),
    ("api_key_shape", "glpat-" + "Ab1" * 8, "docs/gl.md"),
    ("bearer_token", "Authorization: Bearer " + "dGhpc2lzYXZlcnlsb25nYmFzZTY0dG9rZW4",
     "docs/rest.md"),
    ("credential_assignment", "password=Sup3r-Secret-Value",
     "configs/app.ini"),
    ("credential_assignment", "api_key = \"0123456789abcdefXYZ\"",
     "settings.toml"),
    ("credential_assignment", "export DB_PASSWORD='hunter2-hunter2'",
     "run.sh"),
    ("secret_path", "key", ".env"),
    ("secret_path", "key", "id_rsa"),
    ("secret_path", "key", "server.pem"),
    ("private_net", "submit to 10.0.1.47 via slurm", "docs/cluster.md"),
    ("private_net", "ssh 192.168.1.10", "scripts/setup.sh"),
    ("private_net", "portal at https://172.16.4.1", "README.md"),
])
def test_each_rule_refuses_a_planted_sample(rule, text, relpath, tmp_path):
    _plant(tmp_path, relpath, text)
    found = secret_guard.run_tree(tmp_path)
    hits = [v for v in found if f"rule {rule}" in v]
    assert hits, f"rule {rule} silent on its planted sample in {relpath}"


# -- 2. the stand-ins stay silent ------------------------------------------

@pytest.mark.parametrize("text,relpath", [
    # RFC 5737 / RFC 3849 documentation addresses are the stand-ins.
    ("reach the mirror at 192.0.2.10 or 198.51.100.7 or 203.0.113.9",
     "docs/cluster.md"),
    # placeholder assignments are not credentials
    ("password=${DB_PASSWORD}\napi_key = <your-key-here>\n"
     "token: changeme\ndb_password = placeholder",
     "configs/app.ini"),
    # a placeholder value in a shell file
    ("export API_KEY=\"your-key-here\"", "run.sh"),
    # a JWT-shaped word alone is not a token; fragment, no dot-separated body
    ("the header part eyJhbGciOiJIUzI1NiJ9 names no token", "docs/api.md"),
    # prose identifiers: short, no digit/letter mix at high entropy
    ("parameter conformational_analysis_beta yields the table",
     "docs/methodology.md"),
    # private addresses inside tests/ or examples are test data
    ("ssh 10.0.0.5 in the test fixture", "tests/test_example.py"),
])
def test_the_stand_ins_stay_silent(text, relpath, tmp_path):
    _plant(tmp_path, relpath, text)
    assert secret_guard.run_tree(tmp_path) == [], relpath


def test_high_entropy_token_refused_in_a_non_test_file(tmp_path):
    token = "Zk9mP2qX7rT4bN8cW1dV5yH3jL6fQ"  # 29 near-unique chars, mixed
    _plant(tmp_path, "src/generated.py", f'key = "{token}"\n')
    found = secret_guard.run_tree(tmp_path)
    assert any("rule high_entropy" in v for v in found)
    # the same string in a test file is test data and stays silent
    _plant(tmp_path, "tests/test_generated.py", f'key = "{token}"\n')
    assert not [v for v in secret_guard.run_tree(tmp_path)
                if "high_entropy" in v and "tests/" in v], found


def test_private_key_in_a_binary_file_is_found(tmp_path):
    blob = b"\x89PNG\r\n" + _PEM.encode() + b"\x00\x01"
    _plant(tmp_path, "images/planted.png", "")
    (tmp_path / "images" / "planted.png").write_bytes(blob)
    found = secret_guard.run_tree(tmp_path)
    assert any("rule private_key" in v for v in found)


def test_allow_marker_exempts_its_line_only(tmp_path):
    token = "Zk9mP2qX7rT4bN8cW1dV5yH3jL6fQ"
    _plant(tmp_path, "docs/synthetic.md",
           f"shape: sk-{token[:20]}  # secret-guard: allow (synthetic)\n"
           f"shape: sk-{token[:20]}\n")
    found = secret_guard.run_tree(tmp_path)
    assert found, "the unmarked line must stay refused"
    assert all("docs/synthetic.md:1" not in v for v in found)


def test_fixtures_are_not_a_whole_directory_exemption(tmp_path):
    # A real-shaped secret in tests/fixtures/x.txt without the inline
    # marker IS reported: the tests/fixtures/ path must not silently pass.
    token = "Zk9mP2qX7rT4bN8cW1dV5yH3jL6fQ"
    _plant(tmp_path, "tests/fixtures/x.txt",
           "OPENAI_API_KEY=" + "sk-proj-" + token[:24] + "-plant\n")
    found = secret_guard.run_tree(tmp_path)
    assert any("api_key_shape" in v for v in found), found
    # the same line carrying the documented inline marker is allowed
    _plant(tmp_path, "tests/fixtures/x.txt",
           "OPENAI_API_KEY=" + "sk-proj-" + token[:24] + "-plant"
           "  # secret-guard: allow (synthetic fixture)\n")
    found = secret_guard.run_tree(tmp_path)
    assert not [v for v in found if "tests/fixtures/x.txt" in v], found


def test_a_bare_allow_marker_is_itself_a_finding(tmp_path):
    token = "Zk9mP2qX7rT4bN8cW1dV5yH3jL6fQ"
    # bare marker (no reason): still a finding, even though it tries to allow
    _plant(tmp_path, "docs/synthetic.md",
           f"shape: sk-{token[:20]}  # secret-guard: allow\n")
    found = secret_guard.run_tree(tmp_path)
    assert any("bad_marker" in v for v in found), found
    # empty reason is also a bad marker
    _plant(tmp_path, "docs/synthetic.md",
           f"shape: sk-{token[:20]}  # secret-guard: allow ()\n")
    found = secret_guard.run_tree(tmp_path)
    assert any("bad_marker" in v for v in found), found



# -- 3. the tracked tree is clean -------------------------------------------

def _tracked() -> list[str]:
    out = subprocess.run(["git", "ls-files"], cwd=str(_ROOT),
                         capture_output=True, text=True, timeout=120)
    return [line for line in out.stdout.splitlines() if line.strip()]


def test_git_is_available_or_this_checks_nothing():
    files = _tracked()
    assert files, "git ls-files returned nothing; this guard is blind"
    assert "scripts/secret_guard.py" in files


def test_no_tracked_file_carries_a_generic_secret():
    violations = secret_guard.run_all()
    assert not violations, "\n".join(violations)


def test_staged_mode_reports_a_planted_staged_change(tmp_path):
    """--staged must apply the same table: plant a secret in a scratch
    clone's index and expect exit 1; clean the index again afterwards."""
    git_env = dict(
        GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@example.invalid",
        GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@example.invalid")
    env = {**dict(__import__("os").environ), **git_env}
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True, env=env)
    _plant(tmp_path, "configs/app.cfg", "password=hunter2-hunter2\n")
    r = subprocess.run(
        [sys.executable, str(_ROOT / "scripts" / "secret_guard.py"), "--staged"],
        cwd=str(tmp_path), capture_output=True, text=True)
    # nothing is staged yet -> clean exit
    assert r.returncode == 0, r.stdout + r.stderr
    subprocess.run(["git", "add", "configs/app.cfg"], cwd=str(tmp_path),
                   check=True, env=env)
    r = subprocess.run(
        [sys.executable, str(_ROOT / "scripts" / "secret_guard.py"), "--staged"],
        cwd=str(tmp_path), capture_output=True, text=True)
    assert r.returncode == 1, r.stdout + r.stderr
    assert "configs/app.cfg" in r.stdout