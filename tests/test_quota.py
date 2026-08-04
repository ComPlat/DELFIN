"""HOME usage must come from quota accounting on every file system, never du."""

from pathlib import Path

import pytest

from delfin import quota


# --------------------------------------------------------------------------
# Parsers -- fed with real-world output shapes, no cluster required.

def test_lfs_quota_single_line():
    output = (
        'Disk quotas for usr ka_ew7404 (uid 506470):\n'
        '      Filesystem  kbytes   bquota  blimit  bgrace   files   iquota\n'
        '/home/ka/ka_ibcs/ka_ew7404 250652044  524288000 576716800  '
        '-  479667  5000000\n'
    )
    used, soft = quota.parse_lfs_quota(output)
    assert used == 250652044 * 1024
    assert soft == 524288000 * 1024


def test_lfs_quota_mountpoint_wrapped_to_own_line():
    """lfs wraps long mount points, so column indexing would break."""
    output = (
        '/a/very/long/mount/point/that/wrapped\n'
        '         1024  2048  4096  -  10  0  0  -\n'
    )
    used, soft = quota.parse_lfs_quota(output)
    assert used == 1024 * 1024
    assert soft == 2048 * 1024


def test_lfs_quota_over_limit_marker_is_stripped():
    """Values over quota carry a trailing '*' and must still parse."""
    used, soft = quota.parse_lfs_quota('/home 600000*  500000  550000  none  1  0\n')
    assert used == 600000 * 1024
    assert soft == 500000 * 1024


def test_lfs_quota_without_limit_reports_none():
    used, soft = quota.parse_lfs_quota('/home 12345  0  0  -  7  0  0  -\n')
    assert used == 12345 * 1024
    assert soft is None


def test_lfs_quota_empty_output():
    assert quota.parse_lfs_quota('') is None
    assert quota.parse_lfs_quota('some error text\n') is None


def test_posix_quota():
    output = (
        'Disk quotas for user foo (uid 1000):\n'
        '     Filesystem  blocks   quota   limit   grace   files   quota  limit\n'
        '      /dev/sda1  123456  200000  250000            1234       0      0\n'
    )
    used, soft = quota.parse_posix_quota(output)
    assert used == 123456 * 1024
    assert soft == 200000 * 1024


def test_posix_quota_none_when_no_quotas_configured():
    assert quota.parse_posix_quota('') is None
    assert quota.parse_posix_quota('Disk quotas for user foo (uid 1000): none\n') is None


def test_mmlsquota_gpfs():
    output = (
        '                         Block Limits                    |     File Limits\n'
        'Filesystem type  KB      quota    limit  in_doubt  grace | files  quota\n'
        'gpfs1      USR   123456  200000   250000        0   none |  1234      0\n'
    )
    used, soft = quota.parse_mmlsquota(output)
    assert used == 123456 * 1024
    assert soft == 200000 * 1024


def test_mmlsquota_none_on_header_only():
    assert quota.parse_mmlsquota('Filesystem type KB quota limit\n') is None


def test_beegfs_quota_csv_is_bytes_not_kb():
    output = (
        'name,id,size,hard,files,hard\n'
        'foo,1000,1073741824,10737418240,100,0\n'
    )
    used, soft = quota.parse_beegfs_quota(output)
    assert used == 1073741824
    assert soft == 10737418240


def test_beegfs_quota_zero_hard_limit_is_unlimited():
    used, soft = quota.parse_beegfs_quota('name,id,size,hard\nfoo,1000,2048,0\n')
    assert used == 2048
    assert soft is None


# --------------------------------------------------------------------------
# Provider selection

def test_home_usage_prefers_provider_matching_the_filesystem(monkeypatch):
    calls = []

    def fake_run(cmd):
        calls.append(cmd[0])
        if cmd[0] == 'stat':
            return 'lustre\n'
        if cmd[0] == 'lfs':
            return '/home 1024 2048 4096 - 1 0 0 -\n'
        return None

    monkeypatch.setattr(quota, '_run', fake_run)
    usage = quota.home_usage(Path('/home/foo'), user='foo')

    assert usage is not None
    assert usage.source == 'lustre'
    assert usage.used_bytes == 1024 * 1024
    # No other quota tool was consulted once Lustre answered.
    assert 'quota' not in calls and 'mmlsquota' not in calls


def test_home_usage_falls_back_when_hinted_provider_is_absent(monkeypatch):
    """A Lustre site whose lfs is missing must still fall through to quota(1)."""

    def fake_run(cmd):
        if cmd[0] == 'stat':
            return 'lustre\n'
        if cmd[0] == 'quota':
            return (
                'Disk quotas for user foo (uid 1000):\n'
                '  Filesystem blocks quota limit\n'
                '  /dev/sda1  512    1024  2048\n'
            )
        return None

    monkeypatch.setattr(quota, '_run', fake_run)
    usage = quota.home_usage(Path('/home/foo'), user='foo')

    assert usage is not None
    assert usage.source == 'posix'
    assert usage.used_bytes == 512 * 1024


def test_home_usage_returns_none_when_nothing_reports_quota(monkeypatch):
    """Unknown is the correct answer -- a tree walk never is."""
    monkeypatch.setattr(quota, '_run', lambda cmd: None)
    assert quota.home_usage(Path('/home/foo'), user='foo') is None


def test_home_usage_survives_a_provider_raising(monkeypatch):
    def fake_run(cmd):
        if cmd[0] == 'stat':
            return 'ext2/ext3\n'
        if cmd[0] == 'lfs':
            raise OSError('boom')
        if cmd[0] == 'quota':
            return (
                'Disk quotas for user foo (uid 1000):\n'
                '  /dev/sda1  512    1024  2048\n'
            )
        return None

    monkeypatch.setattr(quota, '_run', fake_run)
    usage = quota.home_usage(Path('/home/foo'), user='foo')
    assert usage is not None and usage.source == 'posix'


def test_home_usage_ignores_a_provider_reporting_zero(monkeypatch):
    """Zero usage means "this tool does not know", not "the home is empty"."""

    def fake_run(cmd):
        if cmd[0] == 'stat':
            return 'lustre\n'
        if cmd[0] == 'lfs':
            return '/home 0 0 0 - 0 0 0 -\n'
        if cmd[0] == 'quota':
            return 'Disk quotas for user foo:\n  /dev/sda1  512  1024  2048\n'
        return None

    monkeypatch.setattr(quota, '_run', fake_run)
    usage = quota.home_usage(Path('/home/foo'), user='foo')
    assert usage is not None and usage.source == 'posix'


# --------------------------------------------------------------------------
# The regression that triggered all of this

def test_no_provider_ever_shells_out_to_du(monkeypatch):
    """Cluster ops blocked the account over `du -sb $HOME`; it must not return."""
    seen = []

    def fake_run(cmd):
        seen.append(cmd)
        return None

    monkeypatch.setattr(quota, '_run', fake_run)
    quota.home_usage(Path('/home/foo'), user='foo')

    assert seen, 'expected the providers to be attempted'
    for cmd in seen:
        assert cmd[0] != 'du', f'directory walk reintroduced: {cmd}'


def test_dashboard_module_contains_no_du_call():
    source = Path(quota.__file__).with_name('dashboard') / '__init__.py'
    text = source.read_text(encoding='utf-8')
    assert "'du'" not in text and '"du"' not in text
    assert 'du -s' not in text
