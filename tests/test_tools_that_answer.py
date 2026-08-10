"""Installed is not the same as working, and only one of the two matters.

Every failure users of this dashboard actually hit came *after* a green
banner: a binary that could not start, an xtb reading somebody else's
parameter files, an ordinary xtb quietly pretending to be g-xTB, a program
that segfaults on every call reported as ok. So the check asks the tools a
question whose answer is known, rather than asking the filesystem whether a
file is there.
"""

from __future__ import annotations

import os
import pathlib
import shutil

import pytest

from delfin import qm_health

_needs_xtb = pytest.mark.skipif(not shutil.which('xtb'), reason='xtb not installed')


# ---------------------------------------------------------------------------
# the conditions, read rather than tested
# ---------------------------------------------------------------------------
def test_a_parameter_file_in_the_home_directory_is_reported(tmp_path, monkeypatch):
    """The one that actually happened.

    xtb reads ``param_gfn2-xtb.txt`` from the home directory instead of the
    parameters compiled into it. A damaged one there ends every GFN2 run in
    ``Error termination. Backtrace:`` while GFN-FF, which does not read them,
    goes on working -- so the user sees a tool that half works and a message
    with nothing in it to act on.
    """
    home = tmp_path / 'home'
    home.mkdir()
    (home / 'param_gfn2-xtb.txt').write_text('not a real parameter file\n')
    monkeypatch.setenv('HOME', str(home))

    seen = qm_health.probe_environment()
    assert seen['home_parameters'] == [str(home / 'param_gfn2-xtb.txt')]
    assert any('instead of its own parameters' in w for w in seen['warnings'])

    (home / 'param_gfn2-xtb.txt').unlink()
    assert qm_health.probe_environment()['home_parameters'] == []


def test_a_stack_limit_that_kills_a_hessian_is_reported(monkeypatch):
    """It cannot be probed, so it is read.

    A 122-atom GFN2 Hessian was measured dying with signal 11 and no message
    at a 512 KiB soft limit and finishing at 2 MiB -- while a water single
    point passes at 512 KiB. No cheap probe can tell the difference, so the
    limit is reported rather than tested.
    """
    import resource

    real = resource.getrlimit
    monkeypatch.setattr(resource, 'getrlimit',
                        lambda which: (512 * 1024, real(which)[1]))
    warnings = qm_health.probe_environment()['warnings']
    assert any('stack limit' in w for w in warnings)

    monkeypatch.setattr(resource, 'getrlimit',
                        lambda which: (8 * 1024 * 1024, real(which)[1]))
    assert not any('stack limit' in w
                   for w in qm_health.probe_environment()['warnings'])


def test_an_unset_thread_count_is_reported(monkeypatch):
    """Unset means one thread per core: on the 384-core machine this was
    measured on, ``xtb --version`` alone spent 2.1 CPU seconds."""
    monkeypatch.delenv('OMP_NUM_THREADS', raising=False)
    assert any('thread count' in w
               for w in qm_health.probe_environment()['warnings'])
    monkeypatch.setenv('OMP_NUM_THREADS', '1')
    assert not any('thread count' in w
                   for w in qm_health.probe_environment()['warnings'])


# ---------------------------------------------------------------------------
# the tools themselves
# ---------------------------------------------------------------------------
@_needs_xtb
def test_xtb_is_asked_a_question_whose_answer_is_known():
    health = qm_health.check_tool('xtb', depth='answer')

    assert health.present and health.runnable
    assert health.healthy, health.why
    assert health.level == 'ok'
    assert health.version, 'it did not say which xtb it is'
    # And the evidence is the run, not a claim about it.
    assert health.evidence['energy'] == pytest.approx(-0.3273, abs=1e-3)
    assert '--gfnff' in health.evidence['command']


@_needs_xtb
def test_the_probe_does_not_leave_its_litter_in_the_working_directory(tmp_path,
                                                                     monkeypatch):
    """xtb drops charges, wbo, xtbrestart and xtbtopo.mol where it runs, and
    GFN-FF adds gfnff_topo. A check must not put those in a user's project."""
    monkeypatch.chdir(tmp_path)
    qm_health.check_tool('xtb', depth='answer')
    assert sorted(p.name for p in tmp_path.iterdir()) == []


@_needs_xtb
def test_a_binary_that_cannot_start_is_told_apart_from_one_that_is_missing(
        tmp_path, monkeypatch):
    """A tool copied out of the environment holding its libraries is present,
    executable, and unable to run. That is not the same problem as an absent
    one and does not have the same answer."""
    fake = tmp_path / 'xtb'
    fake.write_text('#!/bin/sh\n'
                    'echo "xtb: error while loading shared libraries: '
                    'libmctc-lib.so.0: cannot open shared object file" >&2\n'
                    'exit 127\n')
    fake.chmod(0o755)
    monkeypatch.setattr(qm_health, 'find_tool',
                        lambda name: {'path': str(fake), 'source': 'qm_tools',
                                      'dangling': ''})

    health = qm_health.check_tool('xtb', depth='answer')
    assert health.present, 'it is there'
    assert not health.runnable and health.level == 'fail'
    assert 'libmctc-lib.so.0' in health.why
    assert health.repair, 'nothing was offered to put it right'


def test_a_dead_link_says_so_instead_of_saying_nothing(monkeypatch):
    """A dangling qm_tools/bin/xtb is what the installer leaves when the
    environment it pointed into is removed. ``is_file()`` follows the link, so
    it read as "xtb was not found" while ls plainly showed an xtb -- which a
    user cannot act on. "It points at a path that is gone" they can."""
    monkeypatch.setattr(
        qm_health, 'find_tool',
        lambda name: {'path': '', 'source': '', 'dangling': '/gone/bin/xtb'})

    health = qm_health.check_tool('xtb', depth='present')
    assert health.level == 'fail'
    assert '/gone/bin/xtb' in health.why
    assert health.repair == 'relink'


def test_a_program_that_crashes_on_every_call_is_not_probed():
    """anmr segfaults on --help. Running it proves nothing either way, so it
    is never run and never called broken for crashing."""
    assert qm_health.PROBES['anmr']['runnable'] is False
    if shutil.which('anmr'):
        health = qm_health.check_tool('anmr', depth='answer')
        assert health.present and health.level == 'warn'
        assert 'segmentation fault' in health.why


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb not installed')
def test_an_ordinary_xtb_cannot_pass_itself_off_as_gxtb(monkeypatch):
    """An ordinary xtb takes --gxtb, warns once and runs GFN2. The energy is
    the only thing that says so: -5.07 Eh against g-xTB's -76.44."""
    monkeypatch.setattr(
        qm_health, 'find_tool',
        lambda name: {'path': shutil.which('xtb'), 'source': 'PATH',
                      'dangling': ''})

    health = qm_health.check_tool('gxtb', depth='answer')
    assert not health.healthy
    assert 'ordinary xtb' in health.why
    assert health.evidence['energy'] == pytest.approx(-5.0703, abs=1e-3)


# ---------------------------------------------------------------------------
# putting it right
# ---------------------------------------------------------------------------
def test_a_repair_is_shown_before_it_is_run():
    """The same contract the installer already keeps: nothing happens that the
    user has not read first."""
    for name, action in (('xtb', 'install'), ('xtb', 'relink'),
                         ('xtb', 'xtbpath')):
        command = qm_health.repair_command(name, action)
        assert command is None or isinstance(command, list)
    assert qm_health.repair_command('xtb', 'no such action') is None
    for action in qm_health.REPAIRS:
        assert qm_health.REPAIRS[action]


def test_nothing_licensed_is_installed_behind_the_users_back():
    """ORCA, Turbomole and Multiwfn are licensed and are fetched by hand."""
    for name in ('orca', 'turbomole', 'multiwfn', 'gaussian'):
        assert name not in qm_health.INSTALLABLE


@_needs_xtb
def test_a_repair_proves_itself_afterwards():
    """A repair that is not re-checked is a claim."""
    answer = qm_health.repair_tool('xtb', 'xtbpath')
    assert 'health' in answer, 'it did not look again'
    assert answer['health']['name'] == 'xtb'
    assert answer['ok'] is True, answer['status']


def test_several_tools_are_checked_at_once():
    rows = qm_health.check_tools(['xtb', 'crest', 'orca'], depth='present')
    assert [r.name for r in rows] == ['xtb', 'crest', 'orca']
    text = qm_health.format_health(rows)
    assert 'xtb' in text
    # It ends with a count, so a glance is enough.
    assert any(word in text.splitlines()[-1]
               for word in ('ok', 'warn', 'fail', 'absent'))


def test_locating_everything_is_cheap_enough_for_a_background_check():
    """It may run while the dashboard starts, so it must not be felt."""
    import time

    started = time.perf_counter()
    qm_health.check_tools(depth='present')
    assert time.perf_counter() - started < 2.0


# ---------------------------------------------------------------------------
# saying what went wrong, rather than that something did
# ---------------------------------------------------------------------------
def test_the_reason_is_reported_and_not_the_announcement():
    """A user got "GFN2-xTB: Error termination. Backtrace:." and nothing else.

    That line says a run ended, not why. It was chosen by scanning the output
    backwards for the last line containing the word "error" -- and for the
    whole family of Fortran runtime errors that is always the terminator,
    because the numbered backtrace frames below it contain no words at all.
    The reason sat two lines higher and went in the bin with the other
    thirteen kilobytes.
    """
    from delfin.dashboard.gfn_optimize import why_it_stopped

    fabina = (
        'At line 34 of file ../src/mctc/io/write/xyz.f90\n'
        'Fortran runtime error: Unit number is negative and unit was not '
        'already opened with OPEN(NEWUNIT=...)\n\n'
        'Error termination. Backtrace:\n'
        '#0  0x5fb4c3935926 in ???\n'
        '#1  0x5fb4c3936aad in ???\n')
    said = why_it_stopped(fabina)
    assert 'Fortran runtime error' in said
    assert 'xyz.f90' in said, 'the file and line are half the diagnosis'
    assert 'Backtrace' not in said

    # xtb's own error block: every detail line, not the last one.
    block = ('[ERROR] Program stopped due to fatal error\n'
             '-1- Found *very* short distance of 0.000E+00 for H2-H3\n'
             '-2- geometry is not sane\n'
             '####################################\n'
             'abnormal termination of xtb\n')
    said = why_it_stopped(block)
    assert 'short distance' in said and 'not sane' in said
    assert '#' not in said, 'a row of hashes is not a reason'

    # A signal, which has no reason line of its own.
    assert 'SIGABRT' in why_it_stopped(
        'Program received signal SIGABRT: Process abort signal.\n'
        'Backtrace for this error:\n#0 0x7f00 in ???\n')

    # Nothing usable: say so, and say where to look.
    empty = why_it_stopped('   |   banner   |\n')
    assert 'without saying why' in empty
    assert 'disk' in empty


def test_a_terminator_is_never_offered_as_a_reason():
    from delfin.dashboard.gfn_optimize import why_it_stopped

    for terminator in ('Error termination. Backtrace:',
                       'Backtrace for this error:',
                       'ERROR STOP',
                       'abnormal termination of xtb',
                       '[ERROR] Program stopped due to fatal error'):
        said = why_it_stopped(f'something happened\n{terminator}\n')
        assert terminator.lower() not in said.lower(), terminator


def test_the_complaint_is_found_even_though_it_arrives_last():
    """The two streams are merged and are not buffered alike: the terminators
    go to stderr and arrive at once, everything xtb printed to stdout is
    flushed when it exits and lands after them. Measured on a damaged
    parameter file: ERROR STOP on line 73 of the captured output and the
    reason on line 102, the last line of all."""
    from delfin.dashboard.gfn_optimize import why_it_stopped

    merged = ('ERROR STOP \n'
              'Error termination. Backtrace:\n'
              '#0  0x7f11 in ???\n'
              '   Cite this work as:\n'
              '   * C. Bannwarth, S. Ehlert and S. Grimme.\n'
              ' no basis found for atom           1  Z=            8\n')
    said = why_it_stopped(merged)
    assert 'no basis found' in said
    assert 'Bannwarth' not in said, 'the citation list is not a reason'


# ---------------------------------------------------------------------------
# installing means installing
# ---------------------------------------------------------------------------
def test_the_installer_builds_its_own_before_adopting_someone_elses():
    """"Install" mostly did not install.

    It looked for whatever xtb stood first on the PATH and made a symlink to
    it, so two accounts on one cluster ran two different builds from the same
    button -- one of them optimising and the other stopping with a Fortran
    format error inside xtb's own optimiser, on the same structure. Pressing
    Install again could not help: it adopted the same broken build again.

    The managed environment now comes first and a system tool is the fallback,
    which keeps the case that made adoption the default: a cluster with no
    network and xtb behind a module.
    """
    from delfin.dashboard.gfn_optimize import install_script

    script = install_script()
    assert script is not None
    text = script.read_text(encoding='utf-8')

    assert 'PREFER_SYSTEM_TOOLS="${PREFER_SYSTEM_TOOLS:-0}"' in text, (
        'adopting whatever is on the PATH is no longer the default'
    )
    assert 'install_managed_or_adopt' in text
    for tool in ('install_xtb', 'install_crest', 'install_dftbplus'):
        body = text.split(f'{tool}() {{')[1].split('\n}')[0]
        assert 'install_managed_or_adopt' in body, tool
        # And when it does adopt, it says so rather than calling it an install.
    adopt = text.split('install_managed_or_adopt() {')[1].split('\n}')[0]
    assert 'adopting the' in adopt
    assert 'cannot vouch for it' in adopt


def test_each_tool_gets_an_environment_of_its_own():
    """One environment for all three installed the broken xtb.

    Asked for xtb, crest and dftbplus together, conda-forge answers **xtb
    6.6.1** -- whose optimiser dies with "Fortran runtime error: Missing comma
    between descriptors" on a doublet optimisation that 6.7.1 finishes
    normally -- because 6.6.1 is the newest xtb compatible with the other two
    at the same time. Asked separately, each gets its own newest: xtb 6.7.1,
    crest 3.0.2, dftbplus 25.1. Three dry runs of the solver, then the real
    installer, which produced 6.7.1.

    So a user could press Install and end up with a broken xtb while the
    person beside them, whose xtb came from somewhere else, had no trouble.
    """
    from delfin.dashboard.gfn_optimize import install_script

    text = install_script().read_text(encoding='utf-8')

    assert 'install_conda_tool()' in text
    # Never the three of them in one solve again.
    assert 'conda-forge xtb crest dftbplus' not in text
    assert '"xtb=${XTB_VERSION}" crest dftbplus' not in text

    # And a floor under xtb, so the build that could not optimise cannot come
    # back -- a floor rather than a pin, so newer ones are still taken.
    assert 'XTB_MINIMUM="${XTB_MINIMUM:-' in text
    assert '"xtb>=${XTB_MINIMUM}"' in text


def test_a_failure_names_the_binary_that_produced_it():
    """Two accounts on one cluster do not necessarily run the same xtb, and
    without the path and version in the message there is nothing to compare."""
    from delfin.dashboard.gfn_optimize import which_xtb_ran

    said = which_xtb_ran('/opt/xtb/bin/xtb', 'xtb version 6.7.1 (edcfbbe)')
    assert '/opt/xtb/bin/xtb' in said
    assert '6.7.1' in said
    assert 'that build' in said, 'it should say where the fault lives'
    assert which_xtb_ran('', '') == ''


def test_an_xtb_that_cannot_optimise_is_stepped_over(monkeypatch):
    """A tool directory is searched before the PATH, so one broken xtb left in
    it beat the working one standing right behind it -- and every optimisation
    failed for a user whose colleague on the same machine had no trouble.

    6.6.1 dies with "Missing comma between descriptors" at optimizer.f90:639
    on every optimisation, GFN2 and GFN-FF alike; 6.7.1 finishes the same job.
    Only a single point without --opt survives on 6.6.1. So the version is a
    fact about whether the tool can do the work, and it is now read before the
    tool is used.
    """
    from delfin.dashboard import gfn_optimize as gfn

    gfn._XTB_JUDGED.clear()
    monkeypatch.setattr(gfn, 'judge_xtb', lambda path: {
        'ok': '/good/' in str(path), 'version': '6.6.1' if '/old/' in str(path) else '6.7.1',
        'why': 'this is xtb 6.6.1, and below 6.7.0 an optimisation dies',
        'path': str(path)})
    monkeypatch.setattr(gfn, '_xtb_candidates',
                        lambda: ['/old/bin/xtb', '/good/bin/xtb'])

    assert gfn.find_xtb() == '/good/bin/xtb', 'the broken one was used'
    note = gfn.unusable_xtb_note()
    assert '/old/bin/xtb' in note and 'not being used' in note

    # None usable: still name a real binary rather than saying there is none.
    monkeypatch.setattr(gfn, 'judge_xtb', lambda path: {
        'ok': False, 'version': '6.6.1', 'why': 'too old', 'path': str(path)})
    assert gfn.find_xtb() == '/old/bin/xtb'


def test_the_judgement_is_made_once_and_remembered(tmp_path):
    """One --version, about ten milliseconds, per binary per session. Not
    asking cost a user an afternoon."""
    import time

    from delfin.dashboard import gfn_optimize as gfn

    fake = tmp_path / 'xtb'
    fake.write_text('#!/bin/sh\necho "   * xtb version 6.7.1 (abc) compiled"\n')
    fake.chmod(0o755)
    gfn._XTB_JUDGED.clear()

    first = gfn.judge_xtb(str(fake))
    assert first['ok'] and first['version'] == '6.7.1'

    started = time.perf_counter()
    for _ in range(200):
        gfn.judge_xtb(str(fake))
    assert time.perf_counter() - started < 0.2, 'it is asking every time'

    # A reinstall at the same path is a different binary and is asked again.
    fake.write_text('#!/bin/sh\necho "   * xtb version 6.6.1 (old) compiled"\n')
    fake.chmod(0o755)
    assert gfn.judge_xtb(str(fake))['ok'] is False


def test_the_shared_environment_goes_once_nothing_uses_it():
    """A binary still on disk can still be reached -- an explicit path saved in
    the settings, a PATH somebody exported -- so the xtb that cannot optimise
    goes on being the one that runs, long after a working one was installed
    beside it. It is removed only after the replacement is in place, and only
    when no link points into it."""
    from delfin.dashboard.gfn_optimize import install_script

    text = install_script().read_text(encoding='utf-8')
    body = text.split('retire_legacy_env() {')[1].split('\n}')[0]
    assert 'conda-meta' in body, 'it must be sure it is a conda environment'
    assert 'still uses it' in body, 'a tool that needs it keeps it'
    assert 'readlink -f' in body
    # And it runs after the replacement, never before.
    tool = text.split('install_conda_tool() {')[1].split('\n}')[0]
    assert tool.index('link_into_bin') < tool.index('retire_legacy_env')


# ---------------------------------------------------------------------------
# installed where it is looked for
# ---------------------------------------------------------------------------
def test_a_tool_is_looked_for_where_the_buttons_put_it(tmp_path, monkeypatch):
    """It was neither, reliably.

    The Settings tab installs into the user's own copy, the Submit tab's
    button and the agent's installer wrote into the packaged one, and the
    resolver answered "the packaged one" unless a setting had been saved. So a
    user pressed Install, the tools landed in ~/.delfin/qm_tools, and every one
    of them was reported missing -- measured on a real machine with dftb+,
    stda and xtb4stda all present on disk and all reported absent.
    """
    from delfin import qm_runtime

    monkeypatch.delenv('DELFIN_QM_TOOLS_ROOT', raising=False)
    monkeypatch.delenv('DELFIN_QM_ROOT', raising=False)

    user_root = tmp_path / 'user'
    (user_root / 'bin').mkdir(parents=True)
    packaged = tmp_path / 'packaged'
    (packaged / 'bin').mkdir(parents=True)
    monkeypatch.setattr(qm_runtime, 'get_user_qm_tools_root', lambda: user_root)
    monkeypatch.setattr(qm_runtime, 'get_packaged_qm_tools_root', lambda: packaged)

    # Nothing installed yet: the packaged copy is as good an answer as any.
    assert qm_runtime.get_qm_tools_root() == packaged

    # Something installed by the dashboard: that is where it is looked for.
    (user_root / 'bin' / 'xtb').write_text('#!/bin/sh\n')
    assert qm_runtime.get_qm_tools_root() == user_root

    # And both are searched either way, so a tool in either is found.
    looked = [str(p) for p in qm_runtime.iter_qm_tools_bin_dirs()]
    assert str(user_root / 'bin') in looked
    assert str(packaged / 'bin') in looked

    # An explicit choice still wins over both.
    monkeypatch.setenv('DELFIN_QM_TOOLS_ROOT', str(tmp_path / 'chosen'))
    assert qm_runtime.get_qm_tools_root() == (tmp_path / 'chosen')


def test_the_install_button_writes_where_the_resolver_reads(tmp_path,
                                                            monkeypatch):
    """Left to itself the installer installs beside its own file, and the copy
    the Submit tab runs lives inside the package."""
    from delfin.dashboard import gfn_optimize as gfn

    monkeypatch.setenv('DELFIN_QM_TOOLS_ROOT', str(tmp_path / 'tools'))
    command = gfn.install_command('xtb')

    assert command[0] == 'env', 'the root has to reach the script'
    joined = ' '.join(command)
    assert f'DELFIN_QM_TOOLS_ROOT={tmp_path / "tools"}' in joined
    assert 'install_qm_tools.sh xtb' in joined
    assert gfn.install_root() == tmp_path / 'tools'

    # And the script listens: it used to install beside itself whatever it was
    # told, which is why naming the directory changed nothing.
    text = gfn.install_script().read_text(encoding='utf-8')
    assert 'ROOT="${DELFIN_QM_TOOLS_ROOT:-${DELFIN_QM_ROOT:-' in text
