"""Job Status tab: view and cancel running/pending jobs."""

import html as _html
import os
import re
import subprocess
from datetime import datetime, timedelta
from pathlib import Path

import ipywidgets as widgets
from IPython.display import clear_output

from .constants import STATUS_COLORS, JOB_TABLE_CSS
from .backend_base import JobInfo
from .input_processing import parse_inp_resources


#: The job table's own styles. Every class starts with ``djs-`` so nothing
#: leaks into other tabs, and the fonts follow the notebook theme.
_DJS_CSS = """<style>
.djs { font-family: var(--jp-ui-font-family, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif);
       color: #111827; font-size: 13px; container-type: inline-size; width: 100%; }
.djs-summary { display: flex; flex-wrap: wrap; gap: 10px; margin: 2px 0 12px 0; }
.djs-stat { flex: 1 1 150px; max-width: 230px; background: #ffffff; border: 1px solid #e5e7eb;
            border-radius: 8px; padding: 8px 14px; }
.djs-stat-value { font-size: 20px; font-weight: 600; line-height: 1.25; white-space: nowrap;
                  overflow: hidden; text-overflow: ellipsis; }
.djs-stat-suffix { font-size: 13px; font-weight: 500; color: #6b7280; }
.djs-stat-label { font-size: 11px; color: #6b7280; text-transform: uppercase; letter-spacing: .04em; white-space: nowrap; }
.djs-wrap { border: 1px solid #e5e7eb; border-radius: 8px; overflow: hidden; background: #ffffff; }
/* One line per job: fixed columns; text that does not fit ends in an
   ellipsis, and the whole of it is in the cell's tooltip. */
.djs-table { width: 100%; border-collapse: collapse; table-layout: fixed; font-size: 13px; line-height: 1.4; }
.djs-table th, .djs-table td { white-space: nowrap; overflow: hidden; text-overflow: ellipsis; font-size: inherit; }
.djs-table th { background: #f9fafb; color: #6b7280; font-size: 11px; font-weight: 600; text-transform: uppercase;
                letter-spacing: .04em; text-align: left; padding: 8px 10px; border-bottom: 1px solid #e5e7eb; }
.djs-table td { padding: 9px 10px; border-bottom: 1px solid #f3f4f6; vertical-align: middle; }
.djs-table tbody tr:last-child td { border-bottom: none; }
.djs-table tbody tr:hover { background: #fafafa; }
/* Name and reason share what the fixed columns leave. */
.djs-c-status { width: 100px; }
.djs-c-job { width: auto; }
.djs-c-partition { width: 124px; }
.djs-c-resources { width: 150px; }
.djs-c-runtime { width: 196px; }
.djs-c-why { width: auto; }
.djs-c-start { width: 204px; }
.djs-name { font-weight: 600; }
.djs-muted { color: #6b7280; }
.djs-mono { font-family: var(--jp-code-font-family, SFMono-Regular, Menlo, Consolas, monospace); font-size: 12px; }
.djs-pill { display: inline-flex; align-items: center; gap: 6px; padding: 2px 9px; border-radius: 999px;
            font-size: 12px; font-weight: 600; }
.djs-dot { width: 7px; height: 7px; border-radius: 50%; background: currentColor; display: inline-block; flex: none; }
.djs-running { background: #e8f5e9; color: #2e7d32; }
.djs-pending { background: #fff4e5; color: #b26a00; }
.djs-failed { background: #fdecea; color: #d32f2f; }
.djs-other { background: #f3f4f6; color: #374151; }
.djs-tag { display: inline-block; padding: 1px 7px; border-radius: 4px; background: #eef2ff; color: #3730a3;
           font-size: 11.5px; margin-right: 4px; }
.djs-progress { display: flex; align-items: center; gap: 8px; min-width: 0; }
.djs-progress-text { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.djs-bar { flex: none; width: 52px; height: 6px; background: #e5e7eb; border-radius: 3px; overflow: hidden; }
.djs-bar > span { display: block; height: 100%; background: #2e7d32; }
.djs-bar.djs-warn > span { background: #ef6c00; }
.djs-reason { font-weight: 600; color: #b26a00; }
.djs-reason.djs-stuck { color: #d32f2f; }
.djs-note { margin-top: 8px; color: #6b7280; font-size: 12px; }
.djs-empty { padding: 28px 16px; text-align: center; color: #6b7280; border: 1px dashed #d1d5db;
             border-radius: 8px; background: #fafafa; }
.djs-empty-title { font-size: 15px; font-weight: 600; color: #374151; margin-bottom: 4px; }
/* Narrow panes give way column by column, least needed first; why a job
   waits and when it starts stay. */
@container (max-width: 1150px) { .djs-c-partition { display: none; } }
@container (max-width: 1080px) { .djs-c-runtime { display: none; } }
@container (max-width: 780px) { .djs-c-resources { display: none; } .djs-c-start { width: 160px; } }
@container (max-width: 560px) { .djs-c-status { width: 34px; } .djs-pill-text { display: none; }
                                th.djs-c-status { color: transparent; }
                                .djs-pill { padding: 5px; } .djs-c-start { width: 118px; } }
</style>"""

#: squeue's state codes, as a label and a colour.
_STATE_STYLE = {
    'R': ('Running', 'djs-running'), 'RUNNING': ('Running', 'djs-running'),
    'PD': ('Waiting', 'djs-pending'), 'PENDING': ('Waiting', 'djs-pending'),
    'CF': ('Starting', 'djs-running'), 'CONFIGURING': ('Starting', 'djs-running'),
    'CG': ('Finishing', 'djs-other'), 'COMPLETING': ('Finishing', 'djs-other'),
    'RQ': ('Requeued', 'djs-pending'), 'RH': ('Requeue held', 'djs-failed'),
    'S': ('Suspended', 'djs-other'), 'ST': ('Stopped', 'djs-other'),
    'PR': ('Preempted', 'djs-other'), 'CA': ('Cancelled', 'djs-other'),
    'F': ('Failed', 'djs-failed'), 'NF': ('Node failed', 'djs-failed'),
    'TO': ('Timed out', 'djs-failed'), 'OOM': ('Out of memory', 'djs-failed'),
    'BF': ('Boot failed', 'djs-failed'), 'SE': ('Special exit', 'djs-failed'),
}


def _fmt_duration(seconds) -> str:
    """'2 d', '1 d 12 h', '3 h 10 min', '45 min', '< 1 min'."""
    try:
        seconds = int(seconds)
    except (TypeError, ValueError):
        return ''
    if seconds < 60:
        return '< 1 min'
    days, rest = divmod(seconds // 60, 24 * 60)
    hours, mins = divmod(rest, 60)
    if days:
        return f'{days} d {hours} h' if hours else f'{days} d'
    if hours:
        return f'{hours} h {mins} min' if mins else f'{hours} h'
    return f'{mins} min'


def _slurm_seconds(text):
    """Seconds of a SLURM time string, None when there are none to read."""
    if not text:
        return None
    try:
        from .backend_slurm import SlurmJobBackend

        return SlurmJobBackend._time_limit_seconds(str(text))
    except Exception:
        return None


def _fmt_memory(text) -> str:
    """squeue's memory (240000M, 240G) as a person reads it: '240 GB'."""
    raw = str(text or '').strip()
    found = re.match(r'^(\d+(?:\.\d+)?)([KMGT]?)', raw)
    if not found:
        return raw
    megabytes = float(found.group(1)) * {'K': 0.001, 'M': 1, 'G': 1000, 'T': 1_000_000}[found.group(2) or 'M']
    if megabytes >= 1000:
        gigabytes = megabytes / 1000
        return f'{gigabytes:.0f} GB' if gigabytes >= 10 or gigabytes.is_integer() else f'{gigabytes:.1f} GB'
    return f'{megabytes:.0f} MB'


def _fmt_when(value) -> str:
    """'Wed 16 Sep, 21:20' for an ISO time, '' when there is none."""
    try:
        when = datetime.fromisoformat(str(value))
    except (TypeError, ValueError):
        return ''
    return f'{when:%a} {when.day} {when:%b}, {when:%H:%M}'


def _relative_start(start: str, now=None) -> str:
    """'in 2 d 7 h' for an ISO time, '' when there is none to read."""
    try:
        when = datetime.fromisoformat(str(start))
    except (TypeError, ValueError):
        return ''
    seconds = (when - (now or datetime.now())).total_seconds()
    if seconds <= 60:
        return 'any moment now'
    return f'in {_fmt_duration(seconds)}'


#: A waiting job's reason as a short label for its cell; the SLURM code and
#: the whole explanation are in the tooltip.
_REASON_LABELS = {
    'Priority': 'Priority', 'Resources': 'Resources', 'None': 'Scheduling',
    'Dependency': 'Dependency', 'DependencyNeverSatisfied': 'Dependency failed',
    'BeginTime': 'Begin time', 'JobHeldUser': 'Held by you', 'JobHeldAdmin': 'Held by admins',
    'ReqNodeNotAvail': 'Nodes unavailable', 'Reservation': 'Reservation',
    'PartitionTimeLimit': 'Time limit too long', 'PartitionNodeLimit': 'Too many nodes',
    'BadConstraints': 'Cannot be met', 'InvalidAccount': 'Invalid account', 'InvalidQOS': 'Invalid QOS',
    'NodeDown': 'Node down', 'Prolog': 'Starting', 'Cleaning': 'Cleaning up',
}

#: Reasons that do not resolve by waiting.
_STUCK_REASONS = {'DependencyNeverSatisfied', 'JobHeldUser', 'JobHeldAdmin', 'PartitionTimeLimit',
                  'PartitionNodeLimit', 'BadConstraints', 'InvalidAccount', 'InvalidQOS'}


def _reason_label(code: str) -> tuple:
    """(label, stuck) for a reason code."""
    from .backend_slurm import SlurmJobBackend

    head = str(code or '').split(',', 1)[0].strip()
    if SlurmJobBackend._is_env_hold(head):
        return 'Launch held', False
    if head.startswith('QOS') and 'Limit' in head:
        return 'Queue limit', False
    if head.startswith('Assoc') and 'Limit' in head:
        return 'Account limit', False
    return _REASON_LABELS.get(head, head or '—'), head in _STUCK_REASONS


def _fmt_when_short(value) -> str:
    """'16 Sep 21:20' for an ISO time, '' when there is none."""
    try:
        when = datetime.fromisoformat(str(value))
    except (TypeError, ValueError):
        return ''
    return f'{when.day} {when:%b} {when:%H:%M}'


def _is_state(job, *codes) -> bool:
    return str(job.status or '').upper() in codes


def _summary_html(jobs, now, standing=None) -> str:
    standing = standing or {}
    running = [j for j in jobs if _is_state(j, 'R', 'RUNNING')]
    waiting = [j for j in jobs if _is_state(j, 'PD', 'PENDING')]
    cpus = sum(int(str((j.extra or {}).get('cpus', '')).strip() or 0)
               for j in running if str((j.extra or {}).get('cpus', '')).strip().isdigit())
    estimates = []
    for job in waiting:
        try:
            estimates.append(datetime.fromisoformat(str((job.extra or {}).get('start_estimate', ''))))
        except ValueError:
            pass
    next_start = _relative_start(min(estimates).isoformat(), now) if estimates else '—'

    def stat(value, label, color='#111827', suffix='', tooltip=''):
        title = f' title="{_html.escape(tooltip, quote=True)}"' if tooltip else ''
        more = f'<span class="djs-stat-suffix"> {_html.escape(suffix)}</span>' if suffix else ''
        return (f'<div class="djs-stat"{title}><div class="djs-stat-value" style="color:{color};">'
                f'{_html.escape(str(value))}{more}</div><div class="djs-stat-label">{label}</div></div>')

    account = standing.get('account')
    of_account = f' (account {account})' if account else ''

    cpu_limit = standing.get('cpu_limit')
    cpu_notes = []
    if cpu_limit:
        cpu_notes.append(f'{cpus} of the {cpu_limit} CPUs you may use at once{of_account}. '
                         'Jobs that would go beyond wait with the reason "Account limit".')
    caps = [f'{standing["max_jobs"]} running' if standing.get('max_jobs') else '',
            f'{standing["max_submit"]} submitted' if standing.get('max_submit') else '']
    if any(caps):
        cpu_notes.append('At most ' + ' and '.join(cap for cap in caps if cap) + ' jobs.')
    cpu_tile = stat(cpus, 'CPUs in use', '#b26a00' if cpu_limit and cpus >= 0.85 * cpu_limit else '#111827',
                    f'/ {cpu_limit}' if cpu_limit else '', ' '.join(cpu_notes))

    fairshare_tile = ''
    fairshare = standing.get('fairshare')
    if fairshare is not None:
        color = '#2e7d32' if fairshare >= 0.5 else '#b26a00' if fairshare >= 0.2 else '#d32f2f'
        tooltip = (f'Fairshare {fairshare:.2f}{of_account}, on a scale from 0 to 1: the higher, the sooner '
                   'your waiting jobs start. 0.5 means you used exactly your share of the cluster.')
        shares, usage = standing.get('norm_shares'), standing.get('effective_usage')
        if shares and usage is not None:
            tooltip += f' Lately you used {usage / shares:.1f}× your share.'
        tooltip += ' It recovers by itself, as past usage counts less and less.'
        fairshare_tile = stat(f'{fairshare:.2f}', 'Fairshare · high is good', color, tooltip=tooltip)

    return ('<div class="djs-summary">'
            + stat(len(running), 'Running', '#2e7d32' if running else '#111827')
            + stat(len(waiting), 'Waiting', '#b26a00' if waiting else '#111827')
            + cpu_tile
            + stat(next_start, 'Next expected start')
            + fairshare_tile
            + '</div>')


def empty_queue_html(state_filter: str = 'all', summary_jobs=None, standing=None, now=None) -> str:
    """What the tab shows when there is nothing to list: the summary, fairshare
    included, stays in view."""
    title = {'running': 'No running jobs', 'pending': 'No waiting jobs'}.get(
        state_filter, 'No jobs in the queue')
    return (_DJS_CSS + '<div class="djs">'
            + _summary_html(list(summary_jobs or []), now or datetime.now(), standing)
            + '<div class="djs-empty">'
            f'<div class="djs-empty-title">{title}</div>'
            '<div class="djs-sub">Submitted jobs appear here: running ones with their progress, '
            'waiting ones with why they wait and when they should start.</div>'
            '</div></div>')


def build_slurm_job_table(jobs, dropdown_options, now=None, summary_jobs=None, standing=None) -> str:
    """The SLURM job table, one line per job: progress of running jobs, and for
    waiting jobs why they wait and when SLURM expects them to start. Whatever
    does not fit its column ends in an ellipsis and is whole in the tooltip."""
    from .backend_slurm import describe_pending_reason

    def esc(value) -> str:
        return _html.escape(str(value), quote=True)

    def cell(column, content, tooltip=''):
        title = f' title="{esc(tooltip)}"' if tooltip else ''
        return f'<td class="djs-c-{column}"{title}>{content}</td>'

    now = now or datetime.now()
    rows = []
    for job in jobs:
        extra = job.extra or {}
        code = str(job.status or '').upper()
        label, css = _STATE_STYLE.get(code, (code or 'Unknown', 'djs-other'))
        pending = code in ('PD', 'PENDING')
        running = code in ('R', 'RUNNING')

        status = cell('status', f'<span class="djs-pill {css}"><span class="djs-dot"></span>'
                                f'<span class="djs-pill-text">{esc(label)}</span></span>',
                      f'{label} (SLURM state {code})' if code else label)
        name = str(job.name)
        job_col = cell('job', f'<span class="djs-name">{esc(name)}</span> '
                              f'<span class="djs-muted djs-mono">#{esc(job.job_id)}</span>',
                       f'{name} (job {job.job_id})')

        partitions = [part.strip() for part in str(extra.get('partition', '')).split(',') if part.strip()]
        partition = cell('partition', ''.join(f'<span class="djs-tag">{esc(part)}</span>' for part in partitions),
                         'Starts in whichever partition frees first: ' + ', '.join(partitions)
                         if len(partitions) > 1 else ', '.join(partitions))

        cpus = str(extra.get('cpus', '')).strip()
        nodes = str(extra.get('nodes', '')).strip()
        pieces = [f'{cpus} CPUs' if cpus else '', _fmt_memory(extra.get('memory', '')),
                  f'{nodes} nodes' if nodes and nodes != '1' else '']
        resources_text = ' · '.join(piece for piece in pieces if piece) or '—'
        resources = cell('resources', esc(resources_text), resources_text)

        limit = _slurm_seconds(extra.get('time_limit') or job.time_limit)
        used = _slurm_seconds(extra.get('time_used'))
        if running and used is not None and limit:
            percent = max(0, min(100, round(100 * used / limit)))
            warn = ' djs-warn' if percent >= 85 else ''
            text = f'{_fmt_duration(used)} / {_fmt_duration(limit)}'
            runtime = cell('runtime', f'<div class="djs-progress"><div class="djs-bar{warn}">'
                                      f'<span style="width:{percent}%"></span></div>'
                                      f'<span class="djs-progress-text">{esc(text)}</span></div>',
                           f'{_fmt_duration(used)} of {_fmt_duration(limit)} used ({percent} %)')
        elif running and used is not None:
            runtime = cell('runtime', esc(_fmt_duration(used)), _fmt_duration(used))
        else:
            text = f'limit {_fmt_duration(limit)}' if limit else ''
            runtime = cell('runtime', f'<span class="djs-muted">{esc(text)}</span>', text)

        reason = str(extra.get('reason', ''))
        start = str(extra.get('start_estimate', '') or job.start_time or '')
        if pending:
            reason_code = str(extra.get('reason_code') or reason).strip().strip('()') or '—'
            short, stuck = _reason_label(reason_code)
            explanation = describe_pending_reason(reason_code)
            detail = explanation if explanation and explanation != reason_code else ''
            tooltip = f'{short}: {detail}' if detail else short
            if short != reason_code:
                tooltip += f' (SLURM reason {reason_code})'
            why = cell('why', f'<span class="djs-reason{" djs-stuck" if stuck else ""}">{esc(short)}</span>'
                              + (f' <span class="djs-muted">— {esc(detail)}</span>' if detail else ''),
                       tooltip)
            if _fmt_when(start):
                relative = _relative_start(start, now)
                start_col = cell('start', f'{esc(_fmt_when_short(start))} '
                                          f'<span class="djs-muted">· {esc(relative)}</span>',
                                 f'Expected start {_fmt_when(start)}, {relative}')
            else:
                start_col = cell('start', '<span class="djs-muted">not estimated yet</span>',
                                 'SLURM has not estimated a start for this job yet')
        else:
            why = cell('why', f'<span class="djs-mono">{esc(reason)}</span>' if reason else '',
                       f'Running on {reason}' if running and reason else reason)
            started = _fmt_when_short(start) if running else ''
            start_col = cell('start', f'<span class="djs-muted">since</span> {esc(started)}' if started else '',
                             f'Started {_fmt_when(start)}' if started else '')

        rows.append(f'<tr>{status}{job_col}{partition}{resources}{runtime}{why}{start_col}</tr>')
        dropdown_options.append((f'{job.job_id} - {job.name}', job.job_id))

    note = ''
    if any(_is_state(job, 'PD', 'PENDING') for job in jobs):
        note = ('<div class="djs-note"><b>About estimated starts:</b> SLURM assumes every running job '
                'uses its full time limit, so jobs usually start earlier than shown. Hover a cell to '
                'read what does not fit.</div>')

    header = ''.join(
        f'<th class="djs-c-{column}">{title}</th>'
        for column, title in (('status', 'Status'), ('job', 'Job'), ('partition', 'Partition'),
                              ('resources', 'Resources'), ('runtime', 'Runtime'),
                              ('why', 'Node / why waiting'), ('start', 'Start'))
    )
    return (
        _DJS_CSS + '<div class="djs">'
        + _summary_html(list(summary_jobs) if summary_jobs is not None else list(jobs), now, standing)
        + f'<div class="djs-wrap"><table class="djs-table"><thead><tr>{header}</tr></thead><tbody>'
        + ''.join(rows)
        + '</tbody></table></div>'
        + note
        + '</div>'
    )


def create_tab(ctx):
    """Create the Job Status tab.

    Returns ``(tab_widget, refs_dict)``.
    """
    # -- widgets --------------------------------------------------------
    # 'auto', not '100%': the widget's own 2px side margins on top of 100%
    # overflow the tab and give it a horizontal scrollbar.
    job_table_html = widgets.HTML(value='<i>Loading...</i>',
                                  layout=widgets.Layout(width='auto', min_width='0'))
    job_start_html = widgets.HTML(value='')

    job_dropdown = widgets.Dropdown(
        options=[],
        description='Select Job:',
        layout=widgets.Layout(width='400px'),
        style={'description_width': '80px'},
    )

    job_status_output = widgets.Output()

    refresh_button = widgets.Button(
        description='REFRESH', button_style='info', icon='refresh',
        layout=widgets.Layout(width='150px'),
    )

    rebuild_button = widgets.Button(
        description='REBUILD QUEUE', button_style='warning',
        layout=widgets.Layout(width='170px'),
    )

    cancel_button = widgets.Button(
        description='CANCEL JOB', button_style='danger', icon='times',
        layout=widgets.Layout(width='150px'),
    )

    # Filter toggle: show all jobs, only running (R), or only pending (PD).
    # Backed by ``list_jobs()`` which already scopes to the current user
    # via ``squeue -u $USER`` — so this filters the user's OWN jobs by
    # state, never adds other users' jobs.
    filter_toggle = widgets.ToggleButtons(
        options=[('All', 'all'), ('Running', 'running'), ('Pending', 'pending')],
        value='all',
        description='Filter:',
        button_style='',
        layout=widgets.Layout(margin='0 0 0 10px'),
        style={'description_width': '50px', 'button_width': '90px'},
    )

    # -- state ----------------------------------------------------------
    state = {
        'job_data': [],
        'cancel_armed_job_id': None,
        'state_filter': 'all',
    }

    def _reset_cancel_arm():
        state['cancel_armed_job_id'] = None
        cancel_button.description = 'CANCEL JOB'

    # -- handlers -------------------------------------------------------
    def _pid_elapsed_seconds(pid):
        """Best-effort elapsed seconds from /proc; fall back to None.

        The start time is read through ``proc_identity`` rather than by
        splitting the whole stat line: field 2 is the command name in
        brackets and may contain spaces, so index 21 of a whole-line
        split lands on a different field for those processes. Measured
        on a binary named ``a b``: 0 instead of 26055418, which is a job
        reported as having run for the machine's entire uptime.
        """
        from delfin.agent.proc_identity import start_ticks

        try:
            # /proc/uptime gives system uptime in seconds
            with open('/proc/uptime', 'r') as f:
                uptime = float(f.read().split()[0])
            ticks = start_ticks(pid)
            if ticks is None:
                return None
            clk_tck = os.sysconf(os.sysconf_names['SC_CLK_TCK'])
            elapsed = max(0, int(uptime - (ticks / clk_tck)))
            return elapsed
        except Exception:
            return None

    def _parse_control_resources(text):
        pal = None
        maxcore = None
        time_limit = None
        m = re.search(r'^\s*PAL\s*=\s*(\d+)', text, flags=re.MULTILINE)
        if m:
            pal = int(m.group(1))
        m = re.search(r'^\s*maxcore\s*=\s*(\d+)', text, flags=re.MULTILINE)
        if m:
            maxcore = int(m.group(1))
        m = re.search(r'^\s*job_timeout_hours\s*=\s*(\d+)', text, flags=re.MULTILINE)
        if m:
            time_limit = f"{int(m.group(1))}:00:00"
        return pal, maxcore, time_limit

    def _parse_orca_resources(text):
        pal, maxcore = parse_inp_resources(text)
        # tab_job_status historically also accepted bare "maxcore N" without "%".
        if maxcore is None:
            m = re.search(r'^\s*maxcore\s*=?\s*(\d+)', text, flags=re.IGNORECASE | re.MULTILINE)
            if m:
                maxcore = int(m.group(1))
        return pal, maxcore

    def _extract_resources_for_root(root_dir):
        """Try to extract PAL/maxcore/time_limit from CONTROL.txt or ORCA inputs."""
        pal = None
        maxcore = None
        time_limit = None

        candidates = [Path(root_dir) / 'CONTROL.txt', Path(root_dir) / 'builder' / 'CONTROL.txt']
        for p in candidates:
            if p.exists():
                try:
                    text = p.read_text()
                    pal, maxcore, time_limit = _parse_control_resources(text)
                    return pal, maxcore, time_limit
                except Exception:
                    pass

        # Fallback: parse an ORCA input in builder or root
        inp_candidates = []
        for d in [Path(root_dir) / 'builder', Path(root_dir)]:
            if d.exists():
                inp_candidates.extend(sorted(d.glob('*.inp')))
        if inp_candidates:
            try:
                text = inp_candidates[-1].read_text()
                pal, maxcore = _parse_orca_resources(text)
            except Exception:
                pass
        return pal, maxcore, time_limit

    def _detect_running_processes():
        """Fallback: detect running ORCA/DELFIN processes when queue is empty."""
        user = os.environ.get('USER') or os.getlogin()
        cmd = ['ps', '-u', user, '-o', 'pid=,pgid=,etimes=,command=']
        try:
            result = subprocess.run(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False
            )
            if result.returncode != 0:
                return []
        except Exception:
            return []

        orca_markers = []
        resolved_orca_base = str(getattr(ctx, 'orca_base', '') or '').strip()
        if resolved_orca_base:
            orca_markers.extend(
                [
                    resolved_orca_base,
                    os.path.join(resolved_orca_base, 'orca'),
                ]
            )
        orca_markers.extend(
            [
                ' orca ',
                '/orca ',
                ' orca"',
                '/orca"',
            ]
        )
        patterns = tuple(dict.fromkeys(orca_markers + [
            'orca_numfreq',
            'orca_esd',
            'delfin.build_up_complex',
            'delfin.guppy_sampling',
            'delfin-build',
        ]))
        detected = []
        calc_root = str(ctx.calc_dir) if getattr(ctx, 'calc_dir', None) else ''
        by_dir = {}
        for line in result.stdout.splitlines():
            line = line.strip()
            if not line:
                continue
            parts = line.split(None, 3)
            if len(parts) < 4:
                continue
            pid, pgid, etimes, command = parts
            if not any(p in command for p in patterns):
                continue

            job_dir = ''
            try:
                job_dir = os.readlink(f'/proc/{pid}/cwd')
            except Exception:
                pass

            # Try to extract an input file or working path from the command
            name = ''
            m = re.search(r'(/[^\\s]+\\.(?:inp|cisinp\\.tmp|goat|xyz))', command)
            if m:
                name = Path(m.group(1)).name
                if not job_dir:
                    try:
                        job_dir = str(Path(m.group(1)).parent)
                    except Exception:
                        pass

            # If we still don't have a directory, try to extract any absolute path
            if not job_dir:
                m2 = re.search(r'(/[^\\s]+)', command)
                if m2:
                    p = Path(m2.group(1))
                    if p.exists():
                        job_dir = str(p.parent)

            if calc_root:
                if not job_dir:
                    # No directory and we only want calc jobs
                    continue
                if not str(job_dir).startswith(calc_root):
                    # Only show calc jobs to reduce noise
                    continue

            # Derive a human-friendly job name from calc directory
            if calc_root and str(job_dir or '').startswith(calc_root):
                try:
                    rel = Path(job_dir).relative_to(calc_root)
                    # e.g., Ir_8/builder -> Ir_8
                    if len(rel.parts) >= 1:
                        name = rel.parts[0]
                except Exception:
                    pass

            if not name and job_dir:
                name = Path(job_dir).name
            if not name:
                name = 'detected_process'

            # Collapse multiple worker processes by job directory
            key = job_dir
            et = _pid_elapsed_seconds(pid)
            if et is None:
                try:
                    et = int(etimes)
                except Exception:
                    et = 0
            if key not in by_dir:
                by_dir[key] = {
                    'pid': pid,
                    'pgid': pgid,
                    'etimes': et,
                    'name': name,
                    'job_dir': job_dir,
                    'command': command,
                    'count': 1,
                }
            else:
                by_dir[key]['count'] += 1
                # Keep the longest-running process as the representative
                if et > by_dir[key]['etimes']:
                    by_dir[key].update({
                        'pid': pid,
                        'pgid': pgid,
                        'etimes': et,
                        'command': command,
                    })

        for job_dir, info in by_dir.items():
            et = info.get('etimes', 0)
            try:
                hours = et // 3600
                mins = (et % 3600) // 60
                secs = et % 60
                etime_display = f'{hours:02d}:{mins:02d}:{secs:02d}'
            except Exception:
                etime_display = str(et)

            detected.append(JobInfo(
                job_id=f'pgid:{info.get("pgid", info.get("pid"))}',
                name=info.get('name', 'detected_process'),
                mode='detected',
                status='RUNNING',
                submit_time=f'elapsed {etime_display}',
                start_time='',
                pal=0,
                maxcore=0,
                time_limit='-',
                job_dir=info.get('job_dir', ''),
                extra={
                    'command': info.get('command', ''),
                    'proc_count': info.get('count', 1),
                    'pid': info.get('pid'),
                    'pgid': info.get('pgid'),
                },
            ))
        return detected

    def _rebuild_queue_from_detected():
        """Rebuild local queue JSON from detected running processes."""
        if not hasattr(ctx.backend, 'jobs_file'):
            return False, 'Queue rebuild only supported for local backend.'
        detected = _detect_running_processes()
        if not detected:
            return False, 'No running processes detected to rebuild.'

        jobs_file = ctx.backend.jobs_file
        try:
            if jobs_file.exists():
                import json
                data = json.loads(jobs_file.read_text())
            else:
                data = {}
        except Exception:
            data = {}

        # Reset queue to avoid duplicates or stale FAILED entries
        data = {'_next_job_id': 1001, 'jobs': []}

        now = datetime.now()
        seen_dirs = set()
        for job in detected:
            # Normalize to calc root directory for clarity
            job_dir = job.job_dir or ''
            if getattr(ctx, 'calc_dir', None):
                try:
                    rel = Path(job_dir).relative_to(ctx.calc_dir)
                    root_dir = str(Path(ctx.calc_dir) / rel.parts[0])
                except Exception:
                    root_dir = job_dir
            else:
                root_dir = job_dir

            if root_dir in seen_dirs:
                continue
            seen_dirs.add(root_dir)

            pal, maxcore, time_limit = _extract_resources_for_root(root_dir)

            job_id = data['_next_job_id']
            data['_next_job_id'] = job_id + 1
            # Try to back-calculate start_time from elapsed seconds
            start_time = None
            elapsed = 0
            m = re.search(r'elapsed\\s+(\\d+):(\\d+):(\\d+)', job.submit_time or '')
            if m:
                elapsed = int(m.group(1)) * 3600 + int(m.group(2)) * 60 + int(m.group(3))
                start_time = (now - timedelta(seconds=elapsed)).isoformat()

            # Use representative PID to keep status tracking alive
            rep_pid = None
            if isinstance(job.extra, dict):
                rep_pid = job.extra.get('pid')
            if rep_pid is None:
                # Fallback to PGID if PID isn't available
                rep_pid = int(str(job.job_id).split(':')[-1]) if str(job.job_id).startswith('pgid:') else None

            data['jobs'].append({
                'job_id': job_id,
                'name': Path(root_dir).name if root_dir else job.name,
                'mode': 'detected',
                'pid': int(rep_pid) if rep_pid else None,
                'pgid': int(str(job.job_id).split(':')[-1]) if str(job.job_id).startswith('pgid:') else None,
                'status': 'RUNNING',
                'submit_time': now.isoformat(),
                'start_time': start_time,
                'job_dir': root_dir or job.job_dir,
                'time_limit': time_limit or '-',
                'pal': pal or 0,
                'maxcore': maxcore or 0,
                'inp_file': None,
                'override': None,
                'build_mult': None,
                'log_file': None,
            })

        try:
            jobs_file.write_text(__import__('json').dumps(data, indent=2, default=str))
            return True, f'Rebuilt queue with {len(detected)} running job(s).'
        except Exception as e:
            return False, f'Failed to write queue: {e}'

    def refresh_job_list(button=None):
        _reset_cancel_arm()
        with job_status_output:
            clear_output()

        try:
            # The user asked for the queue, so bypass the rate limit.
            all_jobs = ctx.backend.list_jobs(force=True)
            state['job_data'] = all_jobs
            # Fairshare and account limits; the backend throttles these itself.
            try:
                standing = getattr(ctx.backend, 'account_standing', lambda: None)()
            except Exception:
                standing = None
            state['showing_detected'] = False

            # Apply state filter (all / running / pending).
            f = state.get('state_filter', 'all')
            if f == 'running':
                jobs = [
                    j for j in all_jobs
                    if (j.status or '').upper() in ('RUNNING', 'R')
                ]
            elif f == 'pending':
                jobs = [
                    j for j in all_jobs
                    if (j.status or '').upper() in ('PENDING', 'PD')
                ]
            else:
                jobs = list(all_jobs)

            if not jobs:
                # Fallback: detect running processes if local queue is empty
                is_local = not ctx.backend.supports_turbomole  # heuristic
                detected = _detect_running_processes() if is_local else []
                state['showing_detected'] = bool(detected)
                if detected:
                    jobs = detected
                    state['job_data'] = jobs
                    job_start_html.value = (
                        '<div style="margin-top:10px;padding:10px;'
                        'background-color:#e3f2fd;border:1px solid #90caf9;border-radius:4px;">'
                        '<b>Note:</b> Local queue is empty. Showing running processes detected '
                        'on this machine (not from the job queue).</div>'
                    )
                else:
                    job_table_html.value = empty_queue_html(
                        state.get('state_filter', 'all'), summary_jobs=all_jobs, standing=standing)
                    job_start_html.value = ''
                    job_dropdown.options = []
                    return

            dropdown_options = []

            # Build table based on backend type
            is_local = not ctx.backend.supports_turbomole  # heuristic

            table_html = JOB_TABLE_CSS
            if is_local:
                table_html += _build_local_table(jobs, dropdown_options)
            else:
                table_html += _build_slurm_table(jobs, dropdown_options, standing)

            job_table_html.value = table_html

            # The SLURM table carries its own note on start estimates; this
            # box is only for the local backend's detected-process hint.
            if not state.get('showing_detected'):
                job_start_html.value = ''

            if dropdown_options:
                job_dropdown.options = dropdown_options
                job_dropdown.value = None
            else:
                job_dropdown.options = []

            with job_status_output:
                clear_output()
                # Counts use the FILTERED set; show total in parens when
                # a filter is active so the user knows how many were
                # hidden.
                running = sum(
                    1 for j in jobs
                    if (j.status or '').upper() in ('RUNNING', 'R')
                )
                pending = sum(
                    1 for j in jobs
                    if (j.status or '').upper() in ('PENDING', 'PD')
                )
                parts = [f'{len(jobs)} job(s) shown']
                if running:
                    parts.append(f'{running} running')
                if pending:
                    parts.append(f'{pending} in queue')
                if state.get('state_filter', 'all') != 'all':
                    parts.append(f'(of {len(all_jobs)} total)')
                if is_local:
                    # SLURM shows the same counts in the summary above its table.
                    print(', '.join(parts))

        except Exception as e:
            with job_status_output:
                clear_output(wait=True)
                print(f'Error loading jobs: {e}')

    def _build_local_table(jobs, dropdown_options):
        """Build an HTML table for the local backend."""
        tbl = (
            '<table class="job-table"><tr>'
            '<th>JOB ID</th><th>NAME</th><th>MODE</th><th>STATUS</th>'
            '<th>QUEUE</th><th>SUBMITTED</th><th>PAL</th>'
            '<th>TIME LIMIT</th><th>DIRECTORY</th></tr>'
        )

        pending_sorted = [j for j in sorted(jobs, key=lambda j: j.job_id) if j.status == 'PENDING']
        queue_positions = {j.job_id: idx + 1 for idx, j in enumerate(pending_sorted)}

        for job in sorted(jobs, key=lambda j: j.job_id, reverse=True):
            try:
                dt = datetime.fromisoformat(job.submit_time)
                submit_display = dt.strftime('%Y-%m-%d %H:%M')
            except Exception:
                submit_display = job.submit_time

            try:
                dir_display = str(Path(job.job_dir).name)
            except Exception:
                dir_display = job.job_dir

            color = STATUS_COLORS.get(job.status, '#333')
            tbl += (
                f'<tr>'
                f'<td><b>{job.job_id}</b></td>'
                f'<td>{_html.escape(str(job.name))}</td>'
                f'<td>{_html.escape(str(job.mode))}</td>'
                f'<td><b style="color:{color};">{job.status}</b></td>'
                f'<td>{queue_positions.get(job.job_id, "-")}</td>'
                f'<td>{submit_display}</td>'
                f'<td>{job.pal}</td>'
                f'<td>{job.time_limit}</td>'
                f'<td title="{_html.escape(str(job.job_dir))}">{_html.escape(dir_display)}</td>'
                f'</tr>'
            )

            if job.status in ('RUNNING', 'PENDING'):
                dropdown_options.append((f'{job.job_id} - {job.name}', job.job_id))

        tbl += '</table>'
        return tbl

    def _build_slurm_table(jobs, dropdown_options, standing=None):
        """Build an HTML table for the SLURM backend."""
        return build_slurm_job_table(jobs, dropdown_options, summary_jobs=state.get('job_data'),
                                     standing=standing)

    def cancel_selected_job(button):
        with job_status_output:
            clear_output()
            if not job_dropdown.value:
                _reset_cancel_arm()
                print('No job selected.')
                return
            job_id = job_dropdown.value
            if state.get('cancel_armed_job_id') != job_id:
                state['cancel_armed_job_id'] = job_id
                cancel_button.description = 'CONFIRM CANCEL'
                print(f'Click CANCEL JOB again to confirm cancellation of job {job_id}.')
                return
            _reset_cancel_arm()
            success, msg = ctx.backend.cancel_job(job_id)
            print(msg)
            if success:
                refresh_job_list()

    def rebuild_queue(button):
        with job_status_output:
            clear_output()
            _reset_cancel_arm()
            ok, msg = _rebuild_queue_from_detected()
            print(msg)
            if ok:
                refresh_job_list()

    def _on_job_selection_change(change):
        if change.get('name') == 'value':
            _reset_cancel_arm()

    # -- wiring ---------------------------------------------------------
    def _on_filter_change(change):
        if change.get('name') != 'value':
            return
        state['state_filter'] = change['new']
        refresh_job_list()

    refresh_button.on_click(refresh_job_list)
    cancel_button.on_click(cancel_selected_job)
    job_dropdown.observe(_on_job_selection_change, names='value')
    filter_toggle.observe(_on_filter_change, names='value')
    # Only enable rebuild in local backend
    is_local_backend = not ctx.backend.supports_turbomole  # heuristic
    if is_local_backend:
        rebuild_button.on_click(rebuild_queue)
    refresh_job_list()

    # -- layout ---------------------------------------------------------
    header_html = widgets.HTML(
        '<div style="display:flex;flex-wrap:wrap;align-items:baseline;gap:4px 12px;margin:0 0 2px 0;">'
        '<span style="font-size:18px;font-weight:600;color:#111827;">Job Status</span>'
        '<span style="color:#6b7280;font-size:12.5px;">Your own jobs. Waiting jobs show why they wait '
        'and when they are expected to start. Select a job below to cancel it.</span>'
        '</div>'
    )
    tab_widget = widgets.VBox([
        header_html,
        widgets.HBox(
            [filter_toggle, refresh_button],
            layout=widgets.Layout(margin='6px 0 10px 0', align_items='center',
                                  justify_content='space-between', width='100%'),
        ),
        job_table_html, job_start_html,
        widgets.HBox(
            [job_dropdown, cancel_button] + ([rebuild_button] if is_local_backend else []),
            layout=widgets.Layout(margin='12px 0 4px 0', align_items='center'),
        ),
        job_status_output,
    ], layout=widgets.Layout(padding='10px'))

    return tab_widget, {'refresh_job_list': refresh_job_list}
