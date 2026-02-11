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


def create_tab(ctx):
    """Create the Job Status tab.

    Returns ``(tab_widget, refs_dict)``.
    """
    # -- widgets --------------------------------------------------------
    job_table_html = widgets.HTML(value='<i>Loading...</i>')
    job_start_html = widgets.HTML(value='')

    job_dropdown = widgets.Dropdown(
        options=[],
        description='Select Job:',
        layout=widgets.Layout(width='400px'),
        style={'description_width': '80px'},
    )

    job_status_output = widgets.Output()

    refresh_button = widgets.Button(
        description='REFRESH', button_style='info',
        layout=widgets.Layout(width='150px'),
    )

    rebuild_button = widgets.Button(
        description='REBUILD QUEUE', button_style='warning',
        layout=widgets.Layout(width='170px'),
    )

    cancel_button = widgets.Button(
        description='CANCEL JOB', button_style='danger',
        layout=widgets.Layout(width='150px'),
    )

    # -- state ----------------------------------------------------------
    state = {'job_data': []}

    # -- handlers -------------------------------------------------------
    def _pid_elapsed_seconds(pid):
        """Best-effort elapsed seconds from /proc; fall back to None."""
        try:
            # /proc/uptime gives system uptime in seconds
            with open('/proc/uptime', 'r') as f:
                uptime = float(f.read().split()[0])
            with open(f'/proc/{pid}/stat', 'r') as f:
                stat = f.read().split()
            # starttime is field 22 (index 21) in clock ticks
            start_ticks = int(stat[21])
            clk_tck = os.sysconf(os.sysconf_names['SC_CLK_TCK'])
            elapsed = max(0, int(uptime - (start_ticks / clk_tck)))
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
        pal = None
        maxcore = None
        # %pal nprocs N
        m = re.search(r'^\s*%pal\b.*?nprocs\s+(\d+)', text, flags=re.IGNORECASE | re.MULTILINE)
        if m:
            pal = int(m.group(1))
        # %maxcore N or maxcore N
        m = re.search(r'^\s*%maxcore\s+(\d+)', text, flags=re.IGNORECASE | re.MULTILINE)
        if m:
            maxcore = int(m.group(1))
        else:
            m = re.search(r'^\s*maxcore\s+(\d+)', text, flags=re.IGNORECASE | re.MULTILINE)
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

        patterns = (
            '/opt/orca',
            'orca_numfreq',
            'orca_esd',
            'delfin.build_up_complex',
            'delfin-build',
        )
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
                if not job_dir.startswith(calc_root):
                    # Only show calc jobs to reduce noise
                    continue

            # Derive a human-friendly job name from calc directory
            if calc_root and job_dir.startswith(calc_root):
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
        with job_status_output:
            clear_output()

        try:
            jobs = ctx.backend.list_jobs()
            state['job_data'] = jobs

            if not jobs:
                # Fallback: detect running processes if local queue is empty
                is_local = not ctx.backend.supports_turbomole  # heuristic
                detected = _detect_running_processes() if is_local else []
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
                    job_table_html.value = '<p><i>No jobs found.</i></p>'
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
                table_html += _build_slurm_table(jobs, dropdown_options)

            job_table_html.value = table_html

            # Pending start times
            pending_times = ctx.backend.get_pending_start_times()
            if pending_times and isinstance(pending_times[0], dict) and 'start' in pending_times[0]:
                start_html = (
                    '<div style="margin-top:10px;padding:10px;'
                    'background-color:#fff3cd;border:1px solid #ffc107;border-radius:4px;">'
                    '<b>Estimated Start Times (pending jobs):</b><br>'
                    '<table style="font-family:monospace;font-size:12px;margin-top:5px;">'
                )
                for pj in pending_times:
                    start_html += (
                        f"<tr><td><b>{pj['id']}</b></td>"
                        f"<td style='padding-left:15px;'>{pj['name']}</td>"
                        f"<td style='padding-left:15px;'>{pj['start']}</td></tr>"
                    )
                start_html += '</table></div>'
                job_start_html.value = start_html
            else:
                job_start_html.value = ''

            if dropdown_options:
                job_dropdown.options = dropdown_options
                job_dropdown.value = None
            else:
                job_dropdown.options = []

            with job_status_output:
                clear_output()
                running = sum(1 for j in jobs if j.status == 'RUNNING')
                pending = sum(1 for j in jobs if j.status == 'PENDING')
                parts = [f'{len(jobs)} job(s) shown']
                if running:
                    parts.append(f'{running} running')
                if pending:
                    parts.append(f'{pending} in queue')
                print(', '.join(parts))

        except Exception as e:
            with job_status_output:
                clear_output()
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

    def _build_slurm_table(jobs, dropdown_options):
        """Build an HTML table for the SLURM backend."""
        tbl = (
            '<table class="job-table"><tr>'
            '<th>JOBID</th><th>PARTITION</th><th>NAME</th><th>USER</th>'
            '<th>ST</th><th>TIME</th><th>NODES</th><th>NODELIST/REASON</th>'
            '</tr>'
        )

        for job in jobs:
            extra = job.extra or {}
            parts = extra.get('raw_line', '').split(None, 7)
            while len(parts) < 8:
                parts.append('')

            tbl += (
                f'<tr>'
                f'<td><b>{parts[0]}</b></td>'
                f'<td>{parts[1]}</td>'
                f'<td>{parts[2]}</td>'
                f'<td>{parts[3]}</td>'
                f'<td>{parts[4]}</td>'
                f'<td>{parts[5]}</td>'
                f'<td>{parts[6]}</td>'
                f'<td>{parts[7]}</td>'
                f'</tr>'
            )
            dropdown_options.append((f'{job.job_id} - {job.name}', job.job_id))

        tbl += '</table>'
        return tbl

    def cancel_selected_job(button):
        with job_status_output:
            clear_output()
            if not job_dropdown.value:
                print('No job selected.')
                return
            job_id = job_dropdown.value
            success, msg = ctx.backend.cancel_job(job_id)
            print(msg)
            if success:
                refresh_job_list()

    def rebuild_queue(button):
        with job_status_output:
            clear_output()
            ok, msg = _rebuild_queue_from_detected()
            print(msg)
            if ok:
                refresh_job_list()

    # -- wiring ---------------------------------------------------------
    refresh_button.on_click(refresh_job_list)
    cancel_button.on_click(cancel_selected_job)
    # Only enable rebuild in local backend
    is_local_backend = not ctx.backend.supports_turbomole  # heuristic
    if is_local_backend:
        rebuild_button.on_click(rebuild_queue)
    refresh_job_list()

    # -- layout ---------------------------------------------------------
    tab_widget = widgets.VBox([
        widgets.HTML('<h3>Job Status</h3>'),
        widgets.HTML('<p>Select a job and click CANCEL to cancel it. Click REFRESH to update.</p>'),
        job_table_html, job_start_html,
        widgets.HBox(
            [job_dropdown, cancel_button, refresh_button] + ([rebuild_button] if is_local_backend else []),
            layout=widgets.Layout(margin='10px 0', align_items='center'),
        ),
        job_status_output,
    ], layout=widgets.Layout(padding='10px'))

    return tab_widget, {'refresh_job_list': refresh_job_list}
