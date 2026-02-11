"""Job Status tab: view and cancel running/pending jobs."""

import html as _html
from datetime import datetime
from pathlib import Path

import ipywidgets as widgets
from IPython.display import clear_output

from .constants import STATUS_COLORS, JOB_TABLE_CSS


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

    cancel_button = widgets.Button(
        description='CANCEL JOB', button_style='danger',
        layout=widgets.Layout(width='150px'),
    )

    # -- state ----------------------------------------------------------
    state = {'job_data': []}

    # -- handlers -------------------------------------------------------
    def refresh_job_list(button=None):
        with job_status_output:
            clear_output()

        try:
            jobs = ctx.backend.list_jobs()
            state['job_data'] = jobs

            if not jobs:
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

    # -- wiring ---------------------------------------------------------
    refresh_button.on_click(refresh_job_list)
    cancel_button.on_click(cancel_selected_job)
    refresh_job_list()

    # -- layout ---------------------------------------------------------
    tab_widget = widgets.VBox([
        widgets.HTML('<h3>Job Status</h3>'),
        widgets.HTML('<p>Select a job and click CANCEL to cancel it. Click REFRESH to update.</p>'),
        job_table_html, job_start_html,
        widgets.HBox(
            [job_dropdown, cancel_button, refresh_button],
            layout=widgets.Layout(margin='10px 0', align_items='center'),
        ),
        job_status_output,
    ], layout=widgets.Layout(padding='10px'))

    return tab_widget, {'refresh_job_list': refresh_job_list}
