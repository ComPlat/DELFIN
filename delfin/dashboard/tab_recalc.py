"""Recalc tab: resubmit an existing DELFIN job with edited CONTROL.txt."""

import ipywidgets as widgets
from IPython.display import clear_output

from delfin.config import validate_control_text

from .helpers import resolve_time_limit, create_time_limit_widgets
from .input_processing import parse_resource_settings


def create_tab(ctx):
    """Create the Recalc tab.

    Returns ``(tab_widget, refs_dict)``.
    """
    # -- widgets --------------------------------------------------------
    recalc_folder_dropdown = widgets.Dropdown(
        options=[],
        description='Job Folder:',
        layout=widgets.Layout(width='500px'),
        style={'description_width': 'initial'},
    )

    recalc_time_toggle, recalc_custom_time = create_time_limit_widgets()
    recalc_time_toggle.value = '24h'

    recalc_control_widget = widgets.Textarea(
        value='',
        description='CONTROL.txt:',
        layout=widgets.Layout(width='500px', height='400px'),
        style={'description_width': 'initial'},
    )

    recalc_button = widgets.Button(
        description='SUBMIT RECALC', button_style='warning',
        layout=widgets.Layout(width='150px'),
    )

    recalc_refresh_btn = widgets.Button(
        description='REFRESH FOLDERS', button_style='info',
        layout=widgets.Layout(width='150px'),
    )

    recalc_output = widgets.Output()

    # -- handlers -------------------------------------------------------
    def refresh_recalc_folders():
        calc_dir = ctx.calc_dir
        if not calc_dir.exists():
            recalc_folder_dropdown.options = []
            recalc_folder_dropdown.value = None
            return
        folders = sorted([p.name for p in calc_dir.iterdir() if p.is_dir()])
        recalc_folder_dropdown.options = folders
        recalc_folder_dropdown.value = None

    def load_recalc_control(change=None):
        with recalc_output:
            clear_output()
        folder = recalc_folder_dropdown.value
        if not folder:
            recalc_control_widget.value = ''
            return
        control_path = ctx.calc_dir / folder / 'CONTROL.txt'
        if control_path.exists():
            recalc_control_widget.value = control_path.read_text()
        else:
            recalc_control_widget.value = ''
            with recalc_output:
                print(f'CONTROL.txt not found in {control_path.parent}')

    def handle_recalc(button):
        with recalc_output:
            clear_output()

            folder = recalc_folder_dropdown.value
            if not folder:
                print('Error: Please select a job folder.')
                return

            job_dir = ctx.calc_dir / folder
            if not job_dir.exists():
                print(f'Error: Job folder does not exist: {job_dir}')
                return

            control_errors = validate_control_text(recalc_control_widget.value)
            if control_errors:
                print('CONTROL.txt validation failed:')
                for err in control_errors:
                    print(f'- {err}')
                return

            control_path = job_dir / 'CONTROL.txt'
            control_path.write_text(recalc_control_widget.value)

            pal, maxcore = parse_resource_settings(recalc_control_widget.value)
            if pal is None or maxcore is None:
                print('Error: PAL/maxcore not found in CONTROL.txt')
                return

            time_limit = resolve_time_limit(recalc_time_toggle, recalc_custom_time, '24:00:00')

            result = ctx.backend.submit_delfin(
                job_dir=job_dir,
                job_name=job_dir.name,
                mode='delfin-recalc',
                time_limit=time_limit,
                pal=pal,
                maxcore=maxcore,
            )

            if result.returncode == 0:
                job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                print('Recalc job successfully submitted!')
                print(f'Job ID: {job_id}')
                print(f'Directory: {job_dir}')
            else:
                print('Error submitting recalc job:')
                print(result.stderr or result.stdout)

    # -- wiring ---------------------------------------------------------
    recalc_folder_dropdown.observe(load_recalc_control, names='value')
    recalc_button.on_click(handle_recalc)
    recalc_refresh_btn.on_click(lambda b: refresh_recalc_folders())
    refresh_recalc_folders()

    # -- layout ---------------------------------------------------------
    spacer = widgets.Label(value='', layout=widgets.Layout(height='10px'))

    tab_widget = widgets.VBox([
        widgets.HTML('<h3>Recalc DELFIN Job</h3>'),
        widgets.HTML('<p>Select a job folder, edit CONTROL.txt, and resubmit.</p>'),
        widgets.HBox([recalc_folder_dropdown, recalc_refresh_btn]),
        recalc_time_toggle,
        recalc_custom_time,
        recalc_control_widget, spacer,
        recalc_button, recalc_output,
    ], layout=widgets.Layout(padding='10px'))

    return tab_widget, {'refresh_recalc_folders': refresh_recalc_folders}
