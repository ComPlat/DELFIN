"""What each dashboard tab is for, as plain data.

Moved out of ``delfin/dashboard/tab_help.py`` on 2026-11-19 so that the
help lookup can be answered without the dashboard. Reading it from the
tab module imported ipywidgets, IPython and traitlets -- a terminal
session, a headless run or an install without the dashboard extras
paid for the whole widget stack to answer "what does this key do", or
could not answer at all. The content is unchanged; only its address is.

This file is data. It must not import anything from delfin.dashboard.
"""

from __future__ import annotations

TABS_HELP: tuple[tuple[str, str, str, str, str], ...] = (
    ("Submit Job", "Configure and start a DELFIN calculation: structure, method, solvent, and the CONTROL.txt of the run.",
     "<ul><li><b>Job Name</b> (Text, tab_submit.py:100): folder name of the run, e.g. Fe_Complex_Ox.</li><li><b>Coordinates</b> (Textarea): paste XYZ coordinates (header auto-removed) or a SMILES string; converter placeholders QUICK|NORMAL|MANTA|ARCHITECTOR|MOLSIMPLIFY|MACE|RXEMBED|MOLASSEMBLER|METALLOGEN (tab_submit.py:95-103).</li><li><b>Converter</b> dropdown [QUICK, NORMAL, MANTA, ARCHITECTOR, MOLSIMPLIFY, MACE, RXEMBED, MOLASSEMBLER, METALLOGEN], default QUICK (tab_submit.py:229).</li><li><b>Charge</b> (IntText, default 0) and <b>Solvent</b> (Text, placeholder 'e.g. water') (tab_submit.py:219-223).</li><li><b>PAL</b> (IntText, default 12) and <b>MaxCore</b> (IntText, default 600) and <b>Time Limit</b> (default 24:00:00) (tab_submit.py:233-241).</li><li><b>BUILD COMPLEX</b> button: opens the metal-complex builder for the pasted SMILES (tab_submit.py:124).</li><li><b>SUBMIT JOB</b> writes the CONTROL.txt and starts the calculation (tab_submit.py:205).</li><li><b>VALIDATE CONTROL</b> checks the CONTROL.txt without submitting; result appears in the validate output (tab_submit.py:209).</li><li>Special submit buttons: <b>SUBMIT ONLY GOAT</b> (Species delta, default -2), <b>SUBMIT DELFIN + CO2</b>, <b>SUBMIT GUPPY</b> (PAL/GOAT/time widgets), <b>SUBMIT FUKUI</b> (Skip cube generation checkbox, PAL 4, MaxCore 3000, Time 04:00:00) (tab_submit.py:245-302).</li></ul>",
     "<b>Typical flow</b>: enter job name &rarr; paste XYZ or SMILES (pick converter; for SMILES a 3D build runs) &rarr; set charge/solvent/PAL/MaxCore/Time Limit &rarr; optionally edit the CONTROL.txt preview &rarr; VALIDATE CONTROL, fix any reported problem &rarr; SUBMIT JOB. The run lands in calc/&lt;job-name&gt;.",
     "<b>Interactions</b>: BUILD COMPLEX hands the structure to the metal-complex builder; a submitted job then appears under Job Status and, once finished, under Calculations. The CONTROL.txt you can edit here is the same file Recalc reads."),
        ("Recalc", "Resume or re-run parts of an existing calculation from an earlier job folder.",
     "<ul><li><b>Job Folder</b> (Text/Dropdown, tab_recalc.py:52): the calc/ folder of the run to resume.</li><li><b>CONTROL.txt</b> (Text, tab_recalc.py:61-62): the CONTROL content for the new attempt (pre-filled from the folder's control).</li><li><b>SUBMIT RECALC</b> starts the resumed run (tab_recalc.py:70).</li><li><b>REFRESH FOLDERS</b> re-reads the list of calc/ folders for the Job-Folder selector (tab_recalc.py:75).</li></ul>",
     "<b>Typical flow</b>: REFRESH FOLDERS &rarr; pick the job folder &rarr; edit CONTROL.txt for the steps you want re-run &rarr; SUBMIT RECALC. The new run writes into its own folder, not over the old one.",
     "<b>Interactions</b>: reads the folders shown in Calculations/Archive; the CONTROL dialect is the one validated in Submit Job (same VALIDATE rules apply)."),
        ("ORCA Builder", "Write and tweak ORCA inputs block by block, then push the result into the CONTROL.txt.",
     "<ul><li><b>Job Name</b> (Text, placeholder 'e.g. water_opt', tab_orca_builder.py:255) and <b>Coordinates</b> (Textarea; COPY COORDINATES copies the current input, tab_orca_builder.py:272-282).</li><li>Numbering helpers: <b>CHECK NUMBERING</b>, <b>APPLY NUMBERING FIX</b>, <b>BACK TO EDITOR</b> (tab_orca_builder.py:287-303).</li><li><b>Charge</b> (IntText, default 0), <b>Multiplicity</b> (tab_orca_builder.py:308-314).</li><li><b>Method</b> dropdown (default PBE0), <b>Job Type</b> [SP, OPT, FREQ, OPT FREQ, ...] (default OPT), <b>Basis Set</b> (default def2-SVP), <b>Dispersion</b> (default D4), <b>RI Approx</b> [RIJCOSX, ...] (default RIJCOSX), <b>Aux Basis</b> (default def2/J), <b>AutoAux</b> checkbox (tab_orca_builder.py:316-348).</li><li><b>Solvation Type</b> [None, CPCM, SMD] and <b>Solvent</b> dropdown (water, acetonitrile, dmso, dmf, methanol, ...), default acetonitrile (tab_orca_builder.py:348-354).</li></ul>",
     "<b>Typical flow</b>: name + coordinates &rarr; choose method/job type/basis/solvation &rarr; the builder renders the ORCA input blocks; apply a numbering fix if CHECK NUMBERING complains &rarr; copy the resulting block into the CONTROL.txt (or submit directly from Submit Job).",
     "<b>Interactions</b>: feeds the same CONTROL dialect Submit Job validates; available on every backend (pure input editor, nothing is submitted from here)."),
        ("TURBOMOLE Builder", "Interactive TURBOMOLE define session plus job submission. <b>Only on SLURM backends that support TURBOMOLE</b> — elsewhere the tab is listed as unavailable with that reason (__init__.py:437-441).",
     "<ul><li><b>Job Name</b>, <b>Coordinates</b> (XYZ or TURBOMOLE $coord format).</li><li><b>TM Module</b> dropdown [ridft, dscf, ricc2, escf, grad, egrad, rdgrad, rigrad, jobex, jobex-noopt, aoforce, NumForce, freeh, ccsdf12, pnoccsd], <b>Parallel</b> [SMP, MPI], <b>CPUs</b>, <b>Mem/CPU</b> (MB), <b>Time Limit</b>.</li><li>Interactive <b>define</b> terminal: Start define, output view, response box + Send, Stop define; the CONTROL editor appears when define ends (tab_turbomole_builder.py:140-284).</li><li><b>Save to Folder</b>, <b>SUBMIT JOB</b>, <b>Load from</b> a prior session and <b>Replay</b> the saved define answers (tab_turbomole_builder.py:610-679); bwHPC wiki and TURBOMOLE manual links (tab_turbomole_builder.py:736-747).</li></ul>",
     "<b>Typical flow</b>: name + coordinates &rarr; Start define &rarr; answer define's prompts in the Send box (every answer is stored, so Load from + Replay restores the session) &rarr; edit the CONTROL from define &rarr; Save to Folder &rarr; SUBMIT JOB (handle_tm_submit, tab_turbomole_builder.py:552-606 calls ctx.backend.submit_turbomole with module, walltime, CPUs, mem/CPU and parallel architecture).",
     "<b>Interactions</b>: standalone (empty refs, tab_turbomole_builder.py:751); the SLURM backend provides the chem/turbomole module loading (backend_slurm.py, turbomole module strings around line 415-418)."),
        ("ChemDarwin", "ChemDarwin seed-generation module integration; the tab appears only when the ChemDarwin tab module is installed (probed lazily at first open, __init__.py:578-600, because the import pulls the stk chain). Default off ('default_visible': False, __init__.py:444).",
     "<ul><li><b>Job name</b> input (required, tab_chemdarwin.py:1265-1267).</li><li><b>Seed input/output</b> widgets (2D + 3D views, classes chemdarwin-seed-out, tab_chemdarwin.py:1275-1277).</li><li><b>Iterations</b> control and custom iteration jump + <b>page navigation</b> (PAGE_SIZE = 50, tab_chemdarwin.py:1250, 1311-1327).</li></ul>",
     "<b>Typical flow</b>: enter a job name &rarr; provide/import a seed structure &rarr; run ChemDarwin generation &rarr; browse the generated molecules page by page (50 per page), picking candidates to continue with.",
     "<b>Interactions</b>: candidates can feed further steps (e.g. Swarm or Submit Job); needs RDKit/stk installed — if the import fails the tab stays hidden with its availability reason."),
        ("DELFIN Agent", "Chat with the AI agent: ask about results, methods, and let it run calculations.",
     "<ul>"
     "<li>Chat input and transcript; every model- or command-authored string is HTML-escaped before rendering (tab_agent.py:205).</li>"
     "<li>Chips above the input shape the next message: <b>permission chip</b> (amber styling, the permission ladder with one rung per restriction — tab_agent.py:2509, :2575-2577, :2774-2779; the backend mapping is delfin.agent.permission_profile.PROFILE_TO_MODE, :2811-2814), plus attach and mode chips.</li>"
     "<li>Slash commands: /model, /effort, /mode, /perms, /perm-cycle, /reset, /sync, /pair and more (command tables tab_agent.py:2867, :2923, :3410).</li>"
     "<li>Model prices are read from one table at module import (tab_agent.py:17); reasoning models running above their profile effort get a visible note (:604-605).</li>"
     "<li>Plan mode: a post-plan permission mode pops the pending plan for approval (:911).</li>"
     "</ul>",
     "<b>Typical flow</b>: pick permission rung &rarr; attach files/keys if needed &rarr; type the task &rarr; Send; stop long turns with the stop control; the agent reports tool use, calculations and jobs in the transcript.",
     "<b>Interactions</b>: runs calculations via Submit Job, reads Calculations/Archive, logs everything in Agent Activity; Office mode locks the permission role to the Office folder (:2823); /perms all_free is the unrestricted rung (:2854-2857)."),
    ("Agent Activity", "Live log of what the agent is doing: sessions, tools, costs, worktrees, schedules.",
     "<ul>"
     "<li>Filter dropdowns: Provider, Mode, Verdict, Class — each with an (all) option (tab_agent_activity.py:695-712).</li>"
     "<li><b>Refresh</b> reloads the snapshot (:717-718); the live-state card always renders, even when empty (:64).</li>"
     "<li>Cards: live state, summary, timeline, cost, MCP servers, worktrees, schedules (:689-694), plus aggregate stats (:444).</li>"
     "</ul>",
     "<b>Typical flow</b>: filter by verdict or class &rarr; Refresh &rarr; follow the timeline while the agent works; use it to audit which tool calls ran and what they cost.",
     "<b>Interactions</b>: mirrors what the DELFIN Agent tab does — no controls here, read-only monitoring."),
    ("Job Status", "State of submitted jobs: running, pending and finished, with start/stop controls.",
     "<ul>"
     "<li><b>Select Job</b> (Dropdown, options filled from the job list, tab_job_status.py:410-411).</li>"
     "<li><b>REFRESH</b> re-reads the job list (tab_job_status.py:419).</li>"
     "<li><b>REBUILD QUEUE</b> rebuilds the scheduler queue after a crash or manual change (tab_job_status.py:424).</li>"
     "<li><b>CANCEL JOB</b> cancels the selected job (danger button, tab_job_status.py:429).</li>"
     "<li><b>Filter</b> [All, Running, Pending] — default All (tab_job_status.py:433-440).</li>"
     "<li>Job table and start/stop info as HTML views (tab_job_status.py:405-408).</li>"
     "</ul>",
     "<b>Typical flow</b>: REFRESH &rarr; filter (e.g. Running) &rarr; select a job to see its state in the table &rarr; CANCEL JOB only if it is stuck; REBUILD QUEUE when the queue display and reality disagree.",
     "<b>Interactions</b>: shows what Submit Job / Recalc / Builder tabs submitted; the job events stream feeds the Agent Activity log (job_events.py)."),
        ("Calculations", "Browse the calc/ folder: inputs, outputs and parsed results, with file management.",
     "<ul>"
     "<li>Folder navigation: <b>&#11014; Up</b>, <b>&#127968; home</b>, <b>&#128260; refresh</b> (tab_calculations_browser.py:374-382); the path prefix is shown as HTML (tab_calculations_browser.py:354).</li>"
     "<li>File operations: <b>&#128465; Delete</b>, <b>Duplicate</b>, <b>+ Folder</b>, <b>Rename</b>, <b>Move</b> — each with its own prompt/status HTML (tab_calculations_browser.py:386-522, 602-650).</li>"
     "<li>Archive round-trip: <b>&#8594; Archive</b> and <b>&#8592; Calculations</b> move folders between calc/ and archive/ (tab_calculations_browser.py:444-452), plus an SSH <b>Transfer</b> view with running-transfer list (tab_calculations_browser.py:459, 687-701).</li>"
     "<li>Structure viewers (3D molecule, trajectory) are built into the browser.</li>"
     "</ul>",
     "<b>Typical flow</b>: navigate to a run folder &rarr; open inputs/outputs or a structure in the viewer &rarr; tidy up with Rename/Duplicate/+ Folder or move finished runs to the Archive; transfers to other machines go through SSH Transfer.",
     "<b>Interactions</b>: Calculations, Archive and Office are three builds of one shared browser module (tab_calculations_browser.py:138-150) — structure IDs are namespaced per build so a viewer opened in one tab cannot hijack another's."),
        ("Archive", "Finished and archived calculation data, served from the archive/ folder.",
     "<ul>"
     "<li>Same browser widgets as Calculations — navigation, Delete/Duplicate/+ Folder/Rename/Move (tab_archive_statistics.py:8-11 builds it with calc_dir replaced by ctx.archive_dir).</li>"
     "<li>Round-trip buttons <b>&#8594; Archive</b> / <b>&#8592; Calculations</b> move folders between the two stores (tab_calculations_browser.py:444-452).</li>"
     "</ul>",
     "<b>Typical flow</b>: in Calculations, select finished runs &rarr; &#8594; Archive; here browse or restore them with &#8592; Calculations when a re-run or Recalc needs them again.",
     "<b>Interactions</b>: archive/ is read-only for the agent tools — file changes here are yours; the Office tab is the third build of the same browser."),
        ("Home", "File browser over your home directory, as a mirror of the Calculations browser.",
     "<ul>"
     "<li>Same browser widgets as Calculations (navigation Up/home/refresh, Delete/Duplicate/+ Folder/Rename/Move); Home builds it with calc_dir = Path.home() (tab_home.py:26-35).</li>"
     "<li>Move targets keep the calculations root, so items can be pushed into calc/ from home (tab_home.py:28-29).</li>"
     "<li>The scan is bounded and dotfiles are hidden (tab_home.py:31-32) — home is large, so deep unbounded scans are avoided.</li>"
     "</ul>",
     "<b>Typical flow</b>: browse ~ &rarr; collect structures/files you need &rarr; Move them into the calculations root &rarr; work with them in Submit Job or Calculations.",
     "<b>Interactions</b>: third build of the shared browser (tab_calculations_browser.py:138-150), besides Calculations and Office."),
        ("Office", "Office-mode file browser for documents and helpers, an exact Calculations-browser mirror.",
     "<ul>"
     "<li>Same browser widgets as Calculations; the browser heading is renamed to Office (tab_office.py:12-14).</li>"
     "</ul>",
     "<b>Typical flow</b>: keep Office-mode documents (templates, PDFs, spreadsheets) in a known folder and manage them here with Rename/Move/Duplicate.",
     "<b>Interactions</b>: shares the browser module with Calculations, Archive and Home (tab_calculations_browser.py:138-150)."),
        ("Literature", "Literature search and storage: Google-Scholar style results with DOI/PDF links.",
     "<ul>"
     "<li>Search runs in a background thread (delfin-lit-search, tab_literature.py:28) so the tab stays responsive.</li>"
     "<li>Each result card shows title link, people, venue, year, citations, and DOI / PDF links when present (tab_literature.py:42-56; DOI renders as doi.org link :53-56).</li>"
     "<li>Same browser widgets below the search: navigation, Delete/Duplicate/+ Folder/Rename/Move (tab_literature.py:234-246) for the stored references.</li>"
     "</ul>",
     "<b>Typical flow</b>: search a topic &rarr; open the DOI/PDF of a hit &rarr; store what matters in the Literature folder via the browser widgets.",
     "<b>Interactions</b>: stored references back method choices in calculations; no submit path from here."),
        ("Remote Archive", "Browse and transfer data from a remote archive server. Available only after enabling the remote archive in Settings.",
     "<ul>"
     "<li>Browser navigation Up / home / refresh / Open (tab_remote_archive.py:163-178).</li>"
     "<li><b>&#8594; Calculations</b> transfer-back button (tab_remote_archive.py:214-215).</li>"
     "<li><b>Download</b> button with status view (tab_remote_archive.py:362-367) — download jobs are created through create_download_job (tab_remote_archive.py:43).</li>"
     "<li>3D visualisation controls exist but are disabled by default (visualize_enabled=False, tab_remote_archive.py:125); a high-res PNG download button sits on the 3D view (tab_remote_archive.py:549, 695).</li>"
     "</ul>",
     "<b>Typical flow</b>: enable the remote archive in Settings &rarr; browse the remote tree &rarr; Download selected data or transfer it straight &#8594; Calculations.",
     "<b>Interactions</b>: endpoint configured in Settings; the local counterpart is the Archive tab."),
        ("Tools", "Tools &amp; platform panel: environment, tool status, installs, and run browsing.",
     "<ul>"
     "<li>One combined HTML view of environment, tool status, capabilities, applications, keys and recent runs (render_all_html joins six renderers, tab_tools.py:304-311).</li>"
     "<li><b>Refresh</b> re-renders the view (tab_tools.py:323).</li>"
     "<li><b>Install all open-source</b> installs every auto-installable open-source tool — ORCA, Turbomole and Multiwfn are never installed automatically (tooltip, tab_tools.py:324-328).</li>"
     "<li><b>Install selected</b> installs the tools checked in the per-tool checkboxes (tab_tools.py:329-333).</li>"
     "<li><b>Run:</b> dropdown shows recent runs with status, file list and log tail (tab_tools.py:336-340, render_run_detail_html :247).</li>"
     "</ul>",
     "<b>Typical flow</b>: check the platform view &rarr; install missing open-source tools (selected or all) &rarr; pick a run to inspect its files and log tail without leaving the tab.",
     "<b>Interactions</b>: installs land where qm_runtime resolves tools; the Settings tab manages the same tool set for QM engines explicitly."),
        ("Settings", "Global DELFIN configuration: clusters, backends, QM tool installs, secrets and defaults.",
     "<ul>"
     "<li><b>Show Values</b> toggle reveals/hides sensitive values (credentials) in the view (tab_settings.py:85-91, 1645).</li>"
     "<li><b>Reload</b> re-reads the configuration; <b>Validate Setup</b> runs the runtime validation (tab_settings.py:91-99).</li>"
     "<li>QM tool management: <b>Prepare qm_tools</b>, <b>Install qm_tools</b>, <b>Update qm_tools</b> plus per-tool install buttons xtb, crest, dftb+, xtb4stda, std2, mopac (tab_settings.py:100-123).</li>"
     "<li>Cluster setup for bwUniCluster: <b>Setup cluster</b>, verification and full install buttons (tab_settings.py:262-272).</li>"
     "<li>The tab reuses the same API keys/credentials as the DELFIN agent (tab_settings.py:3956, 4074).</li>"
     "</ul>",
     "<b>Typical flow</b>: set up the cluster once (Setup cluster &rarr; verify) &rarr; install the QM tools you need &rarr; Validate Setup &rarr; Reload. Remote Archive and backend availability elsewhere in the dashboard follow this configuration.",
     "<b>Interactions</b>: gates the TURBOMOLE Builder (backend availability), Remote Archive (endpoint) and every submit path (backend choice resolved here, tab_settings.py:37)."),
        ("Ketcher", "Draw and edit molecules in the browser (2D structure editor), registered tab (tab_ketcher.py:60).",
     "<ul>"
     "<li>The Ketcher panel (ketcher_panel.build, height 72vh, tab_ketcher.py:28): drawing canvas, file list, Open/Open drawing and Refresh files buttons.</li>"
     "<li>Drawings opened from the Archive or Office tabs land here: the refs dict (open_drawing, refresh_files) is shared via ctx.ketcher_refs and updated in place (:40-52), so the two browser mirrors never hold an empty copy.</li>"
     "<li>SMARTS queries are supported via ketcher_smarts.py.</li>"
     "</ul>",
     "<b>Typical flow</b>: draw or open a structure &rarr; save/export as coordinates &rarr; use them in Submit Job or a Builder tab.",
     "<b>Interactions</b>: receives structures opened in Archive/Office; feeds Submit Job and the builders. Registered tab, present on every backend."),
    ("Swarm", "Metal-complex builder: guided start of a PSO swarm search over candidate structures, plus the live board (tab_swarm.py:544). Read-only until you click.",
     "<ul>"
     "<li><b>Refresh</b> (tab_swarm.py:563) and <b>Stop swarm job</b> — Stop is two-step armed (:571).</li>"
     "<li>Spec editor: <b>Load</b> (:586), <b>New from example</b> (:588), <b>Save</b> (:591), spec <b>Textarea</b> (:603); the saved path is published as ctx.swarm_spec_path and read by the board and Swarm Watch (:729).</li>"
     "<li><b>Check</b> validates the spec (:612, check_start :404), <b>Start swarm</b> starts the job (:616, start_checked :457).</li>"
     "<li>Board and job line refresh every SWARM_REFRESH_S while the tab is on screen (:547-548).</li>"
     "</ul>",
     "<b>Typical flow</b>: New from example or Load &rarr; edit the spec &rarr; Check &rarr; Start swarm &rarr; watch on the board and in Swarm Watch.",
     "<b>Interactions</b>: needs a scheduler (throttled scheduler_client when started from here, :550-551); results are readable in Calculations; the roster is shared with Swarm Watch (:213-224)."),
    ("Swarm Watch", "Live monitoring of a running swarm or wave, also without a roster (tab_swarm.py:947, registered at :1489-1490).",
     "<ul>"
     "<li><b>Refresh</b> (&#8635;, reads the sessions now, tab_swarm.py:1320).</li>"
     "<li>Message box + <b>Send</b> (Enter sends, Shift+Enter newline, :1221, :1337-1342) to steer or annotate the running swarm.</li>"
     "<li>Watched rows and feed render as HTML (:887, :1103); the tab borrows the Agent tab's stylesheet with delfin-swv- class namespace (:1125-1200).</li>"
     "</ul>",
     "<b>Typical flow</b>: start a swarm in the Swarm tab &rarr; open Swarm Watch &rarr; follow rows/feed &rarr; Send messages or stop in the Swarm tab.",
     "<b>Interactions</b>: reads ctx.swarm_spec_path set by Swarm (:729); companion of the Swarm tab."),
    ("Reactions", "Reaction graph of a calculation path, registered tab (tab_reaction_graph.py:1182-1194).",
     "<ul>"
     "<li>ReactionGraphPanel over calc_dir (default ~/calc, :1192-1193): shows the reaction network derived from the calculations in the folder.</li>"
     "<li>Nothing is published onto ctx — the registry keeps the widget and drops the refs (:1185-1190), so no other tab can reach into it.</li>"
     "</ul>",
     "<b>Typical flow</b>: open the tab &rarr; pick a calculation path &rarr; inspect nodes and edges; details in reaction_graph.py / reaction_profile.py.",
     "<b>Interactions</b>: reads the same calc/ folder as Calculations; helps interpret CO2-reduction and redox chains."),
    ("Review", "Review and grade calculation results: structure view with categories, ratings and comments (tab_review.py:1385, registered :1393).",
     "<ul>"
     "<li>Grid of category toggles, each showing its key and the comment of this frame (:26-29); category keys are 1-9, 0, s1-s9 for Shift (category_key :142-143).</li>"
     "<li>An Unsure toggle marks a frame uncertain; per-file and whole-review comments are saved with the rating (:26-29).</li>"
     "<li>Keyboard: removes the word &quot;unsicher&quot; from a comment, H toggles H atoms, R resets the view, Esc leaves a text field (:37-38); marked atoms are ringed in the category colour (:93).</li>"
     "</ul>",
     "<b>Typical flow</b>: pick a calculation &rarr; rate frames by pressing category keys &rarr; comment &rarr; save; ratings are stored with the calculation.",
     "<b>Interactions</b>: reads Calculations data; grades feed later reviews of the same path."),
    ("Pipelines", "Live status of multi-step calculation pipelines (tab_pipeline_status.py:185, PipelineStatusWidget :50).",
     "<ul>"
     "<li>The widget is a pure listener: pipe.on_step(widget.on_step) pushes step events into the display (:53-58).</li>"
     "<li>Returns (widget, refs) with refs {'pipeline_status': widget} (:191-193).</li>"
     "</ul>",
     "<b>Typical flow</b>: start a pipeline (via the agent or a workflow) &rarr; watch each step appear here as it runs.",
     "<b>Interactions</b>: reflects pipeline.py runs; no controls of its own."),
    ("Help", "This tab: feature groups, CONTROL keys, and the dashboard tab guide (tab_help.py:277-293).",
     "<ul>"
     "<li>Feature groups as accordion sections (titles &quot;Group — summary&quot;, :277-280).</li>"
     "<li>The tab guide is the rich _TABS_HELP 5-tuple below (:293): purpose, widgets, typical flow, interactions per tab.</li>"
     "<li>The CONTROL key table is searchable from the header search bar; picking a key routes here pre-filtered (delfin/dashboard/__init__.py:1194-1207).</li>"
     "</ul>",
     "<b>Typical flow</b>: search a key or tab name in the header &rarr; Help opens filtered &rarr; read the accordion section.",
     "<b>Interactions</b>: target of the header-search routing; sources are CONTROL_FIELD_SPECS, ALLOWED_BY_KEY and _TABS_HELP only — no second copies."),
)

__all__ = ["TABS_HELP"]
