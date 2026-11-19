"""Sessions can be renamed; the working directory completes like a shell.

The person (11.11.): renaming sessions, and picking the path the way it is
usual today -- suggestions while typing, Tab completes, arrows choose.
"""
from __future__ import annotations

import inspect

from delfin.dashboard import agent_sessions as AS


def test_a_chosen_name_wins_over_the_first_message():
    state = {"chat_messages": [{"role": "user", "content": "Hi there"}],
             "_messenger_titles": {"sid1": "RX1 review"}}
    assert AS._session_title(state, "sid1") == "RX1 review"
    assert AS._session_title(state, "other") == "Hi there"
    assert AS._session_title(state) == "Hi there"


def test_the_name_is_kept_where_the_save_reads_it():
    src = inspect.getsource(AS.create_tab)
    # the agent tab's save writes _messenger_titles[sid] as the saved title
    assert 'state.setdefault("_messenger_titles", {})' in src
    assert "_store.rename_session(sid, name)" in src
    assert "delfin-session-rename-btn" in src


def test_escape_cancels_a_rename_and_is_bound_once():
    js = AS._RENAME_KEYS_JS
    assert "Escape" in js and ".delfin-live-rename-cancel" in js
    assert "__delfinRenameKeys" in js


def test_the_folder_box_has_keys_like_a_shell():
    js = AS._PATH_KEYS_JS
    for key in ("'Tab'", "'ArrowDown'", "'ArrowUp'", "'Enter'", "'Escape'"):
        assert key in js, key
    assert "__delfinPathKeys" in js
    src = inspect.getsource(AS.create_tab)
    assert "widgets.Text(" in src and "delfin-path-input" in src
    assert "delfin-path-suggest" in src


def test_what_fits_the_typed_path_comes_first():
    src = inspect.getsource(AS.create_tab)
    i = src.index("def _show_suggestions")
    body = src[i:i + 1500]
    assert "fits + [c for c in choices if c not in fits]" in body


def test_every_page_script_is_registered_not_only_the_first():
    import types
    AS._REGISTERED_INIT_JS.clear()
    sent = []
    ctx = types.SimpleNamespace(add_init_js=sent.append)
    for js in (AS._SPLITTER_INIT_JS, AS._RENAME_KEYS_JS, AS._SPLITTER_INIT_JS):
        AS._register_init_js(ctx, js)
    assert sent == [AS._SPLITTER_INIT_JS, AS._RENAME_KEYS_JS]


def test_the_list_runs_its_own_scripts_when_shown():
    html = AS._self_running_html(AS._RENAME_KEYS_JS, AS._PATH_KEYS_JS,
                                 AS._COMPLETE_JS)
    assert html.count('<script type="text/x-delfin" class="delfin-sessions-js">') == 3
    assert "onerror=" in html and "(0,eval)(t.textContent)" in html
    assert "__delfinComplete" in html and "__delfinPathKeys" in html
    import inspect
    assert "_self_running_html(_RENAME_KEYS_JS, _PATH_KEYS_JS" in inspect.getsource(AS.create_tab)
