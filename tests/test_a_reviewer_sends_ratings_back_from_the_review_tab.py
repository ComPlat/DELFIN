"""A reviewer of a team set reserves the block and sends the ratings back from the Review tab.

The buttons call the same functions as ``delfin review claim | submit`` (a pull request each, to a
fork without write access).  An empty name on a team set becomes the GitHub account gh is logged
in with; a file that is no team set shows no such buttons.
"""
from delfin import review_team as rt
from delfin.dashboard import tab_review as tr

XYZ = '3\nframe0\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n'


def _team_folder(tmp_path, monkeypatch):
    (tmp_path / 'a.xyz').write_text(XYZ, encoding='utf-8')
    monkeypatch.setattr(rt, 'find_assignment', lambda raw: tmp_path / 'assignment.json')
    monkeypatch.setattr(rt, 'load_assignment', lambda found: {'shared': [], 'files_total': 1})
    monkeypatch.setattr(rt, 'share_of', lambda raw, assignment, name: (name, ['a.xyz'], None))
    return tmp_path


def test_an_empty_name_on_a_team_set_becomes_the_github_account(tmp_path, monkeypatch):
    folder = _team_folder(tmp_path, monkeypatch)
    monkeypatch.setattr(rt, '_gh_login', lambda run: 'octo-chemist')
    panel = tr.ReviewPanel()
    assert panel.open(str(folder)) is not None
    assert panel.name_input.value == 'octo-chemist'
    assert panel.team_box.layout.display == ''


def test_the_send_button_opens_a_pull_request_and_shows_its_link(tmp_path, monkeypatch):
    folder = _team_folder(tmp_path, monkeypatch)
    calls = []
    monkeypatch.setattr(rt, 'submit', lambda f, name: calls.append(('submit', f, name))
                        or {'pr': 'https://github.com/ComPlat/DELFIN/pull/1', 'next': ''})
    monkeypatch.setattr(rt, 'claim', lambda f, name: calls.append(('claim', f, name))
                        or {'pr': None, 'next': 'gh is not installed'})
    panel = tr.ReviewPanel()
    panel.name_input.value = 'ana'
    panel.open(str(folder))
    panel.send_btn.click()
    assert calls == [('submit', str(folder), 'ana')]
    assert 'Ratings sent' in panel.team_html.value and 'pull/1' in panel.team_html.value
    panel.claim_btn.click()
    assert calls[-1] == ('claim', str(folder), 'ana')
    assert 'Not sent' in panel.team_html.value and 'gh is not installed' in panel.team_html.value


def test_a_file_that_is_no_team_set_shows_no_send_buttons(tmp_path):
    (tmp_path / 'one.xyz').write_text(XYZ, encoding='utf-8')
    panel = tr.ReviewPanel()
    panel.open(str(tmp_path / 'one.xyz'))
    assert panel.team_box.layout.display == 'none'
    assert panel.send_back('submit') is None
