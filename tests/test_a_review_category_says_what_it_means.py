"""A review category button names the defect in words a chemist reads, not its stored code.

The team set's kategorien.txt holds short German codes (LAENGE_ML, DURCHDRINGUNG, ...).  They stay
what is stored; the button shows a label, its tooltip and the fold-out under the buttons say what the
category means, and the team guide's table carries the same label.
"""
from pathlib import Path

from delfin import review as rv

ROOT = Path(__file__).resolve().parents[1]
SET_CATEGORIES = ROOT / 'review_sets' / 'kalib_2026_10' / 'kategorien.txt'
GUIDE = ROOT / 'docs' / 'REVIEW_TEAM_GUIDE.md'


def _codes():
    return [ln.strip() for ln in SET_CATEGORIES.read_text(encoding='utf-8').splitlines() if ln.strip()]


def test_every_category_of_the_team_set_has_a_label_and_a_meaning():
    for code in _codes():
        assert code in rv.CATEGORY_INFO, code
        label, meaning = rv.CATEGORY_INFO[code]
        assert label and label != code or code == 'H', code
        assert len(meaning) > 15, code


def test_labels_are_distinct_and_fit_a_button():
    labels = [rv.category_label(c) for c in _codes()]
    assert len(set(labels)) == len(labels)
    assert max(len(lab) for lab in labels) <= 16


def test_a_category_without_an_entry_shows_its_own_name():
    assert rv.category_label('wrong isomer') == 'wrong isomer'
    assert rv.category_meaning('wrong isomer') == ''


def test_the_guide_table_carries_the_same_label_for_each_code():
    text = GUIDE.read_text(encoding='utf-8')
    for code in _codes():
        assert f'| {rv.category_label(code)} | `{code}` |' in text, code


def test_the_buttons_show_the_label_and_the_tooltip_the_meaning():
    from delfin.dashboard import tab_review as tr
    panel = tr.ReviewPanel()
    panel._build_categories(_codes())
    by_code = {b._review_category: b for b in panel.category_buttons}
    btn = by_code['DURCHDRINGUNG']
    assert btn.description.endswith('Interpenetration')
    assert 'DURCHDRINGUNG' not in btn.description
    assert 'stuck through each other' in btn.tooltip and 'DURCHDRINGUNG' in btn.tooltip
    shown = panel.cat_help_html.value
    for code in _codes():
        assert rv.category_label(code).replace('&', '&amp;') in shown
    assert 'Where two are close' in shown
