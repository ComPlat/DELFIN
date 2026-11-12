"""The page after Log out is DELFIN's, like the sign-in page.

jupyter_server rendered its own logout.html: an unstyled "Successfully
logged out. Proceed to the login page." The sign-in page had been DELFIN's
for a while. Both are now shipped side by side, share one card style, and
are materialised together into the directory handed to jupyter_server.
"""

from __future__ import annotations

from pathlib import Path


def _rendered(tmp_path, monkeypatch):
    from delfin import cli_voila

    monkeypatch.setattr(cli_voila, "_secure_cache_dir",
                        lambda name: (tmp_path / name).mkdir(exist_ok=True)
                        or tmp_path / name)
    folder = cli_voila._prepare_login_template_dir()
    assert folder is not None
    return Path(folder)


def _env(folder):
    import jinja2

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(str(folder)),
                             autoescape=True, extensions=["jinja2.ext.i18n"])
    env.install_null_translations()
    return env


def test_both_pages_are_materialised_with_the_logo(tmp_path, monkeypatch):
    folder = _rendered(tmp_path, monkeypatch)
    names = {p.name for p in folder.iterdir()}
    assert {"login.html", "logout.html", "_delfin_card_style.html"} <= names
    for name in ("login.html", "logout.html"):
        text = (folder / name).read_text()
        assert "__DELFIN_LOGO_URI__" not in text
        assert ('src="data:image/png;base64,' in text
                or "LOGO_TILE_START" not in text)


def test_the_sign_out_page_renders_as_delfins(tmp_path, monkeypatch):
    env = _env(_rendered(tmp_path, monkeypatch))
    page = env.get_template("logout.html").render(
        base_url="/", login_available=True, message=None)
    assert "<title>Signed out · DELFIN Dashboard</title>" in page
    assert 'class="wordmark">DELFIN<' in page
    assert 'href="/login"' in page and "Sign in again" in page
    assert ".submit {" in page                    # the shared card style
    # jupyter_server passes its notice as info: shown, not as an error.
    told = env.get_template("logout.html").render(
        base_url="/", login_available=True,
        message={"info": "Successfully logged out."})
    assert '<div class="message info">Successfully logged out.</div>' in told
    assert ".message.info {" in told
    no_login = env.get_template("logout.html").render(
        base_url="/x/", login_available=False, message=None)
    assert 'href="/x/"' in no_login and "Back to the dashboard" in no_login


def test_the_sign_in_page_still_has_its_style(tmp_path, monkeypatch):
    env = _env(_rendered(tmp_path, monkeypatch))
    page = env.get_template("login.html").render(
        base_url="/", next="/", login_available=True, token_available=True,
        message=None, xsrf_form_html=lambda: "")
    assert "--brand: #1976d2" in page and ".submit {" in page
    assert 'id="login_submit"' in page
