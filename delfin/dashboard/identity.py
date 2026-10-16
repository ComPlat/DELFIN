"""The dashboard server's identity provider: sessions, logout, refused logins.

``delfin-voila`` starts jupyter-server with
``--ServerApp.identity_provider_class=delfin.dashboard.identity.DelfinIdentityProvider``
(when this module imports; otherwise it says so and starts with
jupyter-server's own provider). It is jupyter-server's
``PasswordIdentityProvider`` -- the token check, the login form and the
cookie handling are jupyter-server's -- with four changes:

* every login cookie carries a server-side session id
  (:mod:`login_session`); a cookie whose session is unknown, closed, idle
  or too old authenticates nothing;
* ``/logout`` (jupyter-server's ``LogoutHandler``, which calls
  ``clear_login_cookie``) closes the session of the cookie it carries;
* a wrong token -- on the login form, in the URL, in the ``Authorization``
  header -- is recorded, rate limited, as the security event
  ``dashboard_login_refused`` (:class:`login_session.RefusedLogins`);
* ``allow_password_change`` defaults to False: with the default, a POST of
  the token and a ``new_password`` to ``/login`` wrote a password hash into
  the Jupyter config directory, a credential that outlived the token and
  applied to every Jupyter server of the account.

The identity provider is the one place every token and cookie check
passes through, and jupyter-server documents it as the extension point
(``ServerApp.identity_provider_class``). Rejected for the refused
logins: wrapping ``LoginHandler`` or tornado's 403 path in
``server_guard`` -- a handler wrapper sees the login form but not a token
in a URL or header, and patches a class it does not own.
"""

from __future__ import annotations

import json
import time
from typing import Any

from jupyter_server.auth.identity import PasswordIdentityProvider
from traitlets import default

from delfin.dashboard import login_session


class DelfinIdentityProvider(PasswordIdentityProvider):
    """jupyter-server's token/password provider with server-side sessions."""

    @default("allow_password_change")
    def _no_password_change(self):
        return False

    # -- the session store, made on first use (the port is final by then) ---
    _delfin_sessions = None
    _delfin_refused = None

    def _sessions(self) -> login_session.Sessions:
        if self._delfin_sessions is None:
            from delfin.dashboard.web_security import login_days

            path = None
            app = getattr(self, "parent", None)
            try:
                path = login_session.store_path(app.runtime_dir, app.port)
                path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
            except Exception:
                path = None
            self._delfin_sessions = login_session.Sessions(
                path, idle=login_session.idle_seconds(),
                lifetime=login_days() * 86400.0, clock=time.time, log=self.log)
        return self._delfin_sessions

    def _refused(self) -> login_session.RefusedLogins:
        if self._delfin_refused is None:
            self._delfin_refused = login_session.RefusedLogins(log=self.log)
        return self._delfin_refused

    def _note_refused(self, how: str, handler: Any) -> None:
        try:
            request = handler.request
            self._refused().refused(how, request.remote_ip, request.path)
        except Exception:
            pass

    # -- the cookie --------------------------------------------------------
    def user_to_cookie(self, user):
        data = json.loads(super().user_to_cookie(user))
        data[login_session.COOKIE_FIELD] = self._sessions().open()
        return json.dumps(data)

    def user_from_cookie(self, cookie_value):
        data = json.loads(cookie_value)
        if not isinstance(data, dict):
            return None
        if not self._sessions().touch(data.get(login_session.COOKIE_FIELD)):
            return None
        return super().user_from_cookie(cookie_value)

    def clear_login_cookie(self, handler):
        try:
            raw = handler.get_secure_cookie(self.get_cookie_name(handler),
                                            **self.get_secure_cookie_kwargs)
            if raw:
                data = json.loads(raw.decode())
                if isinstance(data, dict):
                    self._sessions().close(data.get(login_session.COOKIE_FIELD))
        except Exception:
            pass
        super().clear_login_cookie(handler)

    # -- refused logins ----------------------------------------------------
    def process_login_form(self, handler):
        user = super().process_login_form(handler)
        if user is None and handler.get_argument("password", default=""):
            self._note_refused("wrong token on the login form", handler)
        return user

    async def get_user_token(self, handler):
        user = await super().get_user_token(handler)
        if user is None and getattr(handler, "token", None) \
                and self.get_token(handler):
            self._note_refused("wrong token in the request", handler)
        return user


#: What ``delfin-voila`` passes as ``ServerApp.identity_provider_class``.
CLASS_PATH = f"{__name__}.{DelfinIdentityProvider.__name__}"
