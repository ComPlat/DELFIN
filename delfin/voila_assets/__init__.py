"""Static assets shipped with the DELFIN Voila dashboard.

Holds resources served by ``delfin-voila`` that are not Python code, e.g. the
branded ``login_templates/login.html`` that overrides the generic
jupyter_server login page. The directory is a real package so it can be
located from an installed wheel via ``importlib.resources.files``.
"""
