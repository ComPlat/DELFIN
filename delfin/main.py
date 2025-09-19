# delfin/main.py

def main():
    """
    Small wrapper that loads the actual main() function.
    Adjust the import order below if your primary file has a different name.
    """
    # Try several common filenames:
    for modpath in (
        "delfin",      # delfin/delfin.py (common)
        "cli",         # delfin/cli.py
        "app",         # delfin/app.py
        "__main__",    # delfin/__main__.py
    ):
        try:
            mod = __import__(f"{__package__}.{modpath}", fromlist=["main"])
            run = getattr(mod, "main", None)
            if callable(run):
                return run()
        except ModuleNotFoundError:
            continue

    raise SystemExit(
        "No main() found. Create, for example, delfin/delfin.py with a main() function, "
        "or adjust the entry point in pyproject.toml."
    )
