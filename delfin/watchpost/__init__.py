"""Watchpost: a read-only look-out over the user's own account.

It reports traces left by intrusions and misuse -- SSH keys, persistence
in shell start-up files, odd processes and listeners, changed
credentials. It never changes, deletes, blocks or terminates anything
outside its own report directory. See README-less: the CLI --help is
the documentation.
"""
