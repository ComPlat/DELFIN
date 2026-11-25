# Provenance of the a1/a2 guard copies in this directory

`a1_guard.py` and `a2_guard.py` are verbatim copies of the two arms' guards,
taken from their branches so the offender-count comparison can be re-run
without checking either branch out:

    git show wave/w46/x1-a1:tests/test_the_package_compiles_with_no_syntax_warning.py > a1_guard.py
    git show wave/w46/x1-a2:tests/test_no_module_carries_a_syntax_warning.py            > a2_guard.py

`two_escape.py` loads both and prints what each reports for one file holding
two undefined escapes. It is an instrument, not part of the change.

`two_bad.py` appears in the output but is no repo file: the script writes that
fixture into a fresh temporary directory (`two_escape.py`, the
`tempfile.TemporaryDirectory()` block). Line 2 is the opening line of the
fixture's first docstring — the line CPython blames for an escape written
inside it — and line 5 the opening line of the second.
