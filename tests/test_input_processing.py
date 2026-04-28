from delfin.dashboard.input_processing import parse_inp_resources


def test_parse_inp_resources_inline_pal_block():
    """%pal nprocs N end on a single line must be parsed.

    Regression: previous regex required `nprocs` at line start, which
    missed the inline `%pal nprocs 40 end` form. Misparse caused the
    ORCA Builder to silently fall back to the PAL widget value (12)
    while the .inp told ORCA to use 40 procs -> SLURM under-allocated
    memory -> OOM kill (RSS_paper_85_NEB-TS, job 4089862).
    """
    inp = "!PBE0 def2-SVP\n\n%maxcore 6000\n%pal nprocs 40 end\n"
    pal, maxcore = parse_inp_resources(inp)
    assert pal == 40
    assert maxcore == 6000


def test_parse_inp_resources_multiline_pal_block():
    inp = "!PBE0\n\n%pal\n  nprocs 16\nend\n%maxcore 4000\n"
    pal, maxcore = parse_inp_resources(inp)
    assert pal == 16
    assert maxcore == 4000


def test_parse_inp_resources_missing_returns_none():
    pal, maxcore = parse_inp_resources("!PBE0 def2-SVP\n")
    assert pal is None
    assert maxcore is None


def test_parse_inp_resources_nprocs_equals_syntax():
    """ORCA accepts ``nprocs=N``; the parser must too."""
    inp = "!PBE0\n%pal\n  nprocs=24\nend\n%maxcore 4000\n"
    pal, maxcore = parse_inp_resources(inp)
    assert pal == 24
    assert maxcore == 4000


def test_parse_inp_resources_pal_keyword_shortcut():
    """``! PAL8`` keyword shortcut is the third valid PAL syntax."""
    inp = "! PBE0 def2-SVP PAL8\n%maxcore 3000\n"
    pal, maxcore = parse_inp_resources(inp)
    assert pal == 8
    assert maxcore == 3000


def test_parse_inp_resources_case_insensitive_maxcore():
    """``%MaxCore`` and ``%MAXCORE`` must parse equivalently to ``%maxcore``."""
    inp = "!PBE0\n%MaxCore 5000\n%PAL\n  Nprocs 16\nEnd\n"
    pal, maxcore = parse_inp_resources(inp)
    assert pal == 16
    assert maxcore == 5000


def test_parse_inp_resources_maxcore_equals_syntax():
    inp = "!PBE0\n%maxcore=8000\n%pal nprocs 32 end\n"
    pal, maxcore = parse_inp_resources(inp)
    assert pal == 32
    assert maxcore == 8000


def test_pal_block_replace_handles_inline_and_multiline():
    """Regression: tab_orca_builder rewrites %pal when the user changes the
    widget. The old re.sub only matched the multi-line form, leaving the
    inline form (`%pal nprocs 40 end`) untouched and the inp out of sync
    with the widget. The current regex uses (?is)%pal\\b.*?\\bend\\b so
    both forms collapse to the canonical multi-line template.
    """
    import re
    pattern = r'(?is)%pal\b.*?\bend\b'
    new_pal = '%pal\n  nprocs 24\nend'

    inline = "!PBE0\n%pal nprocs 40 end\n%maxcore 6000\n"
    out = re.sub(pattern, new_pal, inline, count=1)
    assert "nprocs 24" in out
    assert "nprocs 40" not in out

    multiline = "!PBE0\n%pal\n  nprocs 40\nend\n%maxcore 6000\n"
    out = re.sub(pattern, new_pal, multiline, count=1)
    assert "nprocs 24" in out
    assert "nprocs 40" not in out
