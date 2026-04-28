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
