"""Caller-side pass-through (QS arbitration, nacht-s11): the caller of
smiles_to_xyz owns the conversion result tuple (xyz, error) and emits the
conversion_failed finding when error is set.  Driven through the REAL adapter
path (SmilesToXyzAdapter.execute), not the checker directly:

1. A SMILES that crashes conversion (cisplatin `Cl[Pt]Cl(N)N`, the
   AtomValenceException defect) yields StepStatus.FAILED AND
   data["findings"] carrying kind conversion_failed.
2. A successful conversion carries the post-conversion checker findings in
   data["findings"] (here: h_ambiguity for unbracketed donors on a metal).
3. A clean, unambiguous conversion carries no findings.
"""
from delfin.tools.adapters.smiles import SmilesToXyzAdapter


def _run(smiles, tmp_path):
    adapter = SmilesToXyzAdapter()
    work = tmp_path / "sx1-adapter-probe"
    work.mkdir(parents=True, exist_ok=True)
    return adapter.execute(work, smiles=smiles)


def test_conversion_crash_yields_finding(tmp_path):
    res = _run("Cl[Pt]Cl(N)N", tmp_path)
    assert res.status.value == "failed" or res.status.name == "FAILED"
    findings = res.data.get("findings")
    assert findings, f"no findings on conversion failure: {res.data}"
    assert "conversion_failed" in [f["kind"] for f in findings], findings


def test_successful_conversion_carries_checker_findings(tmp_path):
    res = _run("[Pt](Cl)(Cl)(N)N", tmp_path)
    findings = res.data.get("findings")
    assert findings is not None, f"findings key missing: {res.data}"
    assert "h_ambiguity" in [f["kind"] for f in findings], findings


def test_clean_conversion_carries_no_findings(tmp_path):
    res = _run("[Pt](Cl)(Cl)([NH3])[NH3]", tmp_path)
    findings = res.data.get("findings")
    assert findings == [], f"unexpected findings: {res.data}"
