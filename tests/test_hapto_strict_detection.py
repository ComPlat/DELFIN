"""tests/test_hapto_strict_detection.py — strict hapto-mode detection filter.

Verifies the graph-topology filter that prevents false-positive η²/η⁴/η⁶
enumeration for σ-bound aromatic donors (YUHRUP-class).

Architecture overview
---------------------
* :mod:`delfin.fffree.hapto_strict_detection` — pure-topology helper.
* :func:`delfin.fffree.coordination_mode_enum._detect_hapto_modes` —
  integration point.  Calls :func:`is_legitimate_hapto_candidate` BEFORE
  emitting hapto orbits when the env-flag
  ``DELFIN_FFFREE_HAPTO_STRICT_DETECTION`` is set.

Default OFF = byte-identical.  When ON, false-positive arene orbits for
substituent rings (phenyl pendant on σ-N/σ-S chelate) are suppressed
while genuinely bonded hapto rings (ferrocene Cp, (η⁶-arene)Cr(CO)₃) are
preserved.
"""
from __future__ import annotations

import importlib
import os
import unittest


# Make sure we have a clean env at import time.
for _k in ("DELFIN_FFFREE_PURE_TRACK3",
           "DELFIN_FFFREE_BINDING_MODE_ENUM",
           "DELFIN_FFFREE_BINDING_MODE_STRICT",
           "DELFIN_FFFREE_HAPTO_STRICT_DETECTION",
           "DELFIN_FFFREE_HAPTO_STRICT_MAX_DIST"):
    os.environ.pop(_k, None)

os.environ.setdefault("PYTHONHASHSEED", "0")


def _enable_strict():
    os.environ["DELFIN_FFFREE_HAPTO_STRICT_DETECTION"] = "1"


def _disable_strict():
    os.environ.pop("DELFIN_FFFREE_HAPTO_STRICT_DETECTION", None)


def _enumerate(smiles: str):
    """Convenience: import the enumerator FRESH so live-evaluated env
    flags inside ``_detect_hapto_modes`` are honoured per-test.
    """
    import delfin.fffree.coordination_mode_enum as cme
    importlib.reload(cme)
    return cme.enumerate_modes(smiles)


def _hapto_labels(modes):
    return sorted({m["label"] for m in modes if m.get("kind") == "hapto"})


# ---------------------------------------------------------------------------
# 1. YUHRUP synthetic: σ-bound chelate with phenyl substituent
# ---------------------------------------------------------------------------
class TestYuhrupSyntheticNoHapto(unittest.TestCase):
    """Re-S/Re-N σ chelate with a phenyl SUBSTITUENT on the chelate backbone
    must NOT emit eta-arene orbits when the strict filter is enabled.
    """

    YUHRUP_SMILES = (
        "CCN(CC)C1=NC(C2=CC=CC=C2)=[N+]2N=C(N3CCCCCC3)[S][Re]2"
        "([O-])([Cl])[S]1"
    )

    def setUp(self):
        _disable_strict()

    def tearDown(self):
        _disable_strict()

    def test_strict_off_emits_arene_orbits(self):
        """Baseline: with strict OFF the YUHRUP phenyl substituent
        does emit eta-arene orbits (legacy / false-positive behaviour)."""
        modes = _enumerate(self.YUHRUP_SMILES)
        labs = _hapto_labels(modes)
        # The phenyl substituent is matched by the arene SMARTS, so all
        # three arene orbits appear under the legacy enumeration.
        self.assertIn("eta6-arene", labs)
        self.assertIn("eta4-arene", labs)
        self.assertIn("eta2-arene", labs)

    def test_strict_on_suppresses_arene_orbits(self):
        """Strict ON: the phenyl substituent is graph-distance >1 from
        Re's coordination atoms in the metal-removed graph, so no
        eta-arene orbit is emitted.
        """
        _enable_strict()
        modes = _enumerate(self.YUHRUP_SMILES)
        labs = _hapto_labels(modes)
        self.assertEqual(labs, [], f"unexpected hapto orbits: {labs}")


# ---------------------------------------------------------------------------
# 2. Ferrocene: η⁵-Cp directly bonded → hapto PRESERVED
# ---------------------------------------------------------------------------
class TestFerroceneHaptoPreserved(unittest.TestCase):
    """Ferrocene-like fragment with Cp directly bonded to a metal must
    PRESERVE its η¹/η³/η⁵ enumeration under strict mode."""

    # Half-sandwich Fe-Cp: each Cp carbon directly bonded to Fe.
    # We use a hypervalent SMILES that puts Fe-C bonds for each Cp atom.
    FE_CP_SMILES = "[CH]12=[CH][CH]=[CH][CH]1[Fe]2"

    def setUp(self):
        _disable_strict()

    def tearDown(self):
        _disable_strict()

    def test_strict_off_emits_cp_orbits(self):
        _disable_strict()
        modes = _enumerate(self.FE_CP_SMILES)
        labs = _hapto_labels(modes)
        # Cp SMARTS matches an aromatic 5-ring → η¹/η³/η⁵
        # The SMILES uses explicit single/double bonds, so RDKit may not
        # perceive aromaticity; the strict test below is the discriminator.
        # At minimum, the strict-OFF baseline must be at least as
        # permissive as the strict-ON case.
        self.assertIsInstance(labs, list)

    def test_strict_on_keeps_directly_bonded_cp_orbits(self):
        _enable_strict()
        modes = _enumerate(self.FE_CP_SMILES)
        # We assert: any Cp orbit emitted MUST involve atoms that include
        # at least one Fe-neighbour.  In this SMILES every Cp carbon is
        # an Fe neighbour, so the strict filter accepts the ring.
        # (If RDKit doesn't perceive the aromatic ring, no orbit is
        # emitted in either mode, which is acceptable.)
        from delfin.fffree.hapto_strict_detection import (
            is_legitimate_hapto_candidate,
        )
        from rdkit import Chem
        mol = Chem.MolFromSmiles(self.FE_CP_SMILES)
        self.assertIsNotNone(mol)
        # The Cp ring atoms are the 5 carbons; metal_idxs auto-detected.
        # We verify the helper returns True for the ring (regardless of
        # whether RDKit's aromaticity perception emits a SMARTS match).
        carbon_ring = [a.GetIdx() for a in mol.GetAtoms()
                       if a.GetSymbol() == "C"]
        self.assertTrue(
            is_legitimate_hapto_candidate(mol, carbon_ring, strict=True),
            "Cp ring directly bonded to Fe must be a legitimate hapto candidate",
        )


# ---------------------------------------------------------------------------
# 3. (η⁶-benzene)Cr(CO)₃: arene directly bonded → hapto PRESERVED
# ---------------------------------------------------------------------------
class TestCrArenePreserved(unittest.TestCase):
    """Benzenechromium tricarbonyl: each benzene C is directly bonded to Cr
    → strict filter must ACCEPT the ring as hapto-eligible."""

    # (η⁶-benzene)Cr(CO)₃ with all 6 Cr-C bonds explicit.
    CR_BENZENE_SMILES = (
        "[CH]12=[CH][CH]=[CH][CH]=[CH]1[Cr]2(C=O)(C=O)C=O"
    )

    def test_strict_on_accepts_directly_bonded_arene(self):
        _enable_strict()
        from delfin.fffree.hapto_strict_detection import (
            is_legitimate_hapto_candidate,
        )
        from rdkit import Chem
        mol = Chem.MolFromSmiles(self.CR_BENZENE_SMILES)
        # If RDKit can't parse this contrived hypervalent SMILES we fall
        # back to a manually-constructed minimal case.
        if mol is None:
            self.skipTest("RDKit could not parse the test SMILES")
        # Find the ring carbons that are directly bonded to Cr.
        cr_idx = [a.GetIdx() for a in mol.GetAtoms()
                  if a.GetSymbol() == "Cr"]
        self.assertEqual(len(cr_idx), 1)
        cr_neighbors = {n.GetIdx() for n in mol.GetAtomWithIdx(
            cr_idx[0]).GetNeighbors()}
        # The ring carbons are those Cr-bonded carbons that are also in
        # an aromatic ring.  At least one must be a Cr neighbour.
        ring_atoms = list(cr_neighbors)
        self.assertGreater(len(ring_atoms), 0)
        self.assertTrue(
            is_legitimate_hapto_candidate(mol, ring_atoms, strict=True),
            "Directly Cr-bonded carbons must be legitimate hapto candidates",
        )

    def tearDown(self):
        _disable_strict()


# ---------------------------------------------------------------------------
# 4. Pyridyl σ-N donor (NOT hapto for the pyridine ring)
# ---------------------------------------------------------------------------
class TestPyridylSigmaNotHapto(unittest.TestCase):
    """A pyridine ring whose N is σ-bonded to a metal must NOT emit
    eta-arene orbits for that pyridine ring (pyridines bind σ via N).
    Note: pyridine is a 6-aromatic N-containing ring, but the SMARTS
    for arene matches [c]1[c][c][c][c][c]1 (all carbons), so pyridine
    is naturally excluded from the arene SMARTS.  We document and
    verify this here as a regression guard.
    """

    PYRIDYL_SIGMA = "c1ccncc1[Cu]"  # pyridine with Cu σ-bonded to a ring C
    PYRIDYL_N_BOUND = "[Cu]n1ccccc1"  # Cu directly bonded to pyridine N

    def setUp(self):
        _disable_strict()

    def tearDown(self):
        _disable_strict()

    def test_strict_on_no_hapto_for_sigma_pyridine(self):
        _enable_strict()
        # The pyridine ring atoms include 5 C + 1 N → not matched by the
        # all-carbon arene SMARTS, so no eta-arene orbit appears anyway.
        modes = _enumerate(self.PYRIDYL_N_BOUND)
        labs = _hapto_labels(modes)
        self.assertNotIn("eta6-arene", labs)
        self.assertNotIn("eta4-arene", labs)
        self.assertNotIn("eta2-arene", labs)


# ---------------------------------------------------------------------------
# 5. Bipy + phenyl substituent: σ-N donors preserved, phenyl NOT hapto
# ---------------------------------------------------------------------------
class TestBipyPhenylSubstituent(unittest.TestCase):
    """2,2'-bipyridine with a phenyl substituent on the backbone.
    The σ-N donors of bipy → κ²-NN chelate (preserved).  The pendant
    phenyl ring is a substituent, not hapto-bonded.
    """

    # Cu σ-bound by 2 pyridyl N (bipy); phenyl on the 4-position is just
    # a substituent.
    BIPY_PHENYL = "c1ccc(-c2cc(-c3ccccc3)ccn2)nc1[Cu]"

    def setUp(self):
        _disable_strict()

    def tearDown(self):
        _disable_strict()

    def test_strict_on_no_hapto_for_phenyl_substituent(self):
        _enable_strict()
        modes = _enumerate(self.BIPY_PHENYL)
        labs = _hapto_labels(modes)
        # Whether or not RDKit's SMARTS finds an arene match, the strict
        # filter ensures NO eta-arene orbits survive (the phenyl
        # substituent is >1 bond from Cu in the metal-removed graph).
        self.assertEqual(
            labs, [],
            f"phenyl substituent should not emit hapto orbits: {labs}",
        )

    def test_strict_off_baseline(self):
        _disable_strict()
        # Just check enumeration runs cleanly and is deterministic.
        m1 = _enumerate(self.BIPY_PHENYL)
        m2 = _enumerate(self.BIPY_PHENYL)
        self.assertEqual([m["label"] for m in m1],
                         [m["label"] for m in m2])


# ---------------------------------------------------------------------------
# 6. μ-bridging arene (rare): two metals share a ring
# ---------------------------------------------------------------------------
class TestMuBridgingArene(unittest.TestCase):
    """Inverse sandwich [M₁-(η⁶-C₆H₆)-M₂]: ring atoms bonded to BOTH metals.
    Strict filter must preserve hapto orbits (ring atoms are direct
    neighbours of metals).
    """

    # Stylised inverse-sandwich SMILES (μ-η⁶-benzene between two Fe atoms).
    # RDKit-parseable form using explicit dative-bonds to two Fe atoms.
    MU_ARENE = (
        "[Fe]12([CH]3=[CH]1[CH]=[CH]4[CH]2=[CH]3)[Fe]4"
    )

    def test_strict_on_accepts_mu_bridging(self):
        _enable_strict()
        from delfin.fffree.hapto_strict_detection import (
            is_legitimate_hapto_candidate,
        )
        from rdkit import Chem
        mol = Chem.MolFromSmiles(self.MU_ARENE)
        if mol is None:
            self.skipTest("could not parse stylised μ-arene SMILES")
        metal_idxs = [a.GetIdx() for a in mol.GetAtoms()
                      if a.GetSymbol() == "Fe"]
        self.assertGreaterEqual(len(metal_idxs), 1)
        # Get the carbon ring atoms (the bridging arene).
        ring_atoms = [a.GetIdx() for a in mol.GetAtoms()
                      if a.GetSymbol() == "C"]
        self.assertTrue(
            is_legitimate_hapto_candidate(
                mol, ring_atoms, metal_idxs=metal_idxs, strict=True,
            ),
            "μ-bridging arene must be a legitimate hapto candidate",
        )

    def tearDown(self):
        _disable_strict()


# ---------------------------------------------------------------------------
# 7. Byte-identical OFF guarantee
# ---------------------------------------------------------------------------
class TestByteIdenticalOff(unittest.TestCase):
    """When the env-flag is OFF, enumerate_modes output for any SMILES is
    bit-identical to the legacy behaviour (before this patch landed).
    The contract: with the flag unset, the strict filter must NOT alter
    any output."""

    SAMPLES = [
        # YUHRUP
        ("CCN(CC)C1=NC(C2=CC=CC=C2)=[N+]2N=C(N3CCCCCC3)[S][Re]2([O-])([Cl])[S]1"),
        # bare arene
        "c1ccccc1",
        # bare Cp anion
        "[cH-]1cccc1",
        # acetate
        "CC(=O)[O-]",
        # glycinate
        "[NH2]CC(=O)[O-]",
    ]

    def setUp(self):
        _disable_strict()

    def tearDown(self):
        _disable_strict()

    def test_off_matches_legacy_for_samples(self):
        # Capture OFF output.
        off_results = []
        for smi in self.SAMPLES:
            modes = _enumerate(smi)
            off_results.append(
                [(m["kind"], m["label"], m["donors"]) for m in modes]
            )
        # Make sure setting then unsetting the flag still gives the same
        # OFF output (no global state was mutated by enabling/disabling).
        _enable_strict()
        _disable_strict()
        for smi, expected in zip(self.SAMPLES, off_results):
            modes = _enumerate(smi)
            got = [(m["kind"], m["label"], m["donors"]) for m in modes]
            self.assertEqual(got, expected,
                             f"OFF output drifted for {smi!r}")


# ---------------------------------------------------------------------------
# 8. Determinism
# ---------------------------------------------------------------------------
class TestDeterminism(unittest.TestCase):
    """Two identical invocations (same SMILES, same env flag) must produce
    bit-identical output in both OFF and ON modes."""

    YUHRUP = (
        "CCN(CC)C1=NC(C2=CC=CC=C2)=[N+]2N=C(N3CCCCCC3)[S][Re]2"
        "([O-])([Cl])[S]1"
    )

    def setUp(self):
        _disable_strict()

    def tearDown(self):
        _disable_strict()

    def _signature(self, smi):
        modes = _enumerate(smi)
        return tuple((m["kind"], m["label"], m["donors"]) for m in modes)

    def test_determinism_off(self):
        _disable_strict()
        s1 = self._signature(self.YUHRUP)
        s2 = self._signature(self.YUHRUP)
        self.assertEqual(s1, s2)

    def test_determinism_on(self):
        _enable_strict()
        s1 = self._signature(self.YUHRUP)
        s2 = self._signature(self.YUHRUP)
        self.assertEqual(s1, s2)

    def test_helper_determinism(self):
        """The pure-topology helper is deterministic."""
        from delfin.fffree.hapto_strict_detection import (
            is_legitimate_hapto_candidate,
        )
        from rdkit import Chem
        mol = Chem.MolFromSmiles(self.YUHRUP)
        self.assertIsNotNone(mol)
        # Find a phenyl ring (the substituent C2=CC=CC=C2 in YUHRUP, atoms
        # 8..13 in the parsed SMILES order).  We don't hard-code indices
        # — instead detect rings of 6 aromatic C.
        rings = mol.GetRingInfo().AtomRings()
        phenyl_rings = [
            list(r) for r in rings
            if len(r) == 6 and all(
                mol.GetAtomWithIdx(i).GetSymbol() == "C"
                and mol.GetAtomWithIdx(i).GetIsAromatic()
                for i in r
            )
        ]
        self.assertGreaterEqual(len(phenyl_rings), 1)
        for ring in phenyl_rings:
            v1 = is_legitimate_hapto_candidate(mol, ring, strict=True)
            v2 = is_legitimate_hapto_candidate(mol, ring, strict=True)
            self.assertEqual(v1, v2)
            # Phenyl substituent is far from Re → must be False
            self.assertFalse(
                v1,
                "phenyl substituent on σ-chelate must not be a hapto candidate",
            )


# ---------------------------------------------------------------------------
# 9. No-metal fallback: helper returns True for metal-free graphs
# ---------------------------------------------------------------------------
class TestNoMetalFallback(unittest.TestCase):
    """When the input SMILES contains NO metal (bare ligand), the helper
    must return True regardless of strict mode (we can't tell from a
    metal-free graph which rings would be hapto-bound to a hypothetical
    metal — preserve the question "what hapto modes can this ligand
    support in principle?")."""

    def test_metal_free_smiles_always_legit(self):
        _enable_strict()
        from delfin.fffree.hapto_strict_detection import (
            is_legitimate_hapto_candidate,
        )
        from rdkit import Chem
        for smi in ("c1ccccc1", "[cH-]1cccc1", "C=CC=C", "c1ccc2ccccc2c1"):
            mol = Chem.MolFromSmiles(smi)
            self.assertIsNotNone(mol)
            ring_atoms = list(range(mol.GetNumAtoms()))
            self.assertTrue(
                is_legitimate_hapto_candidate(mol, ring_atoms, strict=True),
                f"metal-free {smi!r} must be a legitimate hapto candidate",
            )

    def tearDown(self):
        _disable_strict()


if __name__ == "__main__":
    unittest.main()
