---
title: GFN2-xTB ring-closure artefact on the isopropyl cation
topic: chemistry-methods
scope: delfin-development
provenance: delfin-agent | 2026-10-24 | doc:delfin/tools/adapters/xtb.py
review_by: 2027-12-31
---
Situation: A GFN2-xTB geometry optimization or conformer ranking of the isopropyl cation converges to a three-membered ring-like structure instead of the open secondary carbocation, and the error looks like a convergence failure.
Lesson: The collapsed structure is a known GFN2-xTB artefact of the isopropyl cation, not a failed optimization: do not loosen convergence criteria to chase the open form; check the structure, and prefer a higher-level single point on an xTB conformer ensemble for this species.
Why: Tightening the optimizer will reproduce the same minimum; the artefact is in the semi-empirical Hamiltonian's preferences for this cation, so the fix is method choice, not convergence settings.
