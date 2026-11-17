"""Proton coordinator package.

Analogue of :mod:`delfin.co2` for studying the competing hydrogen-evolution
reaction (HER) at a metal catalyst: a proton donor molecule (never a bare
proton) is placed next to the metal centre, its donor H points at the metal,
and an orientation scan followed by a relaxed M...H distance scan maps the
proton/hydrogen-transfer path.

Donor chemistry contract (see chain_setup.bookkeeping):
- protic donors (TFE, HOAc, phenol, H3O(H2O)n): proton transfer,
  products = [M-H] + neutral-or-anionic conjugate base.
- HAT donors (BIH): hydrogen atom transfer, products = [M-H] + donor
  radical (open-shell partner).
"""
