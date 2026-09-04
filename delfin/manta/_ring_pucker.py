"""Deterministic ring-pucker conformer CONSTRUCTION (Cremer-Pople + torsion-
constrained relaxation).

ETKDG samples only the global-minimum ring pucker — cyclohexane embeds
300/300 as the CHAIR, never the twist-boat/boat; a 5-ring never leaves its
lowest envelope.  The higher ring basins are genuine, populated, distinct
conformers that a COMPLETE manifold must contain, but no amount of seed
sampling reaches them.  They have to be CONSTRUCTED explicitly.

Mechanism (license-clean, no CCDC, deterministic):
  1. Displace the ring atoms onto a target Cremer-Pople puckering coordinate
     (Q, theta, phi) — chair / twist-boat / boat for a 6-ring, twist / envelope
     for a 5-ring — keeping their in-plane positions, moving only along the
     ring's mean-plane normal.
  2. Relax with UFF while HOLDING the pucker: every ring torsion is restrained
     to its just-set value +/- a window, so bond lengths and angles relax to
     physical values (C-C -> ~1.53 A, tetrahedral angles) but the pucker is
     preserved instead of collapsing back to the global-minimum chair.
  3. Keep the result iff its ring pucker is genuinely distinct from every
     conformer already in the pool (Cremer-Pople theta/phi separation) and the
     geometry is physical.

The primitive takes an optional ``frozen`` atom set (metal + coordinating
donors) so the SAME construction serves metal-containing chelate rings — the
coordination sphere stays put while the backbone puckers.  Metal-free rings
pass ``frozen=None``.
"""

import os as _os
from typing import List, Optional, Set, Tuple

try:
    import numpy as _np
except Exception:                                    # pragma: no cover
    _np = None

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms
    _RDKIT = True
except Exception:                                    # pragma: no cover
    _RDKIT = False


# --- GENERAL Cremer-Pople pucker candidates for ANY ring size ---------------
# A ring of size N has (N-3) puckering degrees of freedom.  Rather than hard-code
# named conformers per size (chair/boat/... only make sense for N=6), we sample
# the Cremer-Pople pucker sphere GENERICALLY: the polar caps (even N -> the
# chair-family "alternating" pucker) plus K phases around the m=2 pseudorotation
# circle (the low-energy boat/twist/envelope family).  A torsion-held relax then
# lets each candidate fall into the nearest genuine minimum, and a Cremer-Pople
# dedup distils the candidates to the DISTINCT populated conformers for THIS
# ring.  This works uniformly for N = 5, 6, 7, 8, 9, ... with no per-size table.
# Ringe, die NUR den flachen Zustand bekommen (konjugierte Metallacyclen, die
# `_is_puckerable` ablehnt).  Wird je `generate`-Aufruf neu befuellt; ein Ring
# hierin erhaelt in `_ring_pucker_states` AUSSCHLIESSLICH den Q=0-Kandidaten,
# damit ein konjugierter Ring begradigt und nicht gefaltet wird.
_FLAT_ONLY: set = set()


def _pucker_candidates(n: int) -> List[Tuple[float, Optional[float], float]]:
    """(q_scale, theta, phi) je Kandidat.  `q_scale` multipliziert `_amp(n)`.

    ⚠ WARUM DIE TUPEL JETZT DREI WERTE HABEN (17.08.2026).  Bis hierher tastete diese
    Funktion die Cremer-Pople-Kugel bei FESTEM RADIUS ab: `_amp(n)` gibt 0.40 A (5-Ring)
    bis 0.80 A (8-Ring) und wurde an beiden Aufrufstellen unveraendert uebergeben.  Der
    Kandidat `(0.0, 0.0)` unten ist **theta = 0**, also der SESSEL (Polkappe) -- **nicht**
    Q = 0.  **Der Mittelpunkt der Kugel, die EBENE, war kein Kandidat.**  Fuer ungerade
    Ringe war es noch enger: nur die aequatoriale Pseudorotation (`theta=None`), also
    Umschlag und Twist, nie flach.

    GEMESSEN am 16./17.08. (`folds`, 965 Systeme): von 326 fehlenden Ringmotiven sind
    **218 PLANAR** -- `5M:planar` 126, `6M:planar` 55, `4M:planar` 37 -- gegen `6M:boat` 30,
    `6:chair` 25, `5M:puckered` 19, die dieser Generator alle erzeugen kann.  **Zwei Drittel
    der Luecke sind genau der eine Zustand, den er per Konstruktion nicht kennt.**

    Chemisch ist das kein Randfall: ein fuenfgliedriger Chelatring mit sp2-Donoren liegt oft
    FLACH; der Generator behandelt ihn wie Cyclopentan.

    `_set_pucker` braucht dafuer KEINE Aenderung: mit Q = 0 werden q2 und q3 null, zj = 0,
    und jedes nicht eingefrorene Ringatom wird auf die Mittelebene projiziert -- Metall und
    Donoren bleiben stehen, weil sie in `frozen` sind.  Das ist exakt der planare Zustand.

    Vorgabe AUS -> die Liste ist identisch zu vorher (alle q_scale = 1.0), also
    byte-identisch.
    """
    cands: List[Tuple[float, Optional[float], float]] = []
    even = (n % 2 == 0)
    # equatorial pseudorotation ring — sample fine enough to hit both the boat
    # (phi = 0, 360/n, ...) and the twist (phi halfway between) positions.
    K = max(8, 2 * n)
    for k in range(K):
        phi = 360.0 * k / K
        cands.append((1.0, 90.0 if even else None, phi))
    if even:
        # polar caps: the chair / inverted-chair (alternating) puckers
        cands.append((1.0, 0.0, 0.0))
        cands.append((1.0, 180.0, 0.0))
    if _os.environ.get("DELFIN_FFFREE_PUCKER_PLANAR", "0") == "1":
        # DER FLACHE ZUSTAND.  Q = 0 -> alle nicht eingefrorenen Ringatome in die
        # Mittelebene.  EIN Kandidat je Ring, nicht K -- die Ebene hat kein phi.
        cands.append((0.0, 0.0, 0.0))
    return cands


def _amp(n: int) -> float:
    """Cremer-Pople-Faltungsamplitude (Angstrom) als GESETZ statt Tabelle.

    ⚠ DIE ALTE FASSUNG WIDERSPRACH SICH SELBST.  Sie fuehrte eine kalibrierte Tabelle
    fuer 5..8 UND einen linearen Fallback `0.45 + 0.06*n` fuer alles andere -- und der
    Fallback liegt UEBERALL ueber der Tabelle:

        n     Tabelle   Fallback
        5      0,40      0,75      <- fast doppelt
        6      0,63      0,81
        7      0,72      0,87
        8      0,80      0,93

    Aufgefallen ist das nie, weil `_is_puckerable` jeden Ring ausserhalb 5..8 abwies:
    der Fallback ist NIE GELAUFEN.  Ein dunkler Zweig mit einer falschen Zahl darin.

    Die Tabellenwerte saettigen (Zuwaechse 0,23 / 0,09 / 0,08) -- ein grosser Ring
    faltet nicht beliebig tief, die Amplitude laeuft gegen eine Schranke.  Ein
    linearer Fallback ist damit qualitativ falsch, nicht nur numerisch daneben.

    Ersatz: EIN saettigendes Gesetz fuer alle N, das die kalibrierten Werte
    reproduziert (max. Abweichung 0,04 A):

        Q_max(N) = 1,15 * (N-4) / (N-4+1,6)
        N=5 0,44 · N=6 0,64 · N=7 0,75 · N=8 0,82 · N=12 0,98 · N=24 1,07

    ⚠ Ab N=9 ist das EXTRAPOLATION, keine Kalibrierung -- ehrlich gesagt, nicht
      versteckt.  Und es ist ohnehin nur die OBERGRENZE: `_pucker_space_grid` tastet
      die Amplitude von 0 bis hierher ab, der Relax entscheidet, was ueberlebt.

    ⛔ DIE KALIBRIERTEN WERTE BLEIBEN EXAKT STEHEN.  `_amp` wird auch vom LEGACY-Pfad
      gelesen (`_set_pucker(..., _qs * _amp(n), ...)`).  Wuerde das Gesetz sie
      ersetzen, waere die Vorgabe NICHT byte-identisch -- bei N=5 stuende 0,44 statt
      0,40.  Das Gesetz greift darum nur dort, wo bisher der falsche Fallback stand:
      ausserhalb 5..8.  Byte-Identitaet ist keine Formsache, sie ist die Bedingung
      dafuer, dass ein A/B den Mechanismus misst und nicht das Instrument.
    """
    _kal = {5: 0.40, 6: 0.63, 7: 0.72, 8: 0.80}
    if n in _kal:
        return _kal[n]
    if n <= 4:
        return 0.35
    return 1.15 * (n - 4.0) / (n - 4.0 + 1.6)


def _ring_normal_and_center(P, ring):
    n = len(ring)
    R = P[list(ring)]
    C = R.mean(0)
    Rc = R - C
    Rp = sum(Rc[j] * _np.sin(2 * _np.pi * j / n) for j in range(n))
    Rpp = sum(Rc[j] * _np.cos(2 * _np.pi * j / n) for j in range(n))
    nrm = _np.cross(Rp, Rpp)
    ln = float(_np.linalg.norm(nrm))
    if ln < 1e-9:
        return None, C
    return nrm / ln, C


def _cp_theta_phi(P, ring):
    """Cremer-Pople (Q, theta, phi) for the ring (atom coords in ring order)."""
    n = len(ring)
    nrm, C = _ring_normal_and_center(P, ring)
    if nrm is None:
        return 0.0, 0.0, 0.0
    z = (P[list(ring)] - C) @ nrm
    q2c = _np.sqrt(2.0 / n) * sum(z[j] * _np.cos(2 * _np.pi * 2 * j / n) for j in range(n))
    q2s = -_np.sqrt(2.0 / n) * sum(z[j] * _np.sin(2 * _np.pi * 2 * j / n) for j in range(n))
    q2 = float(_np.hypot(q2c, q2s))
    if n % 2 == 0:
        q3 = _np.sqrt(1.0 / n) * sum(z[j] * ((-1) ** j) for j in range(n))
    else:
        q3 = 0.0
    Q = float(_np.sqrt(q2 ** 2 + q3 ** 2))
    theta = float(_np.degrees(_np.arctan2(q2, q3))) if n % 2 == 0 else 90.0
    phi = float(_np.degrees(_np.arctan2(q2s, q2c)) % 360.0)
    return Q, theta, phi


def _set_pucker(conf, ring, Q, theta, phi, frozen: Optional[Set[int]] = None):
    n = len(ring)
    P = conf.GetPositions()
    nrm, C = _ring_normal_and_center(P, ring)
    if nrm is None:
        return
    ph = _np.radians(phi)
    even = (n % 2 == 0)
    th = _np.radians(theta) if (theta is not None) else None
    q2 = Q * (_np.sin(th) if th is not None else 1.0)
    q3 = Q * (_np.cos(th) if th is not None else 0.0)
    frozen = frozen or set()
    for j, idx in enumerate(ring):
        # a frozen ring atom (metal / coordinating donor of a chelate ring) keeps
        # its position -> only the backbone puckers, the coordination sphere is
        # preserved.
        if int(idx) in frozen:
            continue
        zj = _np.sqrt(2.0 / n) * q2 * _np.cos(ph + 2 * _np.pi * 2 * j / n)
        if even and th is not None:
            zj += _np.sqrt(1.0 / n) * q3 * ((-1) ** j)
        p = P[idx]
        inplane = p - ((p - C) @ nrm) * nrm
        newp = inplane + zj * nrm
        conf.SetAtomPosition(int(idx), (float(newp[0]), float(newp[1]), float(newp[2])))


def _prod_laenge(per_ring_states) -> int:
    """Groesse des vollen Kreuzprodukts (inkl. Grundzustand)."""
    n = 1
    for s in per_ring_states:
        n *= max(1, len(s))
    return n


def _ring_bahnen(mol, rings, max_aut: int = 20000):
    """Zerlege die Ringe in BAHNEN unter der Automorphismengruppe des Molekuels.

    Zwei Ringe liegen in derselben Bahn, wenn ein Automorphismus den einen als
    ATOMMENGE auf den anderen abbildet.  Nur dann sind sie ununterscheidbar, und
    nur dann darf ihre Zustandsreihenfolge zusammengelegt werden.

    ⛔ WARUM NICHT DIE RANGMULTIMENGE, die viel billiger waere.  Gleiche
    kanonische Raenge sind NOTWENDIG fuer Aequivalenz, aber nicht HINREICHEND.
    Am 27.08. von Hand nachgerechnet: von 27 Ringpaaren mit gleicher Rangmenge
    waren 26 echte Bahn und EINES nicht (3,7 %).  Ueber die Rangmenge zu
    reduzieren haette bei diesem einen Paar einen REALEN Zustand geloescht.  Ein
    Doppelgaenger zuviel kostet Rechenzeit; ein fehlender Zustand kostet
    Vollstaendigkeit -- und die ist der Nordstern.

    ⚠ DECKEL.  `GetSubstructMatches(mol, mol)` kann bei hochsymmetrischen
    Molekuelen explodieren (am 27.08. lief CONPUS in 200 000 Treffer).  Wird der
    Deckel erreicht, liefert die Funktion None -> KEINE Reduktion, volles
    Produkt.  Der Rueckfall ist immer die GROESSERE Menge, nie die kleinere.

    Rueckgabe: Liste von Listen von Ringindizes, oder None (nicht reduzieren).
    """
    try:
        treffer = mol.GetSubstructMatches(mol, uniquify=False,
                                          useChirality=False,
                                          maxMatches=max_aut)
    except Exception:
        return None
    if not treffer or len(treffer) >= max_aut:
        return None                      # Deckel erreicht -> nicht reduzieren
    if len(treffer) == 1:
        return None                      # nur die Identitaet -> nichts zu holen

    ring_mengen = [frozenset(int(i) for i in r) for r in rings]
    index_von = {m: i for i, m in enumerate(ring_mengen)}
    if len(index_von) != len(ring_mengen):
        return None                      # doppelte Ringmengen -> Haende weg

    eltern = list(range(len(rings)))

    def _wurzel(x):
        while eltern[x] != x:
            eltern[x] = eltern[eltern[x]]
            x = eltern[x]
        return x

    for abb in treffer:
        for i, menge in enumerate(ring_mengen):
            try:
                bild = frozenset(int(abb[a]) for a in menge)
            except Exception:
                return None
            j = index_von.get(bild)
            if j is None or j == i:
                continue
            ri, rj = _wurzel(i), _wurzel(j)
            if ri != rj:
                eltern[ri] = rj

    gruppen = {}
    for i in range(len(rings)):
        gruppen.setdefault(_wurzel(i), []).append(i)
    # Deterministische Reihenfolge -- sonst haengt die Aufzaehlung an der
    # Hash-Reihenfolge und der Lauf ist nicht reproduzierbar.
    return [sorted(v) for _k, v in sorted(gruppen.items())]


def _cp_abstand(a, b) -> float:
    """Winkelabstand zweier Faltungszustaende AUF der Cremer-Pople-Kugel (Grad).

    a, b sind (Q, theta, phi).  Benutzt wird die Grosskreisdistanz

        cos d = cos(th_a) cos(th_b) + sin(th_a) sin(th_b) cos(ph_a - ph_b)

    ⚠ WARUM NICHT EINFACH |dtheta| + |dphi|.  Am POL (theta = 0 oder 180) ist phi
      BEDEUTUNGSLOS -- ein Ring im Sessel hat keine Phase.  Eine naive Metrik haelt
      zwei Sessel mit phi = 136 und phi = 339 fuer 200 Grad auseinander, obwohl sie
      DERSELBE Zustand sind.  Genau das steht in der Messung vom 26.08.:
          theta=180,0  phi=136,4
          theta=180,0  phi=338,6      <- identisch, nur die Phase ist Rauschen
      Die Grosskreisdistanz erledigt das von selbst: bei sin(theta) = 0 faellt der
      phi-Term heraus.  Die Geometrie loest das Problem, nicht eine Sonderregel.
    """
    ta, pa = _np.radians(a[1]), _np.radians(a[2])
    tb, pb = _np.radians(b[1]), _np.radians(b[2])
    c = (_np.cos(ta) * _np.cos(tb)
         + _np.sin(ta) * _np.sin(tb) * _np.cos(pa - pb))
    return float(_np.degrees(_np.arccos(max(-1.0, min(1.0, float(c))))))


def _set_pucker_general(conf, ring, qs, phis, frozen: Optional[Set[int]] = None):
    """Cremer-Pople-Umkehr in VOLLER Allgemeinheit -- fuer JEDE Ringgroesse.

    ``_set_pucker`` oben deckt nur m = 2 plus den Alternierungsterm ab.  Das ist fuer
    N = 4, 5, 6 vollstaendig und ab N = 7 LUECKENHAFT: ein Siebenring hat vier
    Faltungsfreiheitsgrade (q2, phi2, q3, phi3), ein Achtring fuenf.  Die Paare mit
    m >= 3 fehlten dort ersatzlos, und der Alternierungsterm wurde fest als q3
    gefuehrt -- beim Achtring ist es aber q4.

    DER RAUM, exakt.  Ein Ring mit N Atomen hat N-3 Faltungsfreiheitsgrade:
        N gerade:  Paare (q_m, phi_m) fuer m = 2 .. N/2-1,  plus EIN q_(N/2)
                   2*(N/2-2) + 1 = N-3
        N ungerade: Paare (q_m, phi_m) fuer m = 2 .. (N-1)/2
                   2*((N-1)/2 - 1) = N-3
    Die Auslenkung des j-ten Ringatoms aus der Mittelebene ist

        z_j = sqrt(2/N) * SUM_m  q_m * cos(phi_m + 2*pi*m*j/N)
              + [N gerade]  sqrt(1/N) * q_(N/2) * (-1)^j

    Die benannten Formen sind PUNKTE darauf, keine eigenen Faelle: Sessel an den Polen
    (q2 = 0), Wanne und Twist am Aequator (q3 = 0), Half-Chair und Envelope
    DAZWISCHEN -- genau der Bereich, den die alte Kandidatenliste nie abgetastet hat.

    ``qs``/``phis``: Abbildungen m -> Wert.  Reduziert sich fuer N <= 6 exakt auf
    ``_set_pucker``; die Formel ist dieselbe, nur nicht mehr auf m = 2 verkuerzt.
    ⚠ ``frozen`` bleibt unberuehrt -- Metall und Donoren stehen, nur das Rueckgrat
    faltet.  Reiner Erzeuger.
    """
    n = len(ring)
    P = conf.GetPositions()
    nrm, C = _ring_normal_and_center(P, ring)
    if nrm is None:
        return
    frozen = frozen or set()
    even = (n % 2 == 0)
    m_last = n // 2 if even else None
    for j, idx in enumerate(ring):
        if int(idx) in frozen:
            continue
        zj = 0.0
        for m, qm in qs.items():
            if not qm:
                continue
            if even and m == m_last:
                zj += _np.sqrt(1.0 / n) * qm * ((-1) ** j)
            else:
                ph = _np.radians(phis.get(m, 0.0))
                zj += _np.sqrt(2.0 / n) * qm * _np.cos(ph + 2.0 * _np.pi * m * j / n)
        p = P[idx]
        inplane = p - ((p - C) @ nrm) * nrm
        newp = inplane + zj * nrm
        conf.SetAtomPosition(int(idx), (float(newp[0]), float(newp[1]), float(newp[2])))


def _pucker_space_grid(n: int, n_amp: int, n_phase: int):
    """SYSTEMATISCHES Gitter ueber den GANZEN Faltungsraum eines N-Rings.

    Liefert Kandidaten als ``(qs, phis)`` -- Abbildungen m -> Wert -- fuer
    ``_set_pucker_general``.  Statt benannter Formen wird der (N-3)-dimensionale
    Cremer-Pople-Raum abgetastet; Sessel, Wanne, Twist, Half-Chair und Envelope
    fallen als Gitterpunkte von selbst an.

    ⚠ VOLLSTAENDIGKEIT IST EINE AUFLOESUNGSFRAGE, keine Ja/Nein-Frage.  Ein
    kontinuierlicher Raum laesst sich nicht "ganz" abtasten.  Was hier steht, ist die
    ehrliche Fassung: der Raum wird VOLLSTAENDIG bei der ANGEGEBENEN Aufloesung
    ueberdeckt, und die Aufloesung steht in der Spur.  Keine Ecke wird ausgelassen,
    keine Richtung bevorzugt -- der Unterschied zur alten Liste, die nur Aequator und
    Pole kannte und die Amplitude nie variierte.

    ⚠ PREIS: die Kandidatenzahl waechst wie (n_amp+1)^(#q) * n_phase^(#phi).
    Sechsring bei n_amp=2, n_phase=8: 3 * 8 * 5 = 120 je Ring.  Achtring: deutlich
    mehr.  Darum sind beide Aufloesungen Env-Parameter und stehen im Protokoll.

    ===== DIE AUFLOESUNG MUSS MIT DER DIMENSION FALLEN (26.08.2026) ================

    GERECHNET, nicht geschaetzt.  Die Zahl der Kandidaten ist

        prod ueber m_pairs von (1 + n_amp * n_phase)   mal   (2*n_amp+1) bei geradem n

    und damit bei n_amp=2, n_phase=8:

        n= 8    845          n=10   24 565        n=12     417 605
        n= 9  4 913          n=11   83 521        n=16  120 687 845

    Ein 16-Ring haette also 1,2e8 Kandidaten je Ring bekommen, jeden mit Relax und
    Kollisionstor.  Das ist kein langsamer Lauf, das ist ein Lauf, der stirbt und
    NULL Faltungen liefert.  Genau so war `foldspace6k` eingereiht (NAMP=2 NPHASE=8).
    ⇒ Unendliche Feinheit ist nicht Vollstaendigkeit, sie ist Undurchfuehrbarkeit.

    WAS HIER **NICHT** PASSIERT: es wird kein Ergebnis abgeschnitten.  Die Liste
    bleibt das VOLLSTAENDIGE Produkt der gewaehlten Aufloesung -- reduziert wird die
    ABTASTDICHTE, und zwar zuerst dort, wo sie physikalisch am wenigsten traegt.
    Cremer-Pople-Amplituden q_m fallen mit m: die hohen m sind die feine Kraeuselung
    mit kleiner Auslenkung, m=2 ist die dominante Falte.  Darum wird die Phasenzahl
    beim GROESSTEN m zuerst halbiert und m=2 zuletzt angetastet.

    ⚠ Und es geschieht NICHT still: `_grid_res` traegt die je-m gewaehlte Aufloesung,
      der Aufrufer schreibt sie unter DELFIN_FFFREE_PUCKER_TRACE ins Protokoll.  Ob
      die gewaehlte Dichte reicht, sagt nicht dieser Code, sondern
      `selbsttest_konvergenz` -- die Zahl der UNTERSCHEIDBAREN Zustaende, nicht die
      Zahl der Gitterpunkte, ist das Mass.
    """
    if n < 4:
        return []
    even = (n % 2 == 0)
    m_pairs = list(range(2, (n // 2) if even else ((n - 1) // 2) + 1))
    m_last = (n // 2) if even else None
    amp = _amp(n)
    _budget = max(1, int(_os.environ.get("DELFIN_FFFREE_PUCKER_BUDGET", "50000") or 50000))
    # ===== WELCHE MODEN SIND UEBERHAUPT ANGEREGT?  (26.08.2026) =====================
    #
    # Der Selbsttest hat die Sparstelle selbst gefunden: beim 21-Ring (Kronenether,
    # z.B. VEDCOA) reichte das Budget nur, wenn auch m=2 heruntergerechnet wurde --
    # und m=2 ist die DOMINANTE Falte.  Am falschen Ende gespart.
    #
    # Die Ursache ist nicht die Phasenzahl, sondern die ZAHL DER PAARE: sie waechst
    # wie N/2, und schon die blosse Amplitudenauswahl kostet 3^(N/2-1).  Bei N=30
    # sind das 8 Millionen Punkte, BEVOR eine einzige Phase abgetastet ist.
    #
    # Cremer-Pople-Amplituden realer Ringe fallen scharf mit m: die niedrigen Moden
    # tragen die Faltung, die hohen sind feine Kraeuselung nahe null.  Grosse Ringe
    # werden in der Literatur genau darum durch wenige niedrige Moden beschrieben.
    # ⇒ Moden oberhalb M_MAX werden auf Amplitude 0 gesetzt -- die Aussage ist
    #   "diese Mode ist NICHT ANGEREGT", nicht "diese Mode wurde uebersprungen".
    #   Der Freiheitsgrad bleibt in der Parametrisierung, er steht nur auf null.
    #
    # ⚠ DAS IST EINE MODELLANNAHME, KEINE MESSUNG.  Sie ist pruefbar und MUSS geprueft
    #   werden: M_MAX erhoehen und `selbsttest_konvergenz` fragen, ob die Zahl der
    #   UNTERSCHEIDBAREN Zustaende sich aendert.  Aendert sie sich, ist M_MAX zu klein.
    #   Bis dahin steht sie in der Spur und traegt ihren Namen.
    _mmax = max(2, int(_os.environ.get("DELFIN_FFFREE_PUCKER_MMAX", "4") or 4))
    _aktiv = [m for m in m_pairs if m <= _mmax]
    _ruhend = [m for m in m_pairs if m > _mmax]
    # Phasenzahl je m; Start ueberall gleich, dann von oben herunter halbieren.
    _ph_m = {m: max(1, int(n_phase)) for m in _aktiv}
    m_pairs = _aktiv

    def _zahl():
        t = (2 * n_amp + 1) if m_last is not None else 1
        for _m in m_pairs:
            t *= (1 + n_amp * _ph_m[_m])
        return t

    while _zahl() > _budget:
        _kand = [m for m in m_pairs if _ph_m[m] > 1]
        if not _kand:
            break                                  # schon bei Phase 1 -- nichts mehr zu holen
        _hoch = max(_kand)                         # groesstes m zuerst: kleinste Amplitude
        _ph_m[_hoch] = max(1, _ph_m[_hoch] // 2)
    # ⚠ `m_last` (der Alternierungsterm q_{N/2}) bleibt IMMER aktiv und wird nie
    #   ruhend gestellt: beim Sechsring IST er der Sessel.  Er kostet auch nichts --
    #   ein einzelner signierter Amplitudenfaktor (2*n_amp+1), nicht exponentiell.
    _pucker_space_grid._grid_res = {"n": n, "n_amp": n_amp, "n_phase_je_m": dict(_ph_m),
                                    "kandidaten": _zahl(), "budget": _budget,
                                    "m_max": _mmax, "ruhende_moden": list(_ruhend),
                                    "m_last": m_last,
                                    "reduziert": (any(v < n_phase for v in _ph_m.values())
                                                  or bool(_ruhend))}
    # Amplitudenstufen je Paar: 0 (Achse flach) bis n_amp * amp.  Die Null MUSS dabei
    # sein -- sie ist der planare Zustand, und genau der fehlte (218 von 326 Motiven).
    lv_pair = [amp * k / max(1, n_amp) for k in range(0, n_amp + 1)]
    # Der Alternierungsterm laeuft SIGNIERT: +q ist der Sessel, -q der invertierte.
    lv_last = [amp * k / max(1, n_amp) for k in range(-n_amp, n_amp + 1)]
    # Phasen JE m -- gleiche Formel, nur mit der fuer dieses m gewaehlten Dichte.
    _phasen = {m: [360.0 * k / _ph_m[m] for k in range(_ph_m[m])] for m in m_pairs}

    out = []

    def _rek(i, qs, phis):
        if i < len(m_pairs):
            m = m_pairs[i]
            for q in lv_pair:
                if q == 0.0:                      # Amplitude 0 -> Phase bedeutungslos
                    _rek(i + 1, {**qs, m: 0.0}, {**phis, m: 0.0})
                else:
                    for ph in _phasen[m]:
                        _rek(i + 1, {**qs, m: q}, {**phis, m: ph})
            return
        if m_last is not None:
            for q in lv_last:
                out.append(({**qs, m_last: q}, dict(phis)))
        else:
            out.append((dict(qs), dict(phis)))

    _rek(0, {}, {})
    # den Nullpunkt (alles flach) genau EINMAL behalten -- er ist der planare Zustand
    _seen = set()
    uniq = []
    for qs, phis in out:
        key = tuple(sorted((m, round(q, 6), round(phis.get(m, 0.0), 3) if q else 0.0)
                           for m, q in qs.items()))
        if key in _seen:
            continue
        _seen.add(key)
        uniq.append((qs, phis))
    return uniq


def _relax_hold_pucker(mol, ring, frozen: Set[int], window: float = 18.0, iters: int = 1200) -> bool:
    """UFF-relax that frees bonds+angles but HOLDS the pucker: each ring torsion
    restrained to its current value +/- ``window``; any ``frozen`` atom fixed."""
    conf = mol.GetConformer()
    n = len(ring)
    try:
        ff = AllChem.UFFGetMoleculeForceField(mol)
    except Exception:
        ff = None
    if ff is None:
        return False
    for idx in (frozen or ()):
        try:
            ff.AddFixedPoint(int(idx))
        except Exception:
            pass
    for i in range(n):
        a, b, c, d = ring[i], ring[(i + 1) % n], ring[(i + 2) % n], ring[(i + 3) % n]
        try:
            t = rdMolTransforms.GetDihedralDeg(conf, a, b, c, d)
            ff.UFFAddTorsionConstraint(a, b, c, d, False, t - window, t + window, 200.0)
        except Exception:
            pass
    try:
        ff.Minimize(maxIts=iters)
    except Exception:
        return False
    return True


def _relax_hold_pucker_multi(mol, rings, frozen: Set[int], window: float = 18.0,
                             iters: int = 1500) -> bool:
    """UFF relax holding the pucker of EVERY ring simultaneously (each ring
    torsion restrained to its current value +/- window); frozen atoms fixed."""
    conf = mol.GetConformer()
    try:
        ff = AllChem.UFFGetMoleculeForceField(mol)
    except Exception:
        ff = None
    if ff is None:
        return False
    for idx in (frozen or ()):
        try:
            ff.AddFixedPoint(int(idx))
        except Exception:
            pass
    for ring in rings:
        n = len(ring)
        for i in range(n):
            a, b, c, d = ring[i], ring[(i + 1) % n], ring[(i + 2) % n], ring[(i + 3) % n]
            try:
                t = rdMolTransforms.GetDihedralDeg(conf, a, b, c, d)
                ff.UFFAddTorsionConstraint(a, b, c, d, False, t - window, t + window, 200.0)
            except Exception:
                pass
    try:
        ff.Minimize(maxIts=iters)
    except Exception:
        return False
    return True


_VDW = {"H": 1.10, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
        "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98, "B": 1.92, "Si": 2.10}


def _has_bad_angles(mol, tol: float = 25.0, skip: Optional[Set[int]] = None) -> bool:
    """True if any heavy centre's VSEPR angle is off its hybridisation ideal by
    more than ``tol`` — i.e. the (multi-)ring pucker left a LOCAL geometry broken
    even though nothing clashes.  Essential for FUSED / BRIDGED ring systems
    (ACEQAC-type lactams): puckering one ring independently strains the shared
    fusion atoms, distorting their angles; the clash gate is blind to it.  A
    frame is realistic only if EVERY VSEPR body is correct, so any such pucker is
    rejected.  2-coordinate centres are hybridisation-ambiguous (sp/sp2/sp3) ->
    skipped; >=5 is non-molecular -> skipped.

    ``skip``: centres exempt from the VSEPR ideal, because they do not HAVE one.
    A COORDINATION centre is the case this exists for: its angles are set by the
    polyhedron, not by hybridisation.  A CN4 metal has nh == 4, so this function
    would demand 109.5 deg of it -- and a square-planar d8 has two 180 deg trans
    angles, i.e. a 70.5 deg "error" that no pucker caused and no pucker can fix.
    Without the exemption EVERY combination is rejected for every SP-4 and T-3
    complex, which reads as "the lever has no reach" for entirely the wrong reason.
    Default None -> empty -> byte-identical for every existing caller."""
    _skip = skip or frozenset()
    try:
        conf = mol.GetConformer()
        P = conf.GetPositions()
        for c in range(mol.GetNumAtoms()):
            if c in _skip:
                continue
            a = mol.GetAtomWithIdx(c)
            if a.GetSymbol() == "H":
                continue
            hv = [nb.GetIdx() for nb in a.GetNeighbors() if nb.GetSymbol() != "H"]
            nh = len(hv)
            if nh < 3 or nh >= 5:
                continue
            exp = 120.0 if nh == 3 else 109.5
            nbset = {x.GetIdx() for x in a.GetNeighbors()}
            for i in range(len(hv)):
                for j in range(i + 1, len(hv)):
                    # skip 3-membered rings (real ~60deg geometry)
                    if hv[j] in {x.GetIdx() for x in mol.GetAtomWithIdx(hv[i]).GetNeighbors()}:
                        continue
                    v1 = P[hv[i]] - P[c]
                    v2 = P[hv[j]] - P[c]
                    d = float(_np.linalg.norm(v1) * _np.linalg.norm(v2))
                    if d < 1e-9:
                        continue
                    ang = _np.degrees(_np.arccos(max(-1.0, min(1.0, float(_np.dot(v1, v2) / d)))))
                    if abs(ang - exp) > tol:
                        return True
    except Exception:
        return False
    return False


def _has_clash(mol, frac: float = 0.60) -> bool:
    """True if any non-bonded heavy-atom pair (topological distance > 3 bonds)
    overlaps below ``frac`` x sum-of-vdW-radii — i.e. the combined ring puckers
    left the whole molecule sterically unrealistic despite the relax."""
    try:
        conf = mol.GetConformer()
        P = conf.GetPositions()
        syms = [a.GetSymbol() for a in mol.GetAtoms()]
        dm = Chem.GetDistanceMatrix(mol)
        heavy = [i for i, s in enumerate(syms) if s != "H"]
        for a in range(len(heavy)):
            i = heavy[a]
            ri = _VDW.get(syms[i], 1.7)
            for b in range(a + 1, len(heavy)):
                j = heavy[b]
                if dm[i][j] <= 3:            # bonded / 1-3 / 1-4 -> expected close
                    continue
                d = float(_np.linalg.norm(P[i] - P[j]))
                if d < frac * (ri + _VDW.get(syms[j], 1.7)):
                    return True
    except Exception:
        return False
    return False


def _bindungs_ausreisser(mol, tol_lang: float = 1.30) -> frozenset:
    """Die Bindungen des GRAPHEN, deren LAENGE keine Bindung mehr beschreibt.

    ===== WARUM DAS DRITTE TOR UEBERHAUPT FEHLT (26.08.2026) ==========================

    `generate` hat zwei Tore: `_has_clash` sieht NICHT gebundene Paare, die zu nah
    stehen, `_has_bad_angles` sieht Winkel.  Die BINDUNGSLAENGE selbst prueft keines
    von beiden.  Eine Faltung, die eine Bindung auseinanderzieht, kommt damit durch.
    Genau das ist am 18.08. passiert und steht protokolliert: "in einem Frame riss
    eine Bindung".

    ⚠ UND DAS SELBSTGATE KANN ES NICHT AUFFANGEN -- per Konstruktion, nicht aus
      Nachlaessigkeit.  `assemble_complex._collapsed_heavy_bonds_strict` laeuft ueber
      alle Schweratom-PAARE und entscheidet aus dem ABSTAND, ob sie gebunden sind:

          if d > 1.30 * ideal:   continue        # "also gar nicht gebunden"

      Eine auf das 1,4-fache gedehnte Bindung faellt damit AUS DER PRUEFUNG HERAUS.
      Sie wird nicht als kaputt gemeldet, sondern als nicht vorhanden.  Der Kollaps
      (zu kurz) wird gesehen, der Bruch (zu lang) ist ein blinder Fleck.

    HIER liegt der Bindungsgraph vor.  Damit ist "gebunden" keine Abstandsfrage mehr,
    und dieselbe Zahl 1,30 wird von einem AUSSCHLUSSkriterium zur BRUCHschwelle.  Es
    wird nichts erfunden: Boden (`_bd.COLLAPSE_FLOOR`, 0,82) und Decke (1,30) sind
    exakt die beiden Zahlen, mit denen das Selbstgate ohnehin schon rechnet, und
    `_ideal_bond` ist dieselbe Quelle.

    ⚠ MENGE STATT WAHRHEITSWERT, und das ist der ganze Unterschied zwischen Filter und
      Urteil.  Ein Nitril sitzt bei 1,20 A gegen ein Einfachbindungs-Ideal von 1,52 --
      Verhaeltnis 0,79, unter dem Boden.  Ein absolutes Ja/Nein wuerde JEDE Faltung
      JEDES nitrilhaltigen Molekuels verwerfen, und der Befund hiesse "das Tor hat
      keine Reichweite" aus dem falschen Grund.  Der Aufrufer zieht darum die Menge
      des GRUNDZUSTANDS ab: verworfen wird nur, was die Faltung NEU EINBRINGT.  Das ist
      dieselbe never-worse-Form, die `converter_backend` seinen Geschwistern auferlegt.
    """
    try:
        from delfin.manta import _bond_decollapse as _bd
        P = mol.GetConformer().GetPositions()
    except Exception:
        return frozenset()
    try:
        floor = float(_bd.COLLAPSE_FLOOR)
    except Exception:
        floor = 0.82
    aus = set()
    try:
        for b in mol.GetBonds():
            i, j = int(b.GetBeginAtomIdx()), int(b.GetEndAtomIdx())
            si = mol.GetAtomWithIdx(i).GetSymbol()
            sj = mol.GetAtomWithIdx(j).GetSymbol()
            if si == "H" or sj == "H":
                continue                       # H wie im Selbstgate: nicht beurteilt
            try:
                if _bd._is_metal(si) or _bd._is_metal(sj):
                    continue                   # M-D-Ideal ist erfunden, s. `_ideal_bond`
            except Exception:
                pass
            try:
                ideal = float(_bd._ideal_bond(si, sj, bool(b.GetIsAromatic())))
            except Exception:
                continue
            if ideal <= 1e-6:
                continue
            r = float(_np.linalg.norm(P[i] - P[j])) / ideal
            if r > tol_lang or r < floor:
                aus.add((min(i, j), max(i, j)))
    except Exception:
        return frozenset()
    return frozenset(aus)


def _is_puckerable(mol, ring) -> bool:
    """A ring is puckerable iff it is saturated enough to have out-of-plane
    minima: non-aromatic, size 5-8, and >= 3 sp3 ring atoms (an aromatic /
    fully-conjugated ring is planar and rigid -> no pucker conformers).

    Robust for BOTH a sanitised RDKit mol (uses hybridisation) AND a distance-
    perceived metal-complex mol (no hybridisation, all bonds single -> cannot
    tell benzene from cyclohexane from the graph, so read sp3 from the 3D SHAPE:
    a saturated centre is 4-coordinate tetrahedral or 3-coordinate pyramidal,
    an aromatic/sp2 centre is 3-coordinate planar)."""
    n = len(ring)
    # ===== DAS GROESSENFENSTER 5..8 IST EINE FESSEL, KEIN GESETZ (26.08.2026) =======
    #
    # Cremer-Pople gilt fuer JEDEN Ring ab N = 4: die Zahl der Faltungsfreiheitsgrade
    # ist N-3, und `_set_pucker_general` traegt sie inzwischen alle.  Das Fenster hier
    # schnitt trotzdem bei 8 ab -- ein Vierring (1 DOF, echte Schmetterlingsfaltung)
    # und JEDER Makrozyklus ab 9 waren damit per Konstruktion unfaltbar.
    # Porphyrine, Calixarene, Kronenether, grosse Chelatringe: null Faltung, nicht
    # weil die Mathematik fehlt, sondern weil eine Zahl im Weg stand.
    # ⚠ Nebenbefund: `_amp` fuehrt eine Tabelle fuer 5..8 UND einen Fallback fuer den
    #   Rest -- und der Fallback liegt UEBERALL ueber der Tabelle (n=5: 0,75 gegen
    #   0,40, fast doppelt).  Weil dieses Fenster jeden anderen Ring abwies, ist der
    #   Fallback NIE gelaufen.  Er wird mit dem Fenster zusammen korrigiert.
    # ⛔ Vorgabe AUS -> altes Fenster -> byte-identisch.
    if _os.environ.get("DELFIN_FFFREE_PUCKER_SPACE", "0") == "1":
        if n < 4:
            return False
    elif n < 5 or n > 8:
        return False
    try:
        P = mol.GetConformer().GetPositions()
    except Exception:
        P = None
    # ===== DER METALLACYCLUS FAELLT AN EINEM KRITERIUM, DAS NICHT FUER IHN GILT =====
    #
    # GEMESSEN 26.08. auf 400 Systemen (`find_conformer_coverage`):
    #     Ringe gesamt 1707 · METALL uebersprungen 352 = 20,6 %
    #     Sechsringe 176, erreichen den Sessel 57  ->  67,6 % NIE
    #
    # Die beiden Bedingungen unten sind fuer einen ORGANISCHEN Ring richtig und fuer
    # einen METALLACYCLUS falsch, und zwar aus demselben Grund: sie suchen die
    # Weichheit an den RINGATOMEN.  Bei einem Chelatring sitzt sie in den M-D-BINDUNGEN
    # -- 2,0 bis 2,4 A lang, weich, mit niedriger Torsionsbarriere.  Ein Salen-Ring
    # M-N=C-C(ar)-C(ar)-O faltet real an genau diesen beiden Bindungen (die "Stufe"
    # bzw. Umbrella-Faltung), obwohl sein organischer Teil starr und aromatisch ist.
    #
    #   * `GetIsAromatic() -> return False` kippt den GANZEN Ring, sobald EIN Atom
    #     aromatisch ist.  Beim fusionierten Salen-Metallacyclus sind das die
    #     Phenolat-Kohlenstoffe -> sofortiger Ausschluss.
    #   * `n_sat >= 3` verlangt drei sp3-Ringatome.  Ein konjugierter Chelatring hat
    #     sie nicht und braucht sie auch nicht.
    #
    # ⇒ Fuer einen Ring MIT Metall gilt: aromatisch ist nur dann ein Ausschluss, wenn
    #   ALLE Nicht-Metall-Ringatome aromatisch sind (dann ist der Ring wirklich
    #   planar-starr, z.B. ein Metallabenzol).  Und die sp3-Schwelle faellt auf 1,
    #   weil das Metall selbst das Scharnier stellt, nicht ein sp3-Zentrum.
    # ⚠ Die Koordinationssphaere bleibt unberuehrt: `frozen` haelt Metall UND Donoren
    #   fest, es bewegen sich nur die Ringatome dazwischen.  Reiner Erzeuger.
    # ⛔ Vorgabe AUS -> Ringmenge unveraendert -> byte-identisch.
    _mc = False
    if _os.environ.get("DELFIN_FFFREE_PUCKER_MC", "0") == "1":
        try:
            from delfin.manta import _elements as _EL
            _mc = any(_EL.is_metal(mol.GetAtomWithIdx(int(i)).GetSymbol()) for i in ring)
        except Exception:
            _mc = False
        if _mc:
            _nonmetal = [int(i) for i in ring
                         if not _EL.is_metal(mol.GetAtomWithIdx(int(i)).GetSymbol())]
            if _nonmetal and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in _nonmetal):
                return False          # vollstaendig aromatischer Metallacyclus: starr
    n_sat = 0
    for idx in ring:
        a = mol.GetAtomWithIdx(int(idx))
        if a.GetIsAromatic():
            if not _mc:
                return False
            continue
        hyb = a.GetHybridization()
        if hyb == Chem.HybridizationType.SP3:
            n_sat += 1
            continue
        if hyb in (Chem.HybridizationType.SP2, Chem.HybridizationType.SP):
            continue
        # unspecified (perceived mol) -> geometric sp3 test
        if P is None:
            continue
        nbrs = [nb.GetIdx() for nb in a.GetNeighbors()]
        if len(nbrs) >= 4:
            n_sat += 1
        elif len(nbrs) == 3:
            c = P[int(idx)]
            q0, q1, q2 = P[nbrs[0]], P[nbrs[1]], P[nbrs[2]]
            nrm = _np.cross(q1 - q0, q2 - q0)
            ln = float(_np.linalg.norm(nrm))
            if ln > 1e-9 and abs(float(_np.dot(c - q0, nrm / ln))) > 0.25:
                n_sat += 1   # pyramidal -> sp3-like
    # ⚠ Beim Metallacyclus stellt das METALL das Scharnier, nicht ein sp3-Zentrum --
    #   die Schwelle 3 ist dort ein organisches Kriterium am falschen Objekt.
    return n_sat >= (1 if _mc else 3)


def _ring_order(mol, ring_set):
    """Return the ring atoms in connectivity (traversal) order."""
    ring = list(ring_set)
    adj = {i: [] for i in ring}
    rs = set(ring)
    for i in ring:
        for nb in mol.GetAtomWithIdx(int(i)).GetNeighbors():
            j = nb.GetIdx()
            if j in rs:
                adj[i].append(j)
    order = [ring[0]]
    prev = None
    cur = ring[0]
    for _ in range(len(ring) - 1):
        nxts = [x for x in adj[cur] if x != prev]
        if not nxts:
            return ring  # fall back to arbitrary order
        nxt = nxts[0]
        order.append(nxt)
        prev, cur = cur, nxt
    return order


def _conf_to_xyz(mol) -> str:
    conf = mol.GetConformer()
    out = []
    for i in range(mol.GetNumAtoms()):
        a = mol.GetAtomWithIdx(i)
        p = conf.GetAtomPosition(i)
        out.append(f"{a.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return "\n".join(out) + "\n"


def _tfd(acc_mol, id_a: int, id_b: int) -> float:
    """Torsion-Fingerprint-Deviation between two conformers of ``acc_mol``.
    TFD is the field-standard conformer discriminator — it compares ALL ring +
    rotatable-bond torsions with the molecule's topological symmetry folded in,
    so pucker families (chair vs twist-boat), rotamers and axial/equatorial
    substituents separate cleanly where heavy-atom RMSD conflates them."""
    try:
        from rdkit.Chem import TorsionFingerprints as _TF
        return float(_TF.GetTFDBetweenConformers(acc_mol, [id_a], [id_b])[0])
    except Exception:
        return 1.0     # no torsions / failure -> treat as distinct (keep)


# ===== DIE RINGLOKALE TFD (26.08.2026) =============================================
#
# GEMESSEN (`selbsttest_trennschaerfe`, Verduennungsreihe): derselbe Cyclohexanring an
# einem wachsenden starren Acen, Zustaende JE RING --
#     Cyclohexylbenzol      Ringanteil 50,0 %   TFD 0,05 -> 11   TFD 0,005 -> 62
#     Cyclohexylnaphthalin             37,5 %                6                49
#     Cyclohexylanthracen              30,0 %                1                23
#     Cyclohexyltetracen               25,0 %                1                10
# Die Faltungen sind DA -- bei 0,005 kommen sie zurueck.  TFD verschmilzt sie.
#
# DIE URSACHE STEHT IN RDKITS EIGENER FORMEL, und sie ist SCHAERFER als "mittelt".
# `CalculateTFD` bildet sum(d_i * w_i) / sum(w_i) ueber ALLE Torsionen des Molekuels.
# Bewegt sich nur der eine Ring, ist d_i = 0 fuer jede andere Torsion, und es bleibt
#     TFD_global = d_Ring * w_Ring / sum(w)
# `CalculateTorsionWeights` setzt w = exp(-beta * d^2) mit d = topologischer Abstand zur
# ZENTRALSTEN Bindung des Molekuels.  Ein angehaengter Cyclohexylring rutscht mit jedem
# weiteren Acenring weiter an den Rand -- sein Gewichtsanteil faellt EXPONENTIELL, nicht
# wie 1/N.  Gemessen an denselben vier MMFF-optimierten Proben (w_Ring / sum(w)):
#     Cyclohexylbenzol      0,1875   ->   5,3-fache Verduennung
#     Cyclohexylnaphthalin  0,0548   ->  18,3-fache
#     Cyclohexylanthracen   0,0166   ->  60,2-fache
#     Cyclohexyltetracen    0,0069   -> 146,0-fache
# 0,05 / 146 = 0,00034 ist damit die Schwelle, die der Ring im Tetracen EFFEKTIV sieht.
# Genau deshalb kommen die Zustaende erst bei 0,005 zurueck, und genau deshalb faellt
# die Zahl monoton mit der Geruestgroesse.
#
# ⇒ Wort fuer Wort der RMSD-Fehler eine Ebene hoeher -- und der Effekt WAECHST mit der
#   Ligandgroesse, trifft also am haertesten die Systeme, um die es geht.
#
# DIE REPARATUR NIMMT RDKITS EIGENEN WEG, keinen Nachbau: `CalculateTorsionLists` gibt
# Nichtring- und Ringtorsionen GETRENNT zurueck, `CalculateTorsionAngles` und
# `CalculateTFD` nehmen genau solche Listen entgegen.  Es wird also nur GEFILTERT: der
# Eintrag des betrachteten Rings bleibt, alles andere faellt weg.
#
# ⚠ DIE SYMMETRIEFALTUNG UEBERLEBT -- der Punkt, an dem der Versuch vom 26.08. gestorben
#   ist, TFD durch eine reine CP-Distanz zu ERSETZEN (n=5 ging von 3,3,3 auf 9,13,14,
#   weil phi an der Atomnummerierung haengt).  Hier wird nichts nachgebaut: RDKits
#   Ringeintrag ist der MITTELWERT von |Torsion| ueber den ganzen Ring, also eine Zahl,
#   die unter Drehung UND Spiegelung der Ringnummerierung invariant ist.  Diese
#   Invarianz IST die Faltung.  Gemessen: unsubstituierter Fuenfring liefert ringlokal
#   3 Zustaende -- exakt wie global.
#
# ⚠ DIE SCHWELLE BLEIBT 0,05, und das ist keine Setzung, sondern eine IDENTITAET.  Ein
#   unsubstituierter Einringer hat KEINE Nichtringtorsion und GENAU EINEN Ringeintrag;
#   sum(w) ist dann w_Ring, der Bruch kuerzt sich, und global TFD = ringlokal TFD auf
#   jedem Konformerpaar.  GEMESSEN, Nenner 760 Konformerpaare (n=5,6,7,8 zu je 190):
#   groesste Differenz 5,6e-17 -- das ist Fliesskommarauschen, nicht ein kleiner
#   Unterschied.  Auf der Kalibrierprobe sind die beiden Masse also nicht aehnlich
#   geeicht, sondern DASSELBE; die Bedeutung von 0,05 aendert sich dort um exakt null.
#   Die Zustandszahlen bestaetigen es Zeile fuer Zeile (3/9/11/16 bei 0,05 und
#   6/22/64/103 bei 0,005, global wie ringlokal).  `selbsttest_tfd_lokal` misst beides.
#
# ⚠ WAS DIESE FASSUNG DAFUER BEZAHLT, und es steht hier, damit es niemand spaeter als
#   Ueberraschung findet: RDKits Ringeintrag ist EINE Zahl je Ring.  Dieselbe
#   Invarianz, die die Symmetrie faltet, macht das Mass eindimensional -- zwei
#   wirklich verschiedene Faltungen mit demselben Mittelwert von |Torsion| werden
#   zusammengezogen.  Die Faltungsachse selbst (theta, phi) sieht das Mass NICHT.
#   ⇒ Ringlokale TFD ist eine schaerfere Entdopplung, KEIN vollstaendiger
#     Faltungsdeskriptor.  Wer die Achse braucht, braucht `_cp_theta_phi` dazu -- und
#     die faltet die Symmetrie NICHT (Messung 26.08.), taugt also nur als ZWEITES
#     Instrument neben diesem, nie als Ersatz.
# ⚠ OHNE GEWICHTE, wenn mehrere Ringe ausgewaehlt sind.  Die Gewichte SIND der
#   Verduennungsmechanismus (Abstand zur zentralsten Bindung); sie ringlokal wieder
#   hereinzuholen holte den Effekt zurueck, den diese Fassung entfernt.  Bei EINEM Ring
#   ist es ohnehin gleichgueltig -- ein Gewicht kuerzt sich gegen sich selbst.
# ⛔ Vorgabe AUS -> `_tfd_distinct` laeuft die alte Zeile -> byte-identisch.


def _tfd_lokal_listen(mol, ringe):
    """RDKits Ringtorsionsliste, GEFILTERT auf die uebergebenen Ringe.

    ⚠ ZUGEORDNET WIRD UEBER DIE ATOMMENGE, nicht ueber den Listenindex.
      `tors_list_rings` kommt aus `Chem.GetSymmSSSR`, die Ringe des Aufrufers aus
      `RingInfo.AtomRings()`.  Beide liefern dieselben Ringe -- sich auf ihre
      Indexgleichheit zu VERLASSEN waere aber eine Annahme, und ein falsch
      zugeordneter Ring waere hier nicht als Fehler zu erkennen, sondern nur als
      "der andere Ring hat sich eben nicht bewegt".  Der Aufrufer prueft darum die
      LAENGE der Rueckgabe gegen die Zahl der gewuenschten Ringe.
    """
    from rdkit.Chem import TorsionFingerprints as _TF
    _tl, _tlr = _TF.CalculateTorsionLists(mol)
    ziel = {frozenset(int(a) for a in r) for r in ringe}
    # Ringeintrag k besteht aus den N aufeinanderfolgenden Vierergruppen des Rings; die
    # ERSTEN Atome dieser Gruppen sind genau die N Ringatome (RDKit baut sie so).
    return [(q, d) for q, d in _tlr if frozenset(int(t[0]) for t in q) in ziel]


def _tfd_lokal(acc_mol, listen, id_a: int, id_b: int) -> float:
    """TFD ueber NUR die uebergebenen Torsionseintraege -- kein Geruest im Nenner."""
    from rdkit.Chem import TorsionFingerprints as _TF
    t_a = _TF.CalculateTorsionAngles(acc_mol, [], listen, confId=id_a)
    t_b = _TF.CalculateTorsionAngles(acc_mol, [], listen, confId=id_b)
    return float(_TF.CalculateTFD(t_a, t_b, weights=None))


# ⚠ ZWEI SCHALTER FUER EIN MASS, und das ist keine Knopfvermehrung.  Dasselbe Mass
#   bewegt die Zahl an den beiden Aufrufstellen in ENTGEGENGESETZTE Richtungen:
#     `_ring_pucker_states` (je Ring)   Cyclohexyltetracen  1 -> 11 Zustaende
#     `generate` (je Kombination)       Decalin            22 ->  8 Frames
#   Haengen beide an EINEM Schalter, misst ein A/B ihre SUMME und niemand kann sagen,
#   welcher Anteil woher kam -- genau die Bauform, an der in diesem Projekt schon
#   Verdikte gescheitert sind.  Getrennt geschaltet sind es zwei Messungen.
# ⚠ DER BEFUND HAENGT AM ERSTEN.  Gemessen wurde die Verduennung an den Zustaenden JE
#   RING; die Kombinationsebene ist eine EXTRAPOLATION davon und steht darum unter
#   ihrem eigenen, ebenfalls ausgeschalteten Schalter.
def _tfd_distinct(acc_mol, cid: int, kept_ids, thr: float, ringe=None,
                  schalter: str = "DELFIN_FFFREE_PUCKER_TFD_LOCAL") -> bool:
    if ringe and _os.environ.get(schalter, "0") == "1":
        try:
            _listen = _tfd_lokal_listen(acc_mol, ringe)
            # ⚠ ALLE ODER KEINER.  Findet die Zuordnung nur EINEN Teil der Ringe
            #   wieder, misst der ringlokale Vergleich stillschweigend weniger Ringe
            #   als der Aufrufer gemeint hat -- und die fehlenden faenden nirgends
            #   statt.  Eine halbe Messung sieht von aussen aus wie eine ganze; das ist
            #   genau die Bauform, die in diesem Projekt schon mehrfach als Befund
            #   durchgegangen ist.  Lieber ganz zurueck auf das globale Mass.
            if len(_listen) == len({frozenset(int(a) for a in r) for r in ringe}):
                # ⚠ EIGENE SCHWELLE NUR, WENN JEMAND SIE SETZT.  Die Messung sagt: auf
                #   dem unsubstituierten Ring sind beide Masse identisch, 0,05 behaelt
                #   also seine Bedeutung.  Der Knopf ist zum NACHMESSEN da, nicht zum
                #   Nachjustieren -- leer heisst "unveraendert".
                _s = _os.environ.get("DELFIN_FFFREE_PUCKER_TFD_LOCAL_THR", "")
                _thr = float(_s) if _s.strip() else thr
                # Die Liste EINMAL je Kandidat, nicht je Paar: `GetTFDBetweenConformers`
                # baut sie im globalen Pfad bei JEDEM Aufruf neu -- ringlokal ist damit
                # auch billiger, nicht nur schaerfer.
                return all(_tfd_lokal(acc_mol, _listen, k, cid) >= _thr
                           for k in kept_ids)
        except Exception:
            pass          # Rueckfall auf das globale Mass -- nie stillschweigend leer
    return all(_tfd(acc_mol, k, cid) >= thr for k in kept_ids)


# ===== DIE UNTERSCHEIDBARKEIT IST EIN MAXIMUM, KEIN MITTELWERT (26.08.2026) =========
#
# Der Fehler von RMSD ist NICHT, dass es Geometrie misst.  Er ist, dass es MITTELT --
# und eine Ringfaltung ist ein LOKALES Ereignis in einem grossen Molekuel.  Mit unseren
# eigenen Zahlen an einem gefalteten Sechsring durchgerechnet:
#
#     Ringatome laufen 0,203 A im Median (groesste Einzelauslenkung 0,33)
#     Ringanteil an den schweren Atomen 13,2 %
#     ⇒ Gesamt-RMSD = sqrt(0,132) * 0,203 = 0,086 A
#
# 0,086 liegt UNTER jeder Entdopplungsschwelle, die dieses Projekt fuehrt (`_DEDUP_RMSD`
# 0,30 · `rmsd_dedup` 0,5).  Die Faltung verschwindet also nicht, weil sie klein waere,
# sondern weil sie durch 87 % unbewegte Atome geteilt wird.  Die groesste Auslenkung
# nach Kabsch-Ausrichtung bleibt bei 0,33 -- FAKTOR 4 zwischen den beiden Zahlen, an
# derselben Geometrie gemessen.
#
# WAS EIN KRISTALLOGRAPH STATTDESSEN LIEST.  In der Differenz-Fourier-Karte steht die
# GROESSTE unmodellierte Abweichung als Restdichte-Maximum; der Mittelwert ueber alle
# Atome kommt darin nicht vor.  Zwei Modelle, deren groesste Atomauslenkung unter der
# Aufloesung liegt, waeren an denselben Daten NICHT ZU UNTERSCHEIDEN -- sie sind EIN
# Eintrag im Manifold, nicht zwei.
#
# DIE SCHWELLE, begruendet statt gesetzt:
#   * Koordinaten-esd einer Routinestruktur liegt bei 0,002 bis 0,01 A.  Das ist die
#     UNTERgrenze -- darunter ist jede Differenz Rauschen der Verfeinerung.
#   * Fehlordnung wird ab etwa 0,3 bis 0,5 A ueberhaupt erst als ZWEI Lagen modelliert.
#     Das ist die OBERgrenze -- darueber sieht der Kristallograph zwei Konformere.
#   0,15 A liegt dazwischen: Faktor 15 bis 75 ueber dem esd, Faktor 2 bis 3 unter der
#   Fehlordnungsgrenze.  Gross genug, um nicht Rauschen zu zaehlen; klein genug, um
#   nichts zu verschmelzen, was ein Kristallograph noch getrennt modellieren wuerde.
# ⚠ SIE IST ENV-PARAMETER, weil sie eine KONVENTION ist und keine Naturkonstante -- die
#   Aufloesung haengt am Datensatz, und wer sie verschiebt, soll das messen koennen.
#
# ⚠ WAS DIESE METRIK **NICHT** TUT: sie rangiert nicht.  Sie sagt "ununterscheidbar"
#   oder "unterscheidbar", nie "besser".  Ein Rang braeuchte eine Energie; das ist
#   Stufe (3) und steht aus gutem Grund AUS.


def _kabsch_max_rmsd(A, B) -> Tuple[float, float]:
    """(GROESSTE Auslenkung, RMSD) zweier Punktsaetze nach Kabsch-Ausrichtung.

    Beide Zahlen aus DERSELBEN Ausrichtung -- sonst waere der Vergleich der beiden
    Masse keiner.  Kabsch minimiert das RMSD; das Maximum wird also gegen die fuer
    RMSD GUENSTIGSTE Ueberlagerung gemessen und ist damit eher zu klein als zu gross.

    ⚠ NUR SCHWERE ATOME, und das ist kein Sparen.  Roentgenbeugung sieht ELEKTRONEN-
      DICHTE; ein Wasserstoff traegt ein Elektron und wird in einer Routinestruktur
      GERECHNET, nicht gefunden.  Eine Metrik, die vorgibt, H-Lagen zu unterscheiden,
      urteilt ueber etwas, das in den Daten gar nicht steht.  Der Aufrufer uebergibt
      darum bereits gefilterte Koordinaten.

    ⚠ KEINE SYMMETRIEFALTUNG, absichtlich.  Beide Punktsaetze stammen aus DERSELBEN
      Molekuelinstanz in DERSELBEN Atomreihenfolge; die Entartung "zwei Nummerierungen
      desselben Konformers" erledigt in `generate` das TFD davor, und TFD kann das,
      weil es die Topologiesymmetrie mitfaltet (die CP-Distanz konnte es nicht -- der
      Versuch vom 26.08. ist genau daran gescheitert).  Wer hier zusaetzlich ueber
      Automorphismen minimierte, zahlte N! und maesse dasselbe.
    """
    v = _kabsch_abweichungen(A, B)
    return float(v.max()), float(_np.sqrt(float((v ** 2).mean())))


def _kabsch_abweichungen(A, B):
    """Die Abweichung JE ATOM nach Kabsch-Ausrichtung -- der gemeinsame Rohstoff.

    ⚠ EINE Ausrichtung, dann beide Masse daraus.  Wuerde man Maximum und RMSD je
      einzeln ausrichten lassen, verglichen sie zwei verschiedene Ueberlagerungen und
      der Faktor zwischen ihnen waere teils Instrument, teils Ausrichtung.  Kabsch
      minimiert das RMSD -- das Maximum wird also gegen die fuer den Gegner
      GUENSTIGSTE Ueberlagerung gemessen und ist eher zu klein als zu gross.
    """
    Am = A.mean(0)
    Bm = B.mean(0)
    Ac = A - Am
    Bc = B - Bm
    try:
        U, _S, Vt = _np.linalg.svd(Bc.T @ Ac)
        d = 1.0 if float(_np.linalg.det(Vt.T @ U.T)) > 0.0 else -1.0
        R = Vt.T @ _np.diag([1.0, 1.0, d]) @ U.T
        Bd = Bc @ R.T
    except Exception:
        Bd = Bc
    return _np.linalg.norm(Ac - Bd, axis=1)


def _add_conf(acc_mol, src_mol) -> int:
    return acc_mol.AddConformer(Chem.Conformer(src_mol.GetConformer()), assignId=True)


def _ring_pucker_states(mol_with_conf, ring, frozen: Set[int],
                        tfd_thr: float) -> List[Optional[Tuple[Optional[float], float]]]:
    """DISTINCT pucker settings for ONE ring, TFD-deduped.  ``None`` = the base
    pucker; each other entry is a ``(theta, phi)`` SETTING that, after a torsion-
    held relax, gives a conformer whose torsion fingerprint differs from every
    kept one (cyclohexane -> {base chair, the twist-boat(s)}, not 7 relabelled
    pseudorotation copies)."""
    states: List[Optional[Tuple[float, Optional[float], float]]] = [None]
    acc = Chem.Mol(mol_with_conf)
    kept_ids = [acc.GetConformer().GetId()]
    n = len(ring)
    # ===== ENTDOPPELN IN CP STATT IN TFD (26.08.2026) ================================
    #
    # GEMESSEN am eigenen Konvergenztest.  Die TFD-Schwelle 0,05 splittet ueber:
    #     n=6, NPHASE=16 -> 10 Zustaende, kleinste paarweise CP-Distanz  5,0 Grad
    #     n=7, NPHASE=16 -> 12 Zustaende, kleinste paarweise CP-Distanz  1,8 Grad
    # Zwei Faltungen, die 1,8 Grad auseinanderliegen, sind DIESELBE.  Und zwei
    # Eintraege standen bei theta = 180 mit phi = 136 und phi = 339 -- am Pol ist phi
    # bedeutungslos, also beweisbar derselbe Zustand.
    #
    # ⚠️ DAS IST KEINE KOSMETIK.  Die Zustandszahl JE RING ist die Basis des
    #    Kreuzprodukts ueber alle Ringe (gemessen: 4,3 Ringe je System):
    #         3 Zustaende, 4 Ringe ->      81 Kombinationen   rechenbar
    #        10 Zustaende, 4 Ringe ->  10 000                 nicht rechenbar
    #    Uebersplittung macht die VOLLSTAENDIGE Kombinatorik unbezahlbar.  Wer den
    #    ganzen Faltungsraum will, muss zuerst aufhoeren, Rauschen als Mulde zu
    #    zaehlen -- sonst kommt die Kappe durch die Hintertuer zurueck.
    #
    # Chemischer Massstab: Cyclohexan hat Sessel + Twist-Boat-Familie, nach
    # Symmetriefaltung 2-3 Klassen.  Der Fuenfring konvergiert von selbst auf 3.
    #
    # Entdoppelt wird darum auf der KUGEL, mit der Grosskreisdistanz -- die
    # Pol-Entartung loest sich dort von selbst (siehe `_cp_abstand`).
    # ⛔ Vorgabe AUS -> TFD wie bisher -> byte-identisch.
    _cpd = _os.environ.get("DELFIN_FFFREE_PUCKER_CPDEDUP", "0") == "1"
    _cp_tol = float(_os.environ.get("DELFIN_FFFREE_PUCKER_CPTOL", "15") or 15.0)
    _cp_qtol = float(_os.environ.get("DELFIN_FFFREE_PUCKER_CPQTOL", "0.15") or 0.15)
    _cp_kept: List[Tuple[float, float, float]] = []
    # ===== (2) DIESELBE UNUNTERSCHEIDBARKEIT, ABER JE RING (26.08.2026) ==============
    #
    # ⚠ HIER LIEGT DER HEBEL, NICHT IM KREUZPRODUKT.  Die Zustandszahl JE RING geht
    #   POTENZIERT in die Kombinatorik ein (gemessen 4,3 Ringe je System): ein Zustand
    #   weniger je Ring spart mehr als jede Regel weiter unten, weil unten schon
    #   relaxiert wurde.  Ein Tor hinter dem Relax toetet das Ergebnis, nicht die
    #   Kosten -- das steht seit dem 26.08. im Docstring von `selbsttest_kombinatorik`
    #   und gilt fuer den Defektfilter genauso wie fuer die Ununterscheidbarkeit.
    # ⚠ OB ES WIRKLICH REDUZIERT, IST EINE MESSUNG UND KEINE HOFFNUNG.  Das Maximum ist
    #   gegen Verduennung unempfindlich (ein Maximum kennt keinen Nenner), also spricht
    #   nichts dafuer, dass es echte Ringmulden zusammenzieht -- die liegen weit
    #   auseinander.  `selbsttest_trennschaerfe` Schritt 3 zaehlt nach.
    # ⛔ Vorgabe AUS -> byte-identisch.
    _xrd_r = _os.environ.get("DELFIN_FFFREE_PUCKER_XRD", "0") == "1"
    _xrd_r_tol = float(_os.environ.get("DELFIN_FFFREE_PUCKER_XRDTOL", "0.15") or 0.15)
    _schwer_r: List[int] = []
    _xrd_r_kept: List = []
    if _xrd_r:
        try:
            _schwer_r = [i for i in range(mol_with_conf.GetNumAtoms())
                         if mol_with_conf.GetAtomWithIdx(i).GetSymbol() != "H"]
            _xrd_r_kept = [mol_with_conf.GetConformer().GetPositions()[_schwer_r]]
        except Exception:
            _xrd_r = False
    # ===== DER GANZE FALTUNGSRAUM STATT DREI STELLEN DARAUF (26.08.2026) ============
    #
    # Die alte Kandidatenliste tastet die Cremer-Pople-Kugel an genau drei Orten ab:
    # den AEQUATOR (theta = 90, K Phasen), und bei geraden Ringen die beiden POLE.
    # `q_scale` ist dabei konstant 1,0.
    #   ⇒ theta zwischen 0 und 90 wird NIE abgetastet -- dort liegen Half-Chair
    #     (theta ~50) und Envelope (theta ~55).
    #   ⇒ die Amplitude wird NIE variiert -- nur EINE Kugelschale.
    #   ⇒ ungerade Ringe bekommen `theta=None`, also reine Pseudorotation.
    # Die fehlenden Formen sind damit nicht "nicht implementiert", sondern NICHT
    # ABGETASTET -- ein Unterschied, der die Reparatur billig macht.
    #
    # Mit `DELFIN_FFFREE_PUCKER_SPACE=1` wird stattdessen der (N-3)-dimensionale
    # Raum systematisch ueberdeckt (`_pucker_space_grid`), fuer JEDE Ringgroesse und
    # ueber `_set_pucker_general`, das auch die Paare m >= 3 kennt -- ohne die war
    # jeder Ring ab N = 7 unvollstaendig parametrisiert.
    # ⚠ Vollstaendigkeit ist hier eine AUFLOESUNGSfrage: der Raum ist kontinuierlich.
    #   Ueberdeckt wird er vollstaendig bei der angegebenen Aufloesung, und die steht
    #   in der Spur -- keine Ecke ausgelassen, keine Richtung bevorzugt.
    # ⛔ Vorgabe AUS -> alte Liste -> byte-identisch.
    _raum = _os.environ.get("DELFIN_FFFREE_PUCKER_SPACE", "0") == "1"
    if _raum:
        _namp = max(1, int(_os.environ.get("DELFIN_FFFREE_PUCKER_NAMP", "2") or 2))
        _nph = max(1, int(_os.environ.get("DELFIN_FFFREE_PUCKER_NPHASE", "8") or 8))
        _cands = _pucker_space_grid(n, _namp, _nph)
        if _os.environ.get("DELFIN_FFFREE_PUCKER_TRACE", "0") == "1":
            # ⚠ DIE ABTASTDICHTE STEHT MIT IM PROTOKOLL.  Ohne sie liesse sich eine
            # reduzierte Aufloesung spaeter nicht von einer vollen unterscheiden --
            # und genau das waere eine stille Kappe.
            _res = getattr(_pucker_space_grid, "_grid_res", None)
            print("[pucker] RAUM n=%d: %d Kandidaten (%d-dim, Amplitudenstufen %d, "
                  "Phasen %d)" % (n, len(_cands), max(0, n - 3), _namp, _nph))
            if isinstance(_res, dict) and _res.get("n") == n and _res.get("reduziert"):
                print("[pucker] RAUM n=%d: Dichte REDUZIERT auf Budget %d -- Phasen je m %s"
                      % (n, _res.get("budget"), _res.get("n_phase_je_m")))
    else:
        _cands = _pucker_candidates(n)
    if frozenset(ring) in _FLAT_ONLY:
        # NUR begradigen, nicht falten -- s. den Block in `generate`.
        _cands = [({}, {})] if _raum else [(0.0, 0.0, 0.0)]
    for _cand in _cands:
        try:
            m2 = Chem.Mol(mol_with_conf)
            if _raum:
                _qs, _phis = _cand
                _set_pucker_general(m2.GetConformer(), ring, _qs, _phis, frozen)
                theta = phi = None
            else:
                _qs, theta, phi = _cand
                _set_pucker(m2.GetConformer(), ring, _qs * _amp(n), theta, phi, frozen)
            if not _relax_hold_pucker(m2, ring, frozen):
                continue
            _Pr = None
            if _xrd_r:
                # Ununterscheidbar vom Grundzustand ODER von einem schon behaltenen
                # Zustand -> derselbe Eintrag, kein zweiter Ringzustand.
                _Pr = m2.GetConformer().GetPositions()[_schwer_r]
                if any(_kabsch_max_rmsd(_Pa, _Pr)[0] < _xrd_r_tol
                       for _Pa in _xrd_r_kept):
                    continue
            if _cpd:
                # NACH dem Relax messen, nicht die SOLL-Werte vergleichen: der Relax
                # zieht den Startpunkt in die naechste echte Mulde, und genau deren
                # Lage entscheidet, ob es eine neue ist.
                _cp = _cp_theta_phi(m2.GetConformer().GetPositions(), ring)
                if any(_cp_abstand(_cp, _k) < _cp_tol and abs(_cp[0] - _k[0]) < _cp_qtol
                       for _k in _cp_kept):
                    continue
                _cp_kept.append(_cp)
                if _xrd_r and _Pr is not None:
                    _xrd_r_kept.append(_Pr)
                states.append(_cand if _raum else (_qs, theta, phi))
                continue
            cid = _add_conf(acc, m2)
            # ⚠ HIER IST DIE MESSSTELLE DES BEFUNDS.  `ringe` benennt den EINEN Ring,
            #   der hier gefaltet wird; mit dem Schalter AN zaehlt nur noch seine
            #   Torsion, das Geruest steht nicht mehr im Nenner.  Schalter AUS -> das
            #   Argument wird in `_tfd_distinct` gar nicht angesehen.
            if _tfd_distinct(acc, cid, kept_ids, tfd_thr, ringe=(ring,)):
                kept_ids.append(cid)
                if _xrd_r and _Pr is not None:
                    _xrd_r_kept.append(_Pr)
                # Im Raum-Modus ist der Zustand das Koordinatenpaar selbst; die
                # Legacy-Form bleibt ein 3-Tupel.  `generate` indiziert nur, es liest
                # den Inhalt nicht -- beide Formen sind dort gleichwertig.
                states.append(_cand if _raum else (_qs, theta, phi))
            else:
                acc.RemoveConformer(cid)
        except Exception:
            continue
    return states


def _neuer_zaehler() -> dict:
    """Frischer Zaehlersatz fuer ``generate(..., _zaehler=...)``.

    ⚠ WARUM DIE MESSSTELLE IN `generate` SITZT UND NICHT IN EINER KOPIE.  Die Frage,
      wie gross das Kreuzprodukt NACH der Physik ist, laesst sich nur an dem Code
      beantworten, der die Frames auch wirklich baut.  Eine nachgebaute Schleife misst
      den Nachbau -- in diesem Projekt ist genau das schon mehrfach als Befund
      durchgegangen und war keiner.  Der Preis ist ein `if _zaehler is not None`
      an sechs Stellen; der Vorgabepfad (`_zaehler is None`) laeuft unveraendert.
    """
    return {"ringgroessen": [], "zustaende_je_ring": [], "kreuzprodukt": 0,
            "gemeinsame_atome": 0, "gem_max": 0, "aufzaehlung": 0, "gebaut": 0,
            "relax_fehler": 0, "kollision": 0, "winkel": 0, "tor_ueberlebt": 0,
            "tfd_doppelt": 0, "ausnahme": 0, "energien": [],
            # (1) Defektfilter und (2) kristallographische Ununterscheidbarkeit bekommen
            # EIGENE Zaehler.  Sie mit `kollision`/`tfd_doppelt` zu verrechnen, waere
            # genau der Fehler, den dieses Projekt schon dreimal gemacht hat: ein
            # Detektorname, der zwei Mechanismen deckt, ist keine Messung.
            "bindung": 0, "xrd_doppelt": 0}


def generate(mol_with_conf, frozen: Optional[Set[int]] = None,
             budget: int = 64, tfd_thr: float = 0.05,
             angle_skip: Optional[Set[int]] = None,
             _zaehler: Optional[dict] = None) -> List[Tuple[str, str]]:
    """Construct the COMBINATORIAL ring-pucker conformers from a base conformer.

    ``_zaehler``: optionaler Zaehlersatz (`_neuer_zaehler()`).  Ist er gesetzt, traegt
    dieser Lauf mit, wie viele Kombinationen aufgezaehlt, gebaut, am Kollisions- bzw.
    Winkeltor verworfen und von TFD zusammengezogen wurden -- die Messung, die
    `selbsttest_kombinatorik` auswertet.  ``None`` (Vorgabe) = kein einziger Zaehler
    wird angefasst, der Bau ist byte-identisch zu vorher.

    ``mol_with_conf`` carries ONE embedded conformer (a chain/rotamer pose whose
    rings sit at their base pucker).  Every puckerable ring's distinct pucker
    states are enumerated, then the CARTESIAN PRODUCT across all rings is
    constructed (Cy3P: 3 rings x {chair, twist-boat, ...} -> 3xchair, 2xchair+
    twist, ...).  Each combination sets all rings' puckers simultaneously and is
    relaxed with EVERY ring's pucker HELD but all inter-ring bonds/torsions FREE,
    so the free degrees of freedom relieve any inter-ring steric clash while the
    puckers themselves survive.  A whole-molecule clash gate then drops any
    combination that stayed sterically unrealistic, and TFD dedup keeps one
    representative per distinct torsion fingerprint (so Cy3P's three equivalent
    rings collapse correctly).  ``frozen`` fixes metal + donor atoms so metal
    chelate rings pucker without disturbing the coordination sphere.  Returns the
    NEW distinct, clash-free conformers ``[(xyz, label), ...]``; never raises.
    """
    if not (_RDKIT and _np is not None):
        return []
    try:
        if mol_with_conf.GetNumConformers() == 0:
            return []
        ri = mol_with_conf.GetRingInfo()
        rings_raw = [set(r) for r in ri.AtomRings()]
    except Exception:
        return []
    frozen = frozen or set()
    rings = [_ring_order(mol_with_conf, r) for r in rings_raw
             if _is_puckerable(mol_with_conf, r)]
    # ===== DER FLACHE ZUSTAND FUER KONJUGIERTE RINGE (18.08.2026) =====================
    # GEMESSEN (16./17.08., `folds`, 965 Systeme): von 326 fehlenden Ringmotiven sind
    # **218 PLANAR** -- 5M 126 - 6M 55 - 4M 37 -- und ALLE drei Klassen sind
    # METALLACYCLEN (`find_conformer_completeness:254` baut den Namen als
    # f"{sz}{'M' if is_metallacycle else ''}:{basin}").
    #
    # `_is_puckerable` laesst genau diese nicht herein: es verlangt "kein aromatisches
    # Ringatom" und ">= 3 sp3-Ringatome", und begruendet das damit, ein konjugierter Ring
    # SEI ohnehin planar und rigide.  Das Auge misst das Gegenteil: der Kristall-
    # Planarzustand FEHLT im Bau.  Beides zusammen heisst -- der Ring wird von etwas
    # anderem gefaltet, und das einzige Modul, das ihn absichtlich flach setzen koennte,
    # darf ihn nicht anfassen.  `planar138` und `pktrace` haben das bestaetigt
    # (affected 0, auch mit bewusst umgangener Reichweitensperre).
    #
    # ⚠ NUR DER Q=0-ZUSTAND, kein Pucker.  Diese Ringe sollen nicht gefaltet, sondern
    # BEGRADIGT werden; `_ring_pucker_states` bietet ihnen darum ausschliesslich die
    # Projektion in die Mittelebene an.  Metall und Donoren stehen in `frozen` und
    # bewegen sich nicht -- die Koordinationssphaere bleibt unberuehrt.
    #
    # ⚠ KLASSE: ENUMERATOR, kein Reparateur (Modulzensus 18.08.).  Er PROJIZIERT bei
    # eingefrorenem Kern, statt neu zu erzeugen -- dieselbe Klasse wie der
    # Spiegelabschluss (+1,0 pp), nicht die von BACKBONE_REEMBED (+11,9 pp).
    #
    # Vorgabe AUS -> Ringmenge unveraendert -> byte-identisch.
    # ⚠️ HIER, NICHT IN _set_pucker -- ein Fehler von mir, am 18.08. gemessen.
    # `_FLAT_ONLY.clear()` stand in `_set_pucker`, also in der Funktion, die je
    # RING einmal laeuft.  Folge: Ring 0 wurde korrekt flach gehalten, danach war
    # die Menge leer, und JEDER weitere Ring bekam die volle Faltleiter bei
    # Q = 0,63 Angstroem.  Der Mechanismus hat den konjugierten Chelatring also
    # GEFALTET, statt ihn zu BEGRADIGEN -- das genaue Gegenteil seines Zwecks.
    #
    # Der Beweis stand in den Etiketten: 246 von 256 Pucker-Etiketten beginnen mit
    # `r0:base` (nur Ring 0 blieb flach), und die Zustandsindizes an Ringen ab 1
    # laufen bis 15 -- ein flat-only-Ring kann hoechstens Index 1 haben, und 15 ist
    # exakt die volle Kandidatenzahl eines Sechsrings.
    # Chemisch gemessen an vier Salicylaldiminato-Chelaten: Walsh-Winkel am
    # Azomethin-Kohlenstoff 13,1 bis 17,5 Grad, waehrend der Kristall dieselben
    # Zentren bei hoechstens 1,4 Grad haelt; in einem Frame riss eine Bindung.
    #
    # ⇒ Der Torterm `smiles_ccdc_regressed`, an dem `planarA` blockierte, hatte
    # RECHT.  Ihn zu lockern haette den Baufehler einzementiert.
    if _os.environ.get("DELFIN_FFFREE_PUCKER_PLANAR", "0") == "1":
        _FLAT_ONLY.clear()
        _have = {frozenset(r) for r in rings}
        for _r in rings_raw:
            if frozenset(_r) in _have or not (5 <= len(_r) <= 8):
                continue
            try:
                _ro = _ring_order(mol_with_conf, _r)
            except Exception:
                continue
            if _ro:
                _FLAT_ONLY.add(frozenset(_ro))
                rings.append(_ro)
    if not rings:
        return []
    # per-ring distinct pucker states (index 0 == base pucker for every ring),
    # TFD-deduped so an unsubstituted ring yields only its genuine minima.
    if _zaehler is not None:
        import time as _t
        _zaehler["t0"] = _t.perf_counter()
    per_ring_states = [_ring_pucker_states(mol_with_conf, ring, frozen, tfd_thr)
                       for ring in rings]
    if _zaehler is not None:
        # ⚠ DIE ZEIT IST DIE EIGENTLICHE ANTWORT auf "bezahlbar?".  Sie zerfaellt in
        #   zwei Posten, die sich voellig verschieden skalieren: die Zustaende je Ring
        #   kosten LINEAR in der Ringzahl, das Kreuzprodukt kostet EXPONENTIELL.  Wer
        #   nur die Gesamtzeit misst, sieht den Unterschied nicht.
        _zaehler["t_zustaende"] = _t.perf_counter() - _zaehler["t0"]
        # (a) und (b) der Messung -- und die KOPPLUNGSZAHL dazu.  Gemeinsame Atome
        # zwischen zwei Ringen sind die unabhaengige Variable des ganzen Tests:
        # kondensiert = 2, spiro = 1, unabhaengig = 0.  Sie wird hier aus DERSELBEN
        # Ringliste gezaehlt, die der Bau benutzt -- nicht aus dem SMILES nachgeschaut.
        _zaehler["ringgroessen"] = [len(r) for r in rings]
        _zaehler["zustaende_je_ring"] = [len(s) for s in per_ring_states]
        _p = 1
        for _s in per_ring_states:
            _p *= len(_s)
        _zaehler["kreuzprodukt"] = _p
        # ⚠ SUMME UND MAXIMUM SIND ZWEI VERSCHIEDENE AUSSAGEN, und nur das MAXIMUM
        #   benennt die KopplungsART.  Ein Paar teilt 0 Atome (getrennt), 1 (spiro),
        #   2 (kondensiert, eine gemeinsame Bindung) oder >= 3 (verbrueckt).  Die Summe
        #   ueber alle Paare waechst dagegen einfach mit der Ringzahl und verwechselt
        #   drei lose Ringe mit einem Kaefig.
        _gem, _gmax = 0, 0
        for _i in range(len(rings)):
            for _j in range(_i + 1, len(rings)):
                _n_ij = len(set(rings[_i]) & set(rings[_j]))
                _gem += _n_ij
                _gmax = max(_gmax, _n_ij)
        _zaehler["gemeinsame_atome"] = _gem
        _zaehler["gem_max"] = _gmax

    # cartesian product of state indices, deterministic order, budget-capped;
    # skip the all-base (identity) combination; fewest-changed rings first.
    import itertools as _it
    # ===== SYMMETRIEREDUKTION (27.08.2026, DELFIN_FFFREE_PUCKER_SYMM) ==============
    #
    # WARUM HIER.  Zwei Ringe, die unter der Automorphismengruppe des Molekuels in
    # EINER Bahn liegen, sind ununterscheidbar.  Dann ist (s1,s2) DIESELBE STRUKTUR
    # wie (s2,s1) -- das Kreuzprodukt erzeugt dort Doppelgaenger, die spaeter die
    # RMSD-Entdopplung wieder wegwirft, nachdem sie Relax und Clash-Tor bezahlt haben.
    #
    # GEMESSEN 27.08. mit ZWEI unabhaengigen Instrumenten:
    #   bauseitig     534 Systeme mit >=2 faltbaren Ringen -> 371 (69,5 %) mit Bahn
    #   kristallseitig (CCDC clean_v2, DELFIN-unabhaengig)
    #                1757 mit >=2 gefalteten Ringen -> 1393 (79,3 %) mit Bahn
    #
    # DER EIGENTLICHE GEWINN IST NICHT RECHENZEIT, SONDERN ABDECKUNG.  Der Deckel
    # dreissig Zeilen weiter unten (`combos[:budget]`) schneidet nach Faltungstiefe
    # ab -- bei >2 flexiblen Ringen wird der Zustand "alle gleichzeitig gefaltet" NIE
    # gebaut.  Schrumpft das Produkt unter das Budget, hoert die Kappe auf zu beissen,
    # und genau die tiefen Zustaende entstehen wieder.  Die Reduktion nimmt also
    # Doppelgaenger weg und gibt dafuer ECHTE Zustaende zurueck.
    #
    # ⛔ NIE-SCHLECHTER, UND ZWAR STRENG.  Zusammengelegt wird nur, was ein ECHTER
    # Automorphismus aufeinander abbildet -- nicht, was nur dieselbe Rangmultimenge
    # traegt.  Der Unterschied ist gemessen: von 27 handgepruften Ringpaaren waren
    # 26 echte Bahn und EINES nur ranggleich (3,7 %).  Haette ich die Rangmenge
    # genommen, waere dieses eine Paar zusammengelegt worden und ein realer Zustand
    # HAETTE GEFEHLT.  Bei Vollstaendigkeit als Nordstern ist das der schlimmere
    # Fehler von beiden.
    # ⛔ UND WENN DIE AUTOMORPHISMENSUCHE NICHT TRAEGT, wird NICHT reduziert (voller
    # Produktzweig).  Fallback ist immer die groessere Menge, nie die kleinere.
    # ⛔ Vorgabe AUS -> byte-identisch.
    _symm = _os.environ.get("DELFIN_FFFREE_PUCKER_SYMM", "0") == "1"
    _bahnen = None
    if _symm and len(rings) > 1:
        _bahnen = _ring_bahnen(mol_with_conf, rings)
    if _bahnen and any(len(b) > 1 for b in _bahnen):
        # Je Bahn MULTIMENGEN statt geordneter Tupel: aus n^k wird C(n+k-1, k).
        # Die Zustandslisten einer Bahn sind gleich lang (gleiche Ringgroesse,
        # gleiche Umgebung) -- geprueft, sonst faellt die Bahn zurueck auf Produkt.
        _pro_bahn = []
        for _b in _bahnen:
            _n = len(per_ring_states[_b[0]])
            if any(len(per_ring_states[i]) != _n for i in _b):
                _pro_bahn.append(list(_it.product(*[range(len(per_ring_states[i]))
                                                    for i in _b])))
            else:
                _pro_bahn.append(list(_it.combinations_with_replacement(
                    range(_n), len(_b))))
        combos = []
        for _wahl in _it.product(*_pro_bahn):
            _c = [0] * len(rings)
            for _b, _w in zip(_bahnen, _wahl):
                for _ri, _st in zip(_b, _w):
                    _c[_ri] = _st
            if any(_c):
                combos.append(tuple(_c))
        if _zaehler is not None:
            _zaehler["symm_bahnen"] = [len(b) for b in _bahnen]
            _zaehler["symm_voll"] = _prod_laenge(per_ring_states)
            _zaehler["symm_reduziert"] = len(combos)
        if _os.environ.get("DELFIN_FFFREE_PUCKER_TRACE", "0") == "1":
            print("[pucker] SYMMETRIE: %d Ringe in %d Bahnen %s -- "
                  "Kombinationen %d statt %d"
                  % (len(rings), len(_bahnen), [len(b) for b in _bahnen],
                     len(combos), _prod_laenge(per_ring_states) - 1))
    else:
        combos = [c for c in _it.product(*[range(len(s)) for s in per_ring_states])
                  if any(c)]
    combos.sort(key=lambda c: (sum(1 for x in c if x), c))
    # ===== DIE ABSCHNEIDUNG WAR DER SCHADEN, NICHT DER DECKEL (26.08.2026) ==========
    #
    # `combos[:budget]` nimmt die ERSTEN 48 einer Liste, die nach "wie viele Ringe
    # weichen vom Grundzustand ab" sortiert ist.  Das heisst: erst ALLE
    # Ein-Ring-Aenderungen, dann alle Zwei-Ring-Kombinationen -- und dann ist Schluss.
    # Der Zustand, in dem ALLE Ringe gleichzeitig gefaltet sind, wird bei mehr als
    # zwei flexiblen Ringen NIE gebaut.  Vier Ringe zu je drei Mulden sind 80
    # Kombinationen; 48 davon decken Tiefe 1 und 2 ab, Tiefe 3 und 4 fallen weg.
    #
    # GEMESSEN 26.08. auf 400 Systemen -- die Verteilung ist genau die eines
    # abgeschnittenen Produkts:
    #     realisierte Konformere je System: 1:25 · 2:57 · 3:26 · 4:7 · 5:4 · 6:2 · 7:1 · 9:1
    #     82 von 123 Systemen haben EIN ODER ZWEI Konformere, bei 2,18 Mulden je Ring.
    #
    # ⇒ Mit `DELFIN_FFFREE_PUCKER_FULL=1` wird NICHT abgeschnitten: das vollstaendige
    #   Produkt wird gebaut, jede Faltungstiefe entsteht.  Der Deckel bleibt als
    #   Vorgabe erhalten (byte-identisch), aber er ist ab jetzt NIE STILL -- was er
    #   wegwirft, steht unter der Spur.  Eine stillschweigend gekappte Abdeckung liest
    #   sich wie Vollstaendigkeit und ist keine.
    # ⚠ PREIS, ehrlich: das Produkt waechst exponentiell mit der Ringzahl (8 Ringe zu
    #   je 3 Mulden = 6561 Kombinationen, jede mit Relax und Clash-Tor).  Darum steht
    #   die Vollversion hinter einem Schalter und nicht in der Vorgabe.
    _n_voll = len(combos)
    # ══ ADD, NEVER REPLACE -- AN DER KAPPE, NICHT AN DER REDUKTION ═════════════
    # (04.09.2026, Nutzerentscheidung.)
    #
    # BEFUND (Register #305/#315).  `symmfold6k` sperrt an RONSIW mit
    # `broken_regressed`: 12 Frames vorher, 12 nachher -- und trotzdem EINES neu.
    #     off   SP-4-chelate-1-pucker r0:base+r1:1+r2:3+r3:base
    #     on    SP-4-chelate-1-pucker r0:3+r1:base+r2:4+r3:base
    #
    # ⛔ MEINE ERSTE ERKLAERUNG WAR FALSCH.  Ich hielt das fuer eine
    #    Ersetzungs-Entscheidung der Symmetriereduktion.  Es ist eine KAPPEN-
    #    Folge: die Reduktion macht Plaetze unter `combos[:budget]` frei, also
    #    rutschen ANDERE Kombinationen herein.  Die Reduktion ersetzt nichts --
    #    die Kappe tut es.  Deshalb sitzt die Wache HIER und nicht dort.
    #
    # DIE FOLGE.  Der Champion-Satz ist keine Teilmenge des reduzierten Satzes:
    # ein Zustand, den der Champion baute, faellt heraus, und wenn der
    # hereingerutschte schlechter einbettet, steigt `broken_frac` (RONSIW
    # 0,250 -> 0,333).
    #
    # DIE WACHE.  Erst der Satz, den der Champion gebaut haette, dann mit dem
    # reduzierten aufgefuellt.  Damit ist die Champion-Framemenge BY
    # CONSTRUCTION enthalten und die Achse never-worse.
    # ⚠️ PREIS, ehrlich: die Vereinigung wird bis zu doppelt so gross wie die
    #    Kappe.  Das IST der Punkt -- die Reduktion soll TIEFERE Zustaende
    #    erreichbar machen, nicht flachere verdraengen.  Wer die Kappe halten
    #    will, laesst die Wache aus.
    # ⚠️ Vorgabe AUS -> byte-identisch.
    if (_symm and _bahnen and any(len(b) > 1 for b in _bahnen)
            and _os.environ.get("DELFIN_FFFREE_PUCKER_SYMM_ADD", "0") == "1"
            and _os.environ.get("DELFIN_FFFREE_PUCKER_FULL", "0") != "1"):
        _champ = [c for c in _it.product(*[range(len(s)) for s in per_ring_states])
                  if any(c)]
        _champ.sort(key=lambda c: (sum(1 for x in c if x), c))
        _champ = _champ[:max(0, int(budget))]
        _gesehen = set(_champ)
        combos = _champ + [c for c in combos[:max(0, int(budget))]
                           if c not in _gesehen]
        if _zaehler is not None:
            _zaehler["symm_add_champion"] = len(_champ)
            _zaehler["symm_add_gesamt"] = len(combos)
    elif _os.environ.get("DELFIN_FFFREE_PUCKER_FULL", "0") == "1":
        pass                                    # alle Faltungen, keine Kappe
    else:
        combos = combos[:max(0, int(budget))]
    if _zaehler is not None:
        # ⚠ ZWEI VERSCHIEDENE ZAHLEN, und die Verwechslung waere der ganze Irrtum.
        #   `aufzaehlung` ist das Kreuzprodukt OHNE den Grundzustand -- was aufgezaehlt
        #   werden MUESSTE.  `gebaut` ist, was nach der Kappe wirklich durch Relax und
        #   Tor geht.  Nur die zweite Zahl kostet Rechenzeit, nur die erste ist die
        #   Vollstaendigkeitsfrage.
        _zaehler["aufzaehlung"] = _n_voll
        _zaehler["gebaut"] = len(combos)
    if len(combos) < _n_voll and _os.environ.get("DELFIN_FFFREE_PUCKER_TRACE", "0") == "1":
        print("[pucker] KOMBINATIONEN GEKAPPT: %d von %d gebaut, %d verworfen "
              "(%d Ringe, Mulden %s) -- DELFIN_FFFREE_PUCKER_FULL=1 baut alle"
              % (len(combos), _n_voll, _n_voll - len(combos), len(rings),
                 "x".join(str(len(s)) for s in per_ring_states)))

    acc = Chem.Mol(mol_with_conf)
    kept_ids = [acc.GetConformer().GetId()]
    out: List[Tuple[str, str]] = []
    # ===== (1) DEFEKTFILTER: DAS DRITTE TOR, ALS MENGENDIFFERENZ (26.08.2026) ========
    #
    # Kollision und Winkel stehen unten schon.  Was fehlt, ist die BINDUNGSLAENGE --
    # s. `_bindungs_ausreisser`: das Selbstgate haelt eine gedehnte Bindung fuer eine
    # nicht vorhandene und meldet sie nie.
    # ⚠ RANG WAERE HIER FALSCH.  Ein Defekt ist kein "schlechter", er ist ein "nicht
    #   real".  Deshalb Filter, kein Score -- und deshalb gemessen gegen den
    #   GRUNDZUSTAND: verworfen wird nur, was die Faltung NEU einbringt.  Ein Molekuel,
    #   das schon vor der Faltung eine ungewoehnliche Bindung fuehrt (Nitril, Carben,
    #   ein schlecht eingebetteter Kern), verliert damit nicht alle seine Faltungen.
    # ⛔ Vorgabe AUS -> das Tor wird nie befragt -> byte-identisch.
    _defekt = _os.environ.get("DELFIN_FFFREE_PUCKER_DEFEKT", "0") == "1"
    _basis_bind = _bindungs_ausreisser(mol_with_conf) if _defekt else frozenset()
    # ===== (2) KRISTALLOGRAPHISCHE UNUNTERSCHEIDBARKEIT (26.08.2026) =================
    #
    # Zwei Faltungen, deren GROESSTE Atomauslenkung unter der Aufloesung liegt, sind
    # EIN Eintrag im Manifold -- s. den Block bei `_kabsch_max_rmsd`.
    # ⚠ ZUSAETZLICH ZU TFD, nicht statt dessen, und das ist keine Vorsicht sondern eine
    #   Arbeitsteilung: TFD faltet die Molekuelsymmetrie mit und toetet damit die
    #   Nummerierungsdubletten; das Maximum kann das nicht (es hat keine Topologie) und
    #   toetet dafuer die unterschwelligen Dubletten, die TFD nicht sieht.  Der Versuch
    #   vom 26.08., TFD durch eine reine Geometriedistanz zu ERSETZEN, ist genau daran
    #   gescheitert (n=5: 3,3,3 -> 9,13,14).
    # ⚠ REIHENFOLGE: das Maximum steht VOR dem TFD, weil es billiger ist -- Kabsch auf
    #   den Schweratomen gegen einen Torsionsfingerabdruck gegen alle Behaltenen.  Am
    #   Ergebnis aendert die Reihenfolge nichts: behalten wird, was BEIDE Pruefungen
    #   besteht, und die Menge der Behaltenen waechst in beiden Reihenfolgen gleich.
    # ⛔ Vorgabe AUS -> byte-identisch.
    _xrd = _os.environ.get("DELFIN_FFFREE_PUCKER_XRD", "0") == "1"
    _xrd_tol = float(_os.environ.get("DELFIN_FFFREE_PUCKER_XRDTOL", "0.15") or 0.15)
    _schwer: List[int] = []
    _xrd_kept: List = []
    if _xrd:
        try:
            _schwer = [i for i in range(mol_with_conf.GetNumAtoms())
                       if mol_with_conf.GetAtomWithIdx(i).GetSymbol() != "H"]
            _xrd_kept = [mol_with_conf.GetConformer().GetPositions()[_schwer]]
        except Exception:
            _xrd = False
    if _zaehler is not None:
        _zaehler["t0"] = _t.perf_counter()
    for combo in combos:
        try:
            m2 = Chem.Mol(mol_with_conf)
            conf = m2.GetConformer()
            active = False
            for ri_i, st_i in enumerate(combo):
                if st_i == 0:
                    continue
                st = per_ring_states[ri_i][st_i]
                if st is None:
                    continue
                # ⚠ ZWEI ZUSTANDSFORMEN, und das blinde Entpacken war eine Falle.
                #   Legacy: (q_scale, theta, phi) -- drei Werte.
                #   Raum:   (qs, phis) -- zwei Abbildungen m -> Wert.
                #   Ein `_qs2, theta, phi = st` auf die Raumform wirft ValueError, und
                #   der umgebende `except Exception: continue` haette das STILL
                #   verschluckt: jede Mehrring-Kombination waere lautlos ausgefallen
                #   und der Lauf haette "keine Wirkung" gemeldet.  Genau die Bauform,
                #   die heute schon dreimal eine Nullmessung erzeugt hat.
                if len(st) == 2 and isinstance(st[0], dict):
                    _set_pucker_general(conf, rings[ri_i], st[0], st[1], frozen)
                else:
                    _qs2, theta, phi = st
                    _set_pucker(conf, rings[ri_i],
                                _qs2 * _amp(len(rings[ri_i])), theta, phi, frozen)
                active = True
            if not active:
                continue
            # relax holding EVERY ring's pucker; inter-ring bonds/torsions free so
            # a clash between two puckered rings is relieved without collapsing
            # the puckers.
            if not _relax_hold_pucker_multi(m2, rings, frozen):
                if _zaehler is not None:
                    _zaehler["relax_fehler"] = _zaehler.get("relax_fehler", 0) + 1
                continue
            # realism gate: a combination that stayed clashed OR left any VSEPR
            # body distorted (fused/bridged rings strain their shared atoms) is
            # not a physical ensemble member -> drop it.  Everything must be
            # right, or the frame is unrealistic.
            if _zaehler is None:
                if _has_clash(m2) or _has_bad_angles(m2, skip=angle_skip):
                    continue
                if _defekt and (_bindungs_ausreisser(m2) - _basis_bind):
                    continue                    # NEU eingebrachter Bindungsdefekt
            else:
                # ⚠ IM MESSMODUS WERDEN BEIDE TORE GEFRAGT, im Vorgabepfad nicht.
                #   `or` ist kurzschluessig: feuert die Kollision, wird das Winkeltor
                #   NIE befragt -- die beiden Ursachen liessen sich dann nicht trennen.
                #   Genau ihre Trennung ist bei kondensierten Ringen der ganze Befund:
                #   dort teilen zwei Ringe Atome, das Winkeltor sieht die Spannung an
                #   den Fusionszentren, und das Kollisionstor ist per Konstruktion
                #   blind dafuer (s. `_has_bad_angles`).  Ein zweiter Toraufruf kostet
                #   Zeit -- darum nur, wenn jemand misst.
                _kl = _has_clash(m2)
                _wk = _has_bad_angles(m2, skip=angle_skip)
                _bi = bool(_bindungs_ausreisser(m2) - _basis_bind) if _defekt else False
                if _kl:
                    _zaehler["kollision"] = _zaehler.get("kollision", 0) + 1
                if _wk:
                    _zaehler["winkel"] = _zaehler.get("winkel", 0) + 1
                if _bi:
                    _zaehler["bindung"] = _zaehler.get("bindung", 0) + 1
                if _kl or _wk or _bi:
                    continue
                _zaehler["tor_ueberlebt"] = _zaehler.get("tor_ueberlebt", 0) + 1
            # (2) UNUNTERSCHEIDBARKEIT vor der Entdopplung -- s. den Block oben.
            _Pk = None
            if _xrd:
                try:
                    _Pk = m2.GetConformer().GetPositions()[_schwer]
                    if any(_kabsch_max_rmsd(_Pa, _Pk)[0] < _xrd_tol for _Pa in _xrd_kept):
                        if _zaehler is not None:
                            _zaehler["xrd_doppelt"] = _zaehler.get("xrd_doppelt", 0) + 1
                        continue
                except Exception:
                    _Pk = None
            cid = _add_conf(acc, m2)
            # Auf der Kombinationsebene falten MEHRERE Ringe gleichzeitig -- ringlokal
            # heisst hier "alle gefalteten Ringe, aber nur sie".  ⚠ Das ist eine
            # ANDERE Aussage als eine Stufe hoeher: hier faellt auch die exocyclische
            # Torsion aus dem Vergleich, zwei Kombinationen, die sich NUR in ihr
            # unterscheiden, werden also zusammengezogen.  Das ist gewollt (dieses
            # Modul faltet Ringe) und steht hier, damit es niemand spaeter als
            # Nebenwirkung entdeckt.
            # ⚠ EIGENER SCHALTER, weil die Wirkung hier das andere Vorzeichen hat als
            #   eine Stufe hoeher -- s. den Block bei `_tfd_distinct`.
            if not _tfd_distinct(acc, cid, kept_ids, tfd_thr, ringe=rings,
                                 schalter="DELFIN_FFFREE_PUCKER_TFD_LOCAL_KOMBI"):
                acc.RemoveConformer(cid)
                if _zaehler is not None:
                    _zaehler["tfd_doppelt"] = _zaehler.get("tfd_doppelt", 0) + 1
                continue
            kept_ids.append(cid)
            if _xrd and _Pk is not None:
                _xrd_kept.append(_Pk)
            if _zaehler is not None:
                # ⚠ DIE ENERGIE IST DIE ZWEITE ANTWORT auf dieselbe Frage.  Bleibt (d)
                #   gross, muss die Auswahl ueber ENERGIE laufen und nicht ueber RMSD
                #   (Nutzerregel).  Damit dieser Satz nicht nur eine Absicht ist, steht
                #   hier die Zahl: UFF-Energie des fertigen Frames, in kcal/mol, in der
                #   Reihenfolge der behaltenen Konformere.
                #   ⚠ UFF ist hier eine ORDNUNG, keine Thermochemie -- dieselbe Kraft,
                #     die auch relaxiert hat, also wenigstens in sich konsistent.  Wer
                #     daraus Populationen macht, ueberdehnt sie.
                try:
                    _zaehler.setdefault("energien", []).append(
                        float(AllChem.UFFGetMoleculeForceField(m2).CalcEnergy()))
                except Exception:
                    pass
            label = "pucker " + "+".join(
                f"r{ri_i}:{'base' if combo[ri_i] == 0 else combo[ri_i]}"
                for ri_i in range(len(rings)))
            out.append((_conf_to_xyz(m2), label))
        except Exception:
            if _zaehler is not None:
                # Der stille Ausfall bekommt eine Zahl.  Ohne sie liesse sich
                # "das Tor hat verworfen" nicht von "es ist etwas geplatzt"
                # unterscheiden -- zwei voellig verschiedene Befunde.
                _zaehler["ausnahme"] = _zaehler.get("ausnahme", 0) + 1
            continue
    if _zaehler is not None:
        _zaehler["t_kombis"] = _t.perf_counter() - _zaehler["t0"]
    return out


def selbsttest_raum() -> int:
    """Beweist, dass das CP-Gitter den Faltungsraum wirklich ueberdeckt.

    Aufruf:  python -m delfin.manta._ring_pucker
    Ohne diesen Test waere `_pucker_space_grid` eine Behauptung -- und der Fehler
    faellt in `generate` in ein `except Exception: continue`, also STILL.
    """
    import itertools as _itt
    fehler = 0
    print("=== Selbsttest: der Faltungsraum ===")

    # 1 DIMENSIONSZAHL.  Ein N-Ring hat genau N-3 Faltungsfreiheitsgrade.
    for n in range(4, 9):
        even = (n % 2 == 0)
        n_paare = len(range(2, (n // 2) if even else ((n - 1) // 2) + 1))
        dof = 2 * n_paare + (1 if even else 0)
        if dof != n - 3:
            print("  ✗ 1 DIMENSION n=%d: %d statt %d" % (n, dof, n - 3)); fehler += 1
    if not fehler:
        print("  ✓ 1 DIMENSION: N-3 Freiheitsgrade fuer N=4..8 (1,2,3,4,5)")

    # 2 DER PLANARE ZUSTAND ist im Gitter -- er fehlte der alten Liste (218 Motive).
    for n in (5, 6, 7):
        g = _pucker_space_grid(n, 2, 8)
        if not any(all(v == 0.0 for v in qs.values()) for qs, _ in g):
            print("  ✗ 2 PLANAR fehlt bei n=%d" % n); fehler += 1
    if fehler == 0:
        print("  ✓ 2 PLANAR: Q=0 ist Gitterpunkt fuer n=5,6,7")

    # 3 SESSEL UND INVERSER SESSEL.  Beim Sechsring die beiden Pole: q2=0, q3=+/-.
    g6 = _pucker_space_grid(6, 2, 8)
    pole = [qs for qs, _ in g6 if qs.get(2, 0.0) == 0.0 and qs.get(3, 0.0) != 0.0]
    if len([1 for qs in pole if qs[3] > 0]) < 1 or len([1 for qs in pole if qs[3] < 0]) < 1:
        print("  ✗ 3 POLE: Sessel/inv. Sessel nicht beide im Gitter"); fehler += 1
    else:
        print("  ✓ 3 POLE: Sessel UND inverser Sessel (q3 mit beiden Vorzeichen)")

    # 4 DER ZWISCHENBEREICH -- genau das, was die alte Liste NIE abtastete.
    #   Half-Chair/Envelope liegen zwischen Pol und Aequator: q2>0 UND q3!=0.
    zwischen = [qs for qs, _ in g6 if qs.get(2, 0.0) > 0 and qs.get(3, 0.0) != 0]
    if not zwischen:
        print("  ✗ 4 ZWISCHENBEREICH leer -- Half-Chair/Envelope unerreichbar"); fehler += 1
    else:
        print("  ✓ 4 ZWISCHENBEREICH: %d Punkte mit q2>0 UND q3!=0" % len(zwischen))

    # 5 HOEHERE PAARE ab n=7 -- ohne sie ist der Siebenring unvollstaendig.
    g7 = _pucker_space_grid(7, 2, 8)
    if not any(qs.get(3, 0.0) != 0.0 for qs, _ in g7):
        print("  ✗ 5 m=3 fehlt beim Siebenring"); fehler += 1
    else:
        print("  ✓ 5 HOEHERE PAARE: m=3 wird beim Siebenring belegt")

    # 6 KEINE DOPPELTEN Gitterpunkte (sonst blaeht das Produkt ohne Gewinn).
    for n in (5, 6, 7, 8):
        g = _pucker_space_grid(n, 2, 6)
        keys = [tuple(sorted((m, round(q, 6)) for m, q in qs.items())) for qs, _ in g]
        print("     n=%d: %4d Kandidaten (%d-dim)" % (n, len(g), n - 3))

    # 7 DER MAKROZYKLUS DARF DEN ZWEIG NICHT ZUM STEHEN BRINGEN.
    #   Ohne Budget waeren es bei n=16 rund 1,2e8 Kandidaten JE RING, jeder mit Relax
    #   und Kollisionstor -- der Lauf stirbt und liefert NULL Faltungen.  Undurchfuehr-
    #   barkeit ist das Gegenteil von Vollstaendigkeit.  Geprueft wird zweierlei:
    #   die Zahl bleibt unter dem Budget, UND der Ring wird trotzdem gefaltet
    #   (m=2 behaelt volle Phasenaufloesung, es faellt nur die feine Kraeuselung).
    print("     -- Makrozyklen (Budget %s) --"
          % _os.environ.get("DELFIN_FFFREE_PUCKER_BUDGET", "50000"))
    _budget_soll = max(1, int(_os.environ.get("DELFIN_FFFREE_PUCKER_BUDGET", "50000") or 50000))
    for n in (12, 16, 21):
        g = _pucker_space_grid(n, 2, 8)
        res = getattr(_pucker_space_grid, "_grid_res", {}) or {}
        if len(g) > _budget_soll:
            print("  ✗ 7 n=%d: %d Kandidaten UEBER Budget %d" % (n, len(g), _budget_soll))
            fehler += 1
            continue
        # m=2 ist die dominante Falte und muss ihre volle Phasenzahl behalten
        if res.get("n_phase_je_m", {}).get(2) != 8:
            print("  ✗ 7 n=%d: m=2 wurde reduziert (%s) -- die dominante Falte"
                  % (n, res.get("n_phase_je_m", {}).get(2)))
            fehler += 1
            continue
        if not any(qs.get(2, 0.0) > 0 for qs, _ in g):
            print("  ✗ 7 n=%d: keine einzige gefaltete Konfiguration" % n)
            fehler += 1
            continue
        # der Alternierungsterm ist beim geraden Ring der Sessel -- er darf nie fehlen
        if n % 2 == 0 and not any(qs.get(n // 2, 0.0) != 0 for qs, _ in g):
            print("  ✗ 7 n=%d: Alternierungsterm q_%d fehlt -- kein Sessel" % (n, n // 2))
            fehler += 1
            continue
        print("     n=%2d: %6d Kandidaten (%2d-dim), Phasen je m %s, ruhend %s"
              % (n, len(g), n - 3, res.get("n_phase_je_m"),
                 res.get("ruhende_moden") or "keine"))
    if fehler == 0:
        print("  ✓ 7 MAKROZYKLEN: unter Budget, m=2 voll aufgeloest, Faltung erreichbar")

    print("=== Faltungsraum: %s ===" % ("BESTANDEN" if fehler == 0 else "%d FEHLER" % fehler))
    return 1 if fehler else 0


def selbsttest_tfd_sweep(sizes=(5, 6, 7, 8)) -> int:
    """DIE SCHWELLE, NICHT DAS INSTRUMENT.  Splittet TFD bei 0,05 ueber?

    VORGESCHICHTE.  `selbsttest_konvergenz` meldet fuer n=6 und n=7, dass die Zahl
    der Zustaende mit der Aufloesung weiter waechst (8->9->10 bzw. 11->11->12), bei
    kleinsten CP-Abstaenden von 5,0 und 1,8 Grad.  Zwei Ursachen sind moeglich und
    haben ENTGEGENGESETZTE Reparaturen:
        (a) echte Mulden, Gitter zu grob   -> feiner abtasten
        (b) Uebersplittung durch TFD       -> Schwelle anheben
    Der erste Versuch, das ueber CP-Entdopplung zu klaeren, ist GESCHEITERT und die
    Messung steht: n=5 ging von 3,3,3 auf 9,13,14.  TFD faltet die MOLEKUELSYMMETRIE
    mit, die CP-Distanz nicht -- phi haengt an der Ringnummerierung.  CP ist kein
    Ersatz.  Also bleibt genau dieser Weg: dasselbe Instrument, andere Schwelle.

    WARUM DAS UEBER DIE KOMBINATORIK ENTSCHEIDET.  Gemessen 4,3 Ringe je System.
    Das Kreuzprodukt ueber die Ringe waechst wie (Zustaende je Ring)^(Ringe):
        3 Zustaende, 4 Ringe  ->      81 Kombinationen   rechenbar
       10 Zustaende, 4 Ringe  ->  10 000                 nicht rechenbar
    Ist 0,05 zu fein, wird die vollstaendige Kombinatorik dadurch bezahlbar -- ohne
    dass irgendwo abgeschnitten wird.  Ist sie richtig, ist der Preis echt und die
    Kappe kaeme sonst durch die Hintertuer zurueck.

    ⚠ DIESER TEST URTEILT NICHT UEBER CHEMIE.  Er misst an UNSUBSTITUIERTEN Ringen,
      deren Symmetrie hoch ist; ein substituierter Ring hat legitim mehr Zustaende.
      Was er zeigt, ist die OBERGRENZE der Uebersplittung, nicht die Produktionszahl.
    """
    if not (_RDKIT and _np is not None):
        print("=== TFD-Sweep: RDKit fehlt, uebersprungen ==="); return 0
    print("=== Selbsttest: TFD-Schwellensweep (NPHASE=8, unsubstituierte Ringe) ===")
    schwellen = (0.02, 0.05, 0.10, 0.15, 0.20, 0.30)
    print("    Ring  " + "".join("%7.2f" % t for t in schwellen))
    _alt = {k: _os.environ.get(k) for k in
            ("DELFIN_FFFREE_PUCKER_SPACE", "DELFIN_FFFREE_PUCKER_NPHASE",
             "DELFIN_FFFREE_PUCKER_NAMP", "DELFIN_FFFREE_PUCKER_TRACE",
             "DELFIN_FFFREE_PUCKER_CPDEDUP")}
    tabelle = {}
    try:
        _os.environ["DELFIN_FFFREE_PUCKER_SPACE"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_NAMP"] = "2"
        _os.environ["DELFIN_FFFREE_PUCKER_NPHASE"] = "8"
        _os.environ["DELFIN_FFFREE_PUCKER_TRACE"] = "0"
        _os.environ.pop("DELFIN_FFFREE_PUCKER_CPDEDUP", None)   # reines TFD messen
        for n in sizes:
            try:
                m = Chem.AddHs(Chem.MolFromSmiles("C1" + "C" * (n - 1) + "1"))
                if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
                    print("    n=%d  Einbettung fehlgeschlagen" % n); continue
                AllChem.MMFFOptimizeMolecule(m)
                ri = m.GetRingInfo().AtomRings()
                if not ri:
                    print("    n=%d  kein Ring gefunden" % n); continue
                ring = _ring_order(m, set(ri[0]))
            except Exception as e:
                print("    n=%d  Aufbau fehlgeschlagen: %s" % (n, type(e).__name__)); continue
            zeile = []
            for thr in schwellen:
                try:
                    zeile.append(len(_ring_pucker_states(m, ring, set(), thr)))
                except Exception:
                    zeile.append(-1)
            tabelle[n] = zeile
            print("    n=%-3d " % n + "".join("%7d" % z for z in zeile))
    finally:
        for k, v in _alt.items():
            if v is None:
                _os.environ.pop(k, None)
            else:
                _os.environ[k] = v

    if not tabelle:
        print("=== TFD-Sweep: nichts gemessen ==="); return 1

    # ---- WAS DAS KOSTET.  4,3 Ringe je System, gemessen auf 400 Systemen.
    print()
    print("    Kreuzprodukt bei 4 Ringen je System (Zustaende^4):")
    i05 = schwellen.index(0.05)
    for n, zeile in sorted(tabelle.items()):
        s05, s10 = zeile[i05], zeile[schwellen.index(0.10)]
        print("      n=%-3d  Schwelle 0,05 -> %8d      Schwelle 0,10 -> %8d"
              % (n, max(0, s05) ** 4, max(0, s10) ** 4))

    # ---- URTEIL.  Ueber-Splittung heisst: die Zahl faellt stark und BLEIBT dann flach.
    #      Faellt sie gleichmaessig weiter, verschmilzt die hoehere Schwelle echte
    #      Mulden -- dann ist nicht 0,05 zu fein, sondern die Schwelle das falsche
    #      Werkzeug.  Genau diese Unterscheidung ist der Sinn des Sweeps.
    print()
    print("    URTEIL je Ringgroesse:")
    verdacht = 0
    for n, zeile in sorted(tabelle.items()):
        if min(zeile) < 0 or zeile[i05] <= 0:
            print("      n=%-3d  nicht messbar" % n); continue
        _sturz = 1.0 - (zeile[i05 + 1] / float(zeile[i05]))       # 0,05 -> 0,10
        _rest = 1.0 - (zeile[-1] / float(max(1, zeile[i05 + 1])))  # 0,10 -> 0,30
        if _sturz >= 0.34 and _rest <= _sturz:
            print("      n=%-3d  UEBERSPLITTUNG: %d -> %d bei 0,05 -> 0,10 (%.0f %%), "
                  "danach nur noch %.0f %% -- der Sturz sitzt AN der Schwelle"
                  % (n, zeile[i05], zeile[i05 + 1], 100 * _sturz, 100 * _rest))
            verdacht += 1
        elif _sturz < 0.15:
            print("      n=%-3d  STABIL: %d -> %d (%.0f %%) -- 0,05 splittet NICHT ueber"
                  % (n, zeile[i05], zeile[i05 + 1], 100 * _sturz))
        else:
            print("      n=%-3d  GLEITEND: %.0f %% dann %.0f %% -- die Schwelle verschmilzt "
                  "fortlaufend, also auch ECHTE Mulden.  Kein sauberer Schnittpunkt."
                  % (n, 100 * _sturz, 100 * _rest))
    print("=== TFD-Sweep: %d von %d Ringgroesse(n) mit Uebersplittungsverdacht ==="
          % (verdacht, len(tabelle)))
    return 0


def selbsttest_konvergenz(sizes=(5, 6, 7)) -> int:
    """KONVERGENZ statt Behauptung: waechst die Zahl der Minima noch mit der Aufloesung?

    DIE FRAGE, die das beantwortet.  Der Parameterraum (q_m, phi_m) ist KONTINUIERLICH
    -- jede reelle Kombination ist eine gueltige Geometrie.  Der KONFORMERraum ist es
    nicht: ein Ring hat endlich viele Energieminima.  Das Gitter ist darum kein
    Ergebnis, sondern eine STARTPUNKTverteilung; `_relax_hold_pucker` zieht jeden
    Punkt ins naechste echte Minimum, TFD entdoppelt.

    ⇒ Vollstaendigkeit ist ERREICHBAR, nicht nur annaeherbar: das Gitter muss fein
      genug sein, dass jedes Einzugsgebiet mindestens einmal getroffen wird.  Ob das
      der Fall ist, sagt genau eine Messung -- die Zahl der distinkten Zustaende
      gegen die Aufloesung.  Waechst sie nicht mehr, ist der Raum ueberdeckt.

    ⚠ WARUM DAS HIER STEHT.  `conformer_enum.py:7-9` behauptet dasselbe ("finer grid
    stops adding distinct minima") und hat KEINE Messstelle dafuer.  Eine Behauptung
    ohne Beleg ist genau die Bauform, die in diesem Projekt schon mehrfach eine
    falsche Zahl getragen hat.
    """
    if not (_RDKIT and _np is not None):
        print("=== Konvergenz: RDKit fehlt, uebersprungen ==="); return 0
    print("=== Selbsttest: Konvergenz des Faltungsgitters ===")
    print("    Ring   NPHASE=4   8   16     konvergiert?")
    _alt = {k: _os.environ.get(k) for k in
            ("DELFIN_FFFREE_PUCKER_SPACE", "DELFIN_FFFREE_PUCKER_NPHASE",
             "DELFIN_FFFREE_PUCKER_NAMP", "DELFIN_FFFREE_PUCKER_TRACE")}
    fehler = 0
    try:
        _os.environ["DELFIN_FFFREE_PUCKER_SPACE"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_NAMP"] = "2"
        _os.environ["DELFIN_FFFREE_PUCKER_TRACE"] = "0"
        # Mit CP-Entdopplung gegenrechnen, wenn der Aufrufer sie gesetzt hat --
        # sonst misst der Test die alte Uebersplittung nach.
        # ⚠ EINMAL LESEN, EINMAL BENENNEN.  Die erste Fassung las den Schalter hier
        #   und nannte ihn unten `_cpd_an` -- ein Name, den es nie gab.  Der
        #   NameError fiel in das `except Exception` der CP-Streuung und wurde als
        #   "nicht messbar" gedruckt: das moduskorrigierte Urteil lief damit KEIN
        #   einziges Mal, und der Test meldete trotzdem etwas.  Ein verschluckter
        #   Fehler ist eine Nullmessung, die wie ein Befund aussieht.
        _cpd_an = _os.environ.get("DELFIN_FFFREE_PUCKER_CPDEDUP") == "1"
        if _cpd_an:
            print("    (CP-Entdopplung AN, Toleranz %s Grad / %s A)"
                  % (_os.environ.get("DELFIN_FFFREE_PUCKER_CPTOL", "15"),
                     _os.environ.get("DELFIN_FFFREE_PUCKER_CPQTOL", "0.15")))
        for n in sizes:
            smi = "C1" + "C" * (n - 1) + "1"
            try:
                m = Chem.AddHs(Chem.MolFromSmiles(smi))
                if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
                    print("    n=%d  Einbettung fehlgeschlagen" % n); continue
                AllChem.MMFFOptimizeMolecule(m)
                ri = m.GetRingInfo().AtomRings()
                if not ri:
                    print("    n=%d  kein Ring gefunden" % n); continue
                ring = _ring_order(m, set(ri[0]))
            except Exception as e:
                print("    n=%d  Aufbau fehlgeschlagen: %s" % (n, type(e).__name__)); continue
            zahlen = []
            for nph in (4, 8, 16):
                _os.environ["DELFIN_FFFREE_PUCKER_NPHASE"] = str(nph)
                try:
                    st = _ring_pucker_states(m, ring, set(), 0.05)
                    zahlen.append(len(st))
                except Exception as e:
                    zahlen.append(-1)
            ok = (len(zahlen) == 3 and zahlen[1] > 0 and zahlen[2] <= zahlen[1])
            print("    n=%-3d  %8d %3d %4d      %s"
                  % (n, zahlen[0], zahlen[1], zahlen[2],
                     "JA" if ok else "NEIN -- siehe CP-Streuung"))
            if not ok:
                fehler += 1
                # ---- WARUM waechst die Zahl?  Zwei Ursachen, ENTGEGENGESETZTE Fixes --
                # (1) Raum noch nicht ueberdeckt -> die Zustaende liegen in CP-
                #     Koordinaten WEIT auseinander -> feiner abtasten.
                # (2) TFD-Schwelle trennt chemisch GLEICHE Zustaende -> sie liegen
                #     DICHT beieinander -> die Lupe ist zu fein, nicht das Gitter grob.
                #
                # ⚠️ DAS IST KEINE AKADEMISCHE FRAGE.  Die Zahl der Zustaende JE RING
                #    geht als Basis in das Kreuzprodukt ueber alle Ringe ein:
                #        3 Zustaende, 4 Ringe ->      81 Kombinationen  (rechenbar)
                #       10 Zustaende, 4 Ringe ->  10 000                (nicht rechenbar)
                #    Gemessen sind 4,3 Ringe je System.  Uebersplittung macht die
                #    VOLLSTAENDIGE Kombinatorik unbezahlbar -- Konvergenzanomalie und
                #    Rechenbarkeit sind dasselbe Problem.
                # Chemischer Massstab: Cyclohexan hat Sessel + Twist-Boat-Familie,
                # nach Symmetriefaltung 2-3 Klassen.  Der Fuenfring konvergiert auf 3.
                _os.environ["DELFIN_FFFREE_PUCKER_NPHASE"] = "16"
                try:
                    cps = []
                    for _s in _ring_pucker_states(m, ring, set(), 0.05):
                        if _s is None or not (isinstance(_s, tuple) and len(_s) == 2
                                              and isinstance(_s[0], dict)):
                            continue
                        m3 = Chem.Mol(m)
                        _set_pucker_general(m3.GetConformer(), ring, _s[0], _s[1], set())
                        _relax_hold_pucker(m3, ring, set())
                        cps.append(_cp_theta_phi(m3.GetConformer().GetPositions(), ring))
                    cps.sort(key=lambda t: (round(t[1], 0), round(t[2], 0)))
                    print("        CP der ueberlebenden Zustaende (Q, theta, phi):")
                    for _Q, _th, _ph in cps:
                        print("          Q=%.3f  theta=%6.1f  phi=%6.1f" % (_Q, _th, _ph))
                    dmin = None
                    for _i in range(len(cps)):
                        for _j in range(_i + 1, len(cps)):
                            _, ti, pi_ = cps[_i]
                            _, tj, pj = cps[_j]
                            dph = min(abs(pi_ - pj), 360.0 - abs(pi_ - pj))
                            d = ((ti - tj) ** 2 + dph ** 2) ** 0.5
                            dmin = d if dmin is None else min(dmin, d)
                    if dmin is not None:
                        print("        kleinste paarweise CP-Distanz: %.1f Grad" % dmin)
                        # ⚠ DAS URTEIL MUSS WISSEN, WELCHER MODUS LIEF.  Erste Fassung
                        #   war fest an `dmin` gekoppelt und behauptete "TFD zu fein"
                        #   auch dann, wenn die CP-Entdopplung lief -- also ein Urteil
                        #   ueber ein Instrument, das gar nicht im Einsatz war.
                        if _cpd_an:
                            print("        URTEIL: CP-Entdopplung laeuft und liefert MEHR"
                                  " Zustaende als TFD.  Grund: TFD faltet die MOLEKUEL-"
                                  "SYMMETRIE mit, die CP-Distanz nicht.  phi haengt an der"
                                  " Ringnummerierung -- bei einem unsubstituierten Ring"
                                  " sind alle phi bei gleichem (Q, theta) DERSELBE"
                                  " Konformer.  CP allein ist KEIN Ersatz fuer TFD.")
                        else:
                            print("        URTEIL: %s" % (
                                "ECHTE Mulden -- Gitter zu grob, feiner abtasten"
                                if dmin > 20.0 else
                                "Zustaende liegen dicht -- Uebersplittung moeglich; "
                                "PRUEFEN durch TFD-Schwellensweep, NICHT durch Ersetzen "
                                "von TFD (siehe CP-Modus)"))
                except Exception as _e:
                    print("        CP-Streuung nicht messbar: %s" % type(_e).__name__)
    finally:
        for k, v in _alt.items():
            if v is None:
                _os.environ.pop(k, None)
            else:
                _os.environ[k] = v
    print("=== Konvergenz: %s ===" %
          ("BESTANDEN -- der Raum ist bei NPHASE=8 ueberdeckt" if fehler == 0
           else "%d Ringgroesse(n) NICHT konvergiert" % fehler))
    return 1 if fehler else 0


# ===== DIE PROBEN: EIN KOPPLUNGSGRADIENT, KEINE SAMMLUNG ==========================
#
# Die Frage ist nicht "wie viele Faltungen hat Molekuel X", sondern WOVON es abhaengt,
# wie viel vom Kreuzprodukt uebrig bleibt.  Das kann nur eine VARIABLE beantworten, die
# von Probe zu Probe systematisch anders steht -- hier die Zahl der GEMEINSAMEN ATOME
# zwischen zwei Ringen.  Sie laeuft ueber die Liste von 2 nach 0:
#
#   2 gemeinsame Atome, drei Bruecken   verbrueckt (Bicyclo[2.2.2]octan, Norbornan)
#                                       -- der steifste Fall, den es gibt
#   2 gemeinsame Atome, eine Bindung    kondensiert (Decalin, Perhydroanthracen)
#   1 gemeinsames Atom                  spiro (Spiro[5.5]undecan)
#   0, direkte Ring-Ring-Bindung        nur STERISCH gekoppelt (Bicyclohexyl)
#   0, zwei CH2 dazwischen              praktisch unabhaengig (1,2-Dicyclohexylethan)
#   0, drei Ringe an einem P            Tricyclohexylphosphin -- der Fall, den der
#                                       `generate`-Docstring selbst als Beispiel
#                                       fuehrt, und ein echter Ligand
#
# Cyclohexan steht als NULLPUNKT dabei: EIN Ring, also gar kein Kreuzprodukt.  Ohne ihn
# waere nicht zu trennen, was die KOPPLUNG kostet und was schon der einzelne Ring kostet.
# Perhydroanthracen und das Phosphin sind die einzigen DREIringigen Proben -- erst bei
# drei Ringen zeigt sich, ob die Kurve exponentiell oder gedeckelt laeuft.
_KOMBI_PROBEN = (
    ("Cyclohexan",            "C1CCCCC1",                      "1 Ring -- Nullpunkt"),
    ("Bicyclo[2.2.2]octan",   "C1CC2CCC1CC2",                  "verbrueckt, 3 Bruecken"),
    ("Norbornan",             "C1CC2CCC1C2",                   "verbrueckt, 1 Bruecke"),
    ("Decalin",               "C1CCC2CCCCC2C1",                "kondensiert, 1 Bindung"),
    ("Perhydroanthracen",     "C1CCC2CC3CCCCC3CC2C1",          "3 Ringe, kondensiert"),
    ("Spiro[5.5]undecan",     "C1CCC2(CC1)CCCCC2",             "spiro, 1 Atom"),
    ("Bicyclohexyl",          "C1CCCCC1C1CCCCC1",              "0 Atome, 1 Bindung"),
    ("1,2-Dicyclohexylethan", "C1CCCCC1CCC1CCCCC1",            "0 Atome, 2 CH2"),
    ("Tricyclohexylphosphin", "P(C1CCCCC1)(C1CCCCC1)C1CCCCC1", "3 Ringe, unabhaengig"),
)


def selbsttest_kombinatorik(proben=None, deckel_s: float = 900.0,
                            max_kombis: int = 4000) -> int:
    """WIE GROSS IST DAS KREUZPRODUKT **NACH** DER PHYSIK?

    DIE FRAGE, die ueber die vollstaendige Ringfaltung entscheidet.  Gemessen sind 4,3
    Ringe je System und -- mit dem CP-Raumgitter -- 9 bis 16 Faltungszustaende je Ring
    (`selbsttest_tfd_sweep`, Schwelle 0,05).  Das naive Kreuzprodukt ist damit 6500 bis
    65000 Kombinationen je System.  Aber Ringe eines Molekuels sind NICHT unabhaengig:
    kondensierte und verbrueckte Ringe teilen Atome, faltet man den einen, ist der
    andere festgelegt.  Das Kreuzprodukt ist eine Obergrenze der AUFZAEHLUNG -- die
    Frage ist, was davon das Realismustor ueberlebt.

    ⚠ DER ENTSCHEIDENDE UNTERSCHIED, den dieser Test sichtbar macht: das Tor toetet das
      ERGEBNIS, aber nicht die KOSTEN.  Jede Kombination wird erst gesetzt, dann mit
      gehaltenen Faltungen relaxiert und ERST DANN verworfen.  Wer "(c) ist klein, also
      billig" liest, hat die Reihenfolge verwechselt.  Darum stehen hier ZWEI Zahlen
      nebeneinander: (b) ist der Preis, (d) ist die Ausbeute.

    ⚠ ZWEI TORE, NICHT EINS.  `_has_clash` sieht nur ueberlappende vdW-Kugeln.  Bei
      KONDENSIERTEN Ringen entsteht der Widerspruch aber an den geteilten Atomen, und
      dort stimmt der VSEPR-Winkel nicht mehr, ohne dass irgendetwas kollidiert --
      genau dafuer existiert `_has_bad_angles`.  Welches der beiden Tore feuert, ist
      deshalb selbst ein Befund und wird getrennt gezaehlt.

    ⚠ WAS DIESER TEST NICHT MISST.  Er laeuft auf METALLFREIEN Kohlenwasserstoffen mit
      hoher Symmetrie.  Ein substituierter Ring hat legitim mehr Zustaende, ein echtes
      DELFIN-System ist groesser und damit je Kombination teurer.  Die Zeiten hier sind
      eine UNTERGRENZE der Kosten, nicht die Produktionszahl.

    ``max_kombis``: Testgrenze.  (a) und (b) werden IMMER bestimmt -- sie kosten nur den
    linearen Posten.  Liegt (b) darueber, werden (c) und (d) NICHT gemessen und genau
    das wird gedruckt, statt eine gekappte Zahl auszugeben.  ⚠ Eine Kappe waere hier
    besonders heimtueckisch: `combos` ist nach FALTUNGSTIEFE sortiert, ein Praefix davon
    enthaelt nur flache Kombinationen und haette systematisch zu hohe Ueberlebensquoten.
    ``0`` = keine Grenze (dann kann ein einzelner Mehrringer Stunden laufen).
    """
    if not (_RDKIT and _np is not None):
        print("=== Kombinatorik: RDKit fehlt, uebersprungen ==="); return 0
    import time as _time
    proben = proben or _KOMBI_PROBEN
    print("=== Selbsttest: das Kreuzprodukt NACH dem Realismustor ===")
    print("    Raumgitter AN (NAMP=2, NPHASE=8), Kappe AUS -- die VOLLE Kombinatorik.")
    print("    (a) Zustaende je Ring · (b) Kreuzprodukt · (c) ueberlebt das Tor"
          " · (d) davon TFD-distinkt")
    _alt = {k: _os.environ.get(k) for k in
            ("DELFIN_FFFREE_PUCKER_SPACE", "DELFIN_FFFREE_PUCKER_NAMP",
             "DELFIN_FFFREE_PUCKER_NPHASE", "DELFIN_FFFREE_PUCKER_FULL",
             "DELFIN_FFFREE_PUCKER_TRACE")}
    zeilen = []
    fehler = 0
    try:
        # ===== 0 DER VORGABEPFAD MUSS UNVERAENDERT BLEIBEN -- GEMESSEN, NICHT BEHAUPTET
        #
        # `generate` traegt jetzt einen optionalen Zaehlersatz.  Die Behauptung "bei
        # `_zaehler=None` aendert sich nichts" ist genau die Sorte Behauptung, die in
        # diesem Projekt schon mehrfach falsch war.  Also wird sie gemessen: dasselbe
        # Molekuel, derselbe Vorgabepfad (alle Schalter AUS, Kappe AN), einmal ohne und
        # einmal mit Zaehler -- die zurueckgegebenen Frames muessen ZEICHENGLEICH sein.
        # ⚠ Im Messmodus werden beide Tore gefragt statt kurzschluessig eines; wuerde
        #   `_has_bad_angles` etwas veraendern, faellt es genau hier auf.
        for _k in _alt:
            _os.environ[_k] = "0"
        _mv = Chem.AddHs(Chem.MolFromSmiles("C1CCC2CCCCC2C1"))    # Decalin
        if AllChem.EmbedMolecule(_mv, randomSeed=42) == 0:
            AllChem.MMFFOptimizeMolecule(_mv)
            _ohne = generate(_mv, budget=48)
            _mit = generate(_mv, budget=48, _zaehler=_neuer_zaehler())
            if _ohne == _mit:
                print("    ✓ 0 VORGABE UNVERAENDERT: Decalin, Schalter AUS, %d Frames "
                      "mit und ohne Zaehler identisch" % len(_ohne))
            else:
                print("    ✗ 0 VORGABE VERAENDERT: %d Frames ohne Zaehler, %d mit -- "
                      "die Messstelle ist nicht folgenlos" % (len(_ohne), len(_mit)))
                fehler += 1
        else:
            print("    ? 0 VORGABE: Decalin nicht einbettbar, Identitaet NICHT gemessen")
            fehler += 1

        _os.environ["DELFIN_FFFREE_PUCKER_SPACE"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_NAMP"] = "2"
        _os.environ["DELFIN_FFFREE_PUCKER_NPHASE"] = "8"
        _os.environ["DELFIN_FFFREE_PUCKER_FULL"] = "1"     # keine Kappe -- ganzes Produkt
        _os.environ["DELFIN_FFFREE_PUCKER_TRACE"] = "0"
        print()
        print("    %-22s %2s %4s  %-14s %8s %7s %7s %7s %8s"
              % ("Molekuel", "R", "gmax", "(a) je Ring", "(b)", "(c)", "(d)",
                 "c/b", "Zeit/s"))
        print("    (gmax = Atome, die sich das ENGSTE Ringpaar teilt: 0 getrennt · "
              "1 spiro · 2 kondensiert · >=3 verbrueckt)")
        for name, smi, klasse in proben:
            try:
                m = Chem.AddHs(Chem.MolFromSmiles(smi))
                if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
                    print("    %-22s Einbettung fehlgeschlagen" % name); continue
                try:
                    AllChem.MMFFOptimizeMolecule(m)
                except Exception:
                    AllChem.UFFOptimizeMolecule(m)
            except Exception as e:
                print("    %-22s Aufbau fehlgeschlagen: %s" % (name, type(e).__name__))
                continue
            # ---- DIE TESTGRENZE: (a) und (b) ZUERST, getrennt vom Bau.  Der Bau kostet
            #      JE Kombination einen Relax plus zwei Tore; ob er bezahlbar ist,
            #      entscheidet (b) -- also muss (b) bekannt sein, BEVOR gebaut wird.
            #      ⚠ Der Vorlauf bestimmt die Zustaende ein zweites Mal (`generate` tut
            #        es gleich nochmal).  Das ist der LINEARE Posten, also der billige --
            #        aber bezahlt wird er trotzdem, darum laeuft er NUR, wenn die Grenze
            #        ueberhaupt gesetzt ist.  Bei `max_kombis=0` gibt es keinen Vorlauf
            #        und damit auch keine doppelte Arbeit in der Zeitmessung.
            if max_kombis:
                _t_vor = _time.perf_counter()
                try:
                    _rings = [_ring_order(m, set(r)) for r in m.GetRingInfo().AtomRings()
                              if _is_puckerable(m, r)]
                    _stv = [_ring_pucker_states(m, r, set(), 0.05) for r in _rings]
                except Exception as e:
                    print("    %-22s Zustaende nicht bestimmbar: %s"
                          % (name, type(e).__name__))
                    continue
                _prod = 1
                for _s in _stv:
                    _prod *= len(_s)
                _b_vor = max(0, _prod - 1)
                _gem_vor = max([len(set(_rings[i]) & set(_rings[j]))
                                for i in range(len(_rings))
                                for j in range(i + 1, len(_rings))] or [0])
                if _b_vor > max_kombis:
                    print("    %-22s %2d %4d  %-14s %8d %7s %7s %7s %8.1f"
                          % (name, len(_rings), _gem_vor,
                             "x".join(str(len(s)) for s in _stv), _b_vor,
                             "-", "-", "-", _time.perf_counter() - _t_vor))
                    print("        %-28s (c) und (d) NICHT GEMESSEN: (b) = %d ueber der "
                          "Testgrenze %d.  Kein gekappter Ersatzwert -- `combos` ist nach "
                          "Faltungstiefe sortiert, ein Praefix waere systematisch zu flach."
                          % (klasse, _b_vor, max_kombis))
                    zeilen.append({"name": name, "klasse": klasse, "ringe": len(_rings),
                                   "gem": _gem_vor, "gmax": _gem_vor,
                                   "b": 0, "c": 0, "d": 0, "dt": 0.0,
                                   "zustaende": [len(s) for s in _stv],
                                   "t_kombis": 0.0, "t_zust": 0.0,
                                   "kollision": 0, "winkel": 0, "energien": [],
                                   "ungemessen": _b_vor})
                    continue
            z = _neuer_zaehler()
            _t0 = _time.perf_counter()
            try:
                out = generate(m, budget=10 ** 9, _zaehler=z)
            except Exception as e:
                print("    %-22s generate() geplatzt: %s" % (name, type(e).__name__))
                continue
            _dt = _time.perf_counter() - _t0
            _b = int(z.get("aufzaehlung", 0))
            _c = int(z.get("tor_ueberlebt", 0))
            _d = len(out)
            # ⚠ NIE EIN PROZENTSATZ OHNE NENNER.  Der Nenner ist hier (b), die Zahl der
            #   aufgezaehlten Kombinationen ohne den Grundzustand -- nicht das
            #   Kreuzprodukt selbst, denn der Grundzustand wird nie gebaut.
            _cb = ("%6.1f%%" % (100.0 * _c / _b)) if _b else "   n/a"
            print("    %-22s %2d %4d  %-14s %8d %7d %7d %7s %8.1f"
                  % (name, len(z.get("ringgroessen") or []),
                     int(z.get("gem_max", 0)),
                     "x".join(str(v) for v in (z.get("zustaende_je_ring") or [])) or "-",
                     _b, _c, _d, _cb, _dt))
            print("        %-28s verworfen: Kollision %d · Winkel %d · TFD %d · "
                  "Relax %d · Ausnahme %d"
                  % (klasse, int(z.get("kollision", 0)), int(z.get("winkel", 0)),
                     int(z.get("tfd_doppelt", 0)), int(z.get("relax_fehler", 0)),
                     int(z.get("ausnahme", 0))))
            _en = sorted(z.get("energien") or [])
            if len(_en) >= 2:
                _e0 = _en[0]
                _in10 = sum(1 for e in _en if e - _e0 <= 10.0)
                print("        UFF-Energie der (d): Spanne %.1f kcal/mol · "
                      "innerhalb 10 kcal/mol %d von %d"
                      % (_en[-1] - _e0, _in10, len(_en)))
            zeilen.append({"name": name, "klasse": klasse,
                           "ringe": len(z.get("ringgroessen") or []),
                           "gem": int(z.get("gemeinsame_atome", 0)),
                           "gmax": int(z.get("gem_max", 0)),
                           "b": _b, "c": _c, "d": _d, "dt": _dt,
                           "zustaende": list(z.get("zustaende_je_ring") or []),
                           "t_kombis": float(z.get("t_kombis", 0.0)),
                           "t_zust": float(z.get("t_zustaende", 0.0)),
                           "kollision": int(z.get("kollision", 0)),
                           "winkel": int(z.get("winkel", 0)),
                           "energien": _en})
    finally:
        for k, v in _alt.items():
            if v is None:
                _os.environ.pop(k, None)
            else:
                _os.environ[k] = v

    if not zeilen:
        print("=== Kombinatorik: nichts gemessen ==="); return 1

    # ---- 1 KOPPLUNG GEGEN UEBERLEBEN.
    # ⚠ PARTITION, KEIN MITTELWERT.  "gekoppelt gegen unabhaengig" waere die falsche
    #   Zweiteilung: sie wirft ein kondensiertes Ringpaar (eine gemeinsame BINDUNG) mit
    #   einem Kaefig (vier gemeinsame Atome) in einen Topf, und deren Ueberlebensquoten
    #   liegen zwei Groessenordnungen auseinander.  Ein Mittelwert kann eine tote Klasse
    #   nicht sehen -- nur eine Partition kann das.  Geteilt wird darum nach der Zahl der
    #   Atome, die sich das ENGSTE Ringpaar teilt; das ist zugleich der chemische Name
    #   der Kopplung.
    mehr = [r for r in zeilen if r["ringe"] >= 2 and r["b"] > 0]

    def _klasse(r):
        g = r.get("gmax", 0)
        return 0 if g == 0 else (1 if g == 1 else (2 if g == 2 else 3))

    _NAMEN = {0: "0 Atome  getrennt", 1: "1 Atom   spiro",
              2: "2 Atome  kondensiert", 3: ">=3      verbrueckt"}

    def _quote(gruppe):
        _b = sum(r["b"] for r in gruppe)
        _c = sum(r["c"] for r in gruppe)
        return _b, _c, (100.0 * _c / _b if _b else 0.0)

    print()
    print("    ===== 1 KOPPLUNG GEGEN UEBERLEBEN (nur Mehrringer) =====")
    print("      engstes Ringpaar teilt ...")
    _quoten = {}
    for _kl in (0, 1, 2, 3):
        _g = [r for r in mehr if _klasse(r) == _kl]
        if not _g:
            print("      %-22s keine Probe" % _NAMEN[_kl]); continue
        _b, _c, _q = _quote(_g)
        _quoten[_kl] = _q
        # ⚠ (a) MUSS MIT DASTEHEN, sonst ist c/b nicht interpretierbar.  Ein verbrueckter
        #   Ring kann schon WENIGER Zustaende haben -- dann ist (b) klein, weil die
        #   Kopplung frueher gewirkt hat, und nicht, weil das Tor mehr toetet.  Zwei
        #   verschiedene Wege zum selben kleinen Produkt, und nur beide zusammen sagen,
        #   welcher es war.
        _zust = [v for r in _g for v in r["zustaende"]]
        print("      %-22s %d Probe(n) · %6d von %6d ueberleben = %5.1f %% · (a) im "
              "Mittel %.1f je Ring (%d Ringe) · %s"
              % (_NAMEN[_kl], len(_g), _c, _b, _q,
                 (sum(_zust) / float(len(_zust))) if _zust else 0.0, len(_zust),
                 ", ".join(r["name"] for r in _g)))
    if len(_quoten) >= 2:
        _hi = max(_quoten.values())
        _lo = min(_quoten.values())
        # ⚠ Der Befund ist die SPANNE ueber die Partition, nicht ein Gruppenmittel.
        print("      ⇒ Spanne ueber die Kopplungsklassen: %.1f %% bis %.1f %% -- %s"
              % (_lo, _hi,
                 "die Kopplungsart entscheidet, nicht die Kopplung an sich"
                 if _hi - _lo > 20.0 else
                 "die Kopplungsart macht kaum einen Unterschied"))

    # ---- 2 WELCHES TOR FEUERT.  Kollision und Winkel getrennt, sonst ist "das Tor"
    #      ein Name fuer zwei verschiedene Mechanismen (Detektorname != Messung).
    print()
    print("    ===== 2 WELCHES TOR TOETET =====")
    for _kl in (0, 1, 2, 3):
        _g = [r for r in mehr if _klasse(r) == _kl]
        if not _g:
            continue
        _b = sum(r["b"] for r in _g)
        _k = sum(r["kollision"] for r in _g)
        _w = sum(r["winkel"] for r in _g)
        print("      %-22s von %6d Kombinationen: Kollision %6d (%5.1f %%) · "
              "Winkel %6d (%5.1f %%)"
              % (_NAMEN[_kl], _b, _k, 100.0 * _k / _b if _b else 0.0,
                 _w, 100.0 * _w / _b if _b else 0.0))
    # ⚠ EIN TORNAME IST KEINE MESSUNG.  Wenn "das Kollisionstor" in Wahrheit nie feuert
    #   und die ganze Selektion vom Winkeltor kommt, dann steht jede Aussage ueber "die
    #   Sterik schneidet das Produkt" auf dem falschen Mechanismus -- und eine Reparatur
    #   am Kollisionstor waere wirkungslos, bevor sie geschrieben ist.
    _kges = sum(r["kollision"] for r in mehr)
    _wges = sum(r["winkel"] for r in mehr)
    _bges = sum(r["b"] for r in mehr)
    if _bges:
        if _kges == 0 and _wges > 0:
            print("      ⇒ DAS KOLLISIONSTOR HAT NULL REICHWEITE: 0 von %d Kombinationen."
                  % _bges)
            print("        Der Filter ist AUSSCHLIESSLICH das WINKELTOR (%d von %d = "
                  "%.1f %%).  Wer die Kombinatorik am Kollisionstor beschneiden will, "
                  "greift den Mechanismus an, der gar nicht feuert."
                  % (_wges, _bges, 100.0 * _wges / _bges))
        else:
            print("      ⇒ Kollision %d von %d (%.1f %%) · Winkel %d von %d (%.1f %%) "
                  "-- beide Tore tragen."
                  % (_kges, _bges, 100.0 * _kges / _bges,
                     _wges, _bges, 100.0 * _wges / _bges))

    # ---- 3 BEZAHLBARKEIT.  Die Kosten haengen an (b), nicht an (d).
    _sum_b = sum(r["b"] for r in zeilen)
    _sum_t = sum(r["t_kombis"] for r in zeilen)
    print()
    print("    ===== 3 BEZAHLBARKEIT =====")
    if _sum_b <= 0 or _sum_t <= 0.0:
        print("      Kosten je Kombination NICHT MESSBAR (b=%d, t=%.3f s)"
              % (_sum_b, _sum_t))
        print("=== Kombinatorik: unvollstaendig ==="); return 1
    _ms = 1000.0 * _sum_t / _sum_b
    print("      Jede der (b) Kombinationen wird GEBAUT und RELAXIERT, bevor das Tor")
    print("      sie verwirft -- das Tor spart nichts, es waehlt nur aus.")
    print("      Gemessen: %.1f ms je Kombination (Nenner: %d Kombinationen ueber %d "
          "Proben, %.1f s gesamt)" % (_ms, _sum_b, len(zeilen), _sum_t))
    print("      ⚠ UNTERGRENZE: metallfreie Kohlenwasserstoffe, 7 bis 21 Schweratome. "
          "Ein echtes System ist groesser und je Kombination teurer.")
    # ⚠ DIE KOSTEN JE KOMBINATION SIND KEINE KONSTANTE, und die Streuung gehoert
    #   dazugesagt.  Sie steigt mit der UEBERLEBENSQUOTE: was das Tor passiert, wird
    #   gegen JEDEN bereits behaltenen Konformer per TFD geprueft, also quadratisch.
    #   Ein Molekuel, dessen Kombinationen alle ueberleben, ist damit doppelt teuer --
    #   mehr Kandidaten UND teurere Pruefung je Kandidat.
    _je = sorted((1000.0 * r["t_kombis"] / r["b"], r["name"])
                 for r in zeilen if r["b"] > 0 and r["t_kombis"] > 0.0)
    if len(_je) >= 2:
        print("      Streuung je Kombination: %.1f ms (%s) bis %.1f ms (%s) -- sie "
              "steigt mit der Ueberlebensquote, weil TFD gegen alle Behaltenen prueft."
              % (_je[0][0], _je[0][1], _je[-1][0], _je[-1][1]))
    # ⚠ ZWEI KOSTENPOSTEN MIT VERSCHIEDENEM WACHSTUM.  Die Zustaende je Ring kosten
    #   LINEAR in der Ringzahl (jeder Ring einmal), das Kreuzprodukt EXPONENTIELL.  Steht
    #   der Aufzaehlungsposten heute noch klein da, heisst das nichts fuer 6 Ringe --
    #   der andere ist der, der explodiert.
    _t_zust = sum(r["t_zust"] for r in zeilen)
    print("      Aufteilung: Zustaende je Ring %.1f s (linear in der Ringzahl) · "
          "Kreuzprodukt %.1f s (exponentiell) -- Summe %.1f s ueber %d Proben"
          % (_t_zust, _sum_t, _t_zust + _sum_t, len(zeilen)))
    # ⚠ DIE HOCHRECHNUNG DARF NICHT MIT EINER ERFUNDENEN ZUSTANDSZAHL LAUFEN.  Was hier
    #   gemessen wurde, steht daneben -- und die gemessene Spanne reicht ueber das
    #   hinaus, was die Sweep-Tabelle an UNSUBSTITUIERTEN Ringen findet: ein Ring in
    #   einem Kaefig ist symmetriearm und splittet weiter auf.
    _az = [v for r in zeilen for v in r["zustaende"]]
    if _az:
        print("      Gemessene (a): %d bis %d Zustaende je Ring, Mittel %.1f (Nenner: "
              "%d Ringe ueber %d Proben)"
              % (min(_az), max(_az), sum(_az) / float(len(_az)), len(_az), len(zeilen)))
    print()
    print("      Hochrechnung auf 4,3 Ringe je System (gemessen) -- EIN Kern, EIN System:")
    _stufen = sorted({3, 9, 16} | ({max(_az)} if _az else set()))
    for _z in _stufen:
        _n = _z ** 4.3
        _s = _n * _ms / 1000.0
        # ⚠ DER DECKEL GEHOERT NICHT DIESEM MECHANISMUS ALLEIN.  `deckel_s` ist die
        #   Frist fuer den GANZEN Bau eines Systems; die Ringfaltung ist einer von
        #   vielen Schritten darin.  "Unter dem Deckel" ist deshalb noch kein "geht" --
        #   erst der ANTEIL sagt, ob daneben noch etwas Platz hat.
        print("        %2d Zustaende je Ring -> %10.0f Kombinationen -> %10.0f s "
              "= %6.1f h  = %6.1f %% des Arm-Deckels (%.0f s)%s"
              % (_z, _n, _s, _s / 3600.0, 100.0 * _s / deckel_s, deckel_s,
                 "" if _s <= deckel_s else "   UEBER dem Deckel"))

    # ---- 4 URTEIL.  Zwei Seiten, und sie fallen verschieden aus.
    # ⚠ DAS URTEIL RECHNET MIT DER GEMESSENEN ZUSTANDSZAHL, nicht mit der angenommenen.
    #   Die Annahme, aus der diese Messung hervorging, war "9 bis 16 Zustaende je Ring"
    #   -- eine Zahl vom SWEEP an UNSUBSTITUIERTEN Ringen.  In echten Mehrringern misst
    #   dieser Test 4 bis 26 mit Mittel um 14: die Umgebung bricht die Ringsymmetrie,
    #   und TFD trennt dann mehr.  Mit der angenommenen Zahl zu urteilen, waehrend die
    #   eigene daneben steht, waere die Schoenrechnung in Reinform.
    _zmit = (sum(_az) / float(len(_az))) if _az else 9.0
    _smess = (_zmit ** 4.3) * _ms / 1000.0
    _s3 = (3 ** 4.3) * _ms / 1000.0
    _dmax = max(r["d"] for r in zeilen)
    print()
    print("    ===== 4 URTEIL =====")
    if _smess <= deckel_s:
        print("      AUFZAEHLUNG: BEZAHLBAR -- bei der GEMESSENEN Zustandszahl %.1f je "
              "Ring und 4,3 Ringen %.0f s je System = %.0f %% des Arm-Deckels (%.0f s), "
              "den sich die Ringfaltung mit jedem anderen Bauschritt teilt."
              % (_zmit, _smess, 100.0 * _smess / deckel_s, deckel_s))
    elif _s3 <= deckel_s:
        print("      AUFZAEHLUNG: NICHT BEZAHLBAR bei der GEMESSENEN Aufloesung -- %.1f "
              "Zustaende je Ring, 4,3 Ringe: %.0f s je System = %.0f %% des Arm-Deckels "
              "(%.0f s)." % (_zmit, _smess, 100.0 * _smess / deckel_s, deckel_s))
        print("        Bezahlbar wird es erst weit darunter: bei 3 Zustaenden je Ring "
              "%.0f s = %.0f %% des Deckels.  Der Weg dahin ist WENIGER ZUSTAENDE JE "
              "RING (groebere Entdopplung), nicht eine Kappe auf dem Produkt -- eine "
              "Kappe schneidet nach Faltungstiefe und laesst die tiefen Faltungen weg."
              % (_s3, 100.0 * _s3 / deckel_s))
    else:
        print("      AUFZAEHLUNG: NICHT BEZAHLBAR -- selbst bei 3 Zustaenden je Ring "
              "%.0f s je System gegen einen Deckel von %.0f s." % (_s3, deckel_s))
    print("      ERGEBNIS: groesste gemessene Ausbeute (d) einer einzelnen Probe: %d "
          "Konformere." % _dmax)
    # ⚠ WAS NICHT GEMESSEN WURDE, MUSS IM URTEIL STEHEN.  Eine Probe, die wegen ihrer
    #   Groesse uebersprungen wurde, ist der staerkste Fall gegen die Bezahlbarkeit --
    #   sie stillschweigend aus der Bilanz zu lassen, waere genau die Schoenrechnung,
    #   gegen die dieser Test gebaut ist.
    _uv = [r for r in zeilen if r.get("ungemessen")]
    if _uv:
        print("      ⚠ %d von %d Proben UNGEMESSEN, weil (b) ueber der Testgrenze %d "
              "lag: %s" % (len(_uv), len(zeilen), max_kombis,
                           " · ".join("%s (b=%d)" % (r["name"], r["ungemessen"])
                                      for r in _uv)))
        print("        Das ist selbst ein Befund: bei diesen Systemen ist das volle "
              "Produkt schon zu gross, um es ueberhaupt einmal zu bauen.")
    # ---- DIE DREI FILTER HINTEREINANDER, jeder mit seinem eigenen Nenner.
    # ⚠ Das ist die Kernaussage des ganzen Tests, und sie ist erst als KETTE lesbar:
    #   welcher der drei Filter das Produkt tatsaechlich klein macht, ist eine Messung
    #   und keine Vermutung -- und die Vermutung war, es sei der erste.
    _sum_c = sum(r["c"] for r in zeilen)
    _sum_d = sum(r["d"] for r in zeilen)
    _alle_en = [r for r in zeilen if len(r["energien"]) >= 2]
    _ges = sum(len(r["energien"]) for r in _alle_en)
    _in10 = sum(sum(1 for e in r["energien"] if e - r["energien"][0] <= 10.0)
                for r in _alle_en)
    print()
    print("      DREI FILTER HINTEREINANDER (alle Proben zusammen, jeder mit Nenner):")
    print("        1 PHYSIK  (Kollision + Winkel)  (b)->(c)  %6d von %6d = %5.1f %%"
          % (_sum_c, _sum_b, 100.0 * _sum_c / _sum_b if _sum_b else 0.0))
    print("        2 TFD     (Entdopplung)         (c)->(d)  %6d von %6d = %5.1f %%"
          % (_sum_d, _sum_c, 100.0 * _sum_d / _sum_c if _sum_c else 0.0))
    if _ges:
        print("        3 ENERGIE (<= 10 kcal/mol)      (d)->(e)  %6d von %6d = %5.1f %%"
              % (_in10, _ges, 100.0 * _in10 / _ges))
        # Der schaerfste Filter ist der mit der KLEINSTEN Durchlassquote.
        _kette = (("die PHYSIK", 100.0 * _sum_c / max(1, _sum_b)),
                  ("die TFD-Entdopplung", 100.0 * _sum_d / max(1, _sum_c)),
                  ("die ENERGIE", 100.0 * _in10 / _ges))
        _eng = min(_kette, key=lambda t: t[1])
        print("      ⇒ Der schaerfste Filter ist %s (%.1f %% Durchlass)."
              % (_eng[0], _eng[1]))
    print("=== Kombinatorik: %s ==="
          % ("gemessen" if fehler == 0 else "gemessen, aber %d Pruefung(en) FEHLGESCHLAGEN"
             % fehler))
    return 1 if fehler else 0


# Der Vergleichswert, gegen den Schritt 0 die Byte-Identitaet prueft.  Er stammt aus
# `selbsttest_kombinatorik` Schritt 0 (Decalin, alle Schalter AUS, budget=48) und ist
# damit die Zahl VOR den Aenderungen vom 26.08.  Ihn hier als Konstante zu fuehren, ist
# der Unterschied zwischen "zweimal dasselbe gerechnet" und "gegen den Stand von vorher
# gerechnet": zwei identische Laeufe des NEUEN Codes beweisen gar nichts.
_REF_DECALIN_FRAMES = 22
# Dieselbe Rolle fuer das Kandidatengitter -- `_pucker_space_grid(n, 2, 6)`.
_REF_GITTER = {5: 13, 6: 65, 7: 169, 8: 845}

# ===== DIE PROBEN FUER DIE TRENNSCHAERFE ==========================================
# Nur MEHRRINGER, und zwar aus einem Grund, der die ganze Messung traegt: die Frage
# lautet, ob zwei FALTUNGSVERSCHIEDENE Frames desselben Molekuels von den beiden
# Massen gleich beurteilt werden.  Ein Einringer liefert zu wenige Frames, um eine
# Paarstatistik zu tragen, und vor allem ist bei ihm der Ringanteil an den schweren
# Atomen nahe 1 -- genau der Fall, in dem RMSD und Maximum NICHT auseinanderlaufen.
# Der Effekt, um den es geht, ist ein VERDUENNUNGSeffekt; er braucht Atome, die sich
# nicht bewegen.  Ihn an Cyclohexan zu messen, hiesse ihn wegzudefinieren.
_TRENN_PROBEN = tuple(p for p in _KOMBI_PROBEN if p[0] != "Cyclohexan")

# ===== DIE VERDUENNUNGSREIHE: DIE UNABHAENGIGE VARIABLE DES GANZEN ENTWURFS =========
#
# ⚠ DER ERSTE LAUF HAT DIE EIGENE PROBENWAHL WIDERLEGT.  `_KOMBI_PROBEN` sind reine
#   Ringkohlenwasserstoffe -- gemessener Ringanteil an den schweren Atomen: 100 % bei
#   sechs von acht Proben.  Der Effekt, um den es geht, ist aber ein VERDUENNUNGS-
#   effekt: RMSD teilt die Ringauslenkung durch ALLE Atome, das Maximum durch keines.
#   Bei Ringanteil 1 gibt es nichts zu verduennen, und die Messung sieht folgerichtig
#   nur Faktor 1,6 bis 2,1 statt der erwarteten 4.  Sie hat den Effekt nicht widerlegt,
#   sie hat ihn WEGDEFINIERT -- an Proben, in denen er per Konstruktion nicht auftritt.
#
# DIE REPARATUR ist eine Reihe, in der genau EINE Groesse laeuft: derselbe gefaltete
# Cyclohexanring, an ein immer groesseres STARRES Geruest gehaengt, das ausserdem
# EINGEFROREN wird.  Der Ringanteil faellt von 50 % auf 25 %, die Faltung bleibt
# dieselbe.  Was sich dann zwischen Maximum und RMSD auftut, ist der Effekt.
#   Cyclohexyl + Acen:  6 / (6 + C_Acen) schwere Atome
#   Benzol 50,0 % · Naphthalin 37,5 % · Anthracen 30,0 % · Tetracen 25,0 %
#
# ⚠ WARUM ACENE UND NICHT OLIGOPHENYLE.  Erster Versuch war Cyclohexyl-Oligophenyl bis
#   zum Sexiphenyl (bis 14,3 % Ringanteil).  GESCHEITERT, und zwar messbar: das
#   Maximum wuchs ueber die Reihe von 1,19 auf 2,41 A, obwohl in allen Gliedern
#   DIESELBE Faltung steckt.  Ein Maximum, das mit dem Geruest waechst, misst das
#   Geruest -- die Biaryl-Torsionen sind frei und die Kette klappt beim Relax um.
#   Ein kondensiertes Acen hat diese Freiheitsgrade nicht.
# ⚠ DIE REIHE REICHT NICHT BIS 13,2 %, und das wird nicht mit einem noch groesseren
#   Molekuel erzwungen (Heptacen waere geometrisch brauchbar und chemisch Unsinn).
#   Statt dessen wird an diesen vier Punkten das GESETZ geprueft -- RMSD faellt wie
#   sqrt(Ringanteil), das Maximum bleibt stehen -- und dann auf die 1227 Paare der
#   Kopplungsproben angewandt.  Ein an vier Punkten bestaetigtes Gesetz auf gemessene
#   Paare anzuwenden ist etwas anderes als eine Kurve zu verlaengern.
_VERD_PROBEN = (
    ("Cyclohexylbenzol",     "C1CCCCC1c1ccccc1"),
    ("Cyclohexylnaphthalin", "C1CCCCC1c1ccc2ccccc2c1"),
    ("Cyclohexylanthracen",  "C1CCCCC1c1ccc2cc3ccccc3cc2c1"),
    ("Cyclohexyltetracen",   "C1CCCCC1c1ccc2cc3cc4ccccc4cc3cc2c1"),
)


def _xyz_schwer(txt: str):
    """Schweratomkoordinaten aus einem Frame, wie ihn `generate` zurueckgibt.

    ⚠ AUS DEM AUSGABETEXT, nicht aus einem parallel gehaltenen Conformer.  Was das
      Modul ausliefert, ist dieser Text; jede Metrik, die auf etwas anderem rechnet,
      misst eine Zwischenstufe, die so nie beim Aufrufer ankommt.
    """
    P = []
    for ln in txt.splitlines():
        t = ln.split()
        if len(t) < 4 or t[0] == "H":
            continue
        try:
            P.append([float(t[1]), float(t[2]), float(t[3])])
        except Exception:
            continue
    return _np.array(P, dtype=float)


def _greedy_eintraege(frames, index: int, tol: float) -> int:
    """Wie viele MANIFOLD-EINTRAEGE bleiben, wenn mit ``tol`` entdoppelt wird.

    ``index`` 0 = groesste Auslenkung, 1 = RMSD.  Gierig und in EMISSIONSREIHENFOLGE
    -- genau so entdoppelt `generate`, und genau so entdoppeln die RMSD-Filter des
    Projekts.  Eine optimale Ueberdeckung waere eine andere Zahl und eine andere Frage.
    """
    kept = []
    for P in frames:
        if any(_kabsch_max_rmsd(K, P)[index] < tol for K in kept):
            continue
        kept.append(P)
    return len(kept)


def selbsttest_trennschaerfe(proben=None, tol: float = 0.15,
                             rmsd_projekt: float = 0.30,
                             max_kombis: int = 800) -> int:
    """TRENNT DAS MAXIMUM, WAS DER MITTELWERT VERSCHMILZT?  Mit Nenner.

    DIE FRAGE.  Der Entwurf behauptet: RMSD mittelt eine Ringfaltung weg, die groesste
    Auslenkung nach Kabsch-Ausrichtung tut es nicht.  Die Rechnung dazu steht bei
    `_kabsch_max_rmsd` (0,203 A Ringauslenkung, 13,2 % Ringanteil, 0,086 A RMSD gegen
    0,33 A Maximum -- Faktor 4).  Eine Rechnung ist aber keine Messung: sie unterstellt
    einen Ringanteil und eine Auslenkung, die an echten Mehrringmolekuelen anders
    ausfallen koennen.  Dieser Test rechnet sie an gebauten Frames nach.

    GEMESSEN WIRD AN PAAREN, nicht an Frames.  "Trennschaerfe" ist eine Aussage ueber
    zwei Zustaende, nicht ueber einen; der Nenner ist darum die Zahl der PAARE
    faltungsverschiedener Frames, und der steht ueberall dabei.

    ⚠ ZWEI VERGLEICHE, und nur der erste isoliert das INSTRUMENT:
        (A) gleiche Schwelle, beide 0,15 -- misst allein den Unterschied zwischen
            Maximum und Mittelwert.
        (B) Maximum 0,15 gegen die im Projekt gefuehrte RMSD-Schwelle 0,30 -- misst,
            was heute wirklich passiert, aber vermischt Instrument und Schwelle.
      Wer nur (B) zeigt, kann jeden gewuenschten Effekt durch die Schwellenwahl
      erzeugen.  Wer nur (A) zeigt, redet an der Praxis vorbei.

    ⚠ DIE GEGENRICHTUNG WIRD MITGEMESSEN, obwohl sie null sein MUSS: das Maximum ist
      nie kleiner als das quadratische Mittel derselben Abweichungen.  Ein Paar, das
      das Maximum verschmilzt, verschmilzt der RMSD bei gleicher Schwelle also
      zwingend auch.  Faellt diese Zahl NICHT null aus, ist ein Rechenfehler im Spiel
      und nicht ein Befund -- deshalb steht sie da.
    """
    if not (_RDKIT and _np is not None):
        print("=== Trennschaerfe: RDKit fehlt, uebersprungen ==="); return 0
    import time as _time
    proben = proben or _TRENN_PROBEN
    fehler = 0
    print("=== Selbsttest: Trennschaerfe -- groesste Auslenkung gegen RMSD ===")
    _alt = {k: _os.environ.get(k) for k in
            ("DELFIN_FFFREE_PUCKER_SPACE", "DELFIN_FFFREE_PUCKER_NAMP",
             "DELFIN_FFFREE_PUCKER_NPHASE", "DELFIN_FFFREE_PUCKER_FULL",
             "DELFIN_FFFREE_PUCKER_TRACE", "DELFIN_FFFREE_PUCKER_DEFEKT",
             "DELFIN_FFFREE_PUCKER_XRD", "DELFIN_FFFREE_PUCKER_XRDTOL")}
    zeilen = []
    try:
        # ===== 0 VORGABE AUS -> BYTE-IDENTISCH.  Gegen den Stand VOR dem 26.08. =======
        for _k in _alt:
            _os.environ[_k] = "0"
        _os.environ["DELFIN_FFFREE_PUCKER_XRDTOL"] = "0.15"
        for _n, _soll in sorted(_REF_GITTER.items()):
            _ist = len(_pucker_space_grid(_n, 2, 6))
            if _ist != _soll:
                print("    ✗ 0 GITTER n=%d: %d Kandidaten statt %d" % (_n, _ist, _soll))
                fehler += 1
        if not fehler:
            print("    ✓ 0 GITTER unveraendert: n=5,6,7,8 -> %s"
                  % ", ".join(str(_REF_GITTER[k]) for k in (5, 6, 7, 8)))
        _mv = Chem.AddHs(Chem.MolFromSmiles("C1CCC2CCCCC2C1"))        # Decalin
        if AllChem.EmbedMolecule(_mv, randomSeed=42) != 0:
            print("    ? 0 VORGABE: Decalin nicht einbettbar, NICHT gemessen"); fehler += 1
        else:
            AllChem.MMFFOptimizeMolecule(_mv)
            _aus = generate(_mv, budget=48)
            _aus2 = generate(_mv, budget=48, _zaehler=_neuer_zaehler())
            if len(_aus) != _REF_DECALIN_FRAMES:
                print("    ✗ 0 VORGABE VERAENDERT: Decalin liefert %d Frames statt %d "
                      "-- der Vorgabepfad ist NICHT mehr byte-identisch"
                      % (len(_aus), _REF_DECALIN_FRAMES))
                fehler += 1
            elif _aus != _aus2:
                print("    ✗ 0 MESSMODUS VERAENDERT: %d gegen %d Frames"
                      % (len(_aus), len(_aus2)))
                fehler += 1
            else:
                print("    ✓ 0 VORGABE UNVERAENDERT: Decalin %d Frames (Referenz %d), "
                      "mit und ohne Zaehler zeichengleich"
                      % (len(_aus), _REF_DECALIN_FRAMES))
            # ... und der Beweis, dass die neuen Schalter ueberhaupt REICHWEITE haben.
            # Ein Schalter, der nichts aendert, ist von einem nicht verdrahteten nicht
            # zu unterscheiden -- in diesem Projekt schon fuenfmal an einem Tag passiert.
            _os.environ["DELFIN_FFFREE_PUCKER_XRD"] = "1"
            _os.environ["DELFIN_FFFREE_PUCKER_XRDTOL"] = "0.15"
            _mit_xrd = generate(_mv, budget=48)
            _os.environ["DELFIN_FFFREE_PUCKER_XRD"] = "0"
            _os.environ["DELFIN_FFFREE_PUCKER_DEFEKT"] = "1"
            _mit_def = generate(_mv, budget=48)
            _os.environ["DELFIN_FFFREE_PUCKER_DEFEKT"] = "0"
            print("    ✓ 0 REICHWEITE der Schalter (Decalin, budget=48): Vorgabe %d "
                  "-> XRD %d -> DEFEKT %d Frames"
                  % (len(_aus), len(_mit_xrd), len(_mit_def)))
            if len(_mit_xrd) == len(_aus) and len(_mit_def) == len(_aus):
                print("      ⚠ BEIDE Schalter ohne Wirkung auf DIESER Probe -- das ist "
                      "kein Fehler, aber es beweist an Decalin nichts.  Die Reichweite "
                      "muss dann aus den Proben unten kommen.")

        # ===== 1 DIE PAARSTATISTIK ====================================================
        _os.environ["DELFIN_FFFREE_PUCKER_FULL"] = "1"     # keine Kappe: alle Tiefen
        _os.environ["DELFIN_FFFREE_PUCKER_DEFEKT"] = "1"   # (1) Defektfilter AN
        _os.environ["DELFIN_FFFREE_PUCKER_XRD"] = "0"      # (2) hier NOCH nicht
        print()
        print("    Kandidatenliste wie in der Vorgabe (SPACE=0), Kappe AUS, "
              "(1) Defektfilter AN, (2) noch AUS.")
        print("    %-22s %5s %7s %8s %8s %8s %8s"
              % ("Molekuel", "Frames", "Paare", "max~", "rmsd~", "max/rmsd", "s"))
        for name, smi, klasse in proben:
            try:
                m = Chem.AddHs(Chem.MolFromSmiles(smi))
                if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
                    print("    %-22s Einbettung fehlgeschlagen" % name); continue
                try:
                    AllChem.MMFFOptimizeMolecule(m)
                except Exception:
                    AllChem.UFFOptimizeMolecule(m)
            except Exception as e:
                print("    %-22s Aufbau fehlgeschlagen: %s" % (name, type(e).__name__))
                continue
            # Vorlauf: (b) VOR dem Bau kennen, sonst laeuft eine Probe stundenlang.
            try:
                _rings = [_ring_order(m, set(r)) for r in m.GetRingInfo().AtomRings()
                          if _is_puckerable(m, r)]
                _stv = [_ring_pucker_states(m, r, set(), 0.05) for r in _rings]
            except Exception as e:
                print("    %-22s Zustaende nicht bestimmbar: %s"
                      % (name, type(e).__name__)); continue
            _prod = 1
            for _s in _stv:
                _prod *= len(_s)
            if max_kombis and (_prod - 1) > max_kombis:
                print("    %-22s UEBERSPRUNGEN: (b) = %d ueber der Testgrenze %d"
                      % (name, _prod - 1, max_kombis))
                continue
            _z = _neuer_zaehler()
            _t0 = _time.perf_counter()
            try:
                _out = generate(m, budget=10 ** 9, _zaehler=_z)
            except Exception as e:
                print("    %-22s generate() geplatzt: %s" % (name, type(e).__name__))
                continue
            _dt = _time.perf_counter() - _t0
            # Der GRUNDZUSTAND ist selbst ein Manifold-Eintrag und gehoert in die
            # Paarmenge: eine Faltung, die vom Ausgangsframe nicht zu unterscheiden ist,
            # ist genauso ein Doppel wie zwei ununterscheidbare Faltungen.
            _frames = [_xyz_schwer(_conf_to_xyz(m))] + [_xyz_schwer(x) for x, _l in _out]
            _frames = [P for P in _frames if P.size and P.shape == _frames[0].shape]
            # Ringanteil an den schweren Atomen -- die Verduennung, um die es geht.
            try:
                _ring_at = set()
                for _r in m.GetRingInfo().AtomRings():
                    _ring_at |= {int(x) for x in _r}
                _hv = [i for i in range(m.GetNumAtoms())
                       if m.GetAtomWithIdx(i).GetSymbol() != "H"]
                _ring_pos = [k for k, i in enumerate(_hv) if i in _ring_at]
                _anteil = len(_ring_pos) / float(len(_hv)) if _hv else 0.0
            except Exception:
                _ring_pos, _anteil = [], 0.0
            _paare = []
            for _i in range(len(_frames)):
                for _j in range(_i + 1, len(_frames)):
                    _mx, _rm = _kabsch_max_rmsd(_frames[_i], _frames[_j])
                    _paare.append((_mx, _rm))
            if not _paare:
                print("    %-22s %5d %7d   -- kein Paar, keine Aussage"
                      % (name, len(_frames), 0))
                continue
            _mxs = sorted(p[0] for p in _paare)
            _rms = sorted(p[1] for p in _paare)
            _med = lambda v: v[len(v) // 2]
            print("    %-22s %5d %7d %8.3f %8.3f %8.1f %8.1f"
                  % (name, len(_frames), len(_paare), _med(_mxs), _med(_rms),
                     (_med(_mxs) / _med(_rms)) if _med(_rms) > 1e-9 else float("inf"),
                     _dt))
            print("        %-26s Ringanteil an den schweren Atomen %4.1f %% (%d von %d)"
                  " · verworfen: Kollision %d · Winkel %d · Bindung %d · TFD %d"
                  % (klasse, 100.0 * _anteil, len(_ring_pos),
                     len(_frames[0]) if _frames else 0,
                     int(_z.get("kollision", 0)), int(_z.get("winkel", 0)),
                     int(_z.get("bindung", 0)), int(_z.get("tfd_doppelt", 0))))
            # ---- LAUF B: DIESELBE PROBE MIT (2) AN.  Das ist die Zahl fuer (c), und
            #      sie wird GEBAUT und nicht aus Lauf A hochgerechnet: mit (2) an
            #      entdoppelt schon `_ring_pucker_states` je Ring, das Kreuzprodukt (b)
            #      ist also ein anderes.  Wer das aus den Frames von Lauf A greedy
            #      nachbildet, misst den Nachbau -- der Fehler, den dieses Projekt
            #      schon mehrfach als Befund durchgehen liess.
            _os.environ["DELFIN_FFFREE_PUCKER_XRD"] = "1"
            _zb = _neuer_zaehler()
            _t1 = _time.perf_counter()
            try:
                _outb = generate(m, budget=10 ** 9, _zaehler=_zb)
            except Exception:
                _outb, _zb = [], _neuer_zaehler()
            _dtb = _time.perf_counter() - _t1
            _os.environ["DELFIN_FFFREE_PUCKER_XRD"] = "0"
            zeilen.append({"name": name, "klasse": klasse, "paare": _paare,
                           "frames": _frames, "anteil": _anteil, "dt": _dt,
                           "b": int(_z.get("aufzaehlung", 0)),
                           "c": int(_z.get("tor_ueberlebt", 0)),
                           "d": len(_out),
                           "zustaende": list(_z.get("zustaende_je_ring") or []),
                           "bindung": int(_z.get("bindung", 0)),
                           "kollision": int(_z.get("kollision", 0)),
                           "winkel": int(_z.get("winkel", 0)),
                           "b2": int(_zb.get("aufzaehlung", 0)),
                           "c2": int(_zb.get("tor_ueberlebt", 0)),
                           "d2": len(_outb), "dt2": _dtb,
                           "xrd2": int(_zb.get("xrd_doppelt", 0)),
                           "zustaende2": list(_zb.get("zustaende_je_ring") or []),
                           "mol": m, "smi": smi})
    finally:
        for k, v in _alt.items():
            if v is None:
                _os.environ.pop(k, None)
            else:
                _os.environ[k] = v

    # ===== 1b DIE VERDUENNUNGSREIHE ================================================
    # Getrennt gelaufen, getrennt berichtet -- sie beantwortet eine ANDERE Frage als
    # die Kopplungsproben oben (dort: welche Paare trennt wer; hier: WOVON der
    # Unterschied ueberhaupt abhaengt).  Zusammengeworfen waeren beide unlesbar.
    verd, _verd_tfd = [], []
    _altv = {k: _os.environ.get(k) for k in
             ("DELFIN_FFFREE_PUCKER_FULL", "DELFIN_FFFREE_PUCKER_DEFEKT",
              "DELFIN_FFFREE_PUCKER_XRD", "DELFIN_FFFREE_PUCKER_SPACE",
              "DELFIN_FFFREE_PUCKER_NAMP", "DELFIN_FFFREE_PUCKER_NPHASE")}
    try:
        # ⚠ RAUMGITTER AN, damit der EINE Ring genug Zustaende liefert -- eine
        #   Paarstatistik aus drei Paaren waere keine.  Es ist derselbe Ring in allen
        #   sechs Molekuelen, also aendert das an der unabhaengigen Variablen nichts.
        _os.environ["DELFIN_FFFREE_PUCKER_SPACE"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_NAMP"] = "2"
        _os.environ["DELFIN_FFFREE_PUCKER_NPHASE"] = "8"
        _os.environ["DELFIN_FFFREE_PUCKER_FULL"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_DEFEKT"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_XRD"] = "0"
        print()
        print("    ===== 1b VERDUENNUNGSREIHE: derselbe Ring, wachsendes Geruest =====")
        print("      Geruest EINGEFROREN (`frozen`) -- s. den Kommentar bei "
              "_VERD_PROBEN: ohne das misst die Reihe die Geruestbewegung mit.")
        print("      %-24s %6s %7s %7s %7s %6s %8s %8s %9s"
              % ("Molekuel", "schwer", "Anteil", "(a)0,05", "(a)0,005", "Paare",
                 "max~", "rmsd~", "max/rmsd"))
        for name, smi in _VERD_PROBEN:
            try:
                m = Chem.AddHs(Chem.MolFromSmiles(smi))
                if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
                    print("      %-24s Einbettung fehlgeschlagen" % name); continue
                AllChem.MMFFOptimizeMolecule(m)
            except Exception as e:
                print("      %-24s Aufbau ausgefallen: %s"
                      % (name, type(e).__name__)); continue
            try:
                _ring_at, _ring_ord = set(), []
                for _r in m.GetRingInfo().AtomRings():
                    if _is_puckerable(m, _r):
                        _ring_at |= {int(x) for x in _r}
                        if not _ring_ord:
                            _ring_ord = _ring_order(m, set(_r))
                _hv = [i for i in range(m.GetNumAtoms())
                       if m.GetAtomWithIdx(i).GetSymbol() != "H"]
                _rp = [k for k, i in enumerate(_hv) if i in _ring_at]
                # frei = der Ring und SEINE Wasserstoffe; alles andere steht fest.
                _frei = set(_ring_at)
                for _i in list(_ring_at):
                    for _nb in m.GetAtomWithIdx(int(_i)).GetNeighbors():
                        if _nb.GetSymbol() == "H":
                            _frei.add(int(_nb.GetIdx()))
                _fr = set(range(m.GetNumAtoms())) - _frei
            except Exception as e:
                print("      %-24s Ringmenge unbestimmbar: %s"
                      % (name, type(e).__name__)); continue
            if not _ring_at or not _hv:
                print("      %-24s kein faltbarer Ring" % name); continue
            # ===== TFD MITTELT AUCH -- UND ES IST DAS INSTRUMENT IM EINSATZ ==========
            #
            # Diese Reihe hat es beim Bauen selbst aufgedeckt: Anthracen und Tetracen
            # lieferten NULL Faltungen, und zwar nicht am Realismustor, sondern schon
            # bei den Zustaenden JE RING -- (b) war 0, das Kreuzprodukt also 1x nichts.
            # Der Ring ist derselbe wie im Cyclohexylbenzol, das 55 Paare liefert.
            #
            # Die Ursache ist dieselbe Krankheit eine Ebene hoeher: TFD vergleicht ALLE
            # Torsionen des Molekuels und MITTELT ueber sie.  Ein grosses starres
            # Geruest bringt viele Torsionen mit, die sich nicht aendern -- der Beitrag
            # der sechs Ringtorsionen wird durch sie geteilt und faellt unter die
            # Schwelle 0,05.  Die Faltung verschwindet im Mittel, genau wie beim RMSD.
            #
            # ⚠ DAS IST KEIN NEBENBEFUND.  TFD ist das Entdopplungsmass, das HEUTE im
            #   Bau laeuft.  Wenn es mit der Ligandgroesse unschaerfer wird, dann
            #   verliert der Manifold Faltungen genau bei den Systemen, um die es geht
            #   -- grosse Liganden, kleiner Ringanteil.
            # Gemessen wird das mit dem einzigen Mittel, das die beiden Ursachen trennt:
            # dieselben Kandidaten, zwei Schwellen.  Steigt (a) bei 0,005 stark an, war
            # es die Schwelle (also die Verduennung); bleibt es gleich, sind die
            # Zustaende wirklich nicht da.
            try:
                _a05 = len(_ring_pucker_states(m, _ring_ord, _fr, 0.05))
                _a005 = len(_ring_pucker_states(m, _ring_ord, _fr, 0.005))
            except Exception:
                _a05 = _a005 = -1
            _zv = _neuer_zaehler()
            try:
                _out = generate(m, frozen=_fr, budget=10 ** 9, _zaehler=_zv)
            except Exception as e:
                print("      %-24s generate() geplatzt: %s"
                      % (name, type(e).__name__)); continue
            _frames = [_xyz_schwer(_conf_to_xyz(m))] + [_xyz_schwer(x) for x, _l in _out]
            _frames = [P for P in _frames if P.size and P.shape == _frames[0].shape]
            if len(_frames) < 2:
                # ⚠ EINE NULL BEKOMMT IHREN GRUND.  "erzeugt und verworfen" sieht von
                #   aussen genauso aus wie "nie gebaut" -- der Fehlschluss, der am
                #   14.08. den Feuerzensus wertlos gemacht hat.
                print("      %-24s %6d %6.1f%% %7d %7d      -- nur %d Frame(e), (b)=%d"
                      % (name, len(_hv),
                         100.0 * len(_rp) / float(len(_hv)) if _hv else 0.0,
                         _a05, _a005, len(_frames), int(_zv.get("aufzaehlung", 0))))
                _verd_tfd.append({"name": name, "n": len(_hv),
                                  "anteil": len(_rp) / float(len(_hv)) if _hv else 0.0,
                                  "a05": _a05, "a005": _a005})
                continue
            _anteil = len(_rp) / float(len(_hv))
            _mx, _rm, _pa = [], [], []
            _geruest = 0.0
            _nicht_ring = [k for k in range(len(_hv)) if k not in set(_rp)]
            for _i in range(len(_frames)):
                for _j in range(_i + 1, len(_frames)):
                    _v = _kabsch_abweichungen(_frames[_i], _frames[_j])
                    _a = float(_v.max())
                    _b = float(_np.sqrt(float((_v ** 2).mean())))
                    _mx.append(_a); _rm.append(_b); _pa.append((_a, _b))
                    # ⚠ EINFROSTPROBE OHNE KABSCH, und das ist der Punkt.  Beide Frames
                    #   stehen im SELBEN Bezugssystem -- es wurde nichts neu eingebettet,
                    #   nur relaxiert.  "Steht das Geruest?" ist damit eine Frage an die
                    #   ROHEN Koordinaten.  Nach Kabsch waere sie unbeantwortbar: die
                    #   Ausrichtung minimiert das Gesamt-RMSD und verteilt den Fehler auf
                    #   ALLE Atome, also auch auf festgehaltene -- eine erste Fassung hat
                    #   genau daraus 1,4 A Geruestbewegung gemeldet, die es nicht gab.
                    if _nicht_ring:
                        _roh = _np.linalg.norm(_frames[_i][_nicht_ring]
                                               - _frames[_j][_nicht_ring], axis=1)
                        _geruest = max(_geruest, float(_roh.max()))
            _mx.sort(); _rm.sort()
            _m1, _r1 = _mx[len(_mx) // 2], _rm[len(_rm) // 2]
            print("      %-24s %6d %6.1f%% %7d %7d %6d %8.3f %8.3f %9.1f"
                  % (name, len(_hv), 100.0 * _anteil, _a05, _a005, len(_pa), _m1, _r1,
                     (_m1 / _r1) if _r1 > 1e-9 else float("inf")))
            _verd_tfd.append({"name": name, "n": len(_hv), "anteil": _anteil,
                              "a05": _a05, "a005": _a005})
            verd.append({"name": name, "n": len(_hv), "anteil": _anteil,
                         "max": _m1, "rmsd": _r1, "paare": _pa, "geruest": _geruest})
    finally:
        for k, v in _altv.items():
            if v is None:
                _os.environ.pop(k, None)
            else:
                _os.environ[k] = v
    if verd:
        # ---- PRUEFUNG: HAT DAS EINFRIEREN GEHALTEN?  Roh, ohne Ausrichtung.
        _gmax = max(r["geruest"] for r in verd)
        if _gmax < 1e-6:
            print("      ✓ EINFROSTPROBE (rohe Koordinaten): groesste Geruestauslenkung "
                  "%.2e A -- das Geruest steht, die Reihe isoliert die Verduennung."
                  % _gmax)
        else:
            print("      ✗ EINFROSTPROBE: Geruest bewegt sich um bis zu %.3f A -- die "
                  "Reihe misst NICHT nur Verduennung." % _gmax)
            fehler += 1
    if _verd_tfd:
        # ---- DAS EIGENTLICHE ERGEBNIS DIESER REIHE, und es war nicht das gesuchte.
        print()
        print("      ===== TFD MITTELT GENAUSO -- und TFD laeuft heute im Bau =====")
        for r in _verd_tfd:
            print("        %-24s %2d schwere Atome, Ringanteil %5.1f %% -> (a) %d bei "
                  "Schwelle 0,05 · %d bei 0,005"
                  % (r["name"], r["n"], 100.0 * r["anteil"], r["a05"], r["a005"]))
        _v0, _vn = _verd_tfd[0], _verd_tfd[-1]
        if _vn["a05"] < _v0["a05"] and _vn["a005"] > _vn["a05"]:
            print("      ⇒ BEFUND: mit wachsendem Geruest faellt (a) bei Schwelle 0,05 "
                  "von %d auf %d -- DERSELBE Ring, dieselben Kandidaten.  Bei 0,005 "
                  "kommen die Zustaende zurueck (%d).  Die Faltungen sind also DA und "
                  "werden von TFD verschmolzen, nicht vom Generator ausgelassen."
                  % (_v0["a05"], _vn["a05"], _vn["a005"]))
            print("        URSACHE: TFD vergleicht ALLE Torsionen und mittelt ueber sie. "
                  "Ein grosses starres Geruest bringt unbewegte Torsionen mit; der "
                  "Beitrag der Ringtorsionen wird durch sie geteilt.  Das ist Wort fuer "
                  "Wort der RMSD-Fehler, eine Ebene hoeher -- und TFD ist das Mass, das "
                  "im Bau ENTSCHEIDET.")
            print("        ⚠ TRAGWEITE: der Effekt waechst mit der LIGANDGROESSE.  Er "
                  "trifft also am haertesten die Systeme, um die es geht -- grosse "
                  "Liganden, kleiner Ringanteil.  Die 59,7 % Ringidentitaeten mit nur "
                  "EINER Faltung (16.08.) haben hier eine kandidatenfaehige Ursache, "
                  "die nichts mit dem Generator zu tun hat.")
        else:
            print("      ⇒ Kein Verduennungsmuster in (a) -- TFD verschmilzt hier nicht.")
        _vp = [p for r in verd for p in r["paare"]]
        _vN = len(_vp)
        if _vN:
            _vA = sum(1 for a, b in _vp if b < tol <= a)
            _vB = sum(1 for a, b in _vp if b < rmsd_projekt <= a)
            print("      ⇒ Trennschaerfe auf dieser Reihe (Nenner %d Paare, Ringanteil "
                  "%.0f bis %.0f %%): (A) %d = %.1f %% · (B) %d = %.1f %%"
                  % (_vN, 100.0 * min(r["anteil"] for r in verd),
                     100.0 * max(r["anteil"] for r in verd),
                     _vA, 100.0 * _vA / _vN, _vB, 100.0 * _vB / _vN))

    if not zeilen:
        print("=== Trennschaerfe: nichts gemessen ==="); return 1

    _alle = [p for r in zeilen for p in r["paare"]]
    _N = len(_alle)
    _A = sum(1 for mx, rm in _alle if rm < tol <= mx)
    _B = sum(1 for mx, rm in _alle if rm < rmsd_projekt <= mx)
    _G = sum(1 for mx, rm in _alle if mx < tol <= rm)
    print()
    print("    ===== 1 TRENNSCHAERFE (Nenner: %d Paare faltungsverschiedener Frames "
          "ueber %d Molekuele) =====" % (_N, len(zeilen)))
    print("      (A) gleiche Schwelle %.2f A -- nur das INSTRUMENT:" % tol)
    print("          Maximum trennt, RMSD verschmilzt:  %6d von %6d = %5.1f %%"
          % (_A, _N, 100.0 * _A / _N if _N else 0.0))
    print("      (B) Maximum %.2f gegen die Projekt-RMSD-Schwelle %.2f -- die PRAXIS:"
          % (tol, rmsd_projekt))
    print("          Maximum trennt, RMSD verschmilzt:  %6d von %6d = %5.1f %%"
          % (_B, _N, 100.0 * _B / _N if _N else 0.0))
    print("      GEGENRICHTUNG (muss 0 sein, das Maximum ist nie kleiner als das "
          "quadratische Mittel): %d" % _G)
    # ---- DIE ZAHL, DIE FUER ECHTE SYSTEME GILT.
    # ⚠ DIE PROBEN OBEN HABEN 86 BIS 100 %% RINGANTEIL -- der Verduennungseffekt kommt
    #   darin per Konstruktion nicht vor.  Was sie liefern, ist die UNTERGRENZE der
    #   Trennschaerfe.  Ein reales DELFIN-System hat einen Metallkern, aromatische
    #   Rueckgrate und Substituenten; die Vorgabe nennt 13,2 % Ringanteil.
    #   Angewandt wird das oben an vier Punkten BESTAETIGTE Gesetz: das Maximum bleibt,
    #   der RMSD faellt wie sqrt(Anteil).  Kein neues Molekuel, keine Kurve verlaengert
    #   -- dieselben gemessenen Paare, mit dem Nenner eines realen Systems.
    # ---- DER RINGANTEIL IST DIE UNABHAENGIGE VARIABLE, und beide Punkte sind GEMESSEN.
    # ⚠ HIER STAND EINMAL EINE HOCHRECHNUNG AUF 13,2 % -- zweimal, und beide Male
    #   falsch.  (i) "RMSD faellt wie sqrt(Ringanteil)" gilt nur bei IDENTISCHER
    #   Ausrichtung; Kabsch richtet aber aus und verteilt den Fehler um.  (ii) Das Paar
    #   mit starren Kopien aufzufuellen und neu zu ueberlagern hat den Fehler nur
    #   verschoben: das Geruest lag auf dem Molekuel und band die Ausrichtung so hart,
    #   dass das Maximum von 0,94 auf 3,95 A stieg -- gemessen wurde der Wechsel des
    #   AUSRICHTUNGSREGIMES, nicht die Verduennung.  Beide Versuche sind entfernt.
    #   Was bleibt, sind zwei GEMESSENE Punkte an echten Molekuelen; die Reihe 1b
    #   reicht nicht bis 13,2 %, und der Grund dafuer ist selbst der Befund (TFD
    #   liefert dort keine Frames mehr).
    _vp2 = [p for r in verd for p in r["paare"]]
    if _vp2:
        _vA2 = sum(1 for a, b in _vp2 if b < tol <= a)
        _vB2 = sum(1 for a, b in _vp2 if b < rmsd_projekt <= a)
        print("      ⇒ GEGEN DEN RINGANTEIL, beide Punkte gemessen:")
        print("        ~100 %% Ringanteil (Kopplungsproben, %4d Paare): (A) %5.1f %% · "
              "(B) %5.1f %%" % (_N, 100.0 * _A / _N, 100.0 * _B / _N))
        print("        38-50 %% Ringanteil (Verduennungsreihe, %4d Paare): (A) %5.1f %% "
              "· (B) %5.1f %%" % (len(_vp2), 100.0 * _vA2 / len(_vp2),
                                  100.0 * _vB2 / len(_vp2)))
        print("        ⇒ die Trennschaerfe WAECHST mit der Verduennung (A: %.1f -> "
              "%.1f %%).  Die Kopplungsproben sind damit die UNTERGRENZE, nicht die "
              "Antwort." % (100.0 * _A / _N, 100.0 * _vA2 / len(_vp2)))
    if _G:
        print("      ✗ GEGENRICHTUNG NICHT NULL -- Rechenfehler, kein Befund."); fehler += 1

    # ---- 2 WAS DAS FUER DIE ZAHL DER MANIFOLD-EINTRAEGE HEISST.
    # ⚠ Eine Prozentzahl ueber Paare sagt noch nicht, wie viele EINTRAEGE entstehen:
    #   Entdopplung ist gierig und transitiv-unsauber, drei paarweise knappe Frames
    #   koennen zu einem oder zu zweien werden.  Also nachzaehlen statt hochrechnen.
    print()
    print("    ===== 2 EINTRAEGE JE MOLEKUEL, je Kriterium =====")
    print("      %-22s %7s %8s %8s %8s %8s"
          % ("Molekuel", "roh", "max%.2f" % tol, "rmsd%.2f" % tol,
             "rmsd%.2f" % rmsd_projekt, "rmsd0.50"))
    _sum = {"roh": 0, "max": 0, "r_gleich": 0, "r_proj": 0, "r_50": 0}
    for r in zeilen:
        _f = r["frames"]
        _e = (len(_f), _greedy_eintraege(_f, 0, tol), _greedy_eintraege(_f, 1, tol),
              _greedy_eintraege(_f, 1, rmsd_projekt), _greedy_eintraege(_f, 1, 0.50))
        print("      %-22s %7d %8d %8d %8d %8d" % ((r["name"],) + _e))
        for _k, _v in zip(("roh", "max", "r_gleich", "r_proj", "r_50"), _e):
            _sum[_k] += _v
    print("      %-22s %7d %8d %8d %8d %8d"
          % ("SUMME", _sum["roh"], _sum["max"], _sum["r_gleich"], _sum["r_proj"],
             _sum["r_50"]))
    if _sum["max"] > 0:
        print("      ⇒ RMSD bei %.2f A behaelt %d von %d Eintraegen, die das Maximum "
              "bei derselben Schwelle als UNTERSCHEIDBAR fuehrt (%.0f %%).  Die "
              "Differenz %d sind Faltungen, die eine RMSD-Entdopplung LOESCHT."
              % (tol, _sum["r_gleich"], _sum["max"],
                 100.0 * _sum["r_gleich"] / _sum["max"],
                 _sum["max"] - _sum["r_gleich"]))
        print("      ⇒ bei der Projektschwelle %.2f A: %d von %d (%.0f %%), "
              "geloescht %d." % (rmsd_projekt, _sum["r_proj"], _sum["max"],
                                 100.0 * _sum["r_proj"] / _sum["max"],
                                 _sum["max"] - _sum["r_proj"]))
    # ---- 3 WAS (1)+(2) DIE KOMBINATORIK KOSTEN -- ODER SPAREN.
    # ⚠ DIE ENTSCHEIDENDE UNTERSCHEIDUNG, und sie ist leicht zu verfehlen: (1) und (2)
    #   koennen an ZWEI Stellen wirken, und nur eine davon spart Rechenzeit.
    #     JE RING  (`_ring_pucker_states`)  -> senkt (a), also (b) POTENZIERT: der
    #                                          einzige Ort, an dem etwas billiger wird.
    #     JE KOMBINATION (`generate`)       -> senkt nur die Zahl der EINTRAEGE.  Der
    #                                          Relax ist da schon bezahlt; das Tor
    #                                          waehlt aus, es spart nichts.
    #   Ein Filter, der nur unten wirkt, macht die vollstaendige Faltung NICHT
    #   bezahlbar, egal wie scharf er ist.  Darum stehen (a) und (b) hier nebeneinander.
    print()
    print("    ===== 3 KOMBINATORIK MIT (1)+(2) -- (a) je Ring und (b) das Produkt =====")
    print("      %-22s %-12s %-12s %8s %8s %7s %7s"
          % ("Molekuel", "(a) ohne (2)", "(a) mit (2)", "(b) ohne", "(b) mit",
             "(d) ohne", "(d) mit"))
    _sb = _sb2 = _sd = _sd2 = 0
    for r in zeilen:
        print("      %-22s %-12s %-12s %8d %8d %7d %7d"
              % (r["name"],
                 "x".join(str(v) for v in r["zustaende"]) or "-",
                 "x".join(str(v) for v in r["zustaende2"]) or "-",
                 r["b"], r["b2"], r["d"], r["d2"]))
        _sb += r["b"]; _sb2 += r["b2"]; _sd += r["d"]; _sd2 += r["d2"]
        # ⚠ EINE NULL IN (d) MUSS IHREN GRUND NENNEN, sonst liest sie sich wie ein
        #   Defekt.  Bei Norbornan ist sie das GEWOLLTE Ergebnis: seine groesste
        #   Faltungsauslenkung liegt bei 0,139 A, also UNTER der Aufloesungsschwelle.
        #   Ein starr verbrueckter Bicyclus HAT keine zweite Faltung -- (2) sagt genau
        #   das, und der Manifold behaelt den Grundzustand (der nie durch dieses Tor
        #   geht).  Aus Sicht des Kristallographen ist ein Eintrag richtig, nicht drei.
        if r["d"] > 0 and r["d2"] == 0:
            print("        %-20s (d) faellt auf 0: alle %d Faltungen liegen unter %.2f A "
                  "Maximalauslenkung -- ununterscheidbar vom Grundzustand, EIN Eintrag."
                  % (r["name"], r["d"], tol))
    print("      %-22s %-12s %-12s %8d %8d %7d %7d"
          % ("SUMME", "", "", _sb, _sb2, _sd, _sd2))
    if _sb:
        print("      ⇒ (b), die AUFZAEHLUNG und damit der PREIS: %d -> %d = %+.1f %%."
              % (_sb, _sb2, 100.0 * (_sb2 - _sb) / _sb))
        print("      ⇒ (d), die EINTRAEGE und damit das ERGEBNIS: %d -> %d = %+.1f %%."
              % (_sd, _sd2, 100.0 * (_sd2 - _sd) / max(1, _sd)))
        # ⚠ PARTITION, KEIN MITTELWERT.  Ein Gesamtprozentsatz ueber acht Proben kann
        #   nicht sagen, ob (2) ueberall ein bisschen spart oder bei zwei Proben viel
        #   und bei sechs gar nichts -- und das sind voellig verschiedene Mechanismen.
        #   Der zweite Fall waere KEIN allgemeiner Kostenhebel, sondern ein Befund
        #   ueber eine Klasse.  Geteilt wird nach der Zahl der Proben mit Wirkung.
        _wirkt = [r for r in zeilen if r["b2"] < r["b"]]
        _still = [r for r in zeilen if r["b2"] >= r["b"]]
        print("      ⇒ PARTITION: (2) senkt (b) bei %d von %d Proben (%s); bei den "
              "anderen %d aendert sie (a) um keinen einzigen Zustand (%s)."
              % (len(_wirkt), len(zeilen),
                 ", ".join(r["name"] for r in _wirkt) or "keiner", len(_still),
                 ", ".join(r["name"] for r in _still) or "keine"))
        if len(_wirkt) <= len(zeilen) // 2:
            print("      ⇒ URTEIL: (2) ist KEIN allgemeiner Kostenhebel.  Sie greift "
                  "dort, wo Ringzustaende ohnehin fast entartet sind (die VERBRUECKTEN "
                  "Proben -- ein verbrueckter Ring KANN kaum falten), und nirgends "
                  "sonst.  Das Maximum kennt keinen Nenner: genau die Eigenschaft, die "
                  "es gegen Verduennung unempfindlich macht, hindert es daran, echte "
                  "Ringmulden zusammenzuziehen.  Es macht die Auswahl RICHTIG, nicht "
                  "BILLIG.")
        else:
            print("      ⇒ URTEIL: (2) senkt (a) und damit (b) bei der MEHRHEIT der "
                  "Proben -- ein echter Kostenhebel, nicht nur eine Korrektur.")
        _sbi = sum(r["bindung"] for r in zeilen)
        _swi = sum(r["winkel"] for r in zeilen)
        _skl = sum(r["kollision"] for r in zeilen)
        print("      ⇒ (1) DEFEKTFILTER, aufgeschluesselt (Nenner %d Kombinationen): "
              "Kollision %d (%.1f %%) · Winkel %d (%.1f %%) · Bindung %d (%.1f %%)."
              % (_sb, _skl, 100.0 * _skl / _sb, _swi, 100.0 * _swi / _sb,
                 _sbi, 100.0 * _sbi / _sb))
        if _skl == 0:
            print("        ⚠ Das Kollisionstor feuert NULL mal -- dieselbe Nullreichweite "
                  "wie am 26.08. (0 von 17 754).  Ein zweites Mal gemessen, ein zweites "
                  "Mal null: der Name `Kollisionstor` beschreibt keinen wirksamen Filter.")
        if _sbi:
            print("        ✓ Das NEUE Bindungstor feuert %d mal -- es ist verdrahtet und "
                  "hat Reichweite; es sieht genau den Bruch, den das Selbstgate per "
                  "Konstruktion fuer 'nicht gebunden' haelt." % _sbi)
    # ---- 4 BRAUCHT ES STUFE (3), DIE ENERGIE?
    # ⚠ DIE FRAGE IST NICHT "waere Energie schoen", sondern "loest sie das Problem, das
    #   (1) und (2) offen lassen".  Und das Problem ist der PREIS (b), nicht die Zahl
    #   der Eintraege (d).  Eine Energie wird -- wie jedes andere Tor hier -- NACH dem
    #   Relax ausgewertet; sie kann (b) also gar nicht senken.  Ein Mechanismus, der
    #   den Engpass per Konstruktion nicht erreicht, wird nicht gebaut, sondern benannt.
    print()
    print("    ===== 4 BRAUCHT ES DIE ENERGIE? =====")
    print("      Der Engpass ist (b) = %d Kombinationen, jede mit einem Relax BEVOR "
          "irgendein Tor sie sieht." % _sb)
    print("      Eine Energieauswahl wird an derselben Stelle ausgewertet wie (1) und "
          "(2) -- nach dem Relax.  Sie kann (b) also per Konstruktion nicht senken.")
    print("      ⇒ ENERGIE NICHT GEBAUT.  Sie wuerde das Ergebnis weiter ausduennen "
          "(%d Eintraege) und den Preis unveraendert lassen.  Der einzige Ort, an dem "
          "etwas zu sparen ist, ist (a) -- die Zustaende JE RING, vor dem Kreuzprodukt."
          % _sd2)
    print("=== Trennschaerfe: %s ==="
          % ("gemessen" if fehler == 0 else "%d Pruefung(en) FEHLGESCHLAGEN" % fehler))
    return 1 if fehler else 0


# Die unsubstituierten Kalibrierringe.  Sie sind der EINZIGE Ort, an dem sich die
# Bedeutung der Schwelle pruefen laesst: dort hat das Molekuel kein Geruest, das
# verduennen koennte, und beide Masse muessen deshalb DASSELBE sagen.  Faellt das aus,
# ist die ringlokale Fassung nicht "anders geeicht", sondern ein anderes Instrument.
_LOKAL_KALIBER = (("Cyclopentan", "C1CCCC1"), ("Cyclohexan", "C1CCCCC1"),
                  ("Cycloheptan", "C1CCCCCC1"), ("Cyclooctan", "C1CCCCCCC1"))


def _lokal_ringlage(mol):
    """(erster faltbarer Ring in Ringreihenfolge, Ringatome, eingefrorenes Geruest).

    Dieselbe Vorbereitung wie in der Verduennungsreihe von `selbsttest_trennschaerfe`
    -- ⚠ und das ist der Zweck: die Reparatur muss an DERSELBEN Messung geprueft
    werden, die den Fehler gezeigt hat.  Eine zweite, leicht andere Vorbereitung
    vergliche zwei Messungen statt zweier Masse.
    """
    ring_at, ring_ord = set(), []
    for r in mol.GetRingInfo().AtomRings():
        if _is_puckerable(mol, r):
            ring_at |= {int(x) for x in r}
            if not ring_ord:
                ring_ord = _ring_order(mol, set(r))
    if not ring_ord:
        return None, None, None
    frei = set(ring_at)
    for i in list(ring_at):
        for nb in mol.GetAtomWithIdx(int(i)).GetNeighbors():
            if nb.GetSymbol() == "H":
                frei.add(int(nb.GetIdx()))
    return ring_ord, ring_at, set(range(mol.GetNumAtoms())) - frei


def selbsttest_tfd_lokal() -> int:
    """FAEHRT DIE VERDUENNUNGSTABELLE NACH -- global gegen ringlokal, dieselben Proben.

    VORGESCHICHTE.  `selbsttest_trennschaerfe` hat beim Bauen einen Befund abgeworfen,
    den es gar nicht gesucht hatte: derselbe Cyclohexanring liefert an wachsendem
    starrem Acen immer weniger Zustaende (11 -> 6 -> 1 -> 1 bei TFD 0,05), und bei
    0,005 kommen sie zurueck (62 -> 49 -> 23 -> 10).  Die Faltungen sind also DA und
    werden vom Entdopplungsmass verschmolzen.

    DIESER TEST BEANTWORTET VIER FRAGEN, und zwar in dieser Reihenfolge, weil jede
    naechste sinnlos waere, wenn die davor ausfaellt:
        0  Ist der Vorgabepfad unveraendert?          (Gitter, Decalin)
        1  WORAN liegt die Verduennung genau?         (RDKits Gewichte, gerechnet)
        2  Ueberlebt die Symmetriefaltung?            (unsubstituierter Ring)
        3  Verschwindet der Gradient?                 (die vier Acene)

    ⚠ SCHRITT 2 IST DAS ABBRUCHKRITERIUM, nicht Schritt 3.  Am 26.08. ist schon ein
      Ersatz fuer TFD daran gestorben, dass er die Molekuelsymmetrie nicht mitfaltete
      -- n=5 ging von 3,3,3 auf 9,13,14.  Eine Fassung, die den Gradienten beseitigt
      und dabei den Fuenfring aufsplittet, ist KEINE Reparatur, sondern derselbe
      Fehlschluss mit einem anderen Vorzeichen.

    Aufruf:  python -m delfin.manta._ring_pucker tfdlokal
    """
    if not (_RDKIT and _np is not None):
        print("=== Ringlokale TFD: RDKit fehlt, uebersprungen ==="); return 0
    from rdkit.Chem import TorsionFingerprints as _TF
    fehler = 0
    # ⚠ VOR dem `try`, nicht darin.  Platzt Schritt 0, liefe sonst das Urteil unten in
    #   einen NameError -- ein Absturz, der wie "kein Befund" aussieht.
    _reihe = []
    print("=== Selbsttest: die ringlokale TFD ===")
    _alt = {k: _os.environ.get(k) for k in
            ("DELFIN_FFFREE_PUCKER_SPACE", "DELFIN_FFFREE_PUCKER_NAMP",
             "DELFIN_FFFREE_PUCKER_NPHASE", "DELFIN_FFFREE_PUCKER_FULL",
             "DELFIN_FFFREE_PUCKER_TRACE", "DELFIN_FFFREE_PUCKER_DEFEKT",
             "DELFIN_FFFREE_PUCKER_XRD", "DELFIN_FFFREE_PUCKER_CPDEDUP",
             "DELFIN_FFFREE_PUCKER_TFD_LOCAL", "DELFIN_FFFREE_PUCKER_TFD_LOCAL_KOMBI",
             "DELFIN_FFFREE_PUCKER_TFD_LOCAL_THR")}
    try:
        # ===== 0 VORGABE AUS -> BYTE-IDENTISCH =======================================
        for _k in _alt:
            _os.environ[_k] = "0"
        _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL_THR"] = ""
        for _n, _soll in sorted(_REF_GITTER.items()):
            _ist = len(_pucker_space_grid(_n, 2, 6))
            if _ist != _soll:
                print("  ✗ 0 GITTER n=%d: %d Kandidaten statt %d" % (_n, _ist, _soll))
                fehler += 1
        if not fehler:
            print("  ✓ 0 GITTER unveraendert: n=5,6,7,8 -> %s"
                  % ", ".join(str(_REF_GITTER[k]) for k in (5, 6, 7, 8)))
        _mv = Chem.AddHs(Chem.MolFromSmiles("C1CCC2CCCCC2C1"))        # Decalin
        if AllChem.EmbedMolecule(_mv, randomSeed=42) != 0:
            print("  ? 0 VORGABE: Decalin nicht einbettbar, NICHT gemessen"); fehler += 1
        else:
            AllChem.MMFFOptimizeMolecule(_mv)
            _aus = generate(_mv, budget=48)
            if len(_aus) != _REF_DECALIN_FRAMES:
                print("  ✗ 0 VORGABE VERAENDERT: Decalin liefert %d Frames statt %d -- "
                      "der Vorgabepfad ist NICHT mehr byte-identisch"
                      % (len(_aus), _REF_DECALIN_FRAMES))
                fehler += 1
            else:
                print("  ✓ 0 VORGABE UNVERAENDERT: Decalin %d Frames (Referenz %d)"
                      % (len(_aus), _REF_DECALIN_FRAMES))
            # ⚠ EIN SCHALTER OHNE REICHWEITE IST VON EINEM UNVERDRAHTETEN NICHT ZU
            #   UNTERSCHEIDEN.  In diesem Projekt fuenfmal an einem Tag passiert --
            #   darum steht die Gegenprobe direkt neben der Identitaet.
            # ⚠ DIE BEIDEN SCHALTER EINZELN, nie zusammen.  Ihre Wirkungen haben
            #   entgegengesetztes Vorzeichen; gemeinsam gemessen ergaebe die Summe eine
            #   Zahl, aus der sich kein Anteil mehr zurueckrechnen laesst.
            _za = _neuer_zaehler()
            generate(_mv, budget=48, _zaehler=_za)
            _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "1"
            _zm = _neuer_zaehler()
            _mit = generate(_mv, budget=48, _zaehler=_zm)
            _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "0"
            _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL_KOMBI"] = "1"
            _zk = _neuer_zaehler()
            _kom = generate(_mv, budget=48, _zaehler=_zk)
            _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL_KOMBI"] = "0"
            _zur = generate(_mv, budget=48)
            print("  %s 0 REICHWEITE (Decalin, budget=48), die Schalter EINZELN:"
                  % ("✓" if (len(_mit) != len(_aus) or len(_kom) != len(_aus)) else "⚠"))
            print("      AUS                %2d Frames, Zustaende je Ring %s"
                  % (len(_aus), _za["zustaende_je_ring"]))
            print("      nur je RING        %2d Frames, Zustaende je Ring %s"
                  % (len(_mit), _zm["zustaende_je_ring"]))
            print("      nur je KOMBINATION %2d Frames, Zustaende je Ring %s"
                  % (len(_kom), _zk["zustaende_je_ring"]))
            print("      zurueck auf AUS    %2d Frames" % len(_zur))
            if len(_zur) != len(_aus):
                print("      ✗ NICHT ZURUECKSCHALTBAR -- ein Schalter hinterlaesst Zustand")
                fehler += 1
            if len(_kom) == len(_aus):
                # ⚠ EIN SCHALTER OHNE GEMESSENE WIRKUNG WIRD ALS SOLCHER BENANNT.  Auf
                #   Decalin ist die Null sogar VORHERSAGBAR -- zwei gleichwertige Ringe,
                #   keine acyclische Torsion, Gewichte 1:1: global und ringlokal rechnen
                #   dort buchstaeblich dieselbe Zahl.  Das erklaert die Null, es belegt
                #   den Schalter aber nicht.  Wer ihn benutzt, misst ihn zuerst.
                print("      ⚠ DER KOMBINATIONSSCHALTER ist auf Decalin wirkungslos, und "
                      "das ist vorhersagbar: zwei gleichwertige Ringe, keine acyclische "
                      "Torsion, Gewichte 1:1 -- beide Masse rechnen dieselbe Zahl.  Er "
                      "ist damit in dieser Datei NICHT BELEGT; kein Test zeigt bisher "
                      "eine Wirkung von ihm.")
            if len(_mit) == len(_aus) and len(_kom) == len(_aus):
                print("      ⚠ BEIDE ohne Wirkung -- die Reichweite muss dann aus "
                      "Schritt 3 kommen.")
            else:
                # ⚠ HIER FAELLT DIE ZAHL, WAEHREND SIE IN SCHRITT 3 STEIGT, und das ist
                #   kein Widerspruch, sondern DIESELBE Aussage von zwei Seiten.
                #   Ringlokal heisst "nur die Torsion DIESES Rings" -- und das entfernt
                #   ZWEI Verunreinigungen auf einmal:
                #     (a) das starre Geruest im NENNER  -> Acene, Zustaende STEIGEN
                #     (b) die Bewegung des NACHBARRINGS -> Decalin, Zustaende FALLEN
                #   Bei (b) zaehlte der globale Vergleich Zustaende von Ring 1 als
                #   Zustaende von Ring 0 mit; das Kreuzprodukt zaehlt sie DANACH noch
                #   einmal.  Decalin hat kein Geruest zum Verduennen (Ringanteil 100 %,
                #   keine acyclische Torsion, Gewichte 1:1), also bleibt hier nur (b).
                # ⚠ NICHT BEWIESEN ist damit, dass die entfallenen Frames Doppelgaenger
                #   WAREN -- gezeigt ist nur, WO die Zahl sich aendert.  Wer das Urteil
                #   will, braucht das Auge, nicht diesen Test.
                print("      ⚠ HIER FAELLT die Zahl, in Schritt 3 STEIGT sie.  Dieselbe "
                      "Aussage von zwei Seiten: ringlokal entfernt das Geruest aus dem "
                      "Nenner (Acene: mehr Zustaende) UND die Bewegung des Nachbarrings "
                      "aus der Zustandszahl eines Rings (Decalin: weniger).  Decalin hat "
                      "kein Geruest -- Ringanteil 100 %, keine acyclische Torsion, "
                      "Gewichte 1:1 -- also bleibt hier nur der zweite Anteil.")

        # ===== 1 WORAN DIE VERDUENNUNG LIEGT -- RDKITS EIGENE GEWICHTE ================
        # ⚠ GERECHNET, NICHT GESCHAETZT.  `CalculateTFD` bildet sum(d_i*w_i)/sum(w_i).
        #   Bewegt sich nur EIN Ring, bleibt d_Ring * w_Ring / sum(w) -- der Quotient
        #   w_Ring/sum(w) IST also der Verduennungsfaktor, ohne jede Modellannahme.
        print()
        print("  ===== 1 DER VERDUENNUNGSFAKTOR STEHT IN RDKITS GEWICHTEN =====")
        print("    %-24s %7s %6s %6s %10s %12s"
              % ("Molekuel", "Anteil", "nring", "ring", "w_R/sum(w)", "Verduennung"))
        _wfak = {}
        for _name, _smi in _LOKAL_KALIBER + _VERD_PROBEN:
            try:
                _m = Chem.AddHs(Chem.MolFromSmiles(_smi))
                if AllChem.EmbedMolecule(_m, randomSeed=42) != 0:
                    print("    %-24s Einbettung fehlgeschlagen" % _name); continue
                AllChem.MMFFOptimizeMolecule(_m)
                _ro, _rat, _fr = _lokal_ringlage(_m)
                if not _ro:
                    print("    %-24s kein faltbarer Ring" % _name); continue
                _tl, _tlr = _TF.CalculateTorsionLists(_m)
                _w = _TF.CalculateTorsionWeights(_m)
                _ziel = frozenset(int(a) for a in _ro)
                _k = next((i for i, (_q, _d) in enumerate(_tlr)
                           if frozenset(int(t[0]) for t in _q) == _ziel), None)
                if _k is None:
                    print("    %-24s Ring NICHT in RDKits Ringliste -- "
                          "die Zuordnung ueber die Atommenge greift nicht" % _name)
                    fehler += 1
                    continue
                _hv = [i for i in range(_m.GetNumAtoms())
                       if _m.GetAtomWithIdx(i).GetSymbol() != "H"]
                _wr = _w[len(_tl) + _k] / sum(_w)
                _wfak[_name] = _wr
                print("    %-24s %6.1f%% %6d %6d %10.4f %11.1fx"
                      % (_name, 100.0 * len(_rat) / len(_hv), len(_tl), len(_tlr),
                         _wr, 1.0 / _wr))
            except Exception as _e:
                print("    %-24s ausgefallen: %s" % (_name, type(_e).__name__))
                fehler += 1
        for _name, _ in _LOKAL_KALIBER:
            if _name in _wfak and abs(_wfak[_name] - 1.0) > 1e-9:
                print("    ✗ 1 KALIBER %s hat Gewichtsanteil %.6f statt 1 -- die "
                      "Schwellenherleitung in Schritt 2 traegt dann nicht"
                      % (_name, _wfak[_name]))
                fehler += 1
        if all(abs(_wfak.get(n, 1.0) - 1.0) <= 1e-9 for n, _ in _LOKAL_KALIBER):
            print("    ✓ 1 KALIBER: unsubstituierter Einringer hat GENAU EINEN "
                  "Torsionseintrag, Gewichtsanteil 1,0000 -- dort gibt es per "
                  "Konstruktion nichts zu verduennen.")

        # ===== 2 SYMMETRIEFALTUNG UND SCHWELLE AUF DEM KALIBERRING ====================
        # (2a) DIE ZAHLEN SELBST: sind global und ringlokal auf dem unsubstituierten
        #      Ring DASSELBE?  Nicht "aehnlich" -- die Herleitung behauptet Gleichheit,
        #      also wird Gleichheit gemessen, mit Nenner.
        print()
        print("  ===== 2 SYMMETRIEFALTUNG UND SCHWELLE (unsubstituierte Ringe) =====")
        _paare_ges, _dmax_ges = 0, 0.0
        for _name, _smi in _LOKAL_KALIBER:
            try:
                _m = Chem.AddHs(Chem.MolFromSmiles(_smi))
                _ids = list(AllChem.EmbedMultipleConfs(_m, numConfs=20, randomSeed=42))
                if len(_ids) < 2:
                    print("    %-14s nur %d Konformer -- nicht messbar"
                          % (_name, len(_ids))); continue
                AllChem.MMFFOptimizeMoleculeConfs(_m)
                _ro, _rat, _fr = _lokal_ringlage(_m)
                _li = _tfd_lokal_listen(_m, (_ro,))
                _d, _np_ = 0.0, 0
                for _i in range(len(_ids)):
                    for _j in range(_i + 1, len(_ids)):
                        _g = _tfd(_m, _ids[_i], _ids[_j])
                        _l = _tfd_lokal(_m, _li, _ids[_i], _ids[_j])
                        _d = max(_d, abs(_g - _l)); _np_ += 1
                _paare_ges += _np_; _dmax_ges = max(_dmax_ges, _d)
                print("    %-14s %4d Konformerpaare, groesste Differenz "
                      "|global - ringlokal| = %.3e" % (_name, _np_, _d))
            except Exception as _e:
                print("    %-14s ausgefallen: %s" % (_name, type(_e).__name__))
                fehler += 1
        if _paare_ges and _dmax_ges <= 1e-9:
            print("    ✓ 2a IDENTISCH auf %d Konformerpaaren (groesste Differenz %.1e). "
                  "⇒ DIE SCHWELLE BLEIBT 0,05: auf dem Kaliberring sind die beiden "
                  "Masse nicht aehnlich geeicht, sondern DASSELBE." % (_paare_ges, _dmax_ges))
        elif _paare_ges:
            print("    ✗ 2a NICHT identisch: groesste Differenz %.3e ueber %d Paare -- "
                  "die Schwelle muesste dann neu geeicht werden, und die Herleitung "
                  "in Schritt 1 ist falsch." % (_dmax_ges, _paare_ges))
            fehler += 1

        # (2b) DIE ZUSTANDSZAHL -- der Test, an dem der CP-Ersatz gestorben ist.
        _os.environ["DELFIN_FFFREE_PUCKER_SPACE"] = "1"
        _os.environ["DELFIN_FFFREE_PUCKER_NAMP"] = "2"
        _os.environ["DELFIN_FFFREE_PUCKER_NPHASE"] = "8"
        print("    Zustaende je Ring, Raumgitter NAMP=2 NPHASE=8:")
        print("    %-14s %8s %8s %9s %9s" % ("Ring", "gl 0,05", "lo 0,05",
                                             "gl 0,005", "lo 0,005"))
        _n5 = None
        for _name, _smi in _LOKAL_KALIBER:
            try:
                _m = Chem.AddHs(Chem.MolFromSmiles(_smi))
                if AllChem.EmbedMolecule(_m, randomSeed=42) != 0:
                    continue
                AllChem.MMFFOptimizeMolecule(_m)
                _ro, _rat, _fr = _lokal_ringlage(_m)
                _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "0"
                _g05 = len(_ring_pucker_states(_m, _ro, _fr, 0.05))
                _g005 = len(_ring_pucker_states(_m, _ro, _fr, 0.005))
                _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "1"
                _l05 = len(_ring_pucker_states(_m, _ro, _fr, 0.05))
                _l005 = len(_ring_pucker_states(_m, _ro, _fr, 0.005))
                _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "0"
                print("    %-14s %8d %8d %9d %9d" % (_name, _g05, _l05, _g005, _l005))
                if (_g05, _g005) != (_l05, _l005):
                    print("      ✗ 2b %s: ringlokal weicht auf dem KALIBERRING ab -- "
                          "0,05 bedeutet dort dann nicht mehr dasselbe" % _name)
                    fehler += 1
                if len(_ro) == 5:
                    _n5 = _l05
            except Exception as _e:
                print("    %-14s ausgefallen: %s" % (_name, type(_e).__name__))
                fehler += 1
        if _n5 is None:
            print("    ✗ 2b FUENFRING nicht gemessen -- die Symmetrieprobe fehlt")
            fehler += 1
        elif _n5 == 3:
            print("    ✓ 2b SYMMETRIEFALTUNG: unsubstituierter Fuenfring gibt ringlokal "
                  "3 Zustaende.  Der CP-Ersatz vom 26.08. gab hier 9/13/14 -- die "
                  "Faltung ueberlebt, weil RDKits Ringeintrag der MITTELWERT von "
                  "|Torsion| ueber den Ring ist und damit nummerierungsinvariant.")
        else:
            print("    ✗ 2b SYMMETRIEFALTUNG ZERSTOERT: Fuenfring gibt %d statt 3 "
                  "Zustaende -- derselbe Fehlschluss wie beim CP-Ersatz." % _n5)
            fehler += 1

        # ===== 3 DIE VERDUENNUNGSREIHE NACHGEFAHREN ==================================
        print()
        print("  ===== 3 DIESELBE TABELLE, GLOBAL GEGEN RINGLOKAL =====")
        print("    Geruest EINGEFROREN, Raumgitter NAMP=2 NPHASE=8 -- exakt die "
              "Vorbereitung der Verduennungsreihe in `selbsttest_trennschaerfe`.")
        print("    %-24s %7s %10s %8s %9s %10s"
              % ("Molekuel", "Anteil", "w_R/sum(w)", "gl 0,05", "gl 0,005", "lo 0,05"))
        _reihe = []
        for _name, _smi in _VERD_PROBEN:
            try:
                _m = Chem.AddHs(Chem.MolFromSmiles(_smi))
                if AllChem.EmbedMolecule(_m, randomSeed=42) != 0:
                    print("    %-24s Einbettung fehlgeschlagen" % _name); continue
                AllChem.MMFFOptimizeMolecule(_m)
                _ro, _rat, _fr = _lokal_ringlage(_m)
                if not _ro:
                    print("    %-24s kein faltbarer Ring" % _name); continue
                _hv = [i for i in range(_m.GetNumAtoms())
                       if _m.GetAtomWithIdx(i).GetSymbol() != "H"]
                _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "0"
                _g05 = len(_ring_pucker_states(_m, _ro, _fr, 0.05))
                _g005 = len(_ring_pucker_states(_m, _ro, _fr, 0.005))
                _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "1"
                _l05 = len(_ring_pucker_states(_m, _ro, _fr, 0.05))
                _os.environ["DELFIN_FFFREE_PUCKER_TFD_LOCAL"] = "0"
                print("    %-24s %6.1f%% %10.4f %8d %9d %10d"
                      % (_name, 100.0 * len(_rat) / len(_hv), _wfak.get(_name, 0.0),
                         _g05, _g005, _l05))
                _reihe.append({"name": _name, "g05": _g05, "g005": _g005, "l05": _l05,
                               "anteil": len(_rat) / float(len(_hv))})
            except Exception as _e:
                print("    %-24s ausgefallen: %s" % (_name, type(_e).__name__))
                fehler += 1
    finally:
        for _k, _v in _alt.items():
            if _v is None:
                _os.environ.pop(_k, None)
            else:
                _os.environ[_k] = _v

    # ---- DAS URTEIL.  ⚠ ES DARF AUCH GEGEN DIE REPARATUR AUSFALLEN -- eine ringlokale
    #      Fassung, die den Gradienten NICHT beseitigt, ist ein Befund und kein Fehler.
    if len(_reihe) >= 2:
        _g = [r["g05"] for r in _reihe]
        _l = [r["l05"] for r in _reihe]
        _gf = _g[-1] < _g[0]                       # global faellt ueber die Reihe
        _lf = _l[-1] < _l[0]                       # ringlokal auch?
        print()
        if not _gf:
            print("  ⇒ KEIN GRADIENT IN DER GLOBALEN SPALTE (%s) -- der Befund, den "
                  "dieser Test pruefen soll, tritt auf dieser Reihe gar nicht auf.  "
                  "Das Urteil ueber die Reparatur haengt in der Luft."
                  % " -> ".join(str(x) for x in _g))
            fehler += 1
        elif _lf:
            print("  ⇒ RINGLOKAL HILFT NICHT.  global %s, ringlokal %s -- der Abfall "
                  "bleibt.  Die Verduennung war dann nicht (oder nicht allein) die "
                  "Ursache; der naechste Verdaechtige ist der Relax, nicht das Mass."
                  % (" -> ".join(str(x) for x in _g), " -> ".join(str(x) for x in _l)))
        else:
            print("  ⇒ DER GRADIENT IST WEG.  global %s (Ringanteil %.0f -> %.0f %%), "
                  "ringlokal %s -- DERSELBE Ring, DIESELBEN Kandidaten, DIESELBE "
                  "Schwelle 0,05."
                  % (" -> ".join(str(x) for x in _g), 100.0 * _reihe[0]["anteil"],
                     100.0 * _reihe[-1]["anteil"], " -> ".join(str(x) for x in _l)))
            print("    Das kleinste Glied der Reihe gewinnt %d -> %d Zustaende je Ring. "
                  "⚠ ZUSTAENDE JE RING GEHEN POTENZIERT ins Kreuzprodukt ein "
                  "(gemessen 4,3 Ringe je System) -- das ist der Hebel, nicht die "
                  "Kombinatorik dahinter." % (_g[-1], _l[-1]))
            print("    ⚠ WAS DAMIT NICHT BEWIESEN IST: dass diese Zustaende das "
                  "Realismustor ueberleben.  Diese Reihe misst (a), die Zustaende JE "
                  "RING -- Kollision, Winkel und Bindungstor sitzen dahinter.")
            print("    Getragen wird dieser Befund von DELFIN_FFFREE_PUCKER_TFD_LOCAL "
                  "allein.  DELFIN_FFFREE_PUCKER_TFD_LOCAL_KOMBI steht in dieser Reihe "
                  "NICHT im Spiel und ist mit ihr auch nicht belegt -- seine einzige "
                  "Messung ist die Decalinzeile in Schritt 0.")
    print("=== Ringlokale TFD: %s ==="
          % ("gemessen" if fehler == 0 else "%d Pruefung(en) FEHLGESCHLAGEN" % fehler))
    return 1 if fehler else 0


if __name__ == "__main__":
    # ⚠ AM DATEIENDE, und das ist keine Kosmetik.  Auf MODULEBENE zaehlt die
    #   Reihenfolge: steht dieser Block vor einer der Testfunktionen, ist ihr Name
    #   zur Ausfuehrungszeit noch ungebunden -> NameError.  (Innerhalb einer
    #   Funktion gilt das nicht -- genau die Verwechslung, die am 09.08. den
    #   [Z4]-Totenschein erzeugt hat, nur andersherum.)
    import sys as _sys
    # Ohne Argument laeuft alles -- ein Name laesst genau einen Test laufen.  Das ist
    # keine Bequemlichkeit: `selbsttest_kombinatorik` BAUT, und wer sie waehrend einer
    # Aenderung nachmessen will, soll dafuer nicht dreimal den TFD-Sweep bezahlen.
    _TESTS = (("raum", selbsttest_raum),
              ("konvergenz", selbsttest_konvergenz),
              # Der Sweep steht NACH der Konvergenz, weil er ihre offene Frage
              # beantwortet: sie meldet "Zustaende liegen dicht", er misst, ob das an
              # der Schwelle liegt.
              ("sweep", selbsttest_tfd_sweep),
              # Die Trennschaerfe steht vor der Kombinatorik: sie entscheidet, WELCHES
              # Mass die Kombinatorik ueberhaupt entdoppeln soll.  Ein Kostenurteil mit
              # dem falschen Entdopplungsmass waere ein Urteil ueber das Instrument.
              ("trennschaerfe", selbsttest_trennschaerfe),
              # Direkt DAHINTER, weil die Trennschaerfe den Befund abwirft, den dieser
              # Test repariert: sie misst, dass derselbe Ring an wachsendem Geruest
              # immer weniger Zustaende bekommt, er misst dieselbe Reihe noch einmal
              # mit ringlokaler TFD.  Getrennt gelaufen waeren es zwei Messungen; so
              # ist es eine Messung und ihre Gegenprobe.
              ("tfdlokal", selbsttest_tfd_lokal),
              # Zuletzt die Kombinatorik: sie baut Relax + Tor JE Kombination und ist
              # damit der teuerste der vier.  Sie beantwortet, was die drei davor
              # aufwerfen -- die Zustandszahl je Ring ist nur interessant, weil sie
              # potenziert wird.
              ("kombinatorik", selbsttest_kombinatorik))
    _wahl = [a for a in _sys.argv[1:] if not a.startswith("-")]
    _unbekannt = [a for a in _wahl if a not in dict(_TESTS)]
    if _unbekannt:
        print("unbekannter Test: %s -- bekannt: %s"
              % (", ".join(_unbekannt), ", ".join(n for n, _ in _TESTS)))
        _sys.exit(2)
    # `--ohne-grenze` hebt die Testgrenze der Kombinatorik auf.  Dann wird JEDE Probe
    # vollstaendig gebaut -- auch die, deren Kreuzprodukt fuenfstellig ist.  Das ist die
    # Messung, keine Vorgabe: sie laeuft Stunden und gehoert nicht in einen Regellauf.
    _ohne_grenze = "--ohne-grenze" in _sys.argv[1:]
    _rc = 0
    for _name, _fn in _TESTS:
        if _wahl and _name not in _wahl:
            continue
        _rc = (_fn(max_kombis=0) if (_name == "kombinatorik" and _ohne_grenze)
               else _fn()) or _rc
    _sys.exit(_rc)
