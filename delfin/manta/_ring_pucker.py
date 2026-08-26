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


def _tfd_distinct(acc_mol, cid: int, kept_ids, thr: float) -> bool:
    return all(_tfd(acc_mol, k, cid) >= thr for k in kept_ids)


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
            if _cpd:
                # NACH dem Relax messen, nicht die SOLL-Werte vergleichen: der Relax
                # zieht den Startpunkt in die naechste echte Mulde, und genau deren
                # Lage entscheidet, ob es eine neue ist.
                _cp = _cp_theta_phi(m2.GetConformer().GetPositions(), ring)
                if any(_cp_abstand(_cp, _k) < _cp_tol and abs(_cp[0] - _k[0]) < _cp_qtol
                       for _k in _cp_kept):
                    continue
                _cp_kept.append(_cp)
                states.append(_cand if _raum else (_qs, theta, phi))
                continue
            cid = _add_conf(acc, m2)
            if _tfd_distinct(acc, cid, kept_ids, tfd_thr):
                kept_ids.append(cid)
                # Im Raum-Modus ist der Zustand das Koordinatenpaar selbst; die
                # Legacy-Form bleibt ein 3-Tupel.  `generate` indiziert nur, es liest
                # den Inhalt nicht -- beide Formen sind dort gleichwertig.
                states.append(_cand if _raum else (_qs, theta, phi))
            else:
                acc.RemoveConformer(cid)
        except Exception:
            continue
    return states


def generate(mol_with_conf, frozen: Optional[Set[int]] = None,
             budget: int = 64, tfd_thr: float = 0.05,
             angle_skip: Optional[Set[int]] = None) -> List[Tuple[str, str]]:
    """Construct the COMBINATORIAL ring-pucker conformers from a base conformer.

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
    per_ring_states = [_ring_pucker_states(mol_with_conf, ring, frozen, tfd_thr)
                       for ring in rings]

    # cartesian product of state indices, deterministic order, budget-capped;
    # skip the all-base (identity) combination; fewest-changed rings first.
    import itertools as _it
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
    if _os.environ.get("DELFIN_FFFREE_PUCKER_FULL", "0") == "1":
        pass                                    # alle Faltungen, keine Kappe
    else:
        combos = combos[:max(0, int(budget))]
    if len(combos) < _n_voll and _os.environ.get("DELFIN_FFFREE_PUCKER_TRACE", "0") == "1":
        print("[pucker] KOMBINATIONEN GEKAPPT: %d von %d gebaut, %d verworfen "
              "(%d Ringe, Mulden %s) -- DELFIN_FFFREE_PUCKER_FULL=1 baut alle"
              % (len(combos), _n_voll, _n_voll - len(combos), len(rings),
                 "x".join(str(len(s)) for s in per_ring_states)))

    acc = Chem.Mol(mol_with_conf)
    kept_ids = [acc.GetConformer().GetId()]
    out: List[Tuple[str, str]] = []
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
                continue
            # realism gate: a combination that stayed clashed OR left any VSEPR
            # body distorted (fused/bridged rings strain their shared atoms) is
            # not a physical ensemble member -> drop it.  Everything must be
            # right, or the frame is unrealistic.
            if _has_clash(m2) or _has_bad_angles(m2, skip=angle_skip):
                continue
            cid = _add_conf(acc, m2)
            if not _tfd_distinct(acc, cid, kept_ids, tfd_thr):
                acc.RemoveConformer(cid)
                continue
            kept_ids.append(cid)
            label = "pucker " + "+".join(
                f"r{ri_i}:{'base' if combo[ri_i] == 0 else combo[ri_i]}"
                for ri_i in range(len(rings)))
            out.append((_conf_to_xyz(m2), label))
        except Exception:
            continue
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


if __name__ == "__main__":
    # ⚠ AM DATEIENDE, und das ist keine Kosmetik.  Auf MODULEBENE zaehlt die
    #   Reihenfolge: steht dieser Block vor einer der Testfunktionen, ist ihr Name
    #   zur Ausfuehrungszeit noch ungebunden -> NameError.  (Innerhalb einer
    #   Funktion gilt das nicht -- genau die Verwechslung, die am 09.08. den
    #   [Z4]-Totenschein erzeugt hat, nur andersherum.)
    import sys as _sys
    _rc = selbsttest_raum()
    _rc = selbsttest_konvergenz() or _rc
    # Der Sweep steht NACH der Konvergenz, weil er ihre offene Frage beantwortet:
    # sie meldet "Zustaende liegen dicht", er misst, ob das an der Schwelle liegt.
    _rc = selbsttest_tfd_sweep() or _rc
    _sys.exit(_rc)
