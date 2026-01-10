# occupier_selection.py
# OCCUPIER candidate selection logic extracted from occupier_reports.py

from typing import List, Tuple, Optional, Dict, Any
from decimal import Decimal, ROUND_DOWN


def truncate(x: float, d: int) -> float:
    """Truncate a float to d decimal places."""
    q = Decimal(10) ** -d
    return float(Decimal(str(x)).quantize(q, rounding=ROUND_DOWN))


def fmt_truncate(x: float, d: int) -> str:
    """Format a truncated float with d decimal places."""
    return f"{truncate(x, d):.{d}f}"


class OccupierSelector:
    """Handles OCCUPIER candidate selection logic.

    This class encapsulates the complex selection algorithm that was previously
    embedded in the large generate_summary_report_OCCUPIER function.
    """

    def __init__(self, config: Dict[str, Any], sequence: List[Dict],
                 dev_cache: Dict[int, Optional[float]]):
        self.config = config
        self.sequence = sequence
        self.dev_cache = dev_cache

        # Load configuration parameters
        self._load_parameters()

    def _load_parameters(self):
        """Load selection parameters from config."""
        # Selection mode and precision
        raw_sel = str(self.config.get('occupier_selection', 'tolerance')).lower()
        self.method = raw_sel.split('|')[0].strip()

        if self.method in {'rounding', 'round', 'gerundet', 'runden'}:
            self.method = 'rounding'
        elif self.method in {'tolerance', 'toleranz', 'toleranzband', 'epsilon'}:
            self.method = 'tolerance'
        else:
            self.method = 'truncation'

        try:
            self.precision = int(self.config.get('occupier_precision',
                                               self.config.get('occupier_rounded_value', 6)))
        except (TypeError, ValueError):
            self.precision = 6
        self.precision = max(0, min(10, self.precision))

        try:
            self.epsilon = float(self.config.get('occupier_epsilon', 10.0**(-self.precision)))
            if not (self.epsilon > 0):
                raise ValueError
        except (TypeError, ValueError):
            self.epsilon = 10.0**(-self.precision)

        # Dev max filter
        try:
            self.dev_max = float(self.config.get('dev_max')) if self.config.get('dev_max') is not None else None
        except (TypeError, ValueError):
            self.dev_max = None

        # Selection thresholds
        self.DEV_TINY = 1e-3
        self.DEV_SIMILARITY = float(self.config.get('dev_similarity', 0.15))
        self.DEV_GOOD_MARGIN = 0.30
        self.DEV_HIGH = 0.50
        self.DEV_MATCH_WINDOW = 0.30

        # Window parameters
        self.EPS_AF = float(self.config.get('bs_override_window_h', 0.001))
        self.E_BIAS_H = float(self.config.get('energy_bias_window_h', 0.001))
        self.MIS_BIAS = float(self.config.get('mismatch_bias_window', 0.05))

        # Clean preference parameters
        self.CLEAN_OVERRIDE_H = float(self.config.get('clean_override_window_h', 0.002))
        self.CLEAN_Q_IMPROVE = float(self.config.get('clean_quality_improvement', 0.05))
        self.CLEAN_Q_GOOD = float(self.config.get('clean_quality_good', 0.05))
        self.CLEAN_BIAS_H = float(self.config.get('clean_bias_window_h', 0.001))
        self.QUAL_BIAS_WIN = float(self.config.get('quality_bias_window', 0.05))

    def effective_dev(self, idx: int) -> float:
        """Get effective deviation for an index."""
        dev = self.dev_cache.get(idx)
        if dev is None:
            if (not self.entry_is_bs(idx)) and (self.entry_mult(idx) == 1):
                return 0.0  # closed-shell RKS
            return float("inf")
        return dev

    def entry_record(self, idx: int) -> Dict:
        """Get sequence entry record for index."""
        return next((e for e in self.sequence if e["index"] == idx), {})

    def entry_is_bs(self, idx: int) -> bool:
        """Check if entry is broken symmetry."""
        return bool(self.entry_record(idx).get("BS"))

    def entry_mult(self, idx: int) -> int:
        """Get multiplicity for entry."""
        try:
            return int(self.entry_record(idx).get("m", 0))
        except:
            return 0

    def bs_pair_count(self, idx: int) -> Optional[int]:
        """Get BS pair count for entry."""
        lab = str(self.entry_record(idx).get("BS") or "").strip()
        if not lab or "," not in lab:
            return None
        try:
            _, n = lab.split(",", 1)
            return int(n.strip())
        except Exception:
            return None

    def bs_mismatch(self, idx: int) -> float:
        """Calculate BS pair mismatch."""
        n = self.bs_pair_count(idx)
        if n is None:
            return float("inf")
        return abs(self.effective_dev(idx) - n)

    def is_pseudo_closed_shell(self, idx: int) -> bool:
        """Check if entry is pseudo closed-shell."""
        return (self.entry_is_bs(idx) and
                (self.entry_mult(idx) == 1) and
                (self.effective_dev(idx) <= self.DEV_TINY))

    def within_dev_limit(self, idx: int) -> bool:
        """Check if deviation is within configured limit."""
        if self.dev_max is None:
            return True
        d = self.dev_cache.get(idx)
        return (d is None) or (d <= self.dev_max)

    def energy_key(self, val: float) -> float:
        """Apply energy processing based on method."""
        if self.method == "rounding":
            return round(val, self.precision)
        if self.method == "truncation":
            return truncate(val, self.precision)
        return val  # tolerance -> raw

    def select_candidates(self, fspe_values: List[Optional[float]]) -> Tuple[Optional[int], Optional[float]]:
        """Select the best candidate using OCCUPIER selection logic.

        Returns:
            Tuple of (best_index, best_energy) or (None, None) if no valid candidates
        """
        # Filter valid entries
        valid_all = [(e["index"], f) for e, f in zip(self.sequence, fspe_values) if f is not None]
        valid = [pair for pair in valid_all if self.within_dev_limit(pair[0])] or valid_all

        if not valid:
            return None, None

        energies_by_idx = {i: f for i, f in valid_all}

        # Apply energy-based filtering
        if self.method in ("rounding", "truncation"):
            processed = [(i, f, self.energy_key(f)) for i, f in valid]
            best_key = min(v for _, _, v in processed)
            cands = [(i, f) for i, f, v in processed if v == best_key]
            min_raw = min(f for _, f in valid)
        else:
            min_raw = min(f for _, f in valid)
            eps_eff = self.epsilon + 1e-12
            cands = [(i, f) for i, f in valid if (f - min_raw) <= eps_eff]

        # AF-override: add BS candidates with high contamination
        cands = self._apply_af_override(cands, valid, min_raw)

        # Clean-override: add cleaner candidates in energy window
        cands = self._apply_clean_override(cands, valid, min_raw)

        # Tie-breaking
        return self._resolve_ties(cands, energies_by_idx)

    def _apply_af_override(self, cands: List[Tuple[int, float]],
                          valid: List[Tuple[int, float]],
                          min_raw: float) -> List[Tuple[int, float]]:
        """Apply AF-override logic for highly contaminated candidates."""
        if len(cands) == 1:
            best_idx, _ = cands[0]
            if self.effective_dev(best_idx) >= self.DEV_HIGH:
                extra = []
                for i, f in valid:
                    if i == best_idx or not self.entry_is_bs(i):
                        continue
                    d = self.effective_dev(i)
                    if (d >= self.DEV_HIGH and
                        abs(d - self.effective_dev(best_idx)) <= self.DEV_SIMILARITY and
                        (f - min_raw) <= self.EPS_AF + 1e-12):
                        extra.append((i, f))
                if extra:
                    cands = cands + extra
        return cands

    def _apply_clean_override(self, cands: List[Tuple[int, float]],
                             valid: List[Tuple[int, float]],
                             min_raw: float) -> List[Tuple[int, float]]:
        """Apply clean-override logic to add cleaner candidates."""
        if len(cands) == 1:
            best_idx, _ = cands[0]
            q_best = (abs(self.effective_dev(best_idx) - self.bs_pair_count(best_idx))
                     if self.entry_is_bs(best_idx) else self.effective_dev(best_idx))
            extra = []
            for i, f in valid:
                if i == best_idx:
                    continue
                if (f - min_raw) <= self.CLEAN_OVERRIDE_H + 1e-12:
                    q_i = (abs(self.effective_dev(i) - self.bs_pair_count(i))
                          if self.entry_is_bs(i) else self.effective_dev(i))
                    if ((q_best - q_i) >= self.CLEAN_Q_IMPROVE) or (q_i <= self.CLEAN_Q_GOOD):
                        if (i, f) not in cands:
                            extra.append((i, f))
            if extra:
                cands = cands + extra
        return cands

    def _resolve_ties(self, cands: List[Tuple[int, float]],
                     energies_by_idx: Dict[int, float]) -> Tuple[Optional[int], Optional[float]]:
        """Resolve ties between candidates."""
        if len(cands) == 1:
            return cands[0]
        elif not cands:
            return None, None

        # Reclassify pseudo-CS BS to non-BS
        raw_bs = [(i, f) for i, f in cands if self.entry_is_bs(i)]
        raw_nb = [(i, f) for i, f in cands if not self.entry_is_bs(i)]
        nb_cands = raw_nb + [(i, f) for i, f in raw_bs if self.is_pseudo_closed_shell(i)]
        bs_cands = [(i, f) for i, f in raw_bs if not self.is_pseudo_closed_shell(i)]

        # Apply regime-based selection logic
        pick = self._select_by_regime(nb_cands, bs_cands, cands)

        if pick is None:
            return None, None

        # Apply bias corrections
        return self._apply_bias_corrections(pick, cands, energies_by_idx)

    def _select_by_regime(self, nb_cands: List[Tuple[int, float]],
                         bs_cands: List[Tuple[int, float]],
                         all_cands: List[Tuple[int, float]]) -> Optional[Tuple[int, float]]:
        """Select candidate based on contamination regime."""
        def score_nb(i: int):
            pseudo_flag = 1 if self.is_pseudo_closed_shell(i) else 0
            return (self.effective_dev(i), pseudo_flag, i)

        def score_bs(i: int):
            mis = self.bs_mismatch(i)
            return (mis > self.DEV_MATCH_WINDOW, mis, self.effective_dev(i), i)

        def _min_eff_dev(pairs):
            vals = [self.effective_dev(i) for i, _ in pairs]
            return min(vals) if vals else float("inf")

        min_dev_bs = _min_eff_dev(bs_cands)
        min_dev_nb = _min_eff_dev(nb_cands)

        if bs_cands and nb_cands:
            # A) non-BS clearly cleaner and not high
            if ((min_dev_nb + self.DEV_GOOD_MARGIN < min_dev_bs) and
                (min_dev_nb < self.DEV_HIGH)):
                return min(nb_cands, key=lambda p: score_nb(p[0]))
            # B) similar & both high → best BS by mismatch
            elif ((abs(min_dev_nb - min_dev_bs) <= self.DEV_SIMILARITY) and
                  (min_dev_nb > self.DEV_HIGH) and (min_dev_bs > self.DEV_HIGH)):
                return min(bs_cands, key=lambda p: score_bs(p[0]))
            # C) Fallback
            else:
                def fallback_score(i: int):
                    return (1 if self.entry_is_bs(i) else 0, self.effective_dev(i), i)
                return min(all_cands, key=lambda p: fallback_score(p[0]))
        elif bs_cands:
            return min(bs_cands, key=lambda p: score_bs(p[0]))
        else:
            return min(nb_cands, key=lambda p: score_nb(p[0]))

    def _apply_bias_corrections(self, pick: Tuple[int, float],
                               cands: List[Tuple[int, float]],
                               energies_by_idx: Dict[int, float]) -> Tuple[int, float]:
        """Apply clean-bias and energy-bias corrections."""
        def _qual_metric(i: int) -> float:
            return self.bs_mismatch(i) if self.entry_is_bs(i) else self.effective_dev(i)

        pick_i, _ = pick
        pick_E = energies_by_idx.get(pick_i, float("inf"))
        pick_Q = _qual_metric(pick_i)

        # 1) Clean-Bias: if energy close, take significantly cleaner
        for j, _ in cands:
            if j == pick_i:
                continue
            Ej = energies_by_idx.get(j, float("inf"))
            Qj = _qual_metric(j)
            energy_close = abs(Ej - pick_E) <= self.CLEAN_BIAS_H
            if energy_close and (pick_Q - Qj) >= self.QUAL_BIAS_WIN:
                pick_i, pick_E, pick_Q = j, Ej, Qj

        # 2) Energy-Bias: if quality similar, take lower energy
        for j, _ in cands:
            if j == pick_i:
                continue
            Ej = energies_by_idx.get(j, float("inf"))
            Qj = _qual_metric(j)
            close_in_quality = abs(Qj - pick_Q) <= self.MIS_BIAS
            close_in_energy = abs(Ej - pick_E) <= self.E_BIAS_H
            if close_in_quality and close_in_energy and Ej < pick_E:
                pick_i, pick_E, pick_Q = j, Ej, Qj

        return pick_i, pick_E
