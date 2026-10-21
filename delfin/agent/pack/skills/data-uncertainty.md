---
domains: data
---
# State the uncertainty with every figure
> A number without its spread and population is a claim, not a result.

Use whenever you report a mean, difference, or comparison from data.

1. `summarize(values)` gives n, mean, sample std, median, min, max. Report
   at least mean ± std and n. If n is small (< 5), say "n too small for a
   meaningful spread" instead of quoting std as if it were precise.
2. Convert before comparing, never after: kJ/mol and kcal/mol in one column
   is two populations (`unit_consistency` flags it). Quote the conversion
   factor you applied (1 kcal/mol = 4.184 kJ/mol).
3. A difference between two groups is not an effect until the spreads
   overlap-check says so: report both groups' n, mean and spread, and say
   plainly if the difference is within one spread of zero.
4. A fit or trend claim carries residuals: if you cannot state how far the
   model is from the points, do not claim the trend.
5. Distinguish computed from assumed: a measured spread is quoted with n; an
   estimate from fewer points says so in the same sentence.

Constraints:
- Never report a mean without the n it came from.
- Never drop the low tail because it is inconvenient; name every exclusion.
- If the population splits (e.g. a trend that holds within each subgroup but
  reverses in the aggregate), report the split — the aggregate alone is a
  wrong answer.

Output format: every figure as value ± spread (n = N, rule or method named),
population changes named, comparisons stated with both spreads.
