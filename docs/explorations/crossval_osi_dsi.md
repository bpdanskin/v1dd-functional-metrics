# Split-half cross-validated OSI/DSI vs the naive metric

Our shipped OSI and DSI (`families/drifting_gratings.py`, `dg_metrics_from_trials`) pick the
preferred condition by argmax over the trial-averaged (direction × SF) grid and then measure
selectivity **at that same condition on the same trials**. Selecting and measuring on one
draw biases selectivity **upward**: the argmax rides the noise peak while the orthogonal and
null directions ride their noise troughs, so `(pref − orth)/(pref + orth)` is inflated by
construction. The inflation is worst exactly where it matters least — weakly-driven cells,
where the ranking is mostly noise.

The **split-half cross-validation** (Jun Zhuang's `*_shw_*` columns, `v1dd_physiology`)
breaks the coupling: on each of many iterations, split the 8 trials 4/4, pick the preferred
condition on one half, measure OSI/DSI/gOSI on the *other* half, and average over iterations.
The held-out measurement can't inherit the selection's noise, so the estimate is
(near-)unbiased. Two by-products fall out: the spread of the estimate across splits (how
noisy an 8-trial recording is), and a reliability count — how many splits had the held-out
preferred response clear a z-score threshold above the blank (Jun's `iter_num_2`).

Runs fully locally on the asset. Harness: `code/validation/crossval_osi_dsi.py`
(200 iterations by default; reproduces the shipped naive columns to 1e-7 as a control).
Per-ROI output: `crossval_osi_dsi.csv` (naive vs CV OSI/DSI/gOSI, per-split std,
preferred-direction stability, reliability fractions, join keys). Numbers below are windowed
gratings (`dgw`); full-field (`dgf`) is within 0.02 of every figure.

![Split-half cross-validated OSI/DSI vs naive: (A) OSI distribution collapses toward zero
under cross-validation, (B) most high-OSI cells are non-responsive with true OSI near zero,
(C) OSI/DSI carry a +0.20 selection bias but gOSI carries none, (D) the reliability count
tracks is_responsive](../figures/crossval_osi_dsi.png)

## Result — 39,407 ROIs

**The naive OSI/DSI bias is large, and it is a selection artifact.**
Across all cells, naive OSI averages **0.50** and cross-validated OSI **0.12** (bias
**+0.38**); DSI 0.43 → 0.10 (+0.33). Restricting to the 10,817 responsive cells, the honest
correction is still **+0.20** for both OSI and DSI (OSI 0.64 → 0.44, DSI 0.48 → 0.28). The
artifact is stark at the tail: **46% of all cells have naive OSI > 0.5, and 57% of those are
not responsive** — their true (cross-validated) OSI has a median of **0.16** and a
distribution centred near zero (panel B). The naive metric manufactures an
"orientation-selective" label for the noisy majority.

**gOSI is already unbiased — OSI/DSI are the culprits.**
gOSI (the vector-sum `|Σ r·e^{2iθ}|/Σr`) integrates over all directions instead of
argmax-selecting one, so it has no selection to over-fit: bias **−0.01** for responsive
cells, and naive-vs-CV Spearman **0.97** (panel C). The whole problem is specific to the
argmax-then-ratio construction of OSI and DSI.

**An 8-trial recording is intrinsically noisy for this estimate.**
The per-split standard deviation of CV-OSI is **~0.29** for responsive cells (panel C error
bars) — each 4+4 split is nearly uninformative on its own, and iteration-averaging is what
rescues the mean. The averaged estimate is stable and preserves rank among responsive cells
(naive-vs-CV Spearman 0.88), but the value carries an irreducible 8-trial noise floor that
more iterations cannot remove. Preferred **direction** is also unstable: only 67% of splits
agree with the naive pick even for responsive cells.

**The bias is only weakly SNR-dependent, so no per-cell tuning.**
Across SNR tertiles of responsive cells the OSI bias runs +0.22 / +0.20 / +0.19 (low → high)
— present at every SNR, the responsive/non-responsive split matters far more than SNR. This
mirrors the greedy-RF finding: correct the metric globally, don't tune a knob per cell.

**The reliability count is a re-derivation of `is_responsive`.**
The Jun-style reliability fraction (splits with held-out z ≥ 2) correlates 0.65 with
`dgw_is_responsive`, and a `reliab ≥ 0.5` gate agrees with `is_responsive` on **85%** of
cells (panel D), keeping 20% vs 27%. A separate reliability column, or NaN-gating OSI by
reliability, would largely duplicate a gate the asset already ships.

## What this means for the three ship options

The brief laid out three ways to fold reliability into OSI/DSI. On the data:

- **(a) Replace OSI/DSI with the cross-validated value.** Removes the +0.20 bias and the
  false-selectivity tail in one move, and folds reliability *into the number*: an unreliable
  cell's CV-OSI collapses toward 0 on its own, no separate gate needed — which is exactly the
  "reliability living inside OSI/DSI" the user asked for. It stays one column per metric
  (most compact). Costs: the value is non-deterministic (seed + iteration count) with an
  8-trial noise floor, and it changes every number — but this repo owes no output
  compatibility to `allen_v1dd` or the white paper (`CLAUDE.md`), so that cost is low here.
- **(b) NaN-gate OSI/DSI by reliability.** Most compact in spirit, but the reliability signal
  is ~`is_responsive` (85% agreement), so this mostly reproduces filtering on the
  `is_responsive` column we already have, while making OSI present-for-some-cells-only — the
  cross-cell-comparability problem we deliberately avoided with the greedy-RF global α.
- **(c) Separate reliability column.** Least compact, and again ≈ `is_responsive`. Little
  added information for the extra column.

**Recommendation.** Option (a) for OSI and DSI, keeping gOSI as-is (already unbiased; CV
moves it by 0.01). Replacement is the compact choice that honours "reliability inside the
metric" — the unreliable cells self-report as un-selective rather than needing a companion
gate — and it retires a metric that currently labels a non-responsive plurality as selective.
If the per-cell noise of the CV value is a concern, the fallback is to **keep naive OSI/DSI
but add the CV columns alongside** (additive, like the greedy `_a05` variant), and let
downstream choose; that trades compactness for determinism. The decision is yours — this
exploration is not wired into the pipeline.

## If option (a) is chosen — implementation notes

- Compute per plane inside `dg_metrics_from_trials`, reusing the same trial array; ~15 s for
  the whole asset at 200 iterations, so cost is negligible.
- Fix the seed and iteration count in `config` (e.g. `dg_crossval_iters`, default 200) for
  reproducibility, and document that the value has an 8-trial noise floor the seed pins but
  cannot remove.
- `pref_dir_mean` and `preferred_dir` should come from the naive full-8-trial pick (stable,
  deterministic); cross-validation is about the *selectivity magnitude*, not the reported
  preference.
- Note the change in `docs/comparability.md` (OSI/DSI are now cross-validated, ~0.2 lower than
  the white-paper convention on the same cells) and in the family doc.
