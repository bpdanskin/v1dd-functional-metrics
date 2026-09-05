---
name: fresh-start-metric-changes
description: "Four metrics the fork shipped knowingly wrong for back-compat, two columns it deferred, and the evidence for changing each now."
metadata:
  node_type: memory
  type: project
  modified: 2026-09-03
---

The fork shipped several imperfections deliberately, each justified only by not moving
published columns. This repo owes no such compatibility, so each becomes a fix. **Every one
needs a before/after comparison**, so the reason a number moved is never ambiguous.

## Corrections

1. **`lifetime_sparseness` moves to the condition-mean convention.** The fork computed it
   over *every individual trial* (118 images x 8 trials = 944 numbers per ROI). Vinje &
   Gallant, and de Vries after them, take the mean response per stimulus. Trial variance
   adds spread without adding mean, so the flattened form runs high:

   | | trial-flattened | condition means | de Vries |
   |---|---|---|---|
   | natural_images (118) | 0.9507 | **0.7373** | **0.77** |
   | natural_images_12 (12) | 0.9486 | 0.4023 | — |

   For NI12 the two correlate at **r = -0.005** — the shipped column is effectively
   uninformative about image selectivity. An earlier diagnosis blamed the 2-sample response
   window; that is **not supported**, the convention explains the gap entirely.
2. **Unify the zero-denominator convention.** `_ratio` returns 0 where `_metric_index`
   returns NaN. Use NaN — undefined is not zero.
3. **`ssi_tuning_fit` baseline asymmetry** — it evaluates the fitted curve *including* its
   baseline while the preferred direction feeding it is chosen with the baseline
   subtracted. Make them consistent.
4. **Wrap the array writers.** Provenance and the manifest run *after* them, so a raise in
   either forfeits the whole run's provenance — this cost the first 2026-09-03 run three
   files and five hours, unrebuildable. Record the failure, print loudly, carry on.

## Two columns deferred to "the next schema change" — this is that change

- **`run_corr_dff`** — Pearson correlation of running speed against each ROI's continuous
  dF/F trace (de Vries Fig. 5; their excitatory median ~0.03). Deliberately on **dF/F, not
  events**: a correlation has no denominator, so it does not inherit the sign instability
  that made `run_mod_*` events-only. Needs **no state split**, so it is finite in every
  session — unlike `run_mod_dgf` (all-NaN in 2 of 25) and `run_mod_dgw` (4 of 25).
  Not recoverable post-hoc: no time series ships in the asset.
- **`*_n_trials_at_pref`** for NI / NI12 / NM — the per-ROI denominator needed to attach a
  binomial tail p-value to `frac_responsive_trials`, as is already possible for gratings.
  Cheaper than another trial-level archive (natural movie alone is 3,600 x 9).

## Deliberately not pipeline columns

Post-hoc from the shipped arrays, so they belong in `docs/notebooks/`: binomial p for
gratings, noise correlations (6.88 M within-plane pairs, median +0.095), trial-level
running correlation, RF area. Excluded outright: preferred temporal frequency (undefined —
the grating axes are transposed), CCmax (their 0.25 s window is 1.5 samples at 6 Hz),
decoding and Gabor-wavelet models.

See [[metric-cautions]] for what any of these needs said about it in docs.

## Measured on the reference asset (2026-09-03)

Both P2 changes, shipped convention -> current default, via
`replay_reference.py`'s "intended metric changes" block:

| | median before | median after | ROIs moved | r |
|---|---|---|---|---|
| `dgw_lifetime_sparseness` | 0.6804 | 0.2596 | 39,407 | 0.799 |
| `dgf_lifetime_sparseness` | 0.6874 | 0.2740 | 39,407 | 0.796 |
| `osi`, `dsi` (both types) | — | unchanged | **0** | 1.000 |

**The zero-denominator unification moves nothing on this asset.** No ROI has an exactly
zero `pref + orth` or `pref + null`: invalid ROIs are already NaN, and no valid ROI is
exactly 0.0 at both directions. It is a guard against a case this data does not contain,
not a correction that changes numbers -- say so rather than claiming it fixed something.

The drifting-gratings sparseness drop (0.68 -> 0.26) is much larger than the natural-images
one (0.95 -> 0.74) recorded above, and in the expected direction: gratings have 24
conditions against 118 images, and the `1 - 1/n` normaliser makes sparseness fall with
fewer conditions. Same caution as [[metric-cautions]] -- never compare across stimulus
sets of different size.

## Status 2026-09-03: all four corrections and both columns are implemented

| change | flag | state |
|---|---|---|
| lifetime sparseness over condition means | `lifetime_sparseness_over` | done, measured above |
| zero denominator -> NaN | `zero_denominator_nan` | done; **moves nothing on this data** |
| `ssi_tuning_fit` baseline consistency | `ssi_tuning_fit_includes_baseline` | done; **measured** -- median 0.3291 -> 0.3564, 36,323 of 37,016 ROIs move, r = 0.8407 |
| array writers guarded | — | done in `pipeline._guarded`, covered by a test |
| `run_corr_dff` | — | done in `roi_quality`, on dF/F, no state split |
| `*_n_trials_at_pref` | — | done for `ni` / `ni12` / `nm` |

The wide table is therefore **85 columns**, not 81: `run_corr_dff` plus
`n_trials_at_pref` on each of the three natural-stimulus families.

**Measured 2026-09-03**, once the replay learned to walk planes: `ssi_tuning_fit` moves
on **98 % of ROIs with a fitted curve** (36,323 of 37,016), median 0.3291 -> 0.3564,
r = 0.8407. Driven from the *shipped* fit parameters, so it is exact arithmetic and the
fit's own instability is not in the number. The size is expected -- the baseline was added
to both terms of `(W - F) / (W + F)`, which shifts a ratio rather than cancelling.
