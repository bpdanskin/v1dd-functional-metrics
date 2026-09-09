# Drifting gratings

Direction and spatial-frequency tuning from full-field and windowed drifting gratings. The
two stimulus types share one implementation; surround suppression is what compares them.

## Goal

Characterise each neuron's preferred direction, spatial frequency, and how sharply it is
tuned. The full-field stimulus drives every neuron in the field of view; the windowed
stimulus is a small patch centred on each cortical column's population receptive field, so
it isolates the classical receptive field from its surround.

## Stimulus structure

12 directions × 2 spatial frequencies × 8 trials = 192 grating sweeps, plus a variable
number of blank (mean-luminance) sweeps interleaved. Each sweep lasts ~2 s. Blank sweeps
are **ragged**: 184–187 grating sweeps land in 192 condition slots, so some conditions get 7
trials instead of 8. `frac_responsive_trials` therefore takes 23 distinct values in the
shipped asset, not 9.

The trial array is `(n_rois, 12, 2, n_trials)`, NaN-padded on the trial axis for short
conditions. Every downstream reduction is NaN-aware.

## Selectivity indices

All indices are computed from the condition-mean tuning curve at the preferred spatial
frequency.

**OSI** — orientation selectivity. `(pref − orth) / (pref + orth)`, where `orth` is the
mean of the two orthogonal directions (±90° from preferred). A neuron responding equally to
all orientations scores 0; one responding only at its preferred orientation scores 1.

**DSI** — direction selectivity. `(pref − null) / (pref + null)`, where `null` is the
response at 180° from preferred. OSI compares across orientations; DSI compares the two
directions within the preferred orientation.

**gOSI** — global orientation selectivity. The magnitude of the orientation-tuning vector
normalised by the total response: `|Σ R(θ)·exp(2iθ)| / Σ R(θ)`. Unlike OSI, which reads
only three directions, this uses all twelve and is sensitive to the shape of the whole curve.

**`pref_dir_mean`** — the circular mean direction, treating responses as weights. Computed
with `nan_to_num(nan=0)`, so a NaN response at a direction contributes zero weight rather
than poisoning the mean.

**None of these rectifies the tuning curve**, which is only safe on a non-negative trace.
Deconvolved events are non-negative by construction; on dF/F the sum in `gosi`'s
denominator can approach zero from cancellation, making the ratio unstable. Do not point
`trace_type` at dF/F without understanding this.

## The two preferred conditions

Two argmax definitions are computed, deliberately:

1. **`fillna(-1).argmax`** — what `preferred_dir` and `preferred_sf` in the published table
   use, and what surround suppression keys off. An all-NaN ROI reports condition 0.
2. **NaN-skipping argmax** — what `osi`, `dsi`, and `gosi` use. An all-NaN ROI reports NaN.

With `pref_cond_fillna=False` (the default), both use the NaN-skipping form and agree. The
historical behaviour is `pref_cond_fillna=True`, where they diverge on all-NaN ROIs: a
fabricated preference at condition 0 propagates into surround suppression without touching
any drifting-gratings column. A divergence prints a warning.

## Responsiveness

A trial is significant if its response exceeds the 95th percentile of a single-draw
bootstrap from the spontaneous block (2,500 draws for gratings). `frac_responsive_trials`
is the fraction of trials at the preferred condition that pass.

`is_responsive` requires `frac >= 0.50` **and** `pika_roi_confidence > 0.5`. The 0.50
threshold is much stricter than de Vries' 0.25 — but they had 15 trials per condition;
here, with 8, the false-positive rates are:

| criterion | false-positive rate |
|---|---|
| ≥ 25 % = 2 of 8 | 0.0572 — 10× looser than de Vries |
| ≥ 37.5 % = 3 of 8 | 0.0058 — matched |
| ≥ 50 % = 4 of 8 | 0.00037 — 15× stricter |

So **≥ 37.5 % is the like-for-like comparison**, and comparing `is_responsive` rates at
face value makes V1DD look less responsive than it is.

## The von Mises fit

A two-peak von Mises model is fitted at each spatial frequency:

```
f(θ) = scale₁·exp(k₁·cos(θ − θ₀)) + scale₂·exp(k₂·cos(θ − θ₀ − 180°)) + b
```

Six parameters, all bounded below at zero. The two peaks are locked 180° apart, so the
model describes a direction-tuned neuron whose anti-preferred response can differ in both
height and width.

Trust-region-reflective least squares, retried at 2,000 then 10,000 evaluations. The
default starting guess is `(0.1, 1, 180, 0.01, 1, 0.001)` — fixed, not derived from the
curve.

### Why the fixed guess wins

An expected 2–3× speedup from deriving the start from the curve was measured at **2.1×
the model evaluations** (1,617 against 771 per curve) for a fit of equal quality. Two
further hypotheses about *why* were tested:

| start | x_scale | evals/curve | converged |
|---|---|---|---|
| fixed | default | **771** | 200/200 |
| derived | default | 1,617 | 199/200 |
| fixed | `"jac"` | 7,352 | 190/200 |
| derived | `"jac"` | 9,163 | 188/200 |

Jacobian scaling overflows because **`k` is unbounded above inside an exponential**: large
steps in `k` push `exp(k·cos(...))` to infinity, the residual goes non-finite, and the
trust region collapses. The fixed guess with absolute scaling keeps the search in a safe
region more or less by luck.

The derived guess reaches a lower sum-of-squares on 59 % of curves and a higher one on
33 %, median ratio 0.9993 — a rough objective where no single starting point dominates.
Available as `vonmises_data_p0=True`, off by default.

### Preferred direction from the fit

`vonmises_pref_dir` picks whichever of the two peaks is taller, with the baseline
subtracted: `f(θ₀) − b` vs `f(θ₀ + 180) − b`. This is the direction `ssi_tuning_fit`
evaluates at.

### What `fit_all_sf` controls

By default, the windowed fit runs only at the preferred SF (which is all surround
suppression reads), and the full-field fit runs at every SF. Setting `fit_all_sf=True` fits
every SF for both, which is ~2× slower but fills the `tuning_curves.npz` export
completely. Setting it `False` leaves the unread SF as NaN in the export, which looks like a
failed fit unless you know better.

## How the pipeline runs it

`drifting_gratings_metrics` in `families/drifting_gratings.py`. Called once per plane, once
for full-field and once for windowed. Returns a `DGResult` dataclass carrying the metrics
frame plus the intermediates surround suppression and roi_quality need: trial responses,
running speeds, preferred-condition indices, and fitted parameters.

The windowed variant also records the grating-aperture centre via `window_center`, which
reads the first non-NaN azimuth/elevation from the non-blank sweeps. Two sessions record no
centre at all — see [surround_suppression.md](surround_suppression.md).

Per-trial running speed is computed by prefix-sum over the running trace, padded by
`running_pad_seconds` (0.10 s) on each side of the sweep.

## Columns

| column | meaning |
|---|---|
| `dsi` | direction selectivity index, `(pref − null) / (pref + null)` |
| `frac_responsive_trials` | fraction of preferred-condition trials beating the spontaneous null |
| `gosi` | global orientation selectivity, vector-sum magnitude |
| `is_responsive` | `frac >= 0.50` and `pika_roi_confidence > 0.5` |
| `lifetime_sparseness` | Vinje & Gallant lifetime sparseness over condition means |
| `osi` | orientation selectivity index, `(pref − orth) / (pref + orth)` |
| `preferred_dir` | preferred direction in degrees (0–330, 30° steps) |
| `preferred_sf` | preferred spatial frequency |
| `pref_dir_mean` | circular-mean direction weighted by response |

Full-field and windowed share this column set. The prefix in the wide table
(`dgf_`/`dgw_`) is added by the orchestrator.

Trial-level data is in `tuning_curves.npz`: `dgf_trials` / `dgw_trials`
`(n_rois, 12, 2, n_trials)`, per-trial running speeds `dgf_running` / `dgw_running`
`(12, 2, n_trials)`, fitted parameters `dgf_params` / `dgw_params`
`(n_rois, 2, 6)`, blank-sweep responses, direction and SF lists.

## Sanity checks worth running on a fresh asset

* 12 directions at 30° steps, 2 spatial frequencies. Any other count raises.
* `frac_responsive_trials × n_trials_at_pref` is an integer for 100 % of ROIs.
* `preferred_dir` is one of the 12 stimulus directions (or NaN).
* `is_responsive` is False for every ROI where `pika_roi_confidence <= 0.5`.
* `lifetime_sparseness` medians are in the 0.20–0.35 range (condition-mean convention);
  values above 0.60 suggest the trial-flattened convention — see
  [comparability.md](../comparability.md).
