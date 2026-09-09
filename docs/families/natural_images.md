# Natural images

Image selectivity from two natural-image stimulus sets: the 118-image set and the 12-image
repeat set. Both are handled by one function and share a column layout; they differ in
interpretation and in how `frac_responsive_trials` is computed.

## Goal

Characterise each neuron's preferred natural image, how strongly it responds to it, and how
selective it is across the image set. Two sets run in the same sessions:

* **natural_images** — 118 images, ~8 trials each. The broad survey.
* **natural_images_12** — 12 images drawn from the same 118-image catalog, ~40 trials each.
  Enough repeats for reliable single-cell estimates.

## The response window

`ni_response_frames = 2` — a fixed count of imaging samples after each onset. The original
took the window from an NWB `duration_sec` attribute the current files no longer carry, so
it was **recovered empirically**: scanning against the published table gives a sharp
optimum at 0.33 s (median |diff| in lifetime_sparseness = 1e-16, i.e. exact), while 0.30 s
gives 6e-3 and 0.35 s gives 2e-3.

The reason is discrete. At dt ≈ 0.165 s, imaging samples after an onset sit at
`delta, delta+dt, delta+2dt` with `delta ∈ [0, dt)`. A 0.33 s window is just under
`2·dt = 0.33008`, so it catches **exactly two samples on every trial**; a 0.30 s window
catches two when delta ≤ 0.135 and one otherwise. That varying count rescales each trial
differently, and `lifetime_sparseness` detects it because it is invariant to a *global*
scale but not a per-trial one.

**0.33 s does not generalise.** The margin is 8e-5 s, and `dt` spans 0.16123–0.16671 across
the 25 sessions. A time window catches exactly two samples only where `2·dt` lands just
above 0.33 — the two sessions the value was recovered on. Elsewhere it takes three samples
on up to 4.7 % of trials or one on up to 2.1 %. Expressing the intent as a frame count
removes the dependence on dt entirely.

## Conditions are `image_index`, not `image_order`

`image_order` is the raw presentation slot; `image_index` is the image's identity in the
118-image catalog. `natural_images_12` draws its twelve from that **same** namespace, so its
`pref_img` values are a sparse subset of 0–117 (e.g. 2, 4, 5, ..., 68) rather than 0–11.
Re-ranking them to 0–11 would look tidier and be wrong.

## Responsiveness

**`frac_responsive_trials` is a statistical test here**, unlike natural movie's
`mean(response > 0)`. It is the fraction of preferred-image trials whose response beats a
bootstrapped spontaneous null at p < 0.05 (10,000 single-draw bootstraps). This column
carries bootstrap noise and should be read against a seed control, not against zero.

**`z_score`** — how many standard deviations the preferred-image mean response sits above
the multi-trial spontaneous null. The multi-trial null averages `n_trials` draws per
bootstrap sample. For `natural_images_12` that is 10,000 × 40 = 400,000 window means — the
heaviest single call in the pipeline, and the reason `spontaneous_null` blocks by memory
budget.

## Lifetime sparseness

Computed over condition means (the Vinje & Gallant convention). **Not comparable across
stimulus sets of different size** — the `1 − 1/n` normaliser means the same neurons score
0.74 over 118 images and 0.40 over 12. Never put `ni_` and `ni12_` sparseness on one axis.

The historical column used the trial-flattened form (118 × 8 = 944 values per ROI), which
ran systematically high (median 0.9507 vs 0.7373) and was effectively uninformative about
image selectivity for NI12 (r = −0.005 against the condition-mean form). See
[comparability.md](../comparability.md).

## Reliability

Reported twice:

* **`reliability`** — split-half correlation between repeats, on the trace type the metrics
  use (deconvolved events). With sparse events this can be low even for strongly responsive
  cells.
* **`reliability_dff`** — the same correlation on dF/F, which carries a continuous signal.
  This is what the white paper's Figure 18 reports. Shipping both makes "how reproducible
  are the events every other metric is built on?" a question the asset can answer.

## `n_trials_at_pref`

The number of finite trials at each ROI's preferred image. This is the denominator needed
to attach a binomial tail p-value to `frac_responsive_trials`. `condition_means` is
trial-averaged by design, so the denominators are not otherwise recoverable from the
shipped asset.

## Condition means

The function returns `(metrics, (condition_means, image_ids))`. `condition_means` is
`(n_rois, n_images)` trial-averaged responses in float32, saved to `condition_means.npz`.
`image_ids` matters: for `natural_images_12` the columns are a sparse subset of 0–117,
not 0–11, and the array carries those ids so the images can be identified.

This is the only neuron-by-condition matrix for the natural stimuli in the asset. The
published columns are all reductions of it; keeping it is what makes population-level
analysis (e.g. representational similarity) possible at all.

## How the pipeline runs it

`natural_images_metrics` in `families/natural_images.py`. Called twice per plane: once
with `ns_type="natural_images"` and once with `ns_type="natural_images_12"`. The same
function, parameterised by stimulus type.

When a stimulus is absent, returns `absent_frame` plus `None` for the condition means.

## Columns

| column | meaning |
|---|---|
| `frac_responsive_trials` | fraction of preferred-image trials beating the spontaneous null |
| `lifetime_sparseness` | selectivity over condition means |
| `pref_img` | preferred image index (from the 118-image catalog), −1 if all-NaN |
| `pref_response` | mean response to the preferred image |
| `z_score` | preferred response in SD of the multi-trial spontaneous null |
| `reliability` | split-half correlation on events |
| `reliability_dff` | split-half correlation on dF/F |
| `n_trials_at_pref` | finite trials at the preferred image |

`natural_images` and `natural_images_12` share this column set. The prefix in the wide
table (`ni_`/`ni12_`) is added by the orchestrator.

## Sanity checks worth running on a fresh asset

* `pref_img` values are a subset of the stimulus's `image_index` values (0–117 for NI,
  a sparse twelve for NI12).
* `n_trials_at_pref` is a positive integer, ≤ the maximum trial count for the stimulus.
* `frac_responsive_trials × n_trials_at_pref` is close to an integer (bootstrap noise
  makes exact equality rare, unlike gratings).
* `lifetime_sparseness` medians are ~0.74 for NI and ~0.40 for NI12 (condition-mean
  convention); values above 0.90 suggest the trial-flattened convention.
