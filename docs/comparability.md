# Comparability

This pipeline is a fresh start. It owes no output compatibility to the 2025 `allen_v1dd`
pipeline, the 2019 V1DD white paper, or the assets this project shipped earlier. Where the
numbers line up with earlier work that is worth knowing, and where they deliberately do
not that is worth knowing too — but neither is enforced in code. Nothing here is a
validation gate; the gate is [`pipeline.md`](pipeline.md).

## What changed deliberately, and by how much

Each of these was reachable only once backward compatibility stopped being a requirement.
All are `MetricConfig` flags, so the historical behaviour is one argument away and
`provenance.json` records which one produced the asset.

### Lifetime sparseness now uses condition means

The historical column was computed over **every individual trial response**, flattened
across conditions and trials — 118 images × 8 trials = 944 numbers per ROI. Vinje &
Gallant, and de Vries et al. 2019 after them, take it over the **mean response to each
stimulus** — 118 numbers. Trial-to-trial variance adds spread without adding mean, so the
flattened form runs systematically higher.

| | trial-flattened | condition means | de Vries |
|---|---|---|---|
| natural images (118) | 0.9507 | **0.7373** | **0.77** |
| natural images 12 (12) | 0.9486 | 0.4023 | — |
| drifting gratings, windowed (24) | 0.6804 | **0.2596** | — |
| drifting gratings, full field (24) | 0.6874 | **0.2740** | — |

The condition-mean form is what agrees with published values. An earlier diagnosis blamed
the two-sample response window for the gap; that is **not supported** — the convention
explains it entirely.

Two consequences outlive the change. For `natural_images_12` the historical column was
close to uninformative about image selectivity: it correlates with the condition-mean form
at **r = −0.005**, because at 12 images the flattened form is dominated by trial noise.
And **sparseness is not comparable across stimulus sets of different size** — the `1 − 1/n`
normaliser means the same neurons score 0.74 over 118 images, 0.40 over 12, and 0.26 over
24 grating conditions. Never put two of those on one axis.

### A zero denominator is now NaN everywhere

`osi`/`dsi` returned 0 where the surround-suppression indices returned NaN — two
conventions for "undefined" in one asset. Both are NaN now.

**Measured, this changes nothing on the current data**: 0 of 39,407 ROIs move, r = 1.0000.
Invalid ROIs are already NaN, and no valid ROI is exactly 0.0 at both directions. It is a
guard against a case this dataset does not contain, not a correction.

### `ssi_tuning_fit` evaluates the curve the way it picks the peak

The preferred direction feeding it is chosen with the fitted baseline subtracted, but the
index was then evaluated *including* that baseline. The two are consistent now;
`ssi_tuning_fit_includes_baseline=True` restores the old behaviour.

This one moves nearly everything: median **0.3291 to 0.3564**, with **36,323 of 37,016**
ROIs changing and r = 0.8407. That is expected rather than alarming -- the baseline was
being added to both terms of `(W - F) / (W + F)`, which shifts a ratio instead of
cancelling out of it. The measurement is driven from the *shipped* fit parameters, so it
is exact arithmetic and none of the fit's own instability is in it.

### What did not change: the von Mises starting guess

An expected 2–3× speedup from deriving the fit's initial guess from the curve was
implemented and measured at **2.1× the model evaluations** (1617 against 771 per curve)
for a fit of equal quality — it reaches a lower sum-of-squares on 59 % of curves and a
higher one on 33 %, median ratio 0.9993. Two further hypotheses about *why* were tested
and both refuted, including `x_scale="jac"`, which is 9.5× worse and overflows.
The real conditioning problem is that `k` is unbounded above inside an exponential.
Available as `vonmises_data_p0`, off by default. Full write-up in
[the memory note](../.claude/memory/vonmises-fit-conditioning.md).

## Where this agrees with de Vries et al. 2019

The Brain Observatory is the closest published comparison, on a different dataset with the
same event-extraction algorithm.

* **Event magnitudes agree.** They report a median maximum evoked response of 0.006 AU
  against 0.0004 AU spontaneous; `pref_response` medians here are 0.0074 (natural images)
  and 0.0103 (natural movie), 0.0082 over responsive cells only. Same algorithm, same
  units, same order — the strongest cross-dataset agreement available.
* **Selectivity indices sit in their published range**: full-field medians over responsive
  cells are `dsi` 0.395, `osi` 0.635, `gosi` 0.293; windowed 0.470 / 0.659 / 0.310.
* **Lifetime sparseness agrees once computed their way** — see above.

### One name that means two different things

**Their "reliability" is our `frac_responsive_trials`** — "the percentage of responsive
trials to the cell's preferred stimulus condition", which is what their response-class
clustering is built on. Our `reliability` column is a *correlation between repeats* and
carries no threshold at all. Both are in this asset. Reaching for the wrong one silently
answers a different question.

### A threshold that does not transfer

A responsiveness *fraction* is not comparable across datasets, because the same fraction is
a different statistical test at a different trial count. Each trial has a 5 % chance of
passing by chance, so the false-positive rate is the binomial tail:

| | criterion | false-positive rate |
|---|---|---|
| de Vries 2019 | ≥ 25 % = 4 of 15 | 0.0055 |
| here | ≥ 25 % = 2 of 8 | 0.0572 — 10× looser |
| here | ≥ 37.5 % = 3 of 8 | 0.0058 — matched |
| here | ≥ 50 % = 4 of 8 | 0.00037 — 15× stricter |

So `dg_frac_thresh = 0.50` here is much stricter than their nominally lower 0.25, and
**≥ 0.375 is the like-for-like comparison**. Comparing our `is_responsive` rate against
their 25 % figure without that adjustment makes V1DD look less responsive than it is.

## Where this diverges from the 2019 white paper

* **Receptive-field area is not reported**, and should not be derived from the maps
  casually. At the 9.3° grid, of 7,068 ON fields 3,491 are a single pixel, 95 % are
  fragmented, and only **316 ROIs — 0.8 %** have a compact component of ≥ 6 pixels. V1DD
  never ran the 4.65° sparse noise de Vries had.
* **Surround suppression does not track receptive-field position.** `ssi` against
  RF-to-aperture distance gives r = −0.03, flat from 0° to beyond 37°. The white paper
  listed surround suppression under ongoing analysis and found the same targeting problem
  ("over half" of cells outside the window). `rf_inside_window_on` is therefore a
  *conservative* filter: cells passing it are well targeted, but cells failing it should
  **not** be discarded, because the distance test is not sensitive enough to justify it.
* **`run_mod_*` uses `(R_run − R_stat) / (R_run + R_stat)`**, not the paper's
  `C·(Rmax − Rmin)/Rmax`. The two are monotonically related, and this form matches every
  `ssi_*` column rather than putting a second convention in the same asset.

## An independent derivation agrees on the one thing it overlaps

A separate teaching notebook computed orientation tuning, OSI/gOSI and tuning width from
the same sessions, derived independently of this work and covering full-field gratings
only. It found `is_soma == False` for **1,038 ROIs, all in one session** — exactly the
1,038 ROIs at `pika_roi_confidence ≤ 0.5` in column 4 / volume 1. Two derivations reaching
the same number is the strongest evidence available that those criteria agree row for row.

## Validation against the historical tables is retired

Earlier work compared each family against the shipped `data_frames` reference tables and
against a second bootstrap seed. That comparison was never re-run after 2026-09-01, and it
is **retired by decision rather than left undone**: this pipeline does not owe those tables
its numbers, and several columns now deliberately differ from them.

What replaces it is [the replay gate](pipeline.md#validating-without-the-data), which
checks this code against a shipped asset's own trial arrays — a self-consistency check
rather than a fidelity claim about a pipeline we no longer follow.
