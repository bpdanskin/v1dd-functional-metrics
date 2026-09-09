# ROI quality

Signal quality, spontaneous activity, and locomotion modulation. Every column here
characterises the ROI itself rather than its response to a specific stimulus — baseline
properties that contextualize the per-stimulus metrics.

## Goal

Three questions per cell:

1. **How well is this ROI recorded?** Spectral SNR from the dF/F trace.
2. **What does it do with nothing on the screen?** Spontaneous firing rate.
3. **Does locomotion change its activity?** Running modulation across grating conditions
   and the spontaneous block, plus a continuous correlation.

## Spectral SNR

`spectral_snr` computes per-ROI SNR from the dF/F trace's power spectrum.

* **Signal band:** 0.1–1.5 Hz, where calcium transients live.
* **Noise band:** 2.0–2.1 Hz, where a well-isolated ROI's spectrum is close to flat.

The noise band's mean power per bin, scaled by the number of signal-band bins, estimates
the noise power the signal band would contain if it held nothing but noise. SNR is their
ratio — a power ratio, not decibels (the CCM notebook plots it on a log axis from 1 to
1e5).

The bands and scaling match `estimate_snr_white_noise_model` in the cell-cell correlations
notebook, which ships `snr_by_cell.feather`. So the two assets are directly comparable but
**not bit-identical**: that notebook interpolates every plane onto one reference timebase
before the FFT, while this runs per plane on its own timestamps.

**dF/F only.** On deconvolved events this measures nothing useful: deconvolution has already
removed the noise floor the reference band is meant to sample, so the denominator is
whatever numerical residue is left rather than a noise estimate.

## Spontaneous activity

Spontaneous rates are the mean activity over the spontaneous block (the inter-stimulus
interval where nothing is on the screen), on **deconvolved events**.

| column | what it measures |
|---|---|
| `spont_rate` | mean activity over the whole spontaneous block |
| `spont_rate_run` | mean activity during running frames |
| `spont_rate_stat` | mean activity during stationary frames |

These are gated separately: a session where the animal never ran still has `spont_rate`
and `spont_rate_stat`. 13 of 25 sessions are one-sided (all running or all stationary), so
requiring both states would discard them.

`spont_rate` is **not recoverable** from `spont_rate_run` and `spont_rate_stat`.
`spont_run_frac` is a fraction of *time* on the running trace's own ~59 Hz samples (the
white paper's definition), while the split classifies *imaging frames*; the two are close
but are not the weights that would recombine the state means.

`spont_rate` is the only per-ROI **baseline activity level** in the asset — how much a
neuron does with nothing on the screen. Useful beyond locomotion: as a normaliser for
evoked responses, and for spotting unusually silent or hyperactive cells.

**Why events, not dF/F.** dF/F is defined against a rolling baseline, so its mean over a
long stimulus-free block is ~0 **by construction**. It would not be unstable — it would be
uninformative.

## Locomotion modulation

`(R_run − R_stat) / (R_run + R_stat)`, the same form as every `ssi_*` column.

Three conditions, all on deconvolved events:

* **`run_mod_dgf`** — full-field gratings, pooled over all 192 trials.
* **`run_mod_dgw`** — windowed gratings, pooled over all trials.
* **`run_mod_spont`** — the spontaneous block, where nothing is on the screen. This is the
  control: locomotion modulates cortex whether or not there is a stimulus, so a grating
  index should be read against this rather than against zero.

**Pooled over all trials, not at the preferred condition.** That is forced by the data.
`ssi_running` needs ≥3 running trials at one condition out of 8, and only **6.3 % of ROIs
in 9 of 25 sessions** clear that bar — running is close to all-or-nothing per session, with
many sessions 100 % stationary and two 100 % running. Pooling over all grating trials makes
the same threshold easy wherever the animal ran at all.

### Why events only

A ratio index is not safe on a signed trace. With non-negative events, `R_run + R_stat` is
a sum of magnitudes and vanishes only for a silent cell. On signed dF/F the same
expression breaks in two ways:

* Near-cancelling responses of opposite sign give an unbounded index:
  `(+0.050 − (−0.049)) / (+0.050 + (−0.049))` = 99.
* When both responses are negative the sign inverts: a suppressed cell that is *less*
  suppressed while running scores negative.

The white paper's `C·(Rmax − Rmin)/Rmax` does not rescue the ratio — its denominator can
itself be negative (`max(−0.01, −0.05) = −0.01`), giving −4.0 for that case.

**What no denominator fixes:** when both responses are near zero the ratio is large and
meaningless. Gate on magnitude before trusting a value from a quiet cell — for the gratings,
the raw per-trial responses are in `tuning_curves.npz`; for the spontaneous block, gate on
`spont_rate_run` / `spont_rate_stat`.

## `run_corr_dff`

Pearson correlation of running speed against each ROI's continuous dF/F trace, whole
session. On **dF/F, not events**: a correlation has no denominator, so it does not inherit
the sign instability that made `run_mod_*` events-only. Needs **no state split**, so it is
finite in every session — unlike `run_mod_dgf` (all-NaN in 2 of 25) and `run_mod_dgw`
(4 of 25). De Vries report a median excitatory correlation of ~0.03.

Not recoverable post-hoc: no time series ships in the asset.

## `run_frac` and `spont_run_frac`

Fraction of time spent running at >1 cm/s, reported at full running-trace resolution
(~59 Hz). `run_frac` covers the whole session; `spont_run_frac` covers the spontaneous
block only. The two can differ substantially: a session where the animal ran only during
grating presentations has a moderate `run_frac` and zero `spont_run_frac`.

The white paper gates its locomotion analyses at a running fraction of 0.20.

## How the pipeline runs it

`roi_summary_metrics` in `families/roi_quality.py`. Called once per plane, after both
drifting-gratings calls. Takes both `DGResult` objects (for `run_mod_dgf`/`run_mod_dgw`),
the spontaneous block boundaries, and the running trace.

## Columns

| column | meaning |
|---|---|
| `snr` | dF/F spectral signal-to-noise ratio (power ratio, not dB) |
| `signal_power` | total power in the 0.1–1.5 Hz band |
| `noise_power` | estimated noise power in the signal band |
| `run_frac` | fraction of session spent running |
| `spont_run_frac` | fraction of spontaneous block spent running |
| `spont_rate` | mean event rate during the spontaneous block |
| `spont_rate_run` | mean event rate during running in the spontaneous block |
| `spont_rate_stat` | mean event rate during stationary in the spontaneous block |
| `run_mod_dgf` | running modulation, full-field gratings |
| `run_mod_dgw` | running modulation, windowed gratings |
| `run_mod_spont` | running modulation, spontaneous block |
| `run_corr_dff` | Pearson r of running speed against dF/F, whole session |

## Sanity checks worth running on a fresh asset

* `snr` is positive (power ratio); typical range 1–100,000 on a log axis.
* `run_frac` and `spont_run_frac` are in [0, 1].
* `spont_rate` is ≥ 0 (events are non-negative).
* `run_mod_*` is in [−1, 1] wherever it is finite.
* `run_corr_dff` is in [−1, 1].
* Sessions with `run_frac = 0` have NaN for all running-only columns
  (`spont_rate_run`, `run_mod_*`).
