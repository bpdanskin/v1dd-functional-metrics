---
name: queued-events-vs-dff-study
description: "Queued for after P5/P6: a systematic events-vs-dF/F comparison per metric, and quality checks on L0 event inference at 6 Hz. Includes a discovered inversion between our trace choice and the white paper's stated methods."
metadata:
  node_type: memory
  type: project
  modified: 2026-09-05
---

Requested 2026-09-05, **to start after P5 and P6 have produced a working pipeline and a
full asset.** Two parts: which metrics should be computed from deconvolved events and
which from dF/F, and quality checks on the event inference itself, since V1DD images each
plane at **6 Hz** against the 30 Hz of the comparable Brain Observatory work.

## Start here: our trace choice is inverted relative to the white paper

The V1DD white paper's methods section says, verbatim:

> Except the receptive field mapping, all the other analysis is performed using the df/f
> traces. For the receptive field mapping, events detected by L0-penalized algorithm are
> used to achieve a more accurate estimation of receptive field.

Our `MetricConfig.trace_type` is the **exact opposite**:

| family | ours | the paper |
|---|---|---|
| drifting gratings, full and windowed | `events` | dF/F |
| natural images, images 12, natural movie | `events` | dF/F |
| locally sparse noise (receptive fields) | `dff` | events |

This was inherited from `allen_v1dd`, not chosen here -- the fork recorded "metrics use
deconvolved events, except receptive fields, which use dF/F" as a *fact about the original
code*. So either the paper's methods text does not describe the code that produced its
figures, or the code diverged from it. **Settle this before anything else in the study:**
it decides whether the comparison is "which is better" or "we have been using the other
one all along".

## What we already have to work with

* **`trace_type` is per-family config**, so an A/B is a configuration change, not new code.
  One session all-events against one session all-dF/F is roughly 30 minutes and yields a
  paired comparison across every column.
* **One metric is already shipped both ways**: `reliability` (events) and `reliability_dff`
  (dF/F) exist for all three natural-stimulus families across 39,407 ROIs. Correlating them
  is the cheapest first quantitative handle and needs no new run.
* `roi_quality` already mixes deliberately: `snr` / `signal_power` / `noise_power` and
  `run_corr_dff` come from dF/F; spontaneous rates and `run_mod_*` from events.

## The argument that constrains the answer

Several metrics are **not safe on a signed trace**, and this is already documented rather
than speculative:

* **Any ratio index.** `(a - b) / (a + b)` on signed dF/F is unbounded when the two nearly
  cancel (+0.050 against -0.049 gives 99) and inverts sign when both are negative, so a
  suppressed cell that is less suppressed while running scores negative. This is why
  `run_mod_*` is events-only. The same argument applies to `osi`, `dsi` and every `ssi_*`.
* **`osi` and `gosi` do not rectify the tuning curve.** That is only safe because events
  are non-negative; on dF/F a negative response at the orthogonal direction drives `osi`
  outside [0, 1] and a near-zero signed normaliser makes `gosi` explode. Implementations
  working from dF/F clip at zero first.
* **Lifetime sparseness** assumes non-negative responses; the Vinje & Gallant form is not
  defined for signed data.

So the honest hypothesis is **not** "events or dF/F everywhere" but: ratio and sparseness
metrics need a non-negative trace, while correlations, reliability and SNR do not and may
be better on dF/F. Test that rather than assume it.

## Why 6 Hz makes this sharper than it was for de Vries

Each plane is imaged at **6 Hz** (37 fps across 6 simultaneous planes); the Brain
Observatory ran at 30 Hz -- a 5x difference, which is what the fork's note about "the 5x
rate difference" refers to. A GCaMP6 transient spans only a few samples at 6 Hz, so L0
inference has far less evidence per event and its output depends more strongly on the
assumed decay and the sparsity penalty. Two consequences already visible in the asset:

* the natural-images response window is **2 imaging samples**, so a trial is a mean of two
  numbers;
* the fraction of images with an *exactly zero* response is high -- median 0.63 per ROI --
  which is what makes ratio metrics fragile and sparseness large.

## Quality checks worth building

1. **Event rate and amplitude against ground truth.** Huang et al. 2021 (below) recorded
   spikes and GCaMP6 simultaneously; use it to ask what event rate is plausible, and
   whether ours is depressed by the slower sampling.
2. **Sensitivity to the L0 parameters.** dF/F ships in the NWB, so events can be
   *recomputed* at several sparsity penalties and decay constants and the downstream
   metrics compared. This is the strongest available check and needs no new data.
3. **Zero-inflation per family** -- the fraction of trials, and of ROIs, whose response is
   exactly zero. Report it per family; it predicts where ratio metrics misbehave.
4. **Transient support** -- how many samples a typical detected event spans at 6 Hz,
   against the decay implied by the model. If it is 2-3 samples, event *timing* carries
   little information and only rate and amplitude should be used.
5. **Paired metric agreement** -- for every column, events against dF/F: correlation, rank
   correlation, and the fraction that change sign. Sign changes are the diagnostic that a
   metric is unsafe on a signed trace.

## References, and what each is for

All in `OneDrive - Allen Institute/swdb_material_2026/v1dd_context/`:

* **`Jewell_Witten_2018.pdf`** -- *Exact spike train inference via l0 optimization*, Annals
  of Applied Statistics. **The method that produced our events.** Read for the sparsity
  penalty and the assumed AR decay, both of which are sampling-rate dependent.
* **`elife-51675-v3.pdf`** -- Huang et al. 2021, *Relationship between simultaneously
  recorded spiking activity and fluorescence signal in GCaMP6 transgenic mice*. **Ground
  truth**: simultaneous spikes and fluorescence, so it is the reference for how well any
  inference recovers real spiking and how that degrades with sampling.
* **`de Vries Lecoq Buice 2019.pdf`** and its **supplement** -- the 30 Hz comparison
  dataset and the source of the agreement figures in docs/comparability.md.
* **`V1DD_WhitePaper_v6.pdf`** -- the methods text quoted above, and the imaging rates.

## Deliverable

A `docs/events_vs_dff.md` giving, per metric family, which trace it should use and why,
with the sign-change and zero-inflation evidence behind each call -- plus whatever quality
checks earn a place in the pipeline or in `code/validation/`.

Related: [[metric-cautions]], [[fresh-start-metric-changes]], [[project-status]].
