# Surround suppression

How much the surround inhibits each neuron's response: the windowed grating against the
full-field one, split by locomotion state and estimation method. Plus whether the
windowed-grating aperture actually covered each cell's receptive field.

## Goal

Visual neurons in V1 often respond more strongly to a small stimulus covering their
receptive field than to a large stimulus extending into the surround. `ssi` quantifies
that: `(W − F) / (W + F)`, where W is the windowed (small-patch) response and F is the
full-field response. Positive means the surround suppresses; negative means it facilitates.

## The eight SSI variants

All share the same denominator convention. The **reference condition is always the windowed
stimulus's preferred (direction, SF)**; the full-field response is sampled at *that same*
condition, never at its own preferred one. ROIs whose preferred condition is −1 stay NaN.

| column | what W and F are |
|---|---|
| `ssi` | mean response at the preferred (direction, SF) |
| `ssi_avg` | mean over all conditions |
| `ssi_avg_at_pref_sf` | mean over all directions at the preferred SF |
| `ssi_running` | preferred (direction, SF), running trials only |
| `ssi_running_avg_at_pref_sf` | all directions at preferred SF, running trials only |
| `ssi_stationary` | preferred (direction, SF), stationary trials only |
| `ssi_stationary_avg_at_pref_sf` | all directions at preferred SF, stationary trials only |
| `ssi_tuning_fit` | preferred direction from the fitted von Mises curve |

Running and stationary trials split at exactly 1 cm/s with **strict inequalities on both
sides**: a trial at exactly 1.0 belongs to neither. `ssi_running` and `ssi_stationary`
additionally require at least 3 qualifying trials in *both* stimuli.

### `ssi_tuning_fit`

Instead of taking the mean response at the preferred direction, this evaluates the fitted
von Mises curve at the peak. Both the windowed and full-field fit parameters are needed, and
both must have converged. The preferred direction comes from the windowed fit's own peak
(`vonmises_pref_dir`), with the baseline subtracted.

A historical inconsistency: the peak was selected after subtracting the baseline, but the
index was evaluated *including* it. The baseline adds to both W and F but does not cancel
in a ratio — it shifts `(W − F) / (W + F)` toward zero. The default corrects this; the
historical behaviour is `ssi_tuning_fit_includes_baseline=True`. This moves **36,323 of
37,016 ROIs** with a fitted curve, median 0.3291 to 0.3564, r = 0.84.

## The windowed-grating aperture

The windowed stimulus is presented through a 30° diameter circular aperture (white paper:
"the radius of the window is 15 degrees"), positioned per cortical column to match the
population receptive field.

| column | volume | azimuth | elevation |
|---|---|---|---|
| 1 | all | −8.9° | −12.4° |
| 2 | 1–4 | −19.6° | −10.0° |
| 3 | all | +1.8° | −9.7° |
| 4 | 2–5 | −15.4° | −16.4° |
| 5 | all | +9.9° | −14.4° |

**Two sessions record no centre at all** — column 2 / volume 5 and column 4 / volume 1.
`probe_window_center.py` confirmed the `center_azimuth` / `center_elevation` columns are
absent from their stimulus tables entirely, not present-and-NaN. With
`impute_dgw_center=True` (the default), these are filled from the median of the other
sessions in the same cortical column, and the rows are flagged `dgw_center_inferred=True`.

Column 2 / volume 2 sits 0.2° off the rest of its column. A median tolerates that;
asserting exact equality within a column would fail on real data.

## Receptive-field containment

`ssi` compares a windowed grating response against a full-field one, which only means
"surround suppression" if the window covered the cell's receptive field. A cell whose RF
sat outside the aperture was barely stimulated, and its weak windowed response reads as
suppression when it was a targeting miss.

Two measures, deliberately, because they disagree about which cells to keep:

**`dgw_rf_distance_on` / `dgw_rf_distance_off`** — degrees from the RF centre to the
window centre. The conservative reading. It is also blunt: the centre is an unweighted
centroid on a 9.3° grid, so one marginal pixel moves it ~4.6°.

**`dgw_rf_overlap_on` / `dgw_rf_overlap_off`** — the fraction of the RF's mass falling
inside the aperture, in [0, 1]. More permissive and better behaved: it keeps cells whose
field overlaps the window even though the centroid does not. On this asset, overlap keeps
1,572 cells against 970 at a 0.05 distance cut.

The overlap is weighted by the **post-threshold** receptive-field map. The continuous
pre-threshold map is dominated by noise floor — its mean overlap is 0.086 against the 0.073
a uniform random map would give.

### Neither is a filter

On this asset, overlap correlates with `ssi` at r = +0.07 (n = 6,827) with a
non-monotonic profile. The targeting concern is directionally supported but weak. These
columns are reported so a consumer can judge; gating on them would discard most of the data
on thin evidence.

### How `_window_coverage` works

The aperture area within each stimulus pixel is computed by sub-sampling each pixel on an
8×8 grid (64 sub-cells). A pixel-centre-inside-the-disc test is not good enough: the
pixels are 9.3°, the aperture is 30°, so the disc spans ~3.2 pixels while covering ~8
pixels' worth of area. Whether a centre test counts 5 or 9 depends on alignment — a swing
of ~40 %. Sub-sampling recovers area within ~1 % of πr².

## How the pipeline runs it

`surround_suppression_metrics` in `families/surround_suppression.py`. Called once per
plane, after both drifting-gratings calls. Takes both `DGResult` objects and an optional
containment frame from `window_containment`.

Containment is passed in rather than computed here on purpose: it depends on the
receptive-field map and the aperture position, not on the SSI arithmetic. Computing it
inside would make this family depend on locally sparse noise, so a session missing that
stimulus would take surround suppression down with it. Omit it and the containment columns
are NaN.

## Columns

| column | meaning |
|---|---|
| `ssi` | surround suppression index at the preferred condition |
| `ssi_avg` | SSI averaged over all conditions |
| `ssi_avg_at_pref_sf` | SSI averaged over directions at preferred SF |
| `ssi_running` | SSI at preferred condition, running trials |
| `ssi_running_avg_at_pref_sf` | SSI over directions at preferred SF, running |
| `ssi_stationary` | SSI at preferred condition, stationary trials |
| `ssi_stationary_avg_at_pref_sf` | SSI over directions at preferred SF, stationary |
| `ssi_tuning_fit` | SSI from the fitted von Mises curve |
| `dgw_center_azimuth` | aperture centre, degrees |
| `dgw_center_elevation` | aperture centre, degrees |
| `dgw_center_inferred` | True if the centre was imputed from the column median |
| `dgw_rf_distance_on` | degrees from the ON-subfield centre to the aperture centre |
| `dgw_rf_distance_off` | degrees from the OFF-subfield centre to the aperture centre |
| `dgw_rf_overlap_on` | fraction of ON-subfield mass inside the aperture |
| `dgw_rf_overlap_off` | fraction of OFF-subfield mass inside the aperture |

## Sanity checks worth running on a fresh asset

* `ssi_running` and `ssi_stationary` are NaN wherever the session has fewer than 3 running
  or 3 stationary trials at the preferred condition.
* `dgw_center_inferred` is True for exactly the sessions that record no centre (two in
  this asset), and False everywhere else.
* `dgw_rf_overlap_*` is in [0, 1] wherever it is not NaN.
* `ssi_tuning_fit` is NaN wherever either fit failed (either DGResult has NaN params at the
  preferred SF).
