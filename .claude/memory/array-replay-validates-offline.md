---
name: array-replay-validates-offline
description: "The NWB input is capsule-only, but the shipped asset's arrays plus its metrics table let you validate rewritten metric code offline."
metadata: { node_type: memory, type: project, modified: 2026-09-03 }
---

**There is no NWB data outside Code Ocean.** But the reference asset ships both the
trial-level arrays and the metrics computed from them, so a rewritten metric function can
be checked by replaying the arrays and comparing to the shipped columns — no capsule, no
mounted dataset. This is the main refactor gate; see [[repo-identity-and-the-frozen-fork]]
for the asset path.

What is replayable:

| array | shape | validates |
|---|---|---|
| `dg{w,f}_trials` | (39407, 12, 2, 8) | all drifting-gratings metrics, surround suppression |
| `dg{w,f}_running` | (150, 12, 2, 8) | running modulation, trial-level running correlation |
| `rf_maps` | (39407, 2, 8, 14), pre-threshold | receptive fields |
| `ni_mean` / `ni12_mean` | (39407, 118) / (39407, 12) | natural-image condition means |

## Compare against a noise floor, not a fixed tolerance

**The archives are float32.** A ratio metric like `dsi = (pref - null) / (pref + null)`
loses precision by cancellation in the numerator, so its relative error is about
`eps32 / |dsi|` — unbounded as the metric approaches zero. A flat rtol 1e-6 therefore
fails a tail of ROIs for reasons that have nothing to do with the code, and chasing that
tail is wasted work.

`replay_reference.py` measures the sensitivity instead of guessing: it recomputes every
metric from input perturbed by **half a float32 ULP** and reports agreement-with-reference
beside agreement-with-itself. A disagreement no larger than that floor is storage
precision. Measured on the 2026-09-03 asset — `dgw_pref_dir_mean` 156 ROIs over rtol
against a floor of 156, `dgw_dsi` 168 against 162, `gosi` 6 against 6 — and the
disagreeing ROIs have median `|dsi|` 0.0047 against 0.4153 for agreeing ones, which is the
cancellation signature and confirms the mechanism.

**Compare fitted tuning curves, never the parameters.** Several von Mises parameters are
bounded at zero and park there, so relative error is measured against ~1e-30 and is
meaningless; the curve is also the only thing `ssi_tuning_fit` reads. Even as curves the
fit is genuinely ill-conditioned — 53 % of curve points within 1e-6, 84 % within 1e-4,
97 % within 1e-2 — but the noise floor tracks it at every quantile, so that is a property
of a bounded 6-parameter least squares on 12 points, not of the port.

Do not use rtol 1e-9 as a headline either — at that scale you measure float summation
order, and it once reported 0.1 % agreement for a family whose `preferred_dir` was
exactly identical.

**Three relationships hold exactly and make good assertions:**

* thresholding `rf_maps` at `rf_frac_thresh` (0.25) reproduces `has_rf_on` / `has_rf_off`
  for **100 %** of ROIs;
* `frac_responsive_trials x n` is an integer for **100 %** of ROIs (residual exactly 0),
  where n is the finite trials **at that ROI's preferred condition** -- not the maximum
  across conditions, which is a different and wrong denominator;
* the wide table reproduces column-for-column.

**Not replayable, needs a capsule run:** anything from the continuous trace (dF/F, events),
natural-movie and LSN per-trial responses, and the per-ROI trial denominators for the
natural-stimulus families — `condition_means` is trial-averaged by design.

## Splitting a module by AST drops two things silently

Both bit during the P2 family split of `stimulus_metrics.py`, and both were caught only
because something else checked:

* **`ast.get_source_segment` starts at the `def`/`class` keyword, so decorators are
  excluded.** `@dataclass(frozen=True)` vanished from `MetricConfig`, `DGResult` and
  `WindowCenters`. Compute the start line as
  `min(node.lineno, *(d.lineno for d in node.decorator_list))`. pyflakes caught it as an
  "unused import" of `dataclass`.
* **Module-level subscript assignments have no `ast.Name` target**, so a
  "every top-level name is claimed" check cannot see them. Three
  `OUTPUT_COLUMNS["..."] = OUTPUT_COLUMNS["..."]` alias lines disappeared, leaving five
  of eight families unschedulable. A unit test caught it, not the splitter.

**Lesson: after any AST-driven split, run pyflakes and a test that touches every branch of
the moved code.** The splitter's own assertions checked what it knew to look for.

## Surround suppression replays too (added 2026-09-03)

All eight `ssi_*` columns need nothing the asset does not carry: `dg{w,f}_trials` for the
responses, `dg{w,f}_running` for the running/stationary split, `dg{w,f}_params` for the
fitted curve, and `pref_cond_index` recomputed from the trials. The replay walks planes,
as production does, because running speeds have no ROI axis and key on `plane_key`;
`DGResult` is a dataclass, so stand-ins are cheap to build. All eight reproduce within
their float32 floors, and `dgw_center_*` matches exactly.

**Perturb the input the metric actually reads.** The first attempt built the noise floor
by perturbing trials only, and `ssi_tuning_fit` came out with a floor of *exactly zero* --
because it reads `tuning_params` and never touches the trials. A zero floor against a
non-zero disagreement is the tell that the wrong input is being varied. With both archives
perturbed: 2,959 over tolerance against a floor of 2,993, p99 identical at 8.2e-6.

**Still not replayable:** `roi_quality` and `run_corr_dff` (need the continuous trace),
natural movie (no per-trial array ships), and `dgw_rf_*` (comes from `rf_maps`, which is a
separate exercise from the grating arrays).
