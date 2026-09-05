# How the pipeline is structured

Per-ROI functional metrics from the V1DD NWB sessions. One run reads every session in the
mounted asset, computes each analysis family for every imaging plane, and writes a single
per-ROI table plus the arrays that do not fit a table.

## A run, end to end

`code/run` is what Code Ocean executes on a Reproducible Run. It puts `code/src` on
`PYTHONPATH` and calls `code/run_pipeline.py`, which runs four stages in a fixed order:

1. **Version gate.** Resolve the commit this run is built from, or refuse to start.
2. **Processing.** Every session, every plane, every family; write the asset.
3. **Validation.** Integrity checks, writing artifacts to `/scratch`.
4. **Metadata.** AIND sidecars, written into the asset directory.

The order is not incidental. Metadata runs **last** because `processing.json` records that
validation ran and how it went; written earlier it can only ever describe a validation
that had not happened. Validation is **non-fatal** — a failed check should surface in its
artifacts rather than destroy an asset that took five hours to build.

The version gate runs **first** for the opposite reason: a reproducible run copies `code/`
without `.git`, so the commit cannot be derived inside the capsule, and an asset with no
version cannot be tied to the code that made it. Failing in one second beats discovering
it at the end. The ladder is `$V1DD_CODE_VERSION`, then `code/CODE_VERSION`, then
`git rev-parse`. The value must look like a git object name (7–40 hex): a previous asset
shipped `"version": "= 17cacea..."` because Docker's legacy `ENV <key> <value>` form takes
everything after the first space as the value, and a guard that only asked "is it set?"
let it through.

## Where things go

| | |
|---|---|
| `/data` | read-only mounts; the NWB input asset |
| `/results` | captured by Code Ocean as the output asset |
| `/scratch` | discarded; validation artifacts live here on purpose |

Checking the asset is not part of the asset. Metadata sidecars go **inside** the asset
directory beside the data; executed notebooks and logs stay above it, because they are a
record of the run rather than data.

`V1DD_OUTPUT_TARGET` defaults to `scratch` and only the entry point sets it to `results`,
so opening a notebook interactively cannot accidentally produce something that looks like
a captured asset.

## The response engine

Every metric reduces to one question — *what was this neuron's mean activity in a window
after a stimulus onset?* — asked a few million times. Asking it with a Python loop over
sweeps costs 0.3–1 ms a call, which is 10–30 s for a natural-movie plane and another 15 s
for a 10,000-draw bootstrap, per family, per plane.

`responses.py` answers it with a **prefix sum over time**. Given
`cs[k] = traces[:k].sum(axis=0)`, the mean over any window is `(cs[b] - cs[a]) / (b - a)`:
two gathers and a divide, independent of how many samples the window spans. The thing that
looks like it forces a loop — response windows landing on different numbers of imaging
frames, because onsets are not frame-aligned — is free, since the samples are never
materialised and `b - a` is just an integer vector. Prefix sums are accumulated in float64
regardless of input dtype; a float32 cumsum over 20k frames loses precision visible in the
sixth decimal of a response.

Two subtleties are load-bearing, and both are about matching the original exactly:

**Trial windows are label-closed; bootstrap windows are frame-indexed.** Trial responses
were originally selected with `xarray.sel(time=slice(...))`, which includes *both*
endpoints and so spans a variable, data-dependent number of samples. The bootstrap used a
different primitive — a half-open slice of fixed width `round(w / dt)`. At dt = 0.164 s a
2 s grating window gives 13 samples for a trial and 12 for a bootstrap draw. Both
primitives therefore exist: `window_bounds` for the first, `spontaneous_null` for the
second. A naive `(t >= lo) & (t < hi)` mask drops one sample per trial and shifts every
response slightly, which survives into a metric and reads as an algorithm bug.

The two forms disagree in a way no choice of window length reconciles. A varying sample
count rescales each trial differently, changing the *relative* pattern across conditions;
a fixed count does not. Scale-invariant metrics such as lifetime sparseness can tell them
apart even when the mean response looks similar.

**Trial arrays are NaN-padded.** Conditions are not all presented the same number of
times, so the trial axis is sized to the maximum and short conditions keep NaN in the
tail. Every reduction is nan-aware. A sweep's trial index is its chronological rank within
its condition, which matters because per-trial running speeds are indexed the same way and
the surround-suppression code masks one array with the other.

`responses.py` imports nothing but numpy — arrays in, arrays out — so it can be tested
against a synthetic trace with no data mounted.

## Validating without the data

The NWB input is only mounted in Code Ocean, but a shipped asset carries **both** the
trial-level arrays and the per-ROI metrics derived from them. So the array-to-metric half
of the pipeline can be checked offline:

```bash
python code/validation/replay_reference.py --asset <a shipped asset> --n-fits 200
```

This replays the shipped arrays through the current code and compares against the shipped
columns: seven drifting-gratings metrics for each grating type, and all eight `ssi_*`
columns, which need nothing the asset does not already carry — `dg{w,f}_trials` for the
responses, `dg{w,f}_running` for the state split, `dg{w,f}_params` for the fitted curve.
Surround suppression is walked plane by plane, as production does, because running speeds
have no ROI axis and key on the plane.

Receptive fields replay from `rf_maps`: `has_rf_on/off`, the four ON/OFF centres, and the
four `dgw_rf_*` containment columns, with aperture coverage computed once per distinct
centre. The two sessions recording no centre fall into their own group, where coverage is
undefined and every containment column must come out NaN — an edge case nothing else
exercises against real data.

Not replayable, and reported as skipped: anything needing the continuous trace —
responsiveness, spontaneous rates, `run_corr_dff` — because no time series ships in the
asset; and natural movie, which ships no per-trial array. Natural images are *partly*
replayable from `condition_means` and have not been done yet.

### The receptive-field threshold is a knife edge

The centres reproduce **bit-exactly**, but their noise floor is enormous: perturbing the
stored map by half a float32 ULP moves ~30 % of centres, worst case 60° of azimuth. The
reason is not float precision but the threshold. Each pixel is a fraction of ~44
presentations, so the map takes only 81 distinct values — and `rf_frac_thresh = 0.25` is
exactly 11/44. **23,751 pixels sit precisely on the boundary, across 16.9 % of ROIs**, and
since 49 % of ON fields are a single pixel, one flip relocates or deletes the field.
The comparison is `<`, so those pixels are kept; `<=` would drop all of them. See
[the memory note](../.claude/memory/rf-threshold-is-a-knife-edge.md).

### Why the gate uses a noise floor rather than a tolerance

The archives store trials as **float32**. A ratio metric such as
`dsi = (pref − null) / (pref + null)` loses precision by cancellation in the numerator, so
its relative error is roughly `eps32 / |dsi|` — unbounded as the metric approaches zero. A
flat `rtol = 1e-6` therefore fails a tail of ROIs for reasons that have nothing to do with
the code.

The harness measures that sensitivity instead of guessing at it: it recomputes every
metric from input perturbed by half a float32 ULP and reports agreement-with-reference
beside agreement-with-itself. A disagreement no larger than that floor is storage
precision. On the 2026-09-03 asset the two match closely — `dgw_pref_dir_mean` disagrees
on 156 ROIs against a floor of 156; `dgw_dsi` on 168 against 162 — and the disagreeing
ROIs have median `|dsi|` of 0.0047 against 0.4153 for the agreeing ones, which is the
cancellation signature.

**Perturb the input the metric actually reads.** `ssi_tuning_fit` is computed entirely
from `tuning_params` and never touches the trials, so a floor built by perturbing trials
alone came out exactly 0 — which is the tell that the wrong input was being varied. With
both archives perturbed it agrees closely: 2,959 ROIs over tolerance against a floor of
2,993, p99 identical at 8.2e-6.

**Fitted tuning curves are compared as curves, not parameters.** Several von Mises
parameters are bounded at zero and park there, so their relative error is measured against
~1e-30 and means nothing; the curve is also the only thing `ssi_tuning_fit` reads. Even as
curves the fit is genuinely ill-conditioned — only ~53 % of curve points agree to 1e-6,
84 % to 1e-4, 97 % to 1e-2 — but the noise floor tracks it at every quantile, so that
instability is a property of a bounded 6-parameter least squares on 12 points, not of the
port.

## Unit tests

```bash
python -m pytest                      # from the repo root
python code/validation/run_tests.py <dir>   # same, plus a tests.json summary
```

57 tests, 523 named checks. The JSON summary exists because the run happens where the
data is and the result is read somewhere else: `metadata.py` reads it to record that
checking happened, and by then the capsule log is thousands of lines long.

**The suite must be clean in two shapes** -- a checkout, and a copy with no `.git`, which
is what a reproducible run actually has. Three provenance defects in earlier runs passed
their tests and failed in the capsule because the tests ran in a shape production never
has. The one check that genuinely cannot apply without git -- that the version falls back
to `git rev-parse` -- is `skipif`-marked rather than left to fail. The fork reached the
same place the hard way: its equivalent check failed in the capsule and printed "unit
tests failed" over a clean asset, twice, because its harness could only skip whole files.

### What was not carried over from the fork's suite

| dropped | why |
|---|---|
| `test_reference_tables` | validation against the historical `data_frames` tables is [retired by decision](comparability.md#validation-against-the-historical-tables-is-retired) |
| `test_diff_runs`, `test_preflight` | `compare.py` and `preflight.py` were not carried into this repo |
| `test_v1dd_nwb` sections 5-7 | covered `schema_report` and `checkpoints`, also not carried over |

`test_entrypoint`, `test_import_boundary` and `test_tuning_export` were **rewritten rather
than ported**: the version ladder, the package layering and the array writers all changed
shape, so a line-by-line port would have tested the old design. They are now
`test_entrypoint`, `test_layering` and `test_array_writers`.

### Three relationships that hold exactly

Useful as assertions because any drift in them is unambiguous:

* thresholding `rf_maps` at `rf_frac_thresh` (0.25) reproduces `has_rf_on` / `has_rf_off`
  for **100 %** of ROIs;
* `frac_responsive_trials × n` is an integer for **100 %** of ROIs, where `n` is the
  number of finite trials **at that ROI's preferred condition** (5–8), not the maximum
  across conditions;
* `is_valid` is exactly `pika_roi_confidence > 0.5`, giving 1,038 invalid ROIs, all in
  column 4 / volume 1.

## Configuration

`MetricConfig` in `config.py` carries trace type per family, response windows, bootstrap
counts and thresholds. `REFERENCE_CONFIG` records the historical settings, and any
divergence is reported in the asset's provenance as `differs_from_reference_config`, so a
number that moved can always be traced to the setting that moved it.

Environment variables are listed in the [README](../README.md).

## Local development

```bash
pip install -e code/src[dev]
python -m pytest
python code/run_pipeline.py --check-env
```

`--check-env` resolves the version and paths and exits, which is the cheap way to confirm
a capsule is ready before committing to a multi-hour run.
