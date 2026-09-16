# Receptive fields

Where each neuron responds in visual space. ON and OFF subfield maps are built from the
locally-sparse-noise (LSN) stimulus and reduced to a centre and a presence flag per ROI.

Two methods exist, selected by `config.rf_method`:

* **`"greedy"` (default)** — the greedy pixelwise RF (Millman, ported from `v1dd_physiology`):
  per-pixel stimulus-triggered average of events, a per-pixel bootstrap p-value, and a
  Holm-Šidák correction across all pixels. Two variants ship (strict + sensitive).
* **`"fraction"` (historical)** — the fraction-threshold method described below under
  *The historical fraction method*. This is what `REFERENCE_CONFIG` uses to reproduce the
  old tables and the array-replay validation.

The switch was made after a 25-session comparison: the fraction method fragments (only
**81%** of detected ON fields are a single connected component, mean area 3.3 px vs median
1 px — a long tail of scattered false pixels), while greedy holds **95–97%** single-component
at every depth and SNR, with comparable-or-better detection. See
[explorations/greedy_rf.md](../explorations/greedy_rf.md).

## How the greedy map is built

1. **Sweep responses.** Each LSN frame onset is a sweep; the response is
   `mean(events in [onset, onset + lsn_response_frames*dt])` (2 frames, ~328 ms at 6 Hz).
   Events, no baseline (L0 output is already denoised and non-negative).

2. **Design matrix.** `(2 × n_pixels, n_sweeps)` boolean: the first `n_pixels` rows mark
   ON (bright) pixels, the rest OFF (dark). The template is read from the NWB (this asset
   encodes stimuli as −1/0/1, not the original's 0/127/255).

3. **STA.** For each pixel, `STA = design · sweep_responses` — the summed event response
   over the sweeps where that pixel was active.

4. **Per-pixel bootstrap null.** Resample the sweeps with replacement `rf_greedy_n_boot`
   (5000) times; the p-value is the fraction of shuffled STAs ≥ the observed STA. This is
   the LSN stimulus's own null — no spontaneous block needed.

5. **Holm-Šidák correction.** Across all `2·n_pixels` tests per ROI, control the
   family-wise error rate. A pixel is significant if its corrected p < α. This multiple-
   comparisons control is what prunes the scattered false pixels the fraction method keeps.

6. **Two variants from one p-value array.** `rf_greedy_alpha_strict` (0.01) →
   the canonical columns (`has_rf_on`, …); `rf_greedy_alpha_sens` (0.05) → the parallel
   `_a05` columns, more sensitive at a modest false-positive cost. At 5000 bootstraps the
   sweep **aliases**: α = 0.04..0.01 are identical (only p=0 pixels survive below 0.05), so
   strict α=0.01 stands in for the whole ≤0.04 band. There are effectively two operating
   points, not a continuum.

7. **Centre.** Unweighted centroid of the significant pixels, mapped to degrees by
   interpolation into the stimulus grid (see the scale-bug note below).

α is held **global** (one value per variant) — the greedy tradeoff is stable across depth
and SNR, and a per-cell/per-depth α would make `has_rf` incomparable across cells.

## The historical fraction method (`rf_method="fraction"`)

The map is the fraction of a pixel's presentations that produced a significant response:
each sweep is significant if it exceeds the ROI's spontaneous 95th percentile (10,000
bootstrap draws from the spontaneous block, same window; dF/F historically used a 1 s
subtracted baseline). Pixels below `rf_frac_thresh` (0.25) are zeroed; "has a receptive
field" = at least one pixel survived. Kept because `REFERENCE_CONFIG` reproduces it and the
offline array-replay validates against it. Its two structural weaknesses — the knife-edge
threshold and unusable RF area — are why greedy is now the default:

**The threshold is a knife edge.** Each pixel fraction is a ratio of ~44 presentations, so
the map takes ~81 distinct values (multiples of 1/44); `rf_frac_thresh = 0.25` is 11/44.
~23,751 pixels sit at exactly 0.25 (16.9% of ROIs have one). The comparison is `<`, so a
half-ULP perturbation moves ~30% of RF centres. Greedy has no such threshold — significance
is a corrected p-value, not a fraction cut.

**RF area is near-unusable.** At 9.3° pixels, 95% of fraction ON fields are fragmented and
only 0.8% have a compact ≥6-pixel component; moments of the whole map give a 38° radius,
wider than the monitor. Greedy's MC correction collapses the fragmented tail (mean area
3.3 px → 1.2 px, 95–97% single-component), which is what finally makes RF area a reliable —
if coarse — metric. (The grid is still 9.3°; greedy fixes fragmentation, not resolution.)

## The pixel-to-degree scale has a historical bug

`_rf_pixel_to_degrees` implements two mappings; the original divides the centre-to-centre
range by `n` instead of `n − 1`, compressing the scale by `(n−1)/n` (12.5% altitude, 7.1%
azimuth). `rf_center_scale_bug=False` (the default) is correct; the two differ by exactly
`n/(n−1)`. Applies to both methods.

## How the pipeline runs it

`receptive_field_metrics` in `families/receptive_fields.py`, once per plane, returns a
metrics frame **and** a method-tagged arrays dict. `write_rf_maps` archives to
`receptive_field_maps.npz`:

* greedy: `rf_sta` (STA, float32), `rf_ge` (bootstrap shuffle-counts, uint16 — any α's mask
  is reproducible from these via Holm-Šidák), `strict_mask`, `n_boot`, `alpha_strict`,
  `alpha_sens`, plus `altitudes`, `azimuths`, `seed`.
* fraction: the pre-threshold `rf_maps`, plus the axes and seed.

Surround suppression consumes the **strict** variant (mask for aperture overlap, centres
for RF-to-aperture distance). Low-confidence ROIs (`pika_roi_confidence <= 0.5`) get empty
maps ("excluded", not "no RF").

## Columns

The canonical nine are the **strict** variant; the `_a05` nine are the **sensitive**
variant (identical meaning, α=0.05). Under `rf_method="fraction"` the `_a05` columns are
absent (bools False, centres and areas NaN).

| column | meaning |
|---|---|
| `has_rf_on` / `has_rf_off` | at least one significant ON / OFF pixel (strict) |
| `has_rf_on_or_off` | either subfield present (strict) |
| `azimuth_rf_on` / `altitude_rf_on` | ON subfield centre, degrees (strict) |
| `azimuth_rf_off` / `altitude_rf_off` | OFF subfield centre, degrees (strict) |
| `rf_on_area` / `rf_off_area` | subfield area, deg² (significant-pixel count × pixel area) |
| `*_a05` | the same nine for the sensitive (α=0.05) variant |

Centres and areas are NaN / absent where the corresponding `has_rf_*` is False.

## Sanity checks worth running on a fresh asset

* `has_rf_on_or_off` is exactly `has_rf_on | has_rf_off` (and likewise for `_a05`).
* Sensitive detects ≥ strict per ROI (`has_rf_*_a05 >= has_rf_*`).
* Centres fall within the stimulus grid (±32.55° altitude, ±60.45° azimuth; or the
  compressed ±28.48°/±56.13° under the scale bug).
* Greedy masks reproduce from `rf_ge` + the stored alphas via Holm-Šidák.
* Low-confidence ROIs have no RF.
