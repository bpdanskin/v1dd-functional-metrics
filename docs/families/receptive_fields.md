# Receptive fields

Where each neuron responds in visual space. ON and OFF subfield maps are built from the
locally-sparse-noise stimulus, thresholded, and reduced to a centre and a presence flag per
ROI.

## What makes this family different

Three things set it apart from every other family, all inherited from the original and all
deliberate:

* **dF/F, not deconvolved events**, and the only family with a *subtracted* baseline
  (the 1 s window before onset). Every other family uses events with no baseline.
* **No trial array.** Instead of grouping sweeps by condition and taking means, a design
  matrix records which stimulus pixels were bright (`pixel_on`) and which dark (`pixel_off`)
  on each sweep, and the map is the fraction of a pixel's presentations that produced a
  significant response.
* **No GLM.** The published README describes "a GLM framework"; there is no regression
  anywhere in this code. The design matrix is used purely as a counting indicator.

## How the map is built

1. **Sweep responses.** Each locally-sparse-noise frame onset is a sweep. The response is
   `mean(trace in [onset, onset + 4*dt])` minus `mean(trace in [onset - 1s, onset])`,
   where `dt` is the plane's imaging period. When the trace type is events, baseline
   subtraction is skipped (L0 output is already denoised and non-negative). The baseline
   subtraction matters for dF/F: without it, a neuron with a high sustained rate lights up
   everywhere.

2. **Spontaneous null.** 10,000 bootstrap draws from the spontaneous block, same window,
   same baseline subtraction. The 95th percentile of this distribution is the per-ROI
   threshold: a sweep counts as significant if its response exceeds it.

3. **Design matrix.** `(2 × n_pixels, n_sweeps)` boolean: the first `n_pixels` rows are ON
   (bright), the rest OFF (dark). Gray pixels are neither. The template is read from the
   NWB rather than hard-coded — this asset encodes stimuli as −1/0/1 where the original
   assumed 0/127/255, and hard-coding those would make both matrices all-False.

4. **Fraction map.** For each pixel, the fraction of its presentations that were
   significant: `n_significant / n_presentations`. This gives a continuous
   `(n_rois, 2, 8, 14)` array — ON and OFF subfields on an 8-row by 14-column grid.

5. **Threshold.** Pixel fractions below `rf_frac_thresh` (0.25) are zeroed. "Has a
   receptive field" reduces to "at least one pixel survived."

6. **Centre.** The **unweighted** centroid of the surviving pixel indices, mapped to
   degrees by interpolation into the stimulus grid's known positions. The post-threshold
   fractions are not used as weights.

## Response window smearing at 6 Hz

The response window is `lsn_response_frames` imaging samples (default 4). The
locally-sparse-noise stimulus presents a new pattern every ~250 ms. At 30 Hz, 4 frames =
133 ms — comfortably within one stimulus presentation. At 6 Hz, 4 frames = 660 ms —
**spanning 2.6 consecutive stimulus presentations**. Each sweep response therefore averages
neural activity driven by the target stimulus *and the next two patterns*.

Because consecutive LSN frames are spatially uncorrelated by design, the contamination acts
as noise: it raises the response floor at non-RF pixels and reduces contrast between true
RF pixels and background. This contributes to the 95 % fragmentation rate and the
prevalence of single-pixel fields.

| frames | window (6 Hz) | LSN presentations spanned | notes |
|---|---|---|---|
| 1 | ~165 ms | 0.66 | within one stimulus; may miss calcium peak |
| 2 | ~330 ms | 1.32 | captures most of GCaMP6s transient; ~30% next-stimulus contamination |
| 4 | ~660 ms | 2.64 | current default; substantial smearing |

The 4-frame default was inherited from the 30 Hz pipeline and never adjusted for 6 Hz.
Whether reducing it improves maps is under evaluation — see
`code/validation/rf_trace_comparison.py`, which compares detection rate and contiguity
across 1, 2, and 4 frames with both dF/F and events.

## The threshold is a knife edge

Each pixel fraction is a ratio of ~44 presentations, so the map takes only **81 distinct
values**, all multiples of 1/44. `rf_frac_thresh = 0.25` is exactly 11/44.

| | |
|---|---|
| non-zero pixels across the asset | 8,021,612 |
| pixels sitting at exactly 0.25 | **23,751** |
| ROIs with at least one pixel on the threshold | **6,679 of 39,407 (16.9 %)** |

The comparison is `<` (pixels at exactly 0.25 are kept). Had it been `<=`, all 23,751
would be dropped. Perturbing the stored map by half a float32 ULP — the smallest change
the archive can represent — moves **~30 % of RF centres**, worst case **60 degrees** in
azimuth. The mechanism: a pixel on the boundary flips, and since 49 % of ON fields are a
single pixel, that flip can relocate or delete the entire field.

The centres reproduce **bit-exactly** from the shipped maps. The instability is a property
of the threshold, not of the code.

## RF area is near-unusable at this grid

Pixels are 9.3 degrees across. Of 7,068 ON fields:

| | |
|---|---|
| single-pixel fields | 3,491 (49 %) |
| fragmented (non-contiguous) | 95 % |
| compact component of ≥ 6 pixels | **316 ROIs (0.8 %)** |

Taking moments of the whole thresholded map gives an equivalent radius of 38 degrees —
wider than the monitor. A meaningful area requires extracting the largest connected
component first, but even then almost no ROI has enough pixels. V1DD never ran the 4.65°
sparse noise that de Vries had. RF area is **not reported** and should not be derived from
these maps casually.

## The pixel-to-degree scale has a historical bug

`_rf_pixel_to_degrees` implements two mappings. The original divides the centre-to-centre
*range* by `n` rather than `n - 1`, compressing the scale by `(n−1)/n`: 12.5 % in
altitude (8 rows) and 7.1 % in azimuth (14 columns). The historical table therefore spans
±28.48° and ±56.13° where the screen spans ±32.55° and ±60.45°.

`rf_center_scale_bug=False` (the default) gives the correct mapping. The two differ by
exactly `n/(n−1)`, which makes the correction verifiable rather than merely asserted.

## How the pipeline runs it

`receptive_field_metrics` in `families/receptive_fields.py`. Called once per plane. Returns
a metrics frame **and** the continuous pre-threshold `rf_map` array, which is saved to
`receptive_field_maps.npz` alongside the stimulus grid coordinates and the bootstrap seed.

The pre-threshold map is what ships, not the post-threshold one: graded values before
zeroing are recoverable to post-threshold in one line (`map[map < 0.25] = 0`), but the
reverse is not.

Low-confidence ROIs (`pika_roi_confidence <= 0.5`) have their maps zeroed rather than left
populated. This is documented as "excluded", not "no RF" — a cell the segmentation
distrusts should not contribute a receptive field.

## Columns

| column | meaning |
|---|---|
| `has_rf_on` | at least one ON pixel survived the threshold |
| `has_rf_off` | at least one OFF pixel survived the threshold |
| `has_rf_on_or_off` | either subfield is present |
| `azimuth_rf_on` | ON subfield centre, degrees of visual angle |
| `altitude_rf_on` | ON subfield centre, degrees of visual angle |
| `azimuth_rf_off` | OFF subfield centre, degrees of visual angle |
| `altitude_rf_off` | OFF subfield centre, degrees of visual angle |

Centres are NaN where the corresponding `has_rf_*` is False.

The continuous map is in `receptive_field_maps.npz`, keyed `rf_maps` with shape
`(n_rois, 2, 8, 14)`, plus `altitudes`, `azimuths`, and `seed`.

## Sanity checks worth running on a fresh asset

* Thresholding the shipped `rf_maps` at 0.25 reproduces `has_rf_on` / `has_rf_off` for
  100 % of ROIs — an exact relationship.
* Centres fall within the stimulus grid: altitude within ±32.55°, azimuth within ±60.45°
  (or the compressed ±28.48°/±56.13° under the bug).
* `has_rf_on_or_off` is exactly `has_rf_on | has_rf_off`.
* Low-confidence ROIs have all-zero maps.
