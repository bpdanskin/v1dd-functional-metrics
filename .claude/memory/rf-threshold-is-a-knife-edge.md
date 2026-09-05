---
name: rf-threshold-is-a-knife-edge
description: "rf_frac_thresh = 0.25 sits exactly on an achievable value of a discrete map, so 17% of ROIs have a pixel on the boundary and their centres are unstable to the smallest representable change."
metadata:
  node_type: memory
  type: project
  modified: 2026-09-03
---

Found 2026-09-03 while building the receptive-field replay. **Not a defect in this code** --
the centres reproduce from the shipped maps *bit-exactly* -- but a stability property of
the metric that anyone using RF centres needs to know.

## The map is discrete, and the threshold lands on it

Each pixel of `rf_maps` is the fraction of that pixel's presentations producing a response
above the ROI's bootstrapped spontaneous 95th percentile. With ~44 presentations per pixel
that fraction takes **81 distinct values**, all multiples of 1/44 -- and
`rf_frac_thresh = 0.25` is exactly **11/44**.

Measured on the 2026-09-03 asset:

| | |
|---|---|
| non-zero pixels | 8,021,612 |
| pixels sitting at exactly 0.25 | **23,751** |
| ROIs with at least one pixel on the threshold | **6,679 of 39,407 (16.9 %)** |
| ON fields that are a single pixel | **3,491 of 7,068 (49 %)** |

The comparison is `frac < thresh`, so a pixel at exactly 0.25 is **kept**. Had it been
`<=`, all 23,751 would be dropped and 16.9 % of ROIs would change. That is a real hidden
sensitivity in a line that looks arbitrary.

## What that does to the centres

Perturbing the stored map by half a float32 ULP -- the smallest change the archive can
represent -- moves **~30 % of RF centres** by more than 1e-6 degrees, with a worst case of
**60 degrees in azimuth** and p99 of 12 degrees. Overlap fractions move by up to 0.68.

The mechanism is the two facts above together: a pixel on the boundary flips, and because
half of all fields are a single pixel, that flip can relocate or delete the entire field.
This sharpens the existing caution in [[metric-cautions]] -- "one marginal pixel moves the
centre ~4.6 degrees" -- by giving the reason and a much larger worst case.

## What to consider doing about it

Nothing yet; it is documented rather than changed, and the replay reproduces the shipped
values exactly either way. But if RF centres ever become load-bearing:

* **Move the threshold off the lattice** -- 0.24 or 0.26 rather than 0.25 -- so no
  achievable value sits on the boundary and the metric stops depending on `<` versus `<=`.
* **Weight the centroid by the map** instead of taking it unweighted over surviving pixels,
  which would make a marginal pixel contribute marginally rather than fully.
* Report the pixel count alongside the centre, so a single-pixel field is visible as such.

Related: [[metric-cautions]] on why RF *area* is unusable at this grid, and
[[array-replay-validates-offline]] for how the sensitivity was measured.
