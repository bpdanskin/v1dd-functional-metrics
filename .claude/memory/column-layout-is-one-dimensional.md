---
name: column-layout-is-one-dimensional
description: "The aperture retinotopy identifies the centre column (column 1) and orders columns in azimuth, but does not constrain the second grid axis. Two anatomical frames ship side by side."
metadata:
  node_type: memory
  type: project
  modified: 2026-09-05
---

The NWB records **no usable spatial geometry**: `origin_coords` is `[0,0,0]` and
`grid_spacing` is `[1.0, 1.0]` `"meters"` on every plane. So ROI position in a shared
anatomical frame has to be reconstructed, and the `roi_position` family ships **two
reconstructions rather than one answer** -- see docs/families/roi_position.md.

## What the white paper gives, and what it withholds

The imaging-strategy section states the geometry in words: 800 x 800 x ~800 um volume,
five columns of five stacked volumes, and *"Four columns create a grid covering
800 X 800 um, while the 5th column captured the center."* Two independent cross-checks
pass: 6 planes x 16 um = 80 um volume span matches our 96 um pitch once the gap is
included, and "~500 um" of 2P depth matches our 50-514 um.

**The column-index-to-quadrant mapping is only in Figures 2 and 3, which are images** --
those pages extract as 6 and 7 characters of text. Reading them by eye is still the only
way to confirm the assignment. Until then `*_published` is provisional.

## Column 1 is the centre, and that part is solid

Derived from the windowed-grating aperture positions, which the paper says were set per
column "to align with the population receptive fields of imaged neurons":

| column | azimuth | elevation | distance from the centroid of the others |
|---|---|---|---|
| **1** | -8.9 | -12.4 | **3.1 deg** |
| 2 | -19.6 | -10.0 | 16.8 |
| 3 | +1.8 | -9.7 | 10.9 |
| 4 | -15.4 | -16.4 | 12.2 |
| 5 | +9.9 | -14.4 | 20.6 |

Independently consistent with the paper using column 1 as its representative for the
tiling analysis and calling the rest "the surrounding overlapping columns".

## The second axis is not constrained -- do not pretend otherwise

Fitting the four non-centre columns to the published grid:

| | |
|---|---|
| isotropic | 17.6 um/deg, residual **170.9 um** against 200 um offsets |
| azimuth | **15.0 um/deg** over a 29.5 deg spread |
| elevation | **67.2 um/deg** over a **6.7 deg** spread |
| anisotropy | **4.5** |

15 um/deg is plausible for mouse V1. 67 is just what forcing a 400 um span onto 6.7 deg
produces. **The apertures separate the columns along one cortical axis and barely at all
along the other**, so the retinotopic route gives the centre column and the azimuth
ordering, and nothing trustworthy about the orthogonal axis. `assign_columns` reports the
per-axis scales and the anisotropy for exactly this reason -- a single residual would hide
it.

## The pixel scale is the weakest number, and it does not matter for the comparison

`um_per_pixel` defaults to `800/2/512 = 0.781`, inferred from four columns tiling 800 um.
The paper's ongoing-analysis section says the columns **overlap** and cells are duplicated
between column 1 and its neighbours, which makes 400 um an upper bound on the *spacing*
rather than a field-of-view measurement. So absolute micrometres are approximate --
but **both anatomical frames carry the same factor, so their disagreement is
scale-invariant**, which is what makes comparing them worthwhile at all.

EM coregistration exists for only a small fraction of ROIs, mostly in column 1, so it is a
later sanity check on this estimate rather than a way to derive it.

Related: [[two-anomalous-sessions]] (the two `image_mask` sessions), [[metric-cautions]].
