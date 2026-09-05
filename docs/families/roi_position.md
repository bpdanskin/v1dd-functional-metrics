# ROI position

Where each cell sits. Two very different levels of confidence live in this family, and
they are kept apart deliberately.

**Plane-local position is solid.** It comes straight from the segmentation footprint that
the NWB records per ROI, so `roi_x_px`, `roi_y_px`, `roi_area_px` and `roi_radius_px` are
measurements.

**The anatomical frame is a reconstruction.** The NWB records no usable geometry, so the
columns' positions relative to each other are *inferred* — twice, by two different routes,
and both are shipped so they can be compared rather than reconciled.

## Reading the footprint

The asset stores masks two ways, and the pipeline normalises both to pixel coordinates so
nothing downstream needs to know which it got:

| | sessions | form |
|---|---|---|
| `pixel_mask` | 23 | ragged list of `(x, y, weight)` per ROI |
| `image_mask` | 2 | dense `uint8` per-ROI image, `(n_rois, 512, 512)` |

Two things that could have been ambiguous are not, and both were checked rather than
assumed:

* **The axis convention is explicit.** `pixel_mask` entries are a structured dtype with
  *named* fields `('x', '<u4'), ('y', '<u4'), ('weight', '<f4')`, so there is no guessing
  which coordinate is which. Had they been an unstructured `(n, 3)` array, getting the
  order wrong would have transposed the field of view with nothing downstream noticing.
* **Weighted and unweighted centroids agree.** Every `pixel_mask` weight in the asset is
  exactly 1.0 (one distinct value), and `image_mask` is binary. So the unweighted mean over
  in-mask pixels is not a choice between definitions — it is the only definition available.

Masks are compact somas — 171 to 228 pixels, median 196 — so `roi_area_px` is meaningful
here. That is *not* true of receptive-field area, where the maps are fragmented; see
[receptive_fields.md](receptive_fields.md).

Footprints are read per ROI rather than as a whole column, because a dense mask is
`(n_rois, 512, 512)` — far larger than the footprints it describes. `load_plane`
deliberately drops the mask, so this family reads it separately: seven of the eight other
families do not want it.

An ROI whose footprint is empty gets **NaN**, not `(0, 0)`. A cell with no surviving
pixels has no position, and the origin is a real location.

## The two anatomical frames

### What the file gives us: nothing

Every imaging plane reports `origin_coords` as `[0.0, 0.0, 0.0]` and `grid_spacing` as
`[1.0, 1.0]` in `"meters"`. Both are placeholders. So the NWB supplies neither a pixel
scale nor a shared origin, and `ImagingPlane.location` carries only depth as free text.

Depth is therefore the one axis that *is* recorded: `depth_um` follows
`50 + 96*(volume-1) + 16*plane`, spanning 50–514 µm.

### `*_published` — the white paper's nominal grid

The paper's imaging-strategy section states the design in words: the volume is
800 × 800 × ~800 µm, imaged as five columns of five stacked volumes, and *"Four columns
create a grid covering 800 × 800 µm, while the 5th column captured the center."* Each
volume is 6 planes at 16 µm, spanning 80 µm — which matches our measured 96 µm volume
pitch once the gap is included, and the ~500 µm 2P depth matches our 50–514 µm range.

So the *geometry* is published: one column at the origin, four at (±200, ±200) µm.

**What is not published is which data column sits where.** That is only in Figures 2 and 3,
which are images — the pages carry no extractable text. The assignment used here is
therefore derived (below) and **remains provisional until someone reads those figures**.

### `*_retinotopic` — the layout the data implies

The paper says the windowed-grating aperture was positioned per column *"to align with the
population receptive fields of imaged neurons."* So `dgw_center_azimuth` /
`dgw_center_elevation` is a per-column retinotopic position, and because retinotopy in V1
is monotonic, relative aperture position stands in for relative cortical position.

Measured per column on the 2026-09-03 asset:

| column | azimuth | elevation | distance from the centroid of the others |
|---|---|---|---|
| **1** | −8.9 | −12.4 | **3.1°** |
| 2 | −19.6 | −10.0 | 16.8° |
| 3 | +1.8 | −9.7 | 10.9° |
| 4 | −15.4 | −16.4 | 12.2° |
| 5 | +9.9 | −14.4 | 20.6° |

**Column 1 is the centre**, by a clear margin. That is independently consistent with the
paper, which uses column 1 as its representative for the tiling analysis and calls the rest
"the surrounding overlapping columns."

### The comparison, and why the second axis does not survive it

`assign_columns` matches the four non-centre columns to quadrants by exhaustive search over
all 24 permutations, fitting one isotropic degrees-to-micrometres scale. On the real asset:

| | |
|---|---|
| isotropic fit | **17.6 µm/°**, residual **170.9 µm** against 200 µm offsets |
| azimuth axis | **15.0 µm/°** over a 29.5° spread |
| elevation axis | **67.2 µm/°** over a **6.7°** spread |
| anisotropy (elevation / azimuth) | **4.5** |

The azimuth number is physiologically plausible for mouse V1. The elevation number is not —
67 µm/° is simply what you get forcing a 400 µm cortical span onto a 6.7° retinotopic
spread. **The apertures separate the columns along one cortical axis and barely at all
along the other**, so the retinotopic route can identify the centre column and order the
others in azimuth, but it does not constrain the second grid dimension.

That is why the anisotropy and the per-axis scales are reported in provenance rather than
folded into a single residual: a large ratio is the signal that the second axis is
unconstrained, and it would otherwise hide inside a number that merely looks poor.

**Practical reading.** Trust `*_published` for the grid shape and `*_retinotopic` for the
azimuth ordering and the identity of the centre column. Where they disagree in the
elevation axis, the retinotopic estimate is the weaker one.

## The pixel scale is the least trustworthy number here

`um_per_pixel` defaults to `800 / 2 / 512 = 0.781`, on the reasoning that four columns
tiling an 800 µm grid makes one field of view ~400 µm across 512 pixels. But the paper's
"ongoing analysis" section says the columns **overlap** and cells are duplicated between
column 1 and its neighbours — which makes 400 µm an upper bound on the *spacing*, not a
measurement of the field of view.

Two consequences worth holding onto:

* **Every `*_um` column inherits that uncertainty**, so treat absolute micrometres as
  approximate. Set `um_per_pixel=None` and they come out NaN rather than misleading.
* **The comparison between the two frames does not.** Both are scaled by the same factor,
  so their disagreement is unaffected by getting it wrong.

Because the columns overlap, some neurons appear in more than one column. Nothing in this
family deduplicates them.

## Columns

| column | meaning |
|---|---|
| `roi_x_px`, `roi_y_px` | centroid in field-of-view pixels; NaN if the footprint is empty |
| `roi_area_px` | in-mask pixel count |
| `roi_radius_px` | equivalent-disc radius, `sqrt(area/pi)` |
| `roi_x_um`, `roi_y_um` | plane-local, scaled by `um_per_pixel` |
| `roi_x_um_published`, `roi_y_um_published` | plus the paper's nominal column offset |
| `roi_x_um_retinotopic`, `roi_y_um_retinotopic` | plus the aperture-derived column offset |

The two offset tables, the fitted scales, the anisotropy and the quadrant assignment are
all recorded under `column_layout` in the asset's `provenance.json`.

## Sanity checks worth running on a fresh asset

* Centroids fall inside the field of view — `0 <= x, y < 512`.
* `roi_area_px` medians agree across sessions. They do **not** in the two `image_mask`
  sessions: median 38 pixels in column 4 / volume 1 against 267 in column 2 / volume 5 and
  196 elsewhere. See the note on [the two anomalous sessions](../../.claude/memory/two-anomalous-sessions.md).
* The two anatomical frames agree in azimuth ordering and disagree in elevation — if they
  ever agree in *both*, the fit has probably been over-constrained.
