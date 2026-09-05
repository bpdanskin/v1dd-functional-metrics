---
name: metric-cautions
description: "Cautions that must reach docs/ — sparseness scaling, RF area, the low-confidence session, aperture geometry, and the ROI id that collides."
metadata:
  node_type: memory
  type: project
  modified: 2026-09-03
---

Each of these was learned expensively and is invisible from the code. They belong in the
per-family pages under `docs/`, not in comments.

* **Sparseness is not comparable across stimulus sets of different size.** The `1 - 1/n`
  normaliser depends on `n`: the same neurons score **0.74 over 118 images and 0.40 over
  12**. Never put `ni_` and `ni12_` sparseness on one axis.
* **The RF centre sits on a threshold knife edge** for 17 % of ROIs -- see
  [[rf-threshold-is-a-knife-edge]]. The smallest representable change to the stored map
  moves ~30 % of centres, worst case 60 deg.
* **Receptive-field area is near-unusable at this grid.** Pixels are 9.3 deg, so a fit
  needs the largest connected component first — taking moments of the whole thresholded map
  gives an equivalent radius of 38 deg, wider than the monitor. Of 7,068 ON fields, 3,491
  are a single pixel, **95 % are fragmented**, and only **316 ROIs (0.8 %)** have a compact
  component of 6 or more pixels. V1DD never ran the 4.65 deg sparse noise de Vries had.
* **One session is 67 % low-confidence.** Column 4 / volume 1 has 1,038 of 1,550 ROIs at
  `pika_roi_confidence <= 0.5` — 1,038 of 39,407 overall, all in that one session. An
  independent derivation (`workshop2 - extended version.ipynb`) found the same 1,038 via
  `is_soma == False`. A population average over all ROIs is weighted by one session's
  segmentation quality. Why it is so is a question for whoever produced the filtered NWB
  asset, not something this pipeline can settle.
* **The windowed-grating aperture is per column, 30 deg diameter** (white paper: the radius
  of the window is 15 degrees), positioned per column to match population receptive fields.
  Column 2 / volume 2 sits 0.2 deg off its column — do not test for exact equality within a
  column. **Two sessions record no centre at all** (col 2 / vol 5, col 4 / vol 1).
* **`ssi` does not track RF-to-window distance** — r = -0.03, flat binned means from 0 to
  beyond 37 deg. So `rf_inside_window_on` is a **conservative** filter: cells passing it are
  well targeted, but cells failing it should **not** be discarded, because the distance test
  is not sensitive enough to justify it. Centre distance is the wrong measurement; overlap
  between the RF map and the window disc is the right one.
* **`pref_img` for natural movie is approximate by construction** — the response window
  spans several 1/30 s frames, so read it as "around here in the clip".
* **`roi_unique_id` collides ~2.9x**: 13,555 distinct strings for 39,407 rows. Join on
  `(column, volume, plane, roi)`, or use `roi_key`.
* **Depth is a lattice**: `50 + 96*(volume-1) + 16*plane`, spanning 50-514 um, 30 distinct
  depths. Useful as a sanity check on any session-identity change.

## Working notes still embedded in code (scrub in P5)

`stimulus_metrics.py` and its family modules still carry inline rationale from the fork --
context, decisions and cautions that belong in `docs/`, not in comments. `responses.py`,
`nwb.py` and `paths.py` were scrubbed in P1, with their content moved to
`docs/pipeline.md` and `docs/data_access.md`.

**In P5, once each family page exists, do a final pass over the family modules**: move
anything explanatory into the matching `docs/families/*.md`, drop resolved questions, and
leave only a docstring saying what the function does and how its arguments differ. Do not
delete a note without first checking its content has a home.
