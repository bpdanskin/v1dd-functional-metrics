---
name: two-anomalous-sessions
description: "Column 4/volume 1 and column 2/volume 5 are exceptional in every way we have measured -- mask format, storage format, aperture centre, filtering date. Col4/vol1's plane 0 holds every low-confidence ROI in the asset. Provenance to chase outside the pipeline."
metadata:
  node_type: memory
  type: project
  modified: 2026-09-05
---

Two sessions account for nearly every exception this project has hit. They are the same
two every time, which is why they are worth chasing as one question rather than five.
**None of this is something the pipeline can settle** -- it is a question for whoever
produced the NWB asset.

| | col 4 / vol 1 | col 2 / vol 5 |
|---|---|---|
| session | `409828_2018-11-21_09-22-23` | `409828_2018-12-07_15-05-38` |
| filtering date in the reference run's name | `_filtered_2026-04-16_02-51-11` | `_filtered_2026-08-18_23-10-02` |
| every other session | `_filtered_2026-04-09_*` | — |
| ROI mask column | **`image_mask`** | **`image_mask`** |
| every other session | `pixel_mask` | — |
| storage format | was HDF5, now Zarr | was HDF5, now Zarr |
| windowed-grating aperture centre | **not recorded** -> imputed | **not recorded** -> imputed |
| ROIs | 1,550 | 906 |

So the set {only two with `image_mask`} = {only two that were HDF5} = {only two with no
aperture centre} = {only two with an off-schedule filtering date}. That is not four
coincidences; it is one upstream difference with four visible consequences.

## Col 4 / volume 1: the problem is one plane, not the session

The record used to say "67 % of this session is low-confidence". That understates and
mislocates it. Measured on the 2026-09-03 asset:

| plane | ROIs | low-confidence | median confidence | depth |
|---|---|---|---|---|
| **0** | **1,038** | **1,038 (100 %)** | **0.084** | 50 um |
| 1 | 9 | 0 | 0.915 | 66 um |
| 2 | 34 | 0 | 0.916 | 82 um |
| 3 | 87 | 0 | 0.916 | 98 um |
| 4 | 162 | 0 | 0.916 | 114 um |
| 5 | 220 | 0 | 0.916 | 130 um |

* **Every one of the 1,038 low-confidence ROIs in the entire 39,407-ROI asset is in that
  single plane**, and that plane is 100 % low-confidence at a median confidence of 0.084.
* Its five sibling planes are *starved* -- 9 to 220 ROIs where a normal volume runs
  roughly 400-480 per plane (col 1 / vol 3: 409, 470, 483, 478, 438, 430).
* The counts ramp monotonically with depth, 9 -> 220, which is not how ROI counts behave.

That pattern -- one plane holding 4x the normal count at near-zero confidence while its
neighbours are empty -- reads like a **plane-assignment or segmentation failure that
swept unassigned ROIs into plane 0**, not like poor imaging. Col 2 / vol 5 shows nothing
of the sort: 157/162/187/185/129/86 per plane and **zero** low-confidence ROIs.

`image_mask` mask sizes also differ from each other and from everything else: median
non-zero pixels per ROI is **38** in col 4 / vol 1 plane 0 (min 27) against **267** in
col 2 / vol 5, where `pixel_mask` sessions sit at **196** (171-228). A 38-pixel soma mask
at 512x512 is small enough to question.

## Questions to take upstream

1. Why do these two sessions carry `image_mask` when 23 carry `pixel_mask`? Different
   segmentation tool, different version, or a re-export?
2. What happened in col 4 / vol 1 plane 0? Is it a plane-assignment bug, a rejected-ROI
   dump, or a real acquisition failure? **Should those 1,038 ROIs be in the asset at all?**
3. Why do these two sessions have off-schedule filtering dates -- and why is col 2 / vol 5
   filtered four months after everything else?
4. Why is the windowed-grating aperture position missing for exactly these two? The
   stimulus clearly ran (both have complete SSI and `preferred_dir`), so only the recorded
   position is absent.
5. Is the 38-pixel median mask in col 4 / vol 1 plane 0 the same quantity as the 196-pixel
   median elsewhere, or a different mask definition?

## One record-keeping discrepancy, unresolved

The reference run's own `provenance.json` (2026-09-03) lists all 25 sessions **with**
`_filtered_<date>` suffixes and all 25 as `zarr`. The currently mounted asset has **bare**
names (`409828_2018-11-21_09-22-23.nwb.zarr`). So by the reference run the HDF5-to-Zarr
conversion had already happened, but the suffix was still present *then* and is absent
*now* -- meaning the mount has been re-published since, or the two changes were not
simultaneous. Worth confirming before treating the current mount as the asset the
reference numbers came from: **the ROI sets may not match the 39,407 of the reference.**

Related: [[metric-cautions]], [[repo-identity-and-the-frozen-fork]].
