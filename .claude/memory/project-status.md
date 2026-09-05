---
name: project-status
description: "Where the refactor stands as of 2026-09-05: P0-P4 done plus roi_position, what is verified and how, and the ordered list of what remains."
metadata:
  node_type: memory
  type: project
  modified: 2026-09-05
---

Read this first to pick the work up. Detail lives in the notes it points at, and in
`docs/`.

## Done

| phase | state |
|---|---|
| **P0** skeleton and environment | done -- `code/src` package, `code/run` sets `PYTHONPATH`, version ladder, CCM dropped from the image, `datasets.json`, resources `small` |
| **P1** core ported | done -- five modules into the package, prose moved to `docs/pipeline.md` and `docs/data_access.md` |
| **P2** families and orchestrator | done -- nine output blocks across seven analysis families, `pipeline.py` replaces the 25-cell notebook, four metric corrections, two deferred columns added |
| **P3** output restructure | done, absorbed into P2 -- one wide parquet plus three npz; no per-family CSVs |
| **P4** tests to pytest | done -- 67 tests, clean in a checkout **and** in a `.git`-less copy |
| **roi_position** (added after P4) | built and unit-tested; **never run against real data** |

Wide table is **95 columns** (9 identity + 86 metrics).

## How it is verified

`python -m pytest` (67 tests) and
`python code/validation/replay_reference.py --asset <shipped asset>`, which replays
**34 of the 72 metric columns** the reference asset carries -- drifting gratings, all eight
`ssi_*`, the RF centres and the four `dgw_rf_*` -- against a float32 noise floor. See
[[array-replay-validates-offline]].

Not replayable: responsiveness and `roi_quality` (need the continuous trace), natural
movie (no per-trial array ships). Natural images is **partly** replayable from
`condition_means` and has not been done.

## The first capsule run, 2026-09-05

`code/run --session 1:3` on `f9ca782`. It did what a smoke run is for: **it found a bug
that all 67 tests missed.**

* Session discovery, the aperture pre-pass, the six per-plane family calls and the ROI
  count all worked. **2,708 ROIs across 6 planes -- exactly the reference value**, so the
  mounted asset matches what the replay gate compares against.
* **900.6 s for one session**, so a 25-session run is ~6.2 h against the fork's 5.1 h. The
  extra is mask reading plus `roi_position`. Budget for it before P6.
* `column layout: centre=None` is **correct** for a single-column run -- `assign_columns`
  needs five columns and declines rather than inventing one. It will populate on a full run.
* It then died in `build_wide` with `KeyError: 'roi_position'`.

**The bug.** `roi_position` was added to `FAMILIES`, to `OUTPUT_COLUMNS`, to `PREFIX` and
to the mask loading -- but the one line that actually *calls* the family in `process_plane`
was never applied. A `str.replace` in a batch edit silently did not match, and pyflakes
could not see it because `rpm` is genuinely used elsewhere in the module. So masks were
read for every plane and thrown away, and `Accumulator.tables()` -- which drops families
with no rows -- left no key for `build_wide` to find.

Fixed, plus three guards that would have caught it in seconds:

* `check_families_ran()` runs before `build_wide` and names the family that produced
  nothing, instead of a bare `KeyError` after 15 minutes.
* a test that `process_plane`'s source appends to every family in `FAMILIES`.
* a test that `FAMILIES`, `OUTPUT_COLUMNS` and `PREFIX` agree.

**The lesson, which has now bitten three times:** a silent `str.replace` in a batch edit is
this project's most reliable source of defects -- it also dropped `@dataclass` decorators
and three `OUTPUT_COLUMNS` alias lines during the P2 split. **Assert on every replacement,
or use an editor that fails loudly.**

## What remains, in the order it should happen

1. **Re-run the one-session smoke test**, now that the `roi_position` call exists. The
   first attempt (above) died before the wide table, so nothing downstream of it has ever
   run: the three array writers, provenance, the AIND sidecars, and the validation and
   metadata stages of `run_pipeline.py` are all still unexercised.

   ```bash
   code/run --session 1:3          # ~12 min against ~5 h for the full asset
   ```

   That exercises session discovery, the aperture pre-pass, all nine family calls, the
   wide-table merge, the three array writers, provenance and the AIND sidecars. The run is
   recorded as partial, so it cannot be mistaken for a complete asset. **Do this before
   anything else**: every prior defect in this project's history passed its tests and
   failed in the capsule, and every remaining item is cheaper to judge afterwards.
2. **Exercise `roi_position` on real data** as part of that run. `load_roi_masks` has
   never touched a real NWB, and the two `image_mask` sessions take a different branch --
   see [[two-anomalous-sessions]]. Check centroids land inside 0..512 and that
   `roi_area_px` medians are sane.
3. **Confirm which asset is mounted.** The reference run's provenance lists `_filtered_`
   session names; the current mount has bare ones. If the mount has been re-published the
   ROI set may not match the reference's 39,407, which would make the replay gate compare
   against the wrong target.
4. **P5 documentation** -- six family pages remain (`roi_position` is written), plus the
   working-notes scrub recorded in [[metric-cautions]], example notebooks, and figures
   built from the local 2026-09-03 asset.
5. **P6 capsule run** -- the full asset, against the checklist in
   [[array-replay-validates-offline]].

## Open questions that are not ours to answer

* The two anomalous sessions and the 1,038-ROI plane -- [[two-anomalous-sessions]].
* The true `um_per_pixel`; the NWB records a placeholder -- [[column-layout-is-one-dimensional]].
* Which data column occupies which quadrant, which is only in the white paper's Figures 2
  and 3 and needs reading by eye.

## Optional, not blocking

* Extend the replay to natural images from `condition_means`.
* The von Mises fit is ill-conditioned; bounding `k` is the untried idea --
  [[vonmises-fit-conditioning]].
* The RF threshold sits on an achievable value -- [[rf-threshold-is-a-knife-edge]].
