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

## What remains, in the order it should happen

1. **Run the pipeline end to end, once, on one session.** `pipeline.run()` has never
   executed -- only its helpers are unit-tested, and the four-stage wiring in
   `run_pipeline.py` (processing, validation, metadata) has never run together either.

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
