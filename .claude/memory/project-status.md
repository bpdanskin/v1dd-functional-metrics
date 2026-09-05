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

## The second capsule run, 2026-09-05 -- the pipeline completed

`code/run --session 1:3` on `68255af`. **All eight outputs landed**: the parquet, all
three npz archives, `provenance.json` and the three AIND sidecars. `failed_sessions` and
`failed_outputs` both empty, `complete_asset` correctly **false** (a session filter was
set), eight entries in `differs_from_reference_config`.

**`roi_position` works on real data.** 2,708 ROIs x 95 columns; centroids **100 % inside
0-511**, no NaNs, median area 177 px and median radius 7.5 px -- about **11.7 um
diameter** at the inferred pixel scale, which is a plausible soma. The two anatomical
frames are correctly all-NaN, because one column cannot constrain a layout;
`column_layout.fit` records the reason verbatim: *"need 5 columns with a centre, have 1"*.

### Two tests failed in the capsule and passed locally -- again

`tests.json` reported 68 passed, **2 failed**: `test_a_stamped_version_file_is_read` and
`test_it_refuses_in_a_git_less_copy_with_no_version`. **Both were test defects, not code
defects** -- the pipeline behaved correctly in each case:

* The entry point exports `V1DD_CODE_VERSION` before launching validation and passes the
  environment down, so in a capsule the first test ran with it already set. The env var
  correctly wins over the file, so the test was reading the environment, not the file.
* The second copies `code/` and expects a refusal, but a capsule may have a **stamped**
  `code/CODE_VERSION` -- which correctly supplies a version.

Both now control their own environment (`monkeypatch.delenv`, and blanking the copied
`CODE_VERSION`). Simulating both conditions locally then exposed two *more* tests making
the same assumption, which are fixed the same way: a blank `V1DD_CODE_VERSION` is now
asserted to mean *unset* rather than *malformed*, and the "ships comment-only" check skips
where the file has been stamped, because that is a deployment state rather than a defect.

**The suite is now verified in five shapes** -- clean checkout; env var set; file stamped;
both; and a `.git`-less copy. 70 tests, one skipping where it cannot apply.

This is the third recurrence of the class recorded in [[test-suite-shape]]. The lesson
holds and needs stating more strongly: **a test that reads process state must set that
state itself.**

### Runtime regressed, and it was my design choice

1306.7 s for one session against 900.6 s before `roi_position` existed -- +45 %, which
extrapolates to **~9.1 h** for 25 sessions against the fork's 5.1 h.

Cause: `load_roi_masks` read `pixel_mask` **per ROI**, about 450 store accesses per plane.
The justification for per-ROI reads -- that a dense mask is (n_rois, 512, 512) -- applies
only to `image_mask`. The ragged column is small enough to read whole.

Fixed: the flat array and its cumulative index are read once and sliced in memory, with
the per-ROI path kept as an automatic fallback whenever the flat form does not check out
(length, monotonicity, and ROI count are all verified first). A unit test drives both
paths through a fake ragged column and asserts they agree. **The speedup itself is
unverified against a real NWB** -- confirm on the next capsule run.

## The third capsule run, 2026-09-05 -- `--session 4:1`

The `image_mask` branch and the plane-0 anomaly session, on `d0aafd5`. All eight outputs
landed, 1,550 ROIs x 95 columns, `failed_outputs` empty, and the metadata stage completed
**despite** validation reporting a failure -- which is the intended non-fatal behaviour.

`aperture centres: 0 measured, 0 inferred, 1 unknown` is correct: this session records no
centre and, alone in its column, has no donor to impute from.

**`roi_position` works on the dense-mask path too** -- centroids 100 % inside 0-511, no
NaN -- and it produced a *new* line of evidence on the plane-0 anomaly, now recorded in
[[two-anomalous-sessions]].

### `processing.json` recorded one process instead of two -- again

The exact defect from the fork's first reproducible run, from a new cause that is mine.
`metadata._validation_summary` still expected the **retired validation notebook's**
layout -- `<validation_dir>/checks/validation.json` plus a `tests.json` keyed
`n_pass`/`n_fail` -- while `run_tests.py` writes `<validation_dir>/tests.json` keyed
`n_passed`/`n_failed`. Wrong directory, wrong key names, and a hard dependency on a file
that no longer exists. It fails **silently**: no error, just one process where there
should be two.

Worse, `test_metadata` *passed* throughout, because the ported fixture built the
notebook's layout. The fork's own memory warned about this precise thing: *"test_metadata
built validation artifacts inside the results dir -- a layout production never has -- so
it passed throughout."* The port carried the fiction across.

Fixed: the reader now reads what the writer writes, the fixture builds the real layout,
and a new test drives `run_tests.summarise` output **straight through**
`_validation_summary` so the two halves cannot drift again.

### A test broke on a Unix socket

`test_it_refuses_in_a_git_less_copy_with_no_version` did `copytree(REPO/"code")`, and a
live capsule leaves editor state there -- including `.vscode/code-server-ipc.sock`, which
`copytree` cannot copy. Now it copies only `run_pipeline.py` and `src/v1dd_metrics`,
which is what the test actually needs and is immune to whatever else accumulates.

**Fourth recurrence of the same class.** The rule stands: *a test must construct the shape
production has, and must not assume anything about the directory it runs in.*

### The mask read-speed fix is still unverified

398.5 s for 1,550 ROIs looks better than the previous 1306.7 s for 2,708, but **this
session takes the `image_mask` path, which was never optimised** -- the bulk read applies
to `pixel_mask` only. A `pixel_mask` session (e.g. `--session 1:3`) is still needed to
confirm it.

## What remains, in the order it should happen

1. **Re-run `--session 1:3`** (a `pixel_mask` session) to confirm three things at once:
   the suite is clean in the capsule, `processing.json` now records **two** processes, and
   the bulk mask read actually cuts the time. Everything else is known to work end to end.

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
