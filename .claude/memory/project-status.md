---
name: project-status
description: "Where the refactor stands as of 2026-09-09: P0-P6 done, P5 notebooks and explorer complete. Full 25-session run completed (39,407 ROIs, 5.14 h). Optional: quarto site, comparative dF/F-vs-events analysis."
metadata:
  type: project
  modified: 2026-09-09
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
| **P4** tests to pytest | done -- 76 tests, clean in a checkout, in a `.git`-less copy, and in the capsule |
| **roi_position** (added after P4) | done -- verified in the capsule on both mask formats |

Wide table is **95 columns** (9 identity + 86 metrics).

## How it is verified

`python -m pytest` (76 tests) and
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

### Runtime: what looked like a 45 % regression was a measurement error

**Retracted.** The claim was that `roi_position` cost +45 % (900.6 s -> 1306.7 s) and
that the per-ROI `pixel_mask` read was to blame. **Both halves were wrong**, and the way
they were wrong is the lesson.

The pipeline prints two timers. `... 900.6s` at the end of a session line is
*that session*; `6 planes, 2708 ROIs, 20.5 min` is the whole run, and `wall_seconds` in
provenance is the same quantity. **900.6 s was compared against 1306.7 s, a session timer
against a run timer.** Run 1's comparable number is 1230 s. So:

| run | commit | col/vol | ROIs | wall_seconds |
|---|---|---|---|---|
| 1 | `f9ca782` | 1/3 | 2,708 | 1230 (crashed after this print) |
| 2 | `68255af` | 1/3 | 2,708 | 1306.7 |
| 3 | `d0aafd5` | 4/1 | 1,550 | 737.2 |
| 4 | `05a57c3` | 1/3 | 2,708 | 1324.5 |

Run 1 already read masks per ROI and threw them away, so runs 1 and 2 differ only by the
`roi_position` family itself: **+6 %, not +45 %**.

The second error compounded it. 25 x 22 min = 9.1 h assumed col 1 / vol 3 is typical. It
is the **second largest of the 25 sessions** (2,708 ROIs against a 1,576 mean); only
col 5 / vol 3 is bigger. Runtime is close to linear in ROIs -- 0.489 s/ROI at 2,708 and
0.476 s/ROI at 1,550 -- so the asset projects to **39,407 x ~0.49 s = 5.4 h**, against the
fork's 5.1 h. **There is no runtime problem to solve, and there never was.**

**The lesson:** two numbers printed by the same program are not the same measurement.
Provenance now records `stage_seconds` per family and `mask_reads` per plane, so the next
such question is answered from the asset instead of from arithmetic on log lines.

### The bulk mask read is real but bought nothing measurable

1324.5 s against run 2's 1306.7 s on the identical session: **+1.4 %, i.e. no change.**
Reading the ragged column whole is still the right shape (2 store accesses per plane
rather than ~450) and it is kept, but the 45 % it was meant to recover did not exist. Its
one genuine defect -- that the fallback to per-ROI reads is silent, so a run cannot say
which path it took -- is now fixed by recording `bulk_read` on `RoiMasks` and summarising
it in provenance.

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

## The fourth capsule run, 2026-09-05 -- `--session 1:3` again

`05a57c3`, the same `pixel_mask` session as runs 1 and 2, to close three open questions.
All three closed.

* **The suite is clean in the capsule: 73 passed, 0 failed, 0 skipped.** The five-shape
  verification held; no test read state it had not set.
* **`processing.json` records two processes.** The `_validation_summary` fix works
  against the real writer, and the parameters carry `unit_tests_passed: 73`.
* **The bulk mask read changed nothing measurable** -- see above.

`roi_position` reproduced run 2 exactly: 2,708 ROIs x 95 columns, 409/470/483/478/438/430
per plane, median area 173-185 px across all six planes, centroids 100 % inside 0-511, no
NaN. Both anatomical frames correctly all-NaN on a one-column run.

Two smaller defects found by reading the sidecar rather than the log:

* The Validation process had **`end_date_time: null`**. AIND permits it, but a process
  with no end reads as one that never finished. It now ends at start + the suite's own
  `seconds`.
* Its notes promised *"integrity checks over every row of the asset"*, which this repo
  does not run -- that was the retired notebook's job. The pipeline's row-level guards
  (`check_families_ran`, `check_roi_coverage`) run inside the metrics step and abort it.
  The notes now say what actually happens.

## P6 full run — completed 2026-09-09

`code/run` on `9f9af4a`, full 25-session run. **All checks passed.**

| quantity | expected | actual |
|---|---|---|
| sessions | 25 | 25 |
| planes | 150 | 150 |
| ROIs | 39,407 | 39,407 |
| columns (wide table) | 95 | 95 |
| complete_asset | true | true |
| failed_sessions | [] | [] |
| failed_outputs | [] | [] |
| wall_seconds | ~5.4 h | 5.14 h (18,522 s) |
| mask_reads: pixel_mask | 23 sessions bulk | 138 planes bulk, 0 per_roi_read |
| mask_reads: image_mask | 2 sessions | 12 planes (col2/vol5 + col4/vol1) |
| aperture imputation | 2 inferred | 2 inferred (col2/vol5, col4/vol1) |
| column_layout | populated | populated (um_per_degree=17.58) |

**`stage_seconds` — first measurement of family runtime breakdown:**

| stage | seconds | % of wall |
|---|---|---|
| drifting_gratings_windowed | 6,471 | 34.9% |
| drifting_gratings_full | 6,148 | 33.2% |
| load_traces | 1,332 | 7.2% |
| natural_movie | 360 | 1.9% |
| natural_images_12 | 292 | 1.6% |
| natural_images | 106 | 0.6% |
| roi_summary | 67 | 0.4% |
| receptive_fields | 65 | 0.4% |
| load_masks | 20 | 0.1% |
| surround_suppression | 9 | <0.1% |
| roi_position | 0.5 | <0.1% |

Gratings are **68 %** of wall time, not the ~96 % the fork claimed. The gap is load_traces
(7.2 %) and natural_movie/ni12 (3.5 %), which the fork did not measure separately.

Asset: `results/409828_V1DD_functional_metrics_2026-09-09_01-41-33/`

## P5 — completed 2026-09-09

All seven family doc pages, four infrastructure docs, and two example notebooks
delivered. See [[data-explorer-artifact]] for the interactive explorer.

- `docs/notebooks/asset_access.ipynb` — loading the parquet, NPZ archives, provenance
- `docs/notebooks/example_figures.ipynb` — population selectivity, tuning curves, RF maps,
  depth profiles, retinotopic gradient, image selectivity, running modulation
- `docs/v1dd_explorer.html` — 39,407-ROI 3D Plotly scatter (also published as artifact)

## What remains

The pipeline itself is **done and verified end to end** on both mask formats, both a
normal and an anomalous session, and in five environment shapes. P5 documentation is
complete.

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
