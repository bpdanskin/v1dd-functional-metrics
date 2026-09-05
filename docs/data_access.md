# Reading the NWB sessions

`v1dd_metrics/nwb.py` is the only module that touches `hdmf_zarr` or `pynwb.NWBHDF5IO`.
Everything downstream works on plain arrays, so when the storage layout changes again,
one file moves.

The V1DD analysis code this pipeline descends from read a private Isilon HDF5 tree through
`OPhysClient`/`OPhysSession`, indexing literal paths like
`processing/l0_events_plane3/DfOverF/l0_events`. That tree is gone; the same sessions are
published as NWB.

## Storage format: all Zarr today, both readers retained

**The current asset is entirely NWB-Zarr** (`*.nwb.zarr` directories). An earlier version
mixed those with plain HDF5 `*.nwb` files, and the reader still dispatches on the suffix
via `open_session()`, so both work. That is kept deliberately rather than simplified away:
it costs one branch, and a reader that silently handles only one format fails by returning
a *shorter session list* rather than an error.

Use `find_sessions()` rather than a bare glob for the same reason. Note `*.nwb` does not
match `*.nwb.zarr` (that name ends in `.zarr`), so the two patterns never overlap; a
session present in both formats yields one path, with `prefer` deciding.

### Two traps in session discovery, both of which hung a run

**Never `rglob` into a Zarr store.** A `.nwb.zarr` store holds thousands of chunk
directories. On a mounted asset the crawl is not slow but effectively unbounded: 25 stores
took over 30 seconds to yield the first five results and never finished, so the caller
hung with no error and no output. `_walk_sessions` uses `os.walk` and prunes at the store
boundary. Pruning costs nothing, because nothing inside a Zarr store is ever a session.

**The layout fallback must fire only when the expected layout yields nothing at all.**
Per-format `or` fallbacks look equivalent and are not: on an asset that is uniformly one
format — which is exactly what this one now is — the *other* format's shallow glob is
legitimately empty, so its `or` fires and starts that unbounded crawl. That is what
happened once the mount became all-Zarr.

## Three differences from the old client

* **Traces are transposed.** NWB stores `(n_frames, n_rois)`; the old `get_traces`
  returned `(n_rois, n_time)`. This module keeps NWB's orientation, because that is what
  `responses.prefix_sums` wants. Metric code transposes once, deliberately.
* **One stimulus table, not seven.** NWB concatenates every stimulus into
  `intervals['stimulus_table']` with a union of columns and NaN where a parameter does not
  apply. `stimulus_trials()` slices it back apart. Consequently a blank sweep **cannot** be
  detected as "any NaN in the row" — every drifting-gratings row is NaN in the image and
  frame columns. Pass the parameter columns that define a condition for that stimulus
  explicitly, or `()` for stimuli with no blank sweeps.
* **Running speed is already differentiated.** The old client read cumulative distance and
  took a central difference; NWB ships cm/s directly.

## Session identity

`session_mouse()` returns `(mouse_id, mouse_label)` — `("409828", "M409828")`. The two
forms are kept apart deliberately: the bare number goes inside `roi_unique_id`, the
M-prefixed form is what the `mouse` column and filenames carry. Carrying one string and
prefixing on demand is how you get `MM409828`.

Resolution never silently defaults: `nwb.subject.subject_id` first (present and agreeing
with the directory name in all 25 sessions of this asset), then the leading token of the
session directory name as insurance, then **raise**. A wrong mouse id silently corrupts
every ROI identifier in the output.

## Imaging depth

This asset puts depth in `ImagingPlane.location` as free text (`"50 um"`) rather than in
the structured `origin_coords`. Across the 5×5 grid the values form a lattice —
`50 + 96*(volume-1) + 16*plane`, six planes 16 µm apart within a volume, volumes 96 µm
apart, spanning 50–514 µm.

`plane_depth_um` returns None rather than raising, because no metric depends on depth and
one session with an unexpected string should not stop a run. The pattern is **anchored**:
a number and an optional unit and nothing else. An unanchored search looks equivalent and
is not — on a file writing something like `"VISp layer 2/3"` a loose pattern happily
returns a depth of 2 µm, and a wrong depth is worse than no depth because nothing
downstream can tell it was a parse of the wrong number.

## The locally-sparse-noise template

The one genuinely uncertain piece of the port. The original loaded an 8×14 grid from
`lsn_9deg_28degExclusion_jun_256.npy` and carried a commented-out block labelled
"Incorrect stimulus" pointing at the 16×28 tif — which is what NWB embeds. They are almost
certainly the same stimulus at 2× sampling, so `load_lsn_template` block-reduces and
**asserts every 2×2 block is uniform**. If that assertion fails, the two are not the same
stimulus and receptive fields are not reproducible from this asset; better to find out
there than to ship retinotopy on the wrong scale.

On the current asset `native_shape` is already `[8, 14]`, so the downsample risk has not
materialised. The function returns a dict (`images`, `frames`, `azimuths`, `altitudes`,
`native_shape`, `reduced`, `blocks_uniform`) so a caller can inspect a failure rather than
only catch an exception.

## What this module does not do

It reads what the pipeline needs and nothing else. Describing a whole session — which
tables exist, which columns, what is in them — belongs in `code/validation/`, because that
is a question you ask *about* the data rather than *of* it.
