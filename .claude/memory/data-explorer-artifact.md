---
name: data-explorer-artifact
description: "V1DD Metrics Explorer (docs/v1dd_explorer.html) — 39,407-ROI 3D scatter + click-to-open per-cell side panel (tuning, polar, RF map, preferred NI image). Refreshed from the 2026-09-11 events/RF-window asset."
metadata:
  type: project
  modified: 2026-09-11
---

## Current artifact

URL: `https://claude.ai/code/artifact/abb43a1d-bb9f-419d-b163-53df50cc3c94` (last published pre-2026-09-11; local file is ahead — republish to update)
Local copy: `docs/v1dd_explorer.html` (~14 MB, self-contained except CDN Plotly; approaching the 16 MB artifact limit — the per-cell arrays in detail.json are the main mass, so trim there, e.g. drop RF maps for unused cells or downsample thumbnails, before adding more per-cell data)

Data source refreshed 2026-09-11 from `results/409828_V1DD_functional_metrics_2026-09-11_01-17-17`
(RF now events / 2-frame window, NI 3-frame). Anatomical `x_um`/`y_um` are reused verbatim
from the prior roi_data.json (row order is identical to the parquet on column/volume/plane/roi),
so the validated coregistration transform is NOT recomputed — only metric columns refresh.

### Architecture

Built by `scratchpad/build_explorer_data.py` (roi_data.json + detail.json) then
`scratchpad/assemble.py` concatenates: `part1.html` + `const D =` roi_data.json +
`const DETAIL =` detail.json + `part2.js`. Plotly.js 2.27.0 from cdnjs.
- `roi_data.json` — column-oriented scatter, 39,407 ROIs × 38 fields
- `detail.json` — per-cell panel data: `tuning_b64` (uint8, N×25 = 2sf×12dir mean + blank, normalised per cell, denorm via `tuning_peak`); `rf_b64` (uint8, only the 4,309 has_rf_on_or_off cells × 2×8×14; `rf_index` maps global→row or -1); 118 `ni_thumbs` (JPEG data-URIs, 150px tall) indexed by `ni_pref_img`.
- **Tuning array layout is SF-major** (`sf*12 + dir`), matching the JS reader `TUN[i*25 + sf*NDIR + d]`. The build (`build_explorer_data.py`) MUST `dgw_mean.transpose(0,2,1).reshape(...)` to produce this. A prior version stored it dir-major (`dir*2+sf`) while JS read sf-major, which scrambled the tuning/polar plots so the preferred-direction marker sat off the visible peak (peak==pref_dir only 9% of responsive cells; 99.8% after the fix). If the peak looks off-peak again, check this ordering first.

### Features

- **View toggle**: Anatomical / RF Space
- **Color metrics** grouped dropdown (added run_mod_dgw)
- **Column/volume/plane chip filters**, **threshold sliders** (SNR/OSI/DSI/Conf), wireframe boxes, hide-null
- **Confidence → marker size**; derived `dgw_pref_ori = dgw_pref_dir % 180`
- **Click-to-select cell → right side panel** (NEW 2026-09-11): quality stats, selectivity (OSI/DSI/gOSI/sparseness/SSI/pref dir+SF), running modulation, direction-tuning polar plot, tuning curve (both SFs + blank baseline), preferred natural-image thumbnail, ON/OFF RF heatmaps. Orange highlight ring marks the selected cell.
- **RF panel shows the thresholded field**: pixels < `rf_frac_thresh` (0.25) render as transparent background, suprathreshold pixels keep magnitude over [0.25,max] on a no-white colorscale. This kills the sub-threshold scatter (the continuous pre-threshold map is still what's stored in detail.json; thresholding is display-only in `renderRF`).
- **Panel tuning plots**: fixed per-SF colors shared across polar/direction-curve/orientation-curve (0.04 cpd `#3a96b5` blue, 0.08 cpd `#d68a2e` orange; pref dir/ori marker `#e0484d` red); preferred SF is bolded + labelled "(pref)". The folded-orientation plot is the **surround-suppression** view: windowed (SF color) vs full-field (`#8e5bb5` purple), both folded at the pref SF, captioned with SSI. **SSI facts** (`families/surround_suppression.py`): every variant is `(W−F)/(W+F)` = (windowed − full-field)/(sum). The panel's `ssi` is at the windowed **preferred direction + pref SF** (a single condition, NOT folded); `ssi_avg_at_pref_sf` averages over all 12 directions at pref SF. The plot folds for visual clarity; it is not literally how `ssi` is computed. Full-field tuning (`tuningf_b64`, 24 uint8/cell, same shared peak as windowed) and `ssi_avg_pref_sf` were added to the explorer data for this. A grating-responsiveness badge (green resp-yes / amber resp-no) sits directly above the tuning plots, and non-responsive cells' tuning traces are washed out (opacity 0.35). The preferred-natural-image section shows an amber warning when `ni_frac_responsive_trials < ni_frac_thresh` (0.25) — the image still shows but is flagged as possibly-not-meaningful (mirrors the RF null warning). `ni_frac` was added to the explorer's D columns for this.
- **Notebook SF-index bug (fixed 2026-09-11)**: both `docs/notebooks/*.ipynb` used `sf_idx = int(row.dgw_preferred_sf)` — but `dgw_preferred_sf` is the SF *value* (0.04/0.08), and int() of both is 0, so the example tuning plots always used 0.04 cpd. Fixed to `int(np.argmin(np.abs(spatial_frequencies - dgw_preferred_sf)))`. (Separate from the explorer's base64 read-order bug; the notebooks index the raw npz axes and never had that one.) Notebook outputs are committed from the old 2026-09-09 asset and were not re-executed.
- **Camera (rotation+zoom) persists across re-renders**: tracked in `sceneCamera` via `plotly_relayout` and fed back into every `getTrace` layout (plus `scene.uirevision=view`). `uirevision` alone did NOT hold it because the layout still supplies an explicit camera — the explicit feedback is what works. Reset to default only on **view toggle** (anat↔rf) and the **Reset** button (both set `sceneCamera=null`). Note: because the layout now carries the tracked camera, Plotly's double-click no longer returns to the default angle — Reset is the way back.
- **Selection pins the camera** (`selectCell` snapshots the live camera into `pinnedCam`, sets `camPinned=true`). While pinned, the `plotly_relayout` handler reverts any camera change back to `pinnedCam` and stops updating `sceneCamera`; a fresh `pointerdown` (or view-toggle/Reset) clears the pin. This defeats a **metric-dependent, delayed** gl3d camera reset the user hit: selecting a cell while coloring by `rf_az_on` (a `div`-scale, mostly-null metric) snapped the view to default, while `preferred_direction` did not — and it could NOT be reproduced in the Electron test pane at all (raw restyle/resize/select all held the camera). The pin was validated by **simulating** the reset (forcing camera to default while pinned → watchdog reverted it to the user's camera; relayout count stayed ~5, no loop). Also snapshot the camera (`snapshotCam`, deep copy) because `layout.scene.camera` is a live ref Plotly can mutate. If this class of gl3d fragility keeps recurring, switch the main view to three.js (camera fully decoupled from data ops) — the user is open to it.

### Selection mechanism (important gotcha)

scatter3d does **not** reliably fire `plotly_click` (verified: 0 clicks registered even on
enlarged markers), but `plotly_hover` is dependable. So selection tracks `hoverIdx` from
plotly_hover and commits it on a capture-phase `pointerup` on the chart div. `customdata` on
the trace carries the global ROI index. Do not "fix" this back to plotly_click.

**Selection is Ctrl/⌘+click, not plain click** (changed 2026-09-11 after the user hit a hard
freeze). Two reasons: (1) plain left-click must stay free for orbit/drag; (2) the original
code called a full `Plotly.react` of all 39,407 WebGL points on every selection to redraw the
highlight ring, which thrashed the WebGL context and hung the window. Now: plain click/drag =
orbit; Ctrl/⌘+click (moved <5px, `ctrlKey||metaKey`) = select. The highlight is a **persistent
trace at index 1** updated with `Plotly.restyle(...,[1])` — selection never re-renders the
cloud. Hover text is precomputed once (`HOVER_ALL`), not rebuilt per render. If the freeze ever
returns, the next step is replacing the Plotly scatter3d main view with three.js points +
OrbitControls + raycaster picking (panels can stay Plotly); the user is open to this.

### Coordinate decisions

- Column offsets use azimuth scale (14.95 µm/°) for both X and Y axes
- Left hemisphere confirmed (de Vries 2019: monitor at right eye, craniotomy over contralateral left V1)
- Z axis [530, 30] reversed (depth increases downward)

**Axis labels corrected 2026-09-09** via coregistration sanity check against EM data
(v1169, 564 matched ROIs, affine fit r > 0.997 all axes):
- **X: low = lateral, high = medial** (was medial/lateral — swapped per user's EM convention: left hemisphere, increasing EM_X = lateral→medial)
- **Y: low = posterior, high = anterior** (was anterior/posterior — swapped; EM_Z post→ant increasing maps to +2P_Y)
- **Depth: pia/wm confirmed correct** (EM_Y pia→wm maps to +2P_depth)

**Why:** The user stated EM convention for left hemisphere: if increasing Z = posterior→anterior, then increasing X = lateral→medial. The Z convention is the user's assumption to verify offline — both X and Y corrections depend on it. See [[column-layout-is-one-dimensional]].
