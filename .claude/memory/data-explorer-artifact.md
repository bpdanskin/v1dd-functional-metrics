---
name: data-explorer-artifact
description: "V1DD Metrics Explorer (docs/v1dd_explorer.html) — 39,407-ROI 3D scatter + click-to-open per-cell side panel (tuning, polar, RF map, preferred NI image). Refreshed from the 2026-09-11 events/RF-window asset."
metadata:
  type: project
  modified: 2026-09-11
---

## Current artifact

URL: `https://claude.ai/code/artifact/abb43a1d-bb9f-419d-b163-53df50cc3c94` (**Version 13**,
published 2026-09-15: greedy/CV asset + direction-annotation rebuild (V12) + camera-pin
hardening (V13, see the camera note below).
Local copy: `docs/v1dd_explorer.html` (14.63 MB, self-contained except CDN Plotly; under the
16 MB limit but close — the per-cell arrays in detail.json are the main mass, so trim there
before adding per-cell data. Folding the RF mask into the magnitude array, 0 = not-significant,
saved ~1.8 MB vs shipping a separate mask blob).

Data source refreshed 2026-09-15 from `results/409828_V1DD_functional_metrics_2026-09-14_22-59-50`
(the capsule rerun: **greedy RF** default + RF area + population-RF, **cross-validated OSI/DSI**,
`frac_zero_response` pruned). Anatomical `x_um`/`y_um` are reused verbatim from the prior
roi_data.json (row order identical to the parquet on column/volume/plane/roi), so the
coregistration transform is NOT recomputed — only metric columns refresh. RF panel now reads
the greedy archive (`rf_sta` + `strict_mask`), not the fraction `rf_maps` (which the asset no
longer ships).

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
- **RF panel shows the greedy strict field** (rebuilt 2026-09-15): the per-cell STA magnitude
  (`rf_b64`, per-cell normalised) at the strict-significance pixels only; non-significant pixels
  are 0 → transparent (the mask is folded into the magnitude array in the build, significant
  pixels floored to ≥1 so a faint detection never rounds away). `renderRF` (`rfMaskedMap`)
  treats 0 as null. Header reads "Receptive field (greedy, significant pixels)". Gate is
  `has_rf_on_or_off` (strict); 6,147 cells stored, ~1.24 sig px each. The old fraction-threshold
  (0.25) display is gone with the fraction method.
- **Panel tuning plots**: fixed per-SF colors shared across polar/direction-curve/orientation-curve (0.04 cpd `#3a96b5` blue, 0.08 cpd `#d68a2e` orange; pref dir/ori marker `#e0484d` red); preferred SF is bolded + labelled "(pref)". The folded-orientation plot is the **surround-suppression** view: windowed (SF color) vs full-field (`#8e5bb5` purple), both folded at the pref SF, captioned with SSI. **SSI facts** (`families/surround_suppression.py`): every variant is `(W−F)/(W+F)` = (windowed − full-field)/(sum). The panel's `ssi` is at the windowed **preferred direction + pref SF** (a single condition, NOT folded); `ssi_avg_at_pref_sf` averages over all 12 directions at pref SF. The plot folds for visual clarity; it is not literally how `ssi` is computed. Full-field tuning (`tuningf_b64`, 24 uint8/cell, same shared peak as windowed) and `ssi_avg_pref_sf` were added to the explorer data for this. A grating-responsiveness badge (green resp-yes / amber resp-no) sits directly above the tuning plots, and non-responsive cells' tuning traces are washed out (opacity 0.35). The preferred-natural-image section shows an amber warning when `ni_frac_responsive_trials < ni_frac_thresh` (0.25) — the image still shows but is flagged as possibly-not-meaningful (mirrors the RF null warning). `ni_frac` was added to the explorer's D columns for this.
- **Notebook SF-index bug (fixed 2026-09-11)**: both `docs/notebooks/*.ipynb` used `sf_idx = int(row.dgw_preferred_sf)` — but `dgw_preferred_sf` is the SF *value* (0.04/0.08), and int() of both is 0, so the example tuning plots always used 0.04 cpd. Fixed to `int(np.argmin(np.abs(spatial_frequencies - dgw_preferred_sf)))`. (Separate from the explorer's base64 read-order bug; the notebooks index the raw npz axes and never had that one.) Notebook outputs are committed from the old 2026-09-09 asset and were not re-executed.
- **Camera (rotation+zoom) persists across re-renders**: tracked in `sceneCamera` via `plotly_relayout` and fed back into every `getTrace` layout (plus `scene.uirevision=view`). `uirevision` alone did NOT hold it because the layout still supplies an explicit camera — the explicit feedback is what works. Reset to default only on **view toggle** (anat↔rf) and the **Reset** button (both set `sceneCamera=null`). Note: because the layout now carries the tracked camera, Plotly's double-click no longer returns to the default angle — Reset is the way back.
- **Selection pins the camera** (`selectCell` snapshots the live camera into `pinnedCam`, sets `camPinned=true`). While pinned, the `plotly_relayout` handler reverts any camera change back to `pinnedCam` and stops updating `sceneCamera`; a fresh `pointerdown` (or view-toggle/Reset) clears the pin. **Hardened 2026-09-15 (V13):** the one-shot reassert + relayout-handler revert both MISS a delayed gl3d reset that fires *without* emitting `plotly_relayout` — which is what the user hit after the rebuild (camera snapping to default on ctrl-click-select under every metric, worse than the earlier azimuth-only report). Added `holdPinnedCamera()` — a `setInterval` (90 ms, ~1.6 s window) started in `selectCell` that re-applies `pinnedCam` whenever the live camera drifts from it (camEqual-guarded so it no-ops otherwise), cleared on `pointerdown`. Verified in-pane: forcing the camera to default while pinned snaps back within one tick; a simulated user pointerdown releases the pin so rotation sticks. **CONFIRMED FIXED on live claude.ai (V13) by the user 2026-09-16** — the camera now holds through ctrl-click selection across metrics. So the watchdog is sufficient; the three.js migration is NOT needed for this (keep it in reserve only if a *different* gl3d fragility appears). The native reset still isn't reproducible in the browser pane (only on live gl3d), so future camera changes here must be validated live, not just in-pane. This defeats a **metric-dependent, delayed** gl3d camera reset the user hit: selecting a cell while coloring by `rf_az_on` (a `div`-scale, mostly-null metric) snapped the view to default, while `preferred_direction` did not — and it could NOT be reproduced in the Electron test pane at all (raw restyle/resize/select all held the camera). The pin was validated by **simulating** the reset (forcing camera to default while pinned → watchdog reverted it to the user's camera; relayout count stayed ~5, no loop). Also snapshot the camera (`snapshotCam`, deep copy) because `layout.scene.camera` is a live ref Plotly can mutate. If this class of gl3d fragility keeps recurring, switch the main view to three.js (camera fully decoupled from data ops) — the user is open to it.

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

### Direction annotation — DONE 2026-09-15 (shipped in Version 12)

Implemented in `part2.js` / `part1.html`: `compassTraces()` draws the A/P/M/L rosette as two
scatter3d traces (lines + text) at the pia plane, in data coords so it rotates with the scene
(A/P cyan, M/L pink); `projAnat(x,y)` + `anatomyMode` (the "Snap to anatomy" checkbox `chk-anatomy`)
project the cloud onto the medial/anterior basis so axes become `medial → (µm)` / `anterior → (µm)`;
scope-frame axis titles are `2P X (EM X)` / `2P Y (EM Z)` / `2P Z depth (EM Y)` with the old
lateral/medial·posterior/anterior tick labels removed; column chips carry a `<small>` secondary
label (1 ctr, 2 lat, 3 post, 4 ant, 5 med). The compass basis (`DETAIL.anat`) is computed in
`build_explorer_data.py` from the column centroids (ant = col4−col3, med = col5−col2 orthogonalised);
measured anterior [-0.87,-0.49], medial [0.49,-0.87] in the x_um/y_um frame (~61° off the plot axes).
Verified in-browser: compass renders, snap-to-anatomy rotates + relabels, RF strict-mask panel
shows discrete significant pixels, CV OSI/DSI live in hover/panel. Original spec below.

### Direction annotation — spec (agreed 2026-09-14)

Source: `v1dd_context/.../visual_area_maps/v1dd_cortex_maps.pdf` (5 pages, one per mouse;
p1 = M409828). It draws the 5-column block on the widefield ISI retinotopy. **The block's
imaging fields of view are tilted ~45° to true anatomy** — so the plotted `x_um`/`y_um`
(scan frame; each FOV is an axis-aligned ~394 µm square in the explorer data) are ~45°
rotated (and possibly reflected) from anatomical A/P–M/L. The **column centres**, however,
sit on the anatomical axes.

Anatomical frame per user: **anterior = up, medial = right** (posterior = down, lateral =
left). Measured column→direction map off the vector PDF (bearings within ~10° of cardinal):
- **col4 = anterior, col5 = medial, col3 = posterior, col2 = lateral, col1 = centre.**

Four changes to make:
1. **Column-anchored anatomical compass** (primary). In the plot's own coordinates compute
   `anterior ∝ centroid(col4) − centroid(col3)`, `medial ∝ centroid(col5) − centroid(col2)`,
   orthogonalise medial against anterior; P = −A, L = −M. Draw a small 4-arrow A/P/M/L
   rosette at the block, **rotating with the 3D scene**. Because it is built from column
   identities (not screen geometry), it self-corrects for the 45° AND any reflection/sign.
   Sanity check: A must point toward col4's side, M toward col5's.
2. **"Snap to anatomy" toggle** — button that rotates the scene so anatomy is axis-aligned;
   only THEN may axis titles read A/P and M/L. Default view stays the scope frame.
3. **Remove the anatomical (M/L, A/P) labels from the axis extremes** — inaccurate given the
   45° tilt. (This retires the "Axis labels corrected 2026-09-09" reading below.)
4. **Column labels**: keep the column NUMBER primary; add the anatomical direction as a
   smaller SECONDARY label (e.g. `4` with `anterior` beneath / in parens).
5. **Axis titles carry 2P + EM nomenclature** (for later volume coregistration), NOT anatomy.
   2P↔EM axis map **confirmed 2026-09-15 by querying the CAVE table** (`v1dd_public`,
   `functional_coregistration_manual_2`, server `https://global.em.brain.allentech.org`;
   565 rows deduped on col/vol/plane/roi, 559 joined to the explorer ROIs; `pt_position` in
   voxels 9×9×45 nm). Correlating each 2P axis vs each EM axis: **2P X ↔ EM X r=+0.996,
   2P Y ↔ EM Z r=+0.889, 2P Z/depth ↔ EM Y r=+0.875 — all SAME-sign (no flips).** So titles:
   `2P X (EM X)`, `2P Y (EM Z)`, `2P Z depth (EM Y)`, each same sense (↑2P ⇒ ↑EM; depth and
   EM Y both increase pia→white matter). The r=0.996 also shows **2P and EM are nearly the
   same frame** (a relabeling, not a rotation) — the ~45° tilt is between the block/scope-EM
   plane and the *widefield anatomical* A/P–M/L, which is why anatomy needs the compass (#1)
   and the axes stay in 2P/EM nomenclature. Query form: `client.materialize.query_table(
   "functional_coregistration_manual_2")`, keys col/volume/plane/roi/pt_root_id/pt_position;
   caveclient needs pandas<3 (it downgrades pandas — reinstall the repo's pandas after). Note
   duplicate pt_root_ids across planes; dedupe on (column,volume,plane,roi).

**Supersedes the 2026-09-09 anatomical axis labelling.** That coregistration fit (x_um↔EM_X,
r > 0.997) is still valid as a 2P↔EM correspondence, but its *anatomical* interpretation
("X = lateral↔medial, Y = posterior↔anterior") was the shaky part — the EM tangential axes are
block-aligned (~45° off true anatomy), so those axes are not A/P–M/L. A sign inconsistency
confirms it: the PDF puts col4 = anterior (up), whereas the old "high y_um = anterior" would
make col3 anterior. The compass (PDF-anchored, in data coords) is authoritative and will
reveal the actual data sign convention.

**Axis labels corrected 2026-09-09** (SUPERSEDED — see above; kept for the EM correspondence)
via coregistration sanity check against EM data (v1169, 564 matched ROIs, affine fit
r > 0.997 all axes): X↔EM_X, Y↔EM_Z, depth↔EM_Y. The anatomical direction reading attached to
each is retired in favour of the column-anchored compass. See [[column-layout-is-one-dimensional]].
