---
name: data-explorer-artifact
description: "V1DD Metrics Explorer artifact — 39,407-ROI 3D scatter with wireframe session boxes, hide-null filter, anatomical axis labels. Published and saved to docs/v1dd_explorer.html."
metadata:
  type: project
  modified: 2026-09-09
---

## Current artifact

URL: `https://claude.ai/code/artifact/abb43a1d-bb9f-419d-b163-53df50cc3c94`
Local copy: `docs/v1dd_explorer.html` (6.3 MB, self-contained except CDN Plotly)

### Architecture

Three source files concatenated via Python for UTF-8 correctness:
- `scratchpad/part1.html` — CSS, theme tokens, header controls
- `scratchpad/roi_data.json` — column-oriented JSON, 39,407 ROIs × 27 fields
- `scratchpad/part2.js` — all JS logic

Plotly.js 2.27.0 loaded from cdnjs CDN (artifact CSP).

### Features

- **View toggle**: Anatomical / RF Space
- **16 color metrics** in grouped dropdown (tuning, selectivity, suppression, quality, RF, activity)
- **Column/volume/plane chip filters** with fixed axis ranges
- **Threshold sliders**: SNR, OSI, DSI, Confidence
- **Wireframe boxes**: grey bounding boxes for each of 25 column×volume sessions (anatomical view only)
- **Hide null checkbox**: filters ROIs where current color metric is null (crucial for RF azimuth/elevation, 82% null)
- **Anatomical axis labels**: medial/lateral, anterior/posterior, pia/wm on tick extremes (rotation-invariant)
- **Confidence → marker size**: confMin=0.3, confMax=1.0, sizeMin=1.5, sizeMax=5
- **Derived metric**: dgw_pref_ori = dgw_pref_dir % 180

### Coordinate decisions

- Column offsets use azimuth scale (14.95 µm/°) for both X and Y axes
- Left hemisphere confirmed (de Vries 2019: monitor at right eye, craniotomy over contralateral left V1)
- X increasing = medial→lateral, Y increasing = anterior→posterior (microscope orientation assumed consistent)
- Z axis [530, 30] reversed (depth increases downward)
