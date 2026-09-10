---
name: queued-ni-stimulus-images
description: "Next capsule run: extract the 118+12 natural image stimuli from NWB and save as an NPZ in the asset for quick reference without re-reading NWB."
metadata:
  type: project
  modified: 2026-09-09
---

Collect the actual natural-image stimuli (118 NI + 12 NI12) from the NWB sessions
and save them as an NPZ alongside the asset. The images must share a consistent
index across sessions (verified against `ni_images` and `ni12_images` in
`condition_means.npz`).

**Why:** the `ni_mean` response vectors are indexed by image, but the images
themselves are not shipped — answering "which images drive this neuron?" currently
requires going back to the NWB. A small stimulus archive makes the notebook
self-contained.

**How to apply:** add to the next capsule run. Not blocking P5. Check that
image arrays are stored in the NWB stimulus templates and confirm index alignment
with `condition_means.npz` keys. See [[data-explorer-artifact]] for the asset
structure.
