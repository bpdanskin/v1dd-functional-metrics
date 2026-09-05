---
name: correlations-are-a-separate-pipeline
description: "Cell-cell correlations do not belong in this pipeline — decided 2026-09-03, for shape and size reasons."
metadata:
  node_type: memory
  type: project
  modified: 2026-09-03
---

**Decided 2026-09-03.** Cell-cell correlations stay a separate analysis and become their
own pipeline. Do not propose folding them in without new evidence.

1. **Wrong shape.** Every family here is per-*plane*, one row per ROI. Correlations are
   per-*session* and cross-plane — all 6 planes interpolated onto a reference timebase,
   which is the entire point, since cross-plane pairs carry the connectivity story.
   Folding them in changes the processing loop, not just adds a function.
2. **Wrong size.** Within-session pairs across planes are **37.4 M** (6.88 M within-plane),
   which at 7 stimulus columns is about 1.05 GB upper-triangle. This whole asset is 128 MB.

Compute was never the objection — correlations are matmuls, minutes against a run that is
96 % grating bootstraps.

The reuse inventory is in the frozen fork's `HANDOFF.md`, section "Shared or duplicated?
The inventory before the pipelines split". Headline: **nothing was actually shared** — the
correlations notebook imports nothing from the old `code/utils/`, so every overlap is a
second copy. `functional_similarity.py` belongs to the correlations side, not here.

Also keep out of any reproducible pipeline: the **CAVE coregistration merge**, which needs
an interactive token.
