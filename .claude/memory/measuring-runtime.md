---
name: measuring-runtime
description: "Two runtime claims in this project were wrong: a session timer compared against a run timer, and an extrapolation from the second-largest session. Where the real numbers are and how to read them."
metadata:
  node_type: memory
  type: feedback
  modified: 2026-09-05
---

Wall-clock claims about this pipeline have been wrong twice, both times confidently, and
both times because two numbers were compared that were never the same measurement.

**What the pipeline prints.** Three different quantities, easy to confuse:

* `col1 vol3 <name>   900.6s` — *that session*, printed inside the loop.
* `6 planes, 2708 ROIs, 20.5 min` — the whole session loop, printed after it.
* `wall_seconds` in `provenance.json` — the same quantity as the second, so the printed
  minutes and the recorded seconds agree. It stops **before** `build_wide` and the array
  writers, so it is the metric loop, not the run.

**Error 1.** 900.6 s (session) was compared against 1306.7 s (run) and reported as a
45 % regression caused by reading `pixel_mask` per ROI. The comparable pair is 1230 s and
1306.7 s: **+6 %**, and all of it the new `roi_position` family, since the earlier run
already read the masks and discarded them.

**Error 2.** 25 x 22 min = 9.1 h treated col 1 / vol 3 as typical. It is the
**second largest of 25 sessions** — 2,708 ROIs against a 1,576 mean. Session ROI counts
run 474 to 2,733, a 5.8x spread, so no single session extrapolates.

**How to do it instead.** Runtime is close to linear in ROI count: 0.489 s/ROI at 2,708
and 0.476 s/ROI at 1,550, two sessions that differ in almost every other way. Scale by
`n_rois / 39,407`, not by session count. The full asset projects to **~5.4 h**.

Provenance now records `stage_seconds` (per family, plus `load_traces` and `load_masks`)
and `mask_reads` (planes per mask column, and bulk against per-ROI reads), so the next
such question is answered from the asset. Use those before doing arithmetic on log lines.

**Why to apply it:** a wrong runtime number is not harmless here — it drove an
optimisation that bought nothing and a P6 budget nearly twice the truth.

Related: [[project-status]], [[two-anomalous-sessions]].
