---
name: repo-identity-and-the-frozen-fork
description: "What this repo is, what it replaces, and where the reference asset and prior working notes live."
metadata: { node_type: memory, type: project, modified: 2026-09-03 }
---

`bpdanskin/v1dd-functional-metrics` — a standalone Code Ocean pipeline producing per-ROI
functional metrics from the V1DD NWB sessions. Started 2026-09-03 from the SWDB capsule
template.

It replaces the pipeline that grew inside
`SWDB_2026_Connectomics-functional-metrics`, a fork of the shared AllenSWDB workshop repo.
**That fork is frozen and stays as the reference** — do not develop there, but do read it:
its `HANDOFF.md` and `.claude/memory/` hold the full development record, and its
`results/` holds the shipped assets.

**The reference asset** is
`SWDB_2026_Connectomics-functional-metrics/results/409828_V1DD_stimulus_metrics_2026-09-03_15-55-03`
— complete and clean: 25 sessions, 150 planes, 39,407 ROIs, `complete_asset: true`,
`failed_sessions: []`, built from `eea6957`, 5.1 h. It is the numerical target during the
refactor; see [[array-replay-validates-offline]].

**This is a fresh start.** No output compatibility is owed to the 2025 `allen_v1dd`
pipeline, the 2019 white paper, or the fork's own assets. Comparability is a documentation
topic — see [[fresh-start-metric-changes]] for what that permits.

Related: [[capsule-code-directory-constraint]], [[correlations-are-a-separate-pipeline]].
