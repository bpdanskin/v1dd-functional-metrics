---
name: reference-materials
description: "Where the reference papers and notebooks live, and what each is good for."
metadata:
  node_type: memory
  type: reference
  modified: 2026-09-05
---

`C:/Users/bethanny.danskin/OneDrive - Allen Institute/swdb_material_2026/v1dd_context/`
Not in the repo -- shared from the user's OneDrive.

| file | what it is good for |
|---|---|
| `V1DD_WhitePaper_v6.pdf` | The dataset description. Imaging geometry is in the text (page 5); the column layout is only in Figures 2 and 3, which are images and extract as no text. Methods section states a dF/F-vs-events choice **opposite** to ours -- see [[queued-events-vs-dff-study]]. |
| `Jewell_Witten_2018.pdf` | *Exact spike train inference via l0 optimization*. The algorithm behind our `events` trace. |
| `elife-51675-v3.pdf` | Huang et al. 2021, simultaneous spikes and GCaMP6 fluorescence. Ground truth for judging event inference. |
| `de Vries Lecoq Buice 2019.pdf` + supplement | The 30 Hz Brain Observatory comparison; source of the agreement figures in docs/comparability.md. |
| `workshop2 - extended version.ipynb` | An independent derivation of grating tuning; corroborates the 1,038 low-confidence ROIs via `is_soma == False`. |

PDFs extract with `pypdf`; there is no page renderer in this environment, so figures
cannot be read here. Text extraction of a figure-only page returns the page number and
nothing else -- that is how the column layout turned out to be unavailable.
