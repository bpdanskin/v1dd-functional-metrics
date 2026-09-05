# Outputs

One run writes a single asset directory, `/results/<mouse>_V1DD_functional_metrics_<stamp>/`:

| file | what it holds |
|---|---|
| `stimulus_metrics.parquet` | the per-ROI table — every metric, one row per ROI |
| `tuning_curves.npz` | per-trial grating responses, blank sweeps, von Mises fits, running speed |
| `receptive_field_maps.npz` | pre-threshold ON/OFF subfield maps |
| `condition_means.npz` | neuron-by-image trial-mean responses for the two image sets |
| `provenance.json` | seed, config, sessions, environment, and what differs from the reference config |
| `subject.json`, `data_description.json`, `processing.json` | AIND metadata, written last |

Sidecars sit **inside** the asset directory beside the data; run logs stay above it,
because they are a record of the run rather than data.

## Why there are no per-family CSVs

Earlier assets shipped eight per-family CSVs *and* the wide table. The CSVs repeated the
nine-column identity block eight times and carried nothing the wide table lacks — 56 MB of
CSV encoding the same 72 metric columns the table holds in 15 MB. They are gone. Anything
that read one reads a column subset of `stimulus_metrics.parquet` instead.

## The per-ROI table

Nine identity columns, then the metric blocks below. Families whose published names
already name the stimulus keep an empty prefix.

| family | prefix | n | columns |
|---|---|---|---|
| `roi_summary` | — | 11 | `snr`, `signal_power`, `noise_power`, `run_frac`, `spont_run_frac`, `spont_rate`, `spont_rate_run`, `spont_rate_stat`, `run_mod_dgf`, `run_mod_dgw`, `run_mod_spont` |
| `drifting_gratings_full` | `dgf_` | 9 | `dsi`, `frac_responsive_trials`, `gosi`, `is_responsive`, `lifetime_sparseness`, `osi`, `preferred_dir`, `preferred_sf`, `pref_dir_mean` |
| `drifting_gratings_windowed` | `dgw_` | 9 | the same nine |
| `surround_suppression` | — | 15 | `ssi`, `ssi_avg`, `ssi_avg_at_pref_sf`, `ssi_running`, `ssi_running_avg_at_pref_sf`, `ssi_stationary`, `ssi_stationary_avg_at_pref_sf`, `ssi_tuning_fit`, `dgw_center_azimuth`, `dgw_center_elevation`, `dgw_center_inferred`, `dgw_rf_distance_on`, `dgw_rf_distance_off`, `dgw_rf_overlap_on`, `dgw_rf_overlap_off` |
| `natural_images` | `ni_` | 7 | `frac_responsive_trials`, `lifetime_sparseness`, `pref_img`, `pref_response`, `z_score`, `reliability`, `reliability_dff` |
| `natural_images_12` | `ni12_` | 7 | the same seven |
| `natural_movie` | `nm_` | 7 | the same seven, with `pref_img` a frame index |
| `rf_metrics` | — | 7 | `has_rf_on`, `has_rf_off`, `has_rf_on_or_off`, `azimuth_rf_on`, `altitude_rf_on`, `azimuth_rf_off`, `altitude_rf_off` |
| `roi_position` | — | 10 | `roi_x_px`, `roi_y_px`, `roi_area_px`, `roi_radius_px`, `roi_x_um`, `roi_y_um`, and the two anatomical frames `roi_{x,y}_um_published` / `roi_{x,y}_um_retinotopic` — see [roi_position.md](families/roi_position.md) |

Identity: `roi_unique_id`, `roi_key`, `mouse`, `column`, `volume`, `plane`, `roi`,
`depth_um`, `pika_roi_confidence`.

**Join on `roi_key`, or on `(column, volume, plane, roi)`.** `roi_unique_id` omits the
column and collides — about 13,500 distinct strings for 39,407 rows. `volume` is a string
throughout, because volumes run 1–9 and a–f; a CSV round-trip would re-infer it as an
integer for an all-numeric column.

## The array archives

Each is keyed by `roi_key`, which joins straight to the table.

**`tuning_curves.npz`** — `dg{w,f}_trials` (n_rois, n_dir, n_sf, n_trials) NaN-padded;
`dg{w,f}_blank` (n_rois, max_n_blank) with `dg{w,f}_n_blank` giving each plane's true
width, since sessions ran 5–8 grey sweeps; `dg{w,f}_params` (n_rois, n_sf, 6) von Mises
in scipy parameter order; `dg{w,f}_running` (n_planes, n_dir, n_sf, n_trials) in cm/s with
**no ROI axis**, so it keys on `plane_key` (`roi_key` minus its trailing `_{roi}`); plus
`directions`, `spatial_frequencies`, `trace_type`.

Trial-level rather than trial means because means are one `nanmean` away and the reverse
is not. Only the spatial frequency surround suppression reads is fitted, so the other
column of `*_params` is **NaN by design, not by failed fit** — set `fit_all_sf=True` for a
completeness run.

**`receptive_field_maps.npz`** — `rf_maps` (n_rois, 2, n_rows, n_cols) float32, dim 1 is
ON then OFF. Each value is the fraction of that pixel's presentations producing a response
above the ROI's own bootstrapped spontaneous 95th percentile, **before** the 0.25
threshold. Invalid ROIs are all-zero, so a blank map means *excluded*, not *no receptive
field*. `altitudes`, `azimuths` and `seed` travel with it: without them pixel indices
cannot become degrees and per-pixel significance cannot be reproduced.

**`condition_means.npz`** — `ni_mean` (n_rois, 118) and `ni12_mean` (n_rois, 12) with
their `*_images` ids. Every published natural-image column is a reduction of this matrix,
and a reduction cannot be un-taken. Natural movie is deliberately absent: (n_rois, 3600)
is ~541 MiB even after averaging over repeats.

## Provenance

`provenance.json` records the seed, the full config, every session processed, the
environment, and `differs_from_reference_config` — the settings that would have to change
to reproduce the historical pipeline, as a delta rather than prose so the claim stays
checkable. `complete_asset` is true only when no session failed, no array writer failed,
and no session filter was applied.

`failed_outputs` lists any array archive that could not be written. The writers are
guarded individually because provenance is written *after* them, and an unguarded raise
once cost a five-hour run its record of everything else.
