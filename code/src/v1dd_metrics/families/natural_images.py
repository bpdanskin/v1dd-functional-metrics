"""Natural images, both the 118-image set and the 12-image repeat set.

See docs/families/natural_images.md."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .. import responses as tr
from ..common import _condition_codes, _lifetime_sparseness_chunked, _nanmean
from ..config import DEFAULT_CONFIG, MetricConfig
from ..schema import _reliability_on, absent_frame, roi_frame

def natural_images_metrics(
    plane,
    trials: pd.DataFrame,
    spont: Sequence[float],
    *,
    ns_type: str = "natural_images",
    config: MetricConfig = DEFAULT_CONFIG,
    rng: Optional[np.random.Generator] = None,
    mouse: Optional[str] = None,
) -> pd.DataFrame:
    """Natural-image metrics: one row per ROI.

    Structurally the same as `natural_movie_metrics` — group trials by condition, take
    the condition means, find each neuron's preferred one — with two differences that
    matter:

    * **Conditions are `image_index`, not `image_order`.** `image_order` is the raw
      presentation slot; `image_index` is the image's identity in the 118-image catalog.
      `natural_images_12` draws twelve images from that *same* namespace, so its
      `pref_img` values are a sparse subset of 0..117 (2, 4, 5, ..., 68) rather than
      0..11. Re-ranking them to 0..11 would look tidier and be wrong.
    * **`frac_responsive_trials` is a statistical test here**, unlike natural movie's
      `mean(response > 0)`: the fraction of preferred-image trials whose response beats a
      bootstrapped spontaneous null at p < 0.05. So this column carries bootstrap noise
      and should be read against a seed control, not against zero.

    The response window is `config.ni_response_seconds`. The original took it from an NWB
    `duration_sec` attribute that the current files no longer carry, so it is a recovered
    parameter rather than a known one — see the window probe in the notebook.

    Returns `(metrics, (condition_means, image_ids))`, the same two-value shape
    `receptive_field_metrics` uses. `condition_means` is `(n_rois, n_images)` trial-
    averaged responses, or `None` when the stimulus is absent. `image_ids` matters:
    `natural_images_12` draws from the same 118-image namespace, so its columns are a
    sparse subset of 0..117 rather than 0..11.
    """
    if not len(trials):
        return absent_frame(plane, ns_type, mouse), None
    rng = np.random.default_rng() if rng is None else rng
    traces = plane.traces[config.trace_type[ns_type]]
    window = (0.0, float(config.ni_response_seconds))

    starts = trials["start_time"].to_numpy(dtype=np.float64)
    img = trials["image_index"].to_numpy()
    if np.isnan(img).any():
        raise ValueError(f"{ns_type}: trials contain NaN image_index")
    img = img.astype(int)

    image_ids = np.unique(img)
    code = _condition_codes(img, image_ids)
    n_trials = int(np.bincount(code).max())

    if config.ni_response_frames is not None:
        sweeps = tr.sweep_responses_frames(traces, plane.timestamps, starts,
                                           int(config.ni_response_frames))
    else:
        sweeps = tr.sweep_responses(traces, plane.timestamps, starts, window, None)
    ta = tr.trial_array(sweeps, code, n_trials=n_trials, n_conditions=len(image_ids))

    mean_resp = _nanmean(ta, axis=1)                        # (n_images, n_rois)
    n_rois = plane.n_rois
    roi_ix = np.arange(n_rois)

    all_nan = np.all(~np.isfinite(mean_resp), axis=0)
    pref_idx = np.where(np.isfinite(mean_resp), mean_resp, -np.inf).argmax(axis=0)
    pref_response = np.where(all_nan, np.nan, mean_resp[pref_idx, roi_ix])
    pref_img = np.where(all_nan, -1, image_ids[pref_idx])

    null_single = tr.spontaneous_null(
        traces, plane.timestamps, spont[0], spont[1], window, None,
        n_boot=config.other_n_boot, n_means=1, rng=rng,
        memory_budget_mb=config.memory_budget_mb,
    )
    pref_trials = ta[pref_idx, :, roi_ix]                   # (n_rois, n_trials)
    frac = tr.frac_trials_above_null(pref_trials, null_single, p_thresh=config.sig_p_thresh)

    # The multi-trial null averages n_trials draws per bootstrap sample. For
    # natural_images_12 that is 10,000 x 40 = 400,000 window means -- the heaviest single
    # call in the pipeline, and the reason spontaneous_null blocks by memory budget.
    null_multi = tr.spontaneous_null(
        traces, plane.timestamps, spont[0], spont[1], window, None,
        n_boot=config.other_n_boot, n_means=n_trials, rng=rng,
        memory_budget_mb=config.memory_budget_mb,
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        z_score = (pref_response - null_multi.mean(axis=1)) / null_multi.std(axis=1)

    out = roi_frame(plane, mouse=mouse)
    out["frac_responsive_trials"] = frac
    out["lifetime_sparseness"] = _lifetime_sparseness_chunked(
        ta, over=config.lifetime_sparseness_over)
    out["reliability"] = tr.trial_reliability(ta)
    out["reliability_dff"] = _reliability_on(
        plane, "dff", starts, code, n_trials=n_trials, n_conditions=len(image_ids),
        window=window, frames=config.ni_response_frames)
    out["n_trials_at_pref"] = np.where(
        all_nan, np.nan, np.isfinite(pref_trials).sum(axis=1)).astype(float)
    out["pref_img"] = pref_img
    out["pref_response"] = pref_response
    out["z_score"] = z_score
    # (n_rois, n_images) trial-mean responses, transposed to put ROIs on axis 0 like every
    # other exported array. The published columns are all reductions of this; keeping it is
    # what makes population-level analysis possible at all, since nothing else in the asset
    # carries a neuron-by-condition matrix for the natural stimuli.
    return out, (mean_resp.T.astype(np.float32), np.asarray(image_ids))
