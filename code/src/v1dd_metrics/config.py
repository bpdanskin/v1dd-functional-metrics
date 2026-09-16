"""Tunable settings for every metric family.

``DEFAULT_CONFIG`` computes what we believe is correct; ``REFERENCE_CONFIG``
reproduces the historical behaviour. Divergence is reported in the asset's
provenance as ``differs_from_reference_config``. See docs/pipeline.md."""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Mapping, Optional, Tuple

@dataclass(frozen=True)
class MetricConfig:
    """Knobs, defaulted to computing the right thing.

    Several defaults differ from what reproduces the historical tables —
    `dg_crossval`, `rf_method`, `rf_center_scale_bug`, `pref_cond_fillna`,
    `ni_response_frames`, `impute_dgw_center`, and `lsn_response_frames` — and each is documented
    where it is declared. `fit_all_sf` also differs, but it is a speed knob
    rather than a correction. `REFERENCE_CONFIG` is the historical set, so
    both behaviours are one argument away and which one you asked for is
    written down rather than inferred.
    """

    # --- trace type per family; events everywhere
    trace_type: Mapping[str, str] = field(default_factory=lambda: MappingProxyType({
        "drifting_gratings_full": "events",
        "drifting_gratings_windowed": "events",
        "natural_images": "events",
        "natural_images_12": "events",
        "natural_movie": "events",
        "locally_sparse_noise": "events",
    }))

    # --- response windows
    #: Drifting-grating response window, seconds. The original read `duration_sec` (2.0)
    #: from an NWB attribute; the per-trial rows in this asset say stop-start is 1.985 s.
    #: 2.0 reproduces the published numbers. "per_trial" uses each sweep's own duration.
    dg_response_seconds: Any = 2.0
    #: Natural-images response window, seconds. **Recovered empirically, not read from
    #: the data**: the original took it from an NWB `duration_sec` attribute the current
    #: files no longer carry. Scanning it against the published table gives a sharp
    #: optimum at 0.33 s (median |diff| in lifetime_sparseness = 1e-16, i.e. exact),
    #: while 0.30 s gives 6e-3 and 0.35 s gives 2e-3.
    #:
    #: The reason is discrete. At dt = 0.165 s the samples after an onset sit at
    #: delta, delta+dt, delta+2dt with delta in [0, dt). A 0.33 s window is just under
    #: 2*dt = 0.33008, so it catches **exactly two samples on every trial**; a 0.30 s
    #: window catches two when delta <= 0.135 and one otherwise. That varying count
    #: rescales each trial differently, which `lifetime_sparseness` detects because it
    #: is invariant to a *global* scale but not a per-trial one.
    #:
    #: Used only when `ni_response_frames` is None. Kept because it is what reproduces
    #: the historical tables — see `REFERENCE_CONFIG` — and because the reasoning above
    #: is the evidence for `ni_response_frames = 2`.
    #:
    #: **It does not generalise.** The margin is 8e-5 s, and the pre-flight found dt
    #: spanning 0.16123-0.16671 across the 25 sessions. 0.33 s catches exactly two
    #: samples only where 2*dt lands just above it — dt 0.16504 and 0.16506, which are
    #: precisely the two sessions the value was recovered on. Elsewhere it takes three
    #: samples on up to 4.7 % of trials or one on up to 2.1 %, reintroducing exactly the
    #: per-trial rescaling described above.
    ni_response_seconds: Any = 0.33
    #: Natural-images response window as a FIXED number of imaging samples from each
    #: onset. When set (the default), `ni_response_seconds` is ignored.
    #:
    #: Three samples (~0.5 s at 6 Hz), matching the natural movie window. The historical
    #: value was 2 (recovered from the 0.33 s time window). At 2 frames, L0 event
    #: detection at ~4% per frame per AP gives P(detect in window) ≈ 8%, producing 63%
    #: zero-inflation in the condition means. At 3 frames the detection window rises to
    #: ~12% and the zero-inflation, reliability, and events-vs-dF/F agreement all improve
    #: — the natural movie (already at 3 frames) shows 7% reliability sign disagreement
    #: vs 24% for natural images at 2 frames.
    ni_response_frames: Optional[int] = 3
    #: Natural-movie and LSN windows are counted in *imaging* frames, so they depend on
    #: the plane's own sampling period rather than on the stimulus.
    nm_response_frames: int = 3
    #: LSN response window. The historical value was 4, inherited from the 30 Hz pipeline
    #: where 4 frames = 133 ms (within one 250 ms LSN presentation). At 6 Hz, 4 frames =
    #: 660 ms = 2.6 stimulus presentations of smearing. 2 frames (328 ms, 1.3 presentations)
    #: halves detection rate but produces 68% single-component fields on events vs 43% at
    #: 4 frames on dF/F. See docs/explorations/rf_window_and_trace.md.
    lsn_response_frames: int = 2

    # --- bootstrap
    dg_n_boot: int = 2500
    other_n_boot: int = 10_000
    sig_p_thresh: float = 0.05

    # --- responsiveness thresholds
    dg_frac_thresh: float = 0.50
    ni_frac_thresh: float = 0.25
    rf_frac_thresh: float = 0.25

    # --- receptive field method
    #: ``"greedy"`` (default) is the greedy pixelwise RF (Millman, from v1dd_physiology):
    #: per-pixel STA of events, per-pixel bootstrap p-value (resample sweeps with
    #: replacement), Holm-Šidák correction across all 2*n_pixels tests, detect where the
    #: corrected p < alpha. It dominates the historical fraction-threshold method on both
    #: detection and contiguity at every depth and SNR (see docs/explorations/greedy_rf.md):
    #: fraction gives 81% single-component ON fields, greedy 95-97%. ``"fraction"`` keeps the
    #: historical method (thresholds ``rf_frac_thresh``) and is what REFERENCE_CONFIG uses.
    rf_method: str = "greedy"
    #: Bootstrap resamples for the greedy per-pixel null. At 5000 the alpha sweep aliases:
    #: alpha 0.04..0.01 are identical (only p=0 pixels survive below 0.05), so there are two
    #: operating points — alpha_sens=0.05 and alpha_strict (anywhere <=0.04).
    rf_greedy_n_boot: int = 5000
    #: Default (canonical columns): strict Holm-Šidák alpha -> ~8.6% ON detection, 97%
    #: single-component. The greedy response window reuses ``lsn_response_frames``.
    rf_greedy_alpha_strict: float = 0.01
    #: Sensitive variant (``_a05`` columns): ~11.8% ON detection, 95% single-component —
    #: more yield for downstream at a modest false-positive cost.
    rf_greedy_alpha_sens: float = 0.05

    # --- drifting-grating selectivity
    #: Replace naive OSI/DSI with the **split-half cross-validated** value. The naive
    #: metric picks the preferred condition by argmax and measures selectivity on the same
    #: trials, which biases OSI/DSI **upward** by ~+0.2 on responsive cells and manufactures
    #: apparent selectivity for non-responsive ones (46% of cells have naive OSI>0.5; 57% of
    #: those are not responsive). Cross-validation picks the preferred condition on one trial
    #: half and measures on the other, averaged over `dg_crossval_iters` random splits, which
    #: is (near-)unbiased and folds reliability into the value — an unreliable cell's OSI
    #: collapses toward zero on its own. `gosi`, `preferred_dir`, and `pref_dir_mean` are left
    #: naive: gOSI has no argmax selection so carries no bias, and the preference should stay
    #: the deterministic full-data pick. The value carries an irreducible 8-trial noise floor
    #: the seed pins but cannot remove. `REFERENCE_CONFIG` sets this False (historical naive
    #: OSI/DSI). See docs/explorations/crossval_osi_dsi.md.
    dg_crossval: bool = True
    #: Random trial-splits averaged per ROI. The estimate stabilises within a few hundred;
    #: 200 costs ~15 s over the whole asset.
    dg_crossval_iters: int = 200
    #: Seed for the cross-validation splits, kept separate from the bootstrap `rng` so the
    #: OSI/DSI value is reproducible independent of call order.
    dg_crossval_seed: int = 0

    # --- surround suppression
    running_threshold_cm_s: float = 1.0
    running_pad_seconds: float = 0.10
    #: Radius of the windowed-grating aperture, degrees. **Not recorded anywhere in the
    #: NWB** — there is no size column in the stimulus table. It comes from the V1DD
    #: white paper (Abbasi-Asl et al. 2019), which states a 30 degree diameter twice.
    #: Used only by `window_containment`; no metric that reproduces the historical
    #: tables depends on it.
    dgw_window_radius_deg: float = 15.0

    # --- dF/F spectral SNR (roi_summary)
    snr_signal_band: Tuple[float, float] = (0.1, 1.5)
    snr_noise_band: Tuple[float, float] = (2.0, 2.1)
    ssi_min_trials: int = 3

    # --- expensive extras, absent from every published table
    permutation_test_shuffles: int = 0
    chisq_shuffles: int = 0
    fit_tuning_curves: bool = True
    #: Fit the von Mises curve at EVERY spatial frequency, not just the one surround
    #: suppression reads. False is roughly 2x faster over a full run, because drifting
    #: gratings dominate the runtime and half of every fit was being discarded.
    #:
    #: It changes no published column — `ssi_tuning_fit` is the only consumer and it
    #: reads one SF per ROI — but it is **not** invisible in the asset: `tuning_curves`
    #: exports `dgw_params` / `dgf_params`, and under False the unread SF is NaN. That
    #: reads as a failed fit unless you know better. Set True for a completeness run.
    #:
    #: The original fitted every SF, so `REFERENCE_CONFIG` sets this True and a fast run
    #: therefore shows up in `differs_from_reference_config`. That is deliberate: the
    #: block answers "what did this run do differently", and once the parameters ship,
    #: this qualifies.
    fit_all_sf: bool = False

    #: Fill a session's missing windowed-grating aperture centre from the median of the
    #: other sessions in its cortical column, flagging the rows with
    #: `dgw_center_inferred`. Two of 25 sessions record no centre, and
    #: `probe_window_center.py` confirmed the columns are absent from their stimulus
    #: tables rather than lost by us — see `infer_window_centers`.
    #:
    #: This changes `dgw_center_azimuth` / `dgw_center_elevation` from NaN to a value for
    #: 2,456 ROIs, and therefore the four `dgw_rf_*` containment columns computed from
    #: them. It does not touch any `ssi` column. `REFERENCE_CONFIG` sets it False, since
    #: the original imputed nothing, so a run with it on shows in
    #: `differs_from_reference_config`.
    impute_dgw_center: bool = True

    # --- ROI position
    #: Micrometres per imaging pixel. **Not recorded anywhere in the NWB** --
    #: ``grid_spacing`` is a [1.0, 1.0] 'meters' placeholder on every plane -- and not
    #: stated in the white paper either. Inferred: four columns tile an 800 um grid, so
    #: one field of view is ~400 um across 512 pixels. The paper also says the columns
    #: *overlap*, which makes 400 an upper bound, so treat this as the least trustworthy
    #: number in the asset. Comparing the two anatomical frames is unaffected by it:
    #: both are scaled by the same factor. Set to None to leave every um column NaN.
    um_per_pixel: Optional[float] = 800.0 / 2 / 512

    #: Start the von Mises fit from a guess derived from the curve itself rather than
    #: from the fixed ``p0``. Expected to be faster; measured, it needs **2.1x the model
    #: evaluations** (1617 against 771 per curve) for a fit of equal quality. The cost is
    #: not the starting point but the model: ``k`` is unbounded above inside an
    #: exponential. See docs/families/drifting_gratings.md.
    vonmises_data_p0: bool = False

    #: Whether ``ssi_tuning_fit`` evaluates the fitted curve including its baseline
    #: offset. The preferred direction feeding it is chosen with the baseline
    #: subtracted, so including it here is internally inconsistent.
    ssi_tuning_fit_includes_baseline: bool = False

    #: Whether lifetime sparseness is taken over each condition's mean response
    #: ("conditions", Vinje & Gallant) or over every individual trial ("trials", what
    #: the historical tables shipped). See docs/comparability.md.
    lifetime_sparseness_over: str = "conditions"

    #: Whether a zero denominator yields NaN. The historical `osi`/`dsi` returned 0
    #: there while the surround-suppression indices returned NaN; this makes both NaN.
    zero_denominator_nan: bool = True

    # --- historical compatibility. These default to the CORRECTED behaviour; set both
    # True — or just use REFERENCE_CONFIG — to reproduce the historical tables exactly.
    #: `point_to_alt_azi` in the original divides the centre-to-centre *range* by `n`
    #: rather than `n - 1`, so its degree scale is compressed by `(n-1)/n`: 12.5 % in
    #: altitude (8 rows) and 7.1 % in azimuth (14 columns). The historical tables
    #: therefore span ±28.481° and ±56.132° where the screen actually spans ±32.55° and
    #: ±60.45°. Shipping the compressed scale means anyone who plots retinotopy plots it
    #: wrong, so the default is the true mapping. The two differ by exactly `n/(n-1)`,
    #: which is what makes the correction verifiable rather than merely asserted.
    rf_center_scale_bug: bool = False
    #: The original takes `preferred_dir`/`preferred_sf` from `fillna(-1).argmax`, so an
    #: ROI with no finite response at any condition reports condition **0** rather than
    #: "no preferred condition" — while every other metric in the same function uses a
    #: nan-skipping argmax. False makes the two agree and leaves those ROIs NaN. It
    #: touches only all-NaN rows, but surround suppression keys off the preferred
    #: condition, so a fabricated preference propagates.
    pref_cond_fillna: bool = False

    memory_budget_mb: float = 64.0


DEFAULT_CONFIG = MetricConfig()


REFERENCE_CONFIG = MetricConfig(
    trace_type=MappingProxyType({
        "drifting_gratings_full": "events",
        "drifting_gratings_windowed": "events",
        "natural_images": "events",
        "natural_images_12": "events",
        "natural_movie": "events",
        "locally_sparse_noise": "dff",           # default changed to events
    }),
    dg_crossval=False,             # historical naive OSI/DSI (default now cross-validated)
    rf_method="fraction",          # the historical fraction-threshold RF (default now greedy)
    rf_center_scale_bug=True,      # centres compressed by (n-1)/n
    pref_cond_fillna=True,         # all-NaN ROIs report condition 0
    ni_response_frames=None,       # fall back to the recovered time window
    ni_response_seconds=0.33,
    lsn_response_frames=4,        # default changed to 2
    fit_all_sf=True,               # the original fitted every SF, not just the read one
    impute_dgw_center=False,       # the original imputed no aperture centre
    ssi_tuning_fit_includes_baseline=True,   # inconsistent with its own peak selection
    lifetime_sparseness_over="trials",   # the original flattened every trial
    zero_denominator_nan=False,          # the original's osi/dsi returned 0, not NaN
)
