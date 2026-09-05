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

    Four defaults differ from what reproduces the historical tables —
    `rf_center_scale_bug`, `pref_cond_fillna`, `ni_response_frames` and
    `impute_dgw_center` — and each is documented where it is declared. `fit_all_sf` also
    differs, but it is a speed knob rather than a correction. `REFERENCE_CONFIG` is the historical set, so both behaviours are one
    argument away and which one you asked for is written down rather than inferred.
    """

    # --- trace type per family; events everywhere except receptive fields
    trace_type: Mapping[str, str] = field(default_factory=lambda: MappingProxyType({
        "drifting_gratings_full": "events",
        "drifting_gratings_windowed": "events",
        "natural_images": "events",
        "natural_images_12": "events",
        "natural_movie": "events",
        "locally_sparse_noise": "dff",
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
    #: Two samples, because that is what the 0.33 s window was doing on the sessions it
    #: was tuned against — but expressed in the units the intent actually lives in, so it
    #: cannot drift with the sampling rate. A time window varies with onset phase; a
    #: frame count does not. Since per-trial rescaling changes the relative pattern across
    #: images, scale-invariant metrics like `lifetime_sparseness` can tell the two apart,
    #: which is how the discrepancy was found in the first place.
    ni_response_frames: Optional[int] = 2
    #: Natural-movie and LSN windows are counted in *imaging* frames, so they depend on
    #: the plane's own sampling period rather than on the stimulus.
    nm_response_frames: int = 3
    lsn_response_frames: int = 4

    # --- bootstrap
    dg_n_boot: int = 2500
    other_n_boot: int = 10_000
    sig_p_thresh: float = 0.05

    # --- responsiveness thresholds (different per family in the original)
    #
    # A fraction is not a comparable criterion across datasets, because the same fraction
    # is a different statistical test at a different trial count. Each trial has a 5 %
    # chance of passing `sig_p_thresh` by chance, so the false-positive rate is the
    # binomial tail:
    #
    #   de Vries 2019   >= 25 %  = 4 of 15 trials  ->  0.0055
    #   here            >= 25 %  = 2 of 8          ->  0.0572   (10x LOOSER)
    #   here            >= 37.5 % = 3 of 8         ->  0.0058   (matched)
    #   here            >= 50 %  = 4 of 8          ->  0.00037  (15x STRICTER)
    #
    # So 0.50 here is much stricter than the Brain Observatory's nominally-lower 0.25,
    # and **`>= 0.375` is the like-for-like comparison** — not 0.25. Comparing our
    # `is_responsive` rate against their 25 % figure without that adjustment makes V1DD
    # look less responsive than it is. 0.50 is kept because it reproduces the V1DD white
    # paper's 26 % headline; see [[v1dd-metrics-open-questions]].
    #
    # A second wrinkle: the denominator is not always 8. Blank sweeps are ragged, so
    # 184-187 grating sweeps land in 192 condition slots and some conditions get 7 trials.
    # `frac_responsive_trials` therefore takes 23 distinct values in the shipped asset,
    # not 9 — see [[dg-blank-sweeps-are-ragged]].
    dg_frac_thresh: float = 0.50
    ni_frac_thresh: float = 0.25
    rf_frac_thresh: float = 0.25

    # --- surround suppression
    running_threshold_cm_s: float = 1.0
    running_pad_seconds: float = 0.10
    #: Radius of the windowed-grating aperture, degrees. **Not recorded anywhere in the
    #: NWB** — there is no size column in the stimulus table. It comes from the V1DD
    #: white paper (Abbasi-Asl et al. 2019), which states a 30 degree diameter twice.
    #: Used only by `window_containment`; no metric that reproduces the historical
    #: tables depends on it.
    dgw_window_radius_deg: float = 15.0

    # --- dF/F spectral SNR (roi_summary). Bands from the cell-cell correlations
    # notebook, whose `snr_by_cell.feather` these are meant to be comparable with.
    # Calcium transients occupy the signal band; the narrow band above them samples a
    # spectrum that is close to flat for a well-isolated ROI.
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
    rf_center_scale_bug=True,      # centres compressed by (n-1)/n
    pref_cond_fillna=True,         # all-NaN ROIs report condition 0
    ni_response_frames=None,       # fall back to the recovered time window
    ni_response_seconds=0.33,
    fit_all_sf=True,               # the original fitted every SF, not just the read one
    impute_dgw_center=False,       # the original imputed no aperture centre
    ssi_tuning_fit_includes_baseline=True,   # inconsistent with its own peak selection
    lifetime_sparseness_over="trials",   # the original flattened every trial
    zero_denominator_nan=False,          # the original's osi/dsi returned 0, not NaN
)
