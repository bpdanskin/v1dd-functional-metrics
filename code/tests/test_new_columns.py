"""The two columns added in this refactor, and the helper behind one of them."""

import numpy as np
import pytest

from v1dd_metrics.families.roi_quality import ROI_SUMMARY_COLUMNS, _pearson_columns
from v1dd_metrics.schema import OUTPUT_COLUMNS


def test_pearson_columns_matches_numpy_on_clean_data():
    rng = np.random.default_rng(0)
    y = rng.normal(size=500)
    x = np.stack([y * 2.0 + 1.0,            # perfectly correlated
                  -y,                        # perfectly anticorrelated
                  rng.normal(size=500)], axis=1)
    got = _pearson_columns(x, y)
    want = np.array([np.corrcoef(x[:, i], y)[0, 1] for i in range(x.shape[1])])
    assert np.allclose(got, want, rtol=1e-12)
    assert got[0] == pytest.approx(1.0)
    assert got[1] == pytest.approx(-1.0)


def test_pearson_columns_drops_non_finite_per_roi():
    y = np.arange(10.0)
    x = np.stack([y.copy(), y.copy()], axis=1)
    x[0, 0] = np.nan                        # one ROI loses a frame, the other does not
    got = _pearson_columns(x, y)
    assert got[0] == pytest.approx(1.0), "dropping a frame must not break the correlation"
    assert got[1] == pytest.approx(1.0)


def test_pearson_columns_is_nan_without_variance():
    y = np.arange(10.0)
    flat = np.zeros((10, 1))
    assert np.isnan(_pearson_columns(flat, y)[0]), "a flat trace has no correlation"
    assert np.isnan(_pearson_columns(y[:, None], np.zeros(10))[0])


def test_run_corr_dff_is_published():
    """It exists to be read, so it has to reach the output schema."""
    assert "run_corr_dff" in ROI_SUMMARY_COLUMNS
    assert "run_corr_dff" in OUTPUT_COLUMNS["roi_summary"]


def test_n_trials_at_pref_is_published_for_every_natural_family():
    """The denominator a binomial tail p-value needs, for all three natural stimuli."""
    for fam in ("natural_images", "natural_images_12", "natural_movie"):
        assert "n_trials_at_pref" in OUTPUT_COLUMNS[fam], fam
