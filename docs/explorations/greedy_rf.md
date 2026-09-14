# Greedy pixelwise RF vs the fraction-threshold method

Our historical receptive-field method thresholds, per pixel, the fraction of that pixel's
LSN presentations whose response beats the ROI's spontaneous 95th percentile
(`rf_frac_thresh = 0.25`), with no multiple-comparisons control. It fragments badly, which
is why RF area was never shippable.

The **greedy pixelwise RF** (Dan Millman, from Jun Zhuang's `v1dd_physiology` package —
the upstream of `allen_v1dd` and the white paper) takes a different route: per-pixel
stimulus-triggered average of L0 events, a per-pixel bootstrap p-value (resample sweeps
with replacement, 5000×), and a **Holm-Šidák** correction across all `2·n_pixels` tests.
Detect where the corrected p < α. The multiple-comparisons correction is exactly what our
method lacks.

Both were compared on the true **8×14** grid (the NWB's 16×28 template is a verified 2×2
upsample of an 8×14 stimulus — no finer information), with matched events and a 2-frame
window, so the comparison isolates the significance method.

Harness: `code/validation/greedy_rf_comparison.py` (`--synthetic` for a local ground-truth
check; default runs real sessions in the capsule and sweeps α). Per-ROI output:
`greedy_rf_comparison.csv`.

## Result — 25 sessions, 39,407 ROIs (ON subfield; OFF nearly identical)

| method | detection | single-component | mean area |
|---|---|---|---|
| ours (frac ≥ 0.25) | 6.6% | **81%** | 3.3 px |
| greedy α=0.05 (sensitive) | 11.8% | 95% | 1.2 px |
| greedy α≤0.04 (strict) | 8.6% | **97%** | 1.2 px |

![Greedy vs fraction RF across 25 sessions: (A) detection rate, (B) single-component
fraction, (C) RF-area distribution, (D) contiguity vs SNR](../figures/greedy_rf_comparison.png)

*(A) both greedy variants detect at least as many ON fields as the fraction method. (B)
greedy fields are far more often a single connected component. (C) the fraction method has a
long RF-area tail out to 100+ pixels — the fragments that make its area unusable — while
greedy fields are compact (1–2 px). (D) the fraction method's contiguity collapses to 73% in
the high-SNR tertile, exactly where it detects most; greedy holds 94–98% at every SNR.*

- **De-fragmentation is real and confound-free.** On the cells *both* methods detect,
  single-component rises 84% → 96% and mean largest-component-fraction 0.92 → 0.98. Greedy
  cleans up the same fields, not just a different set.
- **The fraction method's fragmentation is where it detects most.** Stratified by SNR, our
  high-SNR tertile detects 12.5% but is only 73% single-component; greedy holds 94–98%
  single-component at every SNR. Its mean-area 10-px tail (fragments) collapses to ~1.2 px.
- **Detection is comparable-or-better** at both operating points, and the tradeoff is
  **stable across depth and SNR** — so a single global α is justified per variant (no
  per-cell/per-depth tuning, which would make `has_rf` incomparable across cells).
- **The α sweep aliases at 5000 bootstraps**: α = 0.04/0.03/0.02/0.01 are identical (only
  p=0 pixels survive below 0.05). Two operating points, not a continuum; finer tuning would
  need ~10–50× more bootstraps and isn't warranted.
- **Ensembles don't help detection.** Union (ours OR greedy) = 10.6% at lower purity than
  greedy α=0.05 alone (11.8%, 95%). Consensus (both) = 4.5% at 96% — available as a
  high-confidence flag if wanted, but greedy stand-alone dominates.

The synthetic mode (compact planted RFs) confirms the implementation and shows the
sensitivity/specificity knob, but cannot reproduce real-data fragmentation — that is what
the 25-session run demonstrated.

## Decision

Greedy replaces the fraction method as the default (`rf_method="greedy"`), shipping two
variants: **strict** (α=0.01, canonical columns) as the default, and **sensitive** (α=0.05,
`_a05` columns) for downstream analyses that want more yield. The fraction method stays
available behind `rf_method="fraction"` for `REFERENCE_CONFIG` reproduction and array
replay. See [../families/receptive_fields.md](../families/receptive_fields.md).
