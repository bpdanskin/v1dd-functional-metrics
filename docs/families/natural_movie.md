# Natural movie

Responses to a repeated clip, indexed by frame. Structurally similar to natural images but
with important interpretive differences.

## Goal

Characterise each neuron's response across time in a repeated natural movie. Each movie
frame is a "trial" and each full pass through the movie is a "repeat." On this asset:
3,600 frames × 9 repeats.

## Response window and temporal overlap

The response window is `nm_response_frames × dt` — 3 imaging frames, roughly 0.49 s. Movie
frames are presented at 1/30 s apart, so consecutive "trials" overlap heavily and are
strongly autocorrelated.

That overlap has two consequences:

* **`lifetime_sparseness` over 3,600 × 9 such values is not measuring what its name
  suggests.** Sparseness assumes independent conditions; at this overlap, each response
  shares ~93 % of its window with its neighbours.
* **`pref_img` is approximate by construction.** Activity driven by frame *f* lands inside
  the windows of frames *f−15* through *f*. The reported preferred frame can therefore
  **precede** the frame that actually drove the response by up to half a second, and among
  those overlapping windows the argmax is decided partly by how many imaging samples each
  happens to contain. Read `pref_img` as locating a ~0.5 s neighbourhood, not a frame.

## What differs from natural images

Despite sharing a column set and most of the implementation, two things differ
substantively:

### Conditions are contiguous frames

Movie frames run 0, 1, 2, ..., 3599 — contiguous from zero by construction. The code
asserts this, because the original indexes them positionally and by label interchangeably,
which only works if they are the same. Natural images use `image_index` from a catalog,
which can be sparse.

### `frac_responsive_trials` is not a statistical test

Natural images bootstrap a spontaneous null and report the fraction of trials beating it.
Natural movie does **not** — it reports the fraction of repeats whose mean response at the
preferred frame is strictly greater than zero:

```
frac = mean(response > 0, over repeats)
```

No bootstrap is involved. This makes natural movie the one fully deterministic end-to-end
check against a reference table: the same seed, the same data, the same numbers, every
time. It also means `frac_responsive_trials` answers a different question here than it
does for images or gratings — "how often was there any response?" rather than "how often
did the response beat a statistical null?"

## `z_score`

Same as natural images: preferred response in standard deviations of the multi-trial
spontaneous null (10,000 bootstraps, each averaging `n_repeats` draws). The number of draws
is large but the per-draw count (9 repeats) is small, so this is less expensive than NI12's
40-repeat equivalent.

## No condition-means export

Unlike natural images, natural movie does not export a condition-means array. The full
matrix would be `(n_rois, 3,600)` — large and of limited use given the temporal overlap.
All published columns are reductions of the trial array, and the trial array itself is not
shipped (3,600 × 9 × n_rois in float32).

## How the pipeline runs it

`natural_movie_metrics` in `families/natural_movie.py`. Called once per plane. Returns a
metrics frame only (no auxiliary arrays).

## Columns

| column | meaning |
|---|---|
| `frac_responsive_trials` | fraction of repeats with response > 0 at the preferred frame |
| `lifetime_sparseness` | selectivity over condition means (see caveat above) |
| `pref_img` | preferred frame index (0–3599), −1 if all-NaN |
| `pref_response` | mean response to the preferred frame |
| `z_score` | preferred response in SD of the multi-trial spontaneous null |
| `reliability` | split-half correlation on events |
| `reliability_dff` | split-half correlation on dF/F |
| `n_trials_at_pref` | number of finite repeats at the preferred frame |

Same column set as natural images. The prefix in the wide table (`nm_`) is added by the
orchestrator.

## Sanity checks worth running on a fresh asset

* `pref_img` is in [0, 3599] (or −1 for all-NaN ROIs).
* `n_trials_at_pref` is ≤ the number of movie repeats (9 on this asset).
* `frac_responsive_trials` takes only values that are multiples of `1/n_repeats` — no
  bootstrap noise.
* `lifetime_sparseness` is in [0, 1] but interpretively limited by the temporal overlap.
