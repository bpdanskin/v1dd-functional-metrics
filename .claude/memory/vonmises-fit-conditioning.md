---
name: vonmises-fit-conditioning
description: "Why the von Mises tuning fit is expensive: k is unbounded inside an exponential. Three hypotheses tested and refuted, and what to try instead. Covers the data-derived p0 and x_scale='jac'."
metadata:
  node_type: memory
  type: project
  modified: 2026-09-03
---

Explored 2026-09-03 while looking at the drifting-gratings runtime, which is ~96 % of a
full run. **Nothing here was adopted.** The measurements are the point; write them up in
`docs/families/drifting_gratings.md` in P5.

Setup: 200-300 windowed-grating curves from the reference asset, one fit each, model
evaluations counted by wrapping `vonmises_two_peak`. **Count evaluations, not wall time** --
the first timing run put `fixed` first and it absorbed process warm-up, which is how an
unreliable "2.8x slower" figure got recorded here originally.

## Three hypotheses, all refuted

**1. A data-derived `p0` will be faster.** The fork deferred this on the reasoning that
`scale_1 = 0.1` starts ~50x above the median event peak of 2.07e-3, so fits exhaust their
first evaluation budget and retry. Measured, it is **2.1x more expensive** (1617 against
771 evaluations per curve) and converges on 199/200 against 200/200.

**2. It is slow because the derived guess starts on a bound.** Plausible: every parameter
is bounded below at zero, and on sparse events the trough is often exactly zero -- `b`
starts at 0 on 22 % of curves, `scale_2` on 9.8 %, `x0` at 0 deg on 10 %. **Refuted**: a
variant nudged strictly inside every bound needed 2293 evaluations against the plain
derived guess's 2305. Identical.

**3. It is a scaling problem, fixable with `x_scale="jac"`.** Trust-region-reflective
measures its region in absolute parameter units, so parameters of size 1e-4 should be
badly served by the default `x_scale=1.0`. **Refuted, emphatically** -- Jacobian scaling is
far worse:

| start | `x_scale` | evals/curve | converged |
|---|---|---|---|
| fixed | default | **771** | 200/200 |
| derived | default | 1617 | 199/200 |
| fixed | `"jac"` | 7352 | 190/200 |
| derived | `"jac"` | 9163 | 188/200 |

## What the failure actually reveals

`x_scale="jac"` raises `RuntimeWarning: overflow encountered in exp` inside
`vonmises_two_peak`, plus overflows in TRF's own `scale_inv` and cost evaluation.

**The model has `exp(k * cos(theta - x0))` with `k` bounded `(0, inf)` -- unbounded above.**
Scaling by the Jacobian lets the optimizer take large steps in `k`, `exp(k)` overflows, the
residual goes non-finite and the trust region collapses. The fixed `p0` with `k = 1` and
absolute scaling keeps the search in a safe region more or less by luck.

So the conditioning problem is not amplitude scale at all. It is an unbounded parameter
inside an exponential.

## The two starting points find equally good minima

Over 300 curves, the derived guess reaches a lower sum-of-squares on **59.2 %** and a
higher one on **33.3 %**, median SSE ratio **0.9993** -- 0.07 % better on average, at twice
the cost. That spread is the signature of a rough objective where no single starting point
dominates, which is itself an argument against spending more effort on the starting point.

## Worth trying next, in order

1. **Bound `k` above.** 12 directions is 30 deg sampling; a concentration beyond ~50
   describes a peak narrower than the stimulus can resolve. Bounding `k` removes the
   overflow regime entirely and shrinks the search space. Cheapest change, highest
   expected value, and directly testable against SSE.
2. **Fit in normalised units** -- divide each curve by its peak, fit, multiply back. The
   fit then does not care whether it is handed events or dF/F, which is the real fix for
   the scale mismatch below.
3. **Multi-start instead of a better single start**, keeping the lowest SSE. The 59/33
   split says this buys real fit quality; it costs 2-3x, which is the opposite trade from
   the one originally sought, but it is the one the surface supports.
4. **Ask whether the fit is needed.** It has exactly one consumer, `ssi_tuning_fit`, and
   half of every fit was already discarded before speedup 1. A non-parametric estimate of
   the response at the preferred direction may be cheaper and more stable.
5. Record `nfev` and convergence counts in provenance, so this is measurable in production
   rather than only in a benchmark.

## Is the mismatch our events-vs-dF/F choice? No -- it is inherited

`scale_1 = 0.1` really is ~50x above event amplitudes. But the original defaulted to
events too: `DriftingGratings.__init__(..., trace_type: str = "events", ...)`. And
`allen_v1dd/stimulus_analysis/fit_utils.py` line 30 shows the author had already tried an
event-scale guess and reverted:

```python
p0 = (0.1, 1, 180, 0.01, 1, 0.001) # Seems to be a reasonable starting point
# p0 = (1e-4, 0.1, 180, 1e-5, 1e-2, 1e-4)
```

That commented-out line is at event scale, in both amplitude *and* concentration. So two
independent attempts at an event-scale starting point have now been made and both rejected
-- theirs by judgement, ours by measurement. **Do not make a third without first bounding
`k`**, which is the thing neither attempt changed.

Related: [[fresh-start-metric-changes]].
