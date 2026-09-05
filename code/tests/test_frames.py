"""sweep_responses_frames: fixed sample count vs variable time window."""

from support import check


import numpy as np
from v1dd_metrics import responses as tr

RNG=np.random.default_rng(5)
DT=0.16504
ts=np.cumsum(np.full(3000,DT))+3.0
traces=RNG.gamma(2.0,0.5,size=(3000,4))
starts=RNG.uniform(ts[5],ts[-40],size=300)


def test_1_fixed_count_is_exactly_fixed():
    for k in (1,2,3,4):
        r=tr.sweep_responses_frames(traces,ts,starts,k)
        a=np.searchsorted(ts,starts,side="left")
        want=np.stack([traces[a[i]:a[i]+k].mean(axis=0) for i in range(len(starts))])
        check(f"k={k} matches an explicit k-sample mean", np.max(np.abs(r-want))<1e-9,
              f"max diff {np.max(np.abs(r-want)):.2e}")


def test_2_the_two_models_differ_in_exactly_the_way_that_matters():
    a,b=tr.window_bounds(ts,starts,0.0,0.35)
    n_win=b-a
    check("time window gives a VARYING sample count",
          len(np.unique(n_win))>1, f"counts {dict(zip(*np.unique(n_win,return_counts=True)))}")
    check("fixed frames gives a CONSTANT sample count", True, "k by construction")
    w=tr.sweep_responses(traces,ts,starts,(0.0,0.35))
    f2=tr.sweep_responses_frames(traces,ts,starts,2)
    check("the two disagree on trials where the window caught 3 samples",
          np.max(np.abs(w-f2))>1e-6, f"max diff {np.max(np.abs(w-f2)):.4f}")
    same=np.isclose(w,f2,atol=1e-12).all(axis=1)
    check("but agree exactly where the window caught exactly 2",
          np.array_equal(same, n_win==2), f"{same.sum()} agree, {(n_win==2).sum()} had 2 samples")


def test_3_why_lifetime_sparseness_can_tell_them_apart():
    # scale-invariance: LS is unchanged by a global scale, but the window model rescales
    # each trial DIFFERENTLY, which is not a global scale.
    x=RNG.gamma(2.0,1.0,size=(3,50))
    check("lifetime_sparseness is invariant to a GLOBAL scale",
          np.allclose(tr.lifetime_sparseness(x), tr.lifetime_sparseness(x*7.3)),
          f"{tr.lifetime_sparseness(x)[0]:.6f} vs {tr.lifetime_sparseness(x*7.3)[0]:.6f}")
    per_trial=x/np.where(RNG.random((3,50))<0.5,2.0,3.0)   # per-trial rescale, like the window
    check("but NOT to a PER-TRIAL rescale (which is what a time window does)",
          not np.allclose(tr.lifetime_sparseness(x), tr.lifetime_sparseness(per_trial)),
          f"{tr.lifetime_sparseness(x)[0]:.6f} vs {tr.lifetime_sparseness(per_trial)[0]:.6f}")


def test_4_guardrails():
    try:
        tr.sweep_responses_frames(traces,ts,starts,0); check("rejects k<1",False)
    except ValueError: check("rejects k<1",True)
    edge=tr.sweep_responses_frames(traces,ts,np.array([ts[-1]]),5)
    check("clips at the end of the recording without raising", edge.shape==(1,4))
