"""NWB reading: planes, stimulus slicing, spontaneous block, LSN template.

Checked against a mock NWB shaped like the real asset. The original file's
schema_report and checkpoint sections are not here: those validation modules were
not carried into this repo.
"""

from support import check
from v1dd_metrics import nwb as vn
nwbmod = vn


def test_nwb():

    from pathlib import Path
    import json
    import shutil
    import tempfile

    from types import SimpleNamespace

    import numpy as np
    import pandas as pd


    # ------------------------------------------------------------------ mock NWB
    RNG = np.random.default_rng(0)
    N_FRAMES, DT = 6000, 0.16374
    TS = np.cumsum(np.full(N_FRAMES, DT)) + 12.0
    STIM_COLS = nwbmod.STIMULUS_TABLE_COLUMNS


    class Img:
        def __init__(self, d): self.data = d


    class Images:
        def __init__(self, arr): self.images = {str(i): Img(a) for i, a in enumerate(arr)}


    class Table:
        def __init__(self, df): self._df = df
        def to_dataframe(self): return self._df.copy()


    class ImagingPlane:
        """Depth lives in `location` as free text ("50 um"), not in origin_coords."""
        def __init__(self, depth_um): self.location = f"{depth_um} um"


    class Series:
        def __init__(self, data, ts, rois, depth_um=50):
            self.data, self.timestamps, self.rois = data, ts, Table(rois)
            self.rois.table = SimpleNamespace(imaging_plane=ImagingPlane(depth_um))


    class Subject:
        def __init__(self, subject_id="409828"):
            self.subject_id, self.species, self.sex = subject_id, "Mus musculus", "male"


    class Module:
        def __init__(self, d): self.data_interfaces = d
        def __getitem__(self, k): return self.data_interfaces[k]


    def make_stim_table(nm_frames=100):
        rows = []

        def blank_row(**kw):
            r = {c: np.nan for c in STIM_COLS}
            r.update(kw)
            return r

        t = 60.0
        # drifting gratings: 12 dirs x 2 sf x 8 trials, plus 6 blank sweeps
        for name, (caz, cel) in (("drifting_gratings_full", (0.0, 0.0)),
                                 ("drifting_gratings_windowed", (3.2, 2.8))):
            combos = [(d, s) for d in range(0, 360, 30) for s in (0.04, 0.08)] * 8
            RNG.shuffle(combos)
            for d, s in combos:
                rows.append(blank_row(stim_name=name, start_time=t, stop_time=t + 1.985,
                                      temporal_frequency=1.0, spatial_frequency=s,
                                      direction=float(d), center_azimuth=caz,
                                      center_elevation=cel))
                t += 3.0
            for _ in range(6):                       # blank/gray sweeps
                rows.append(blank_row(stim_name=name, start_time=t, stop_time=t + 1.985,
                                      center_azimuth=caz, center_elevation=cel))
                t += 3.0
        # natural images: 118 x 8 ; natural_images_12: 12 x 40
        for name, ids, reps in (("natural_images", range(118), 8),
                                ("natural_images_12", [2, 4, 5, 6, 9, 23, 27, 29, 32, 47, 62, 68], 40)):
            seq = [i for i in ids for _ in range(reps)]
            RNG.shuffle(seq)
            for k, i in enumerate(seq):
                rows.append(blank_row(stim_name=name, start_time=t, stop_time=t + 0.25,
                                      image_order=float(k % len(list(ids))), image_index=float(i)))
                t += 0.317
        # natural movie: frames 0..nm_frames-1, 9 repeats
        for _ in range(9):
            for f in range(nm_frames):
                rows.append(blank_row(stim_name="natural_movie", start_time=t,
                                      stop_time=t + 1 / 30, frame=float(f)))
                t += 1 / 30
        # locally sparse noise
        for f in RNG.integers(0, 1593, size=400):
            rows.append(blank_row(stim_name="locally_sparse_noise", start_time=t,
                                  stop_time=t + 0.25, frame=float(f)))
            t += 0.25
        rows.append(blank_row(stim_name="spontaneous", start_time=t, stop_time=t + 300.0))

        df = pd.DataFrame(rows)[STIM_COLS]
        df["stimulus_condition_id"] = pd.factorize(
            df["stim_name"].astype(str) + "_" + df["direction"].astype(str)
            + "_" + df["spatial_frequency"].astype(str) + "_" + df["frame"].astype(str)
            + "_" + df["image_index"].astype(str))[0].astype(float)
        return df.sort_values("start_time").reset_index(drop=True)


    STIM = make_stim_table()


    def make_epochs(stim):
        rows = []
        for name, g in stim.groupby("stim_name", sort=False):
            s, e = g["start_time"].min(), g["stop_time"].max()
            rows.append({"stim_name": name, "start_time": s, "stop_time": e,
                         "duration": e - s})
        return pd.DataFrame(rows).sort_values("start_time").reset_index(drop=True)


    def make_lsn(shape=(16, 28), uniform=True, n=1593):
        base = RNG.choice([0, 127, 255], size=(n, shape[0] // 2, shape[1] // 2))
        up = np.repeat(np.repeat(base, 2, axis=1), 2, axis=2).astype(np.uint8)
        if not uniform:
            up[0, 0, 0] = 42                     # break one 2x2 block
        return up


    def make_nwb(stim=STIM, lsn_uniform=True, n_planes=2, with_stim_table=True):
        proc = {}
        for p in range(n_planes):
            n_rois = 20 + p
            rois = pd.DataFrame({
                "column": 1, "volume": 3, "plane": p,
                "roi": np.arange(n_rois) * 2,
                "pika_roi_id": [f"x_{i:04d}" for i in range(n_rois)],
                "pika_roi_confidence": np.linspace(0.2, 0.99, n_rois),
                "is_soma": np.linspace(0.2, 0.99, n_rois) > 0.5,
                "pixel_mask": [[[0, 0, 1.0]]] * n_rois,
            })
            rois.index.name = "id"
            d = RNG.gamma(2.0, 0.3, size=(N_FRAMES, n_rois))
            depth = 50 + 16 * p          # the asset's lattice: 6 planes, 16 um apart
            proc[f"plane-{p}"] = Module({"dff": Series(d, TS, rois, depth),
                                         "events": Series(d * 0.5, TS, rois, depth)})
        rs = RNG.gamma(1.0, 1.0, size=40000)
        proc["behavior"] = Module({"running_speed": Series(rs, np.linspace(0, 3600, 40000), None)})

        intervals = {"epochs": Table(make_epochs(stim))}
        if with_stim_table:
            intervals["stimulus_table"] = Table(stim)

        class NWB:
            processing = proc
            subject = Subject()
            stimulus = {"locally_sparse_noise": Images(make_lsn(uniform=lsn_uniform)),
                        "natural_images": Images(RNG.integers(0, 255, (118, 4, 4))),
                        "natural_movie": Images(RNG.integers(0, 255, (100, 4, 4)))}
        NWB.intervals = intervals
        return NWB()


    # ------------------------------------------------------------------ tests
    print("[1] load_plane")
    nwb = make_nwb()
    pl = nwbmod.load_plane(nwb, 0)
    check("PlaneData shapes", pl.traces["dff"].shape == (N_FRAMES, 20), str(pl.traces["dff"].shape))
    check("keeps NWB orientation (n_frames, n_rois)", pl.traces["dff"].shape[0] == len(pl.timestamps))
    check("column/volume/plane parsed", (pl.column, pl.volume, pl.plane) == (1, "3", 0),
          f"{(pl.column, pl.volume, pl.plane)}")
    check("volume normalised to str", isinstance(pl.volume, str))
    check("roi ids come from the table, not positions", pl.roi[1] == 2, str(pl.roi[:4]))
    check("is_valid == pika confidence > 0.5",
          np.array_equal(pl.is_valid, pl.roi_table["pika_roi_confidence"].to_numpy() > 0.5))
    check("dt/fs derived", abs(pl.dt - DT) < 1e-9 and abs(pl.fs - 1 / DT) < 1e-6)
    check("pixel_mask dropped", "pixel_mask" not in pl.roi_table.columns)
    plv = nwbmod.load_plane(nwb, 0, valid_only=True)
    check("valid_only filters traces and table",
          plv.traces["dff"].shape[1] == plv.n_rois == int(pl.is_valid.sum()))

    print("\n[2] stimulus_trials and blank-sweep detection")
    st = nwbmod.load_stimulus_table(nwb)
    dgf, is_blank = nwbmod.stimulus_trials(st, "drifting_gratings_full", nwbmod.DG_PARAM_COLUMNS)
    check("selects one stimulus", set(dgf["stim_name"]) == {"drifting_gratings_full"})
    check("chronological", bool(np.all(np.diff(dgf["start_time"]) > 0)))
    check("blank count correct", int(is_blank.sum()) == 6, f"{is_blank.sum()}")
    naive = dgf.isna().any(axis=1).to_numpy()
    check("the naive 'any NaN' test would flag every row (why we restrict columns)",
          naive.all(), f"{naive.sum()}/{len(naive)}")
    check("restricted test flags far fewer than naive", is_blank.sum() < naive.sum())
    nm, nm_blank = nwbmod.stimulus_trials(st, "natural_movie")
    check("no param columns -> no blanks", nm_blank.sum() == 0)

    print("\n[3] spontaneous_block")
    s0, s1 = nwbmod.spontaneous_block(nwb)
    check("returns a 300 s block", abs((s1 - s0) - 300.0) < 1e-6, f"{s1 - s0:.2f}")
    try:
        nwbmod.spontaneous_block(nwb, which=3)
        check("raises for a missing block index", False)
    except IndexError:
        check("raises for a missing block index", True)

    print("\n[4] load_lsn_template")
    lsn = nwbmod.load_lsn_template(nwb)
    check("16x28 reduces to 8x14", lsn["images"].shape[1:] == (8, 14), str(lsn["images"].shape))
    check("blocks reported uniform", lsn["blocks_uniform"] is True)
    check("reduced flag set", lsn["reduced"] is True)
    check("degrees centred on screen",
          abs(lsn["azimuths"][0] + lsn["azimuths"][-1]) < 1e-9
          and abs(lsn["altitudes"][0] + lsn["altitudes"][-1]) < 1e-9)
    check("azimuth span matches 14 x 9.3 grid",
          abs((lsn["azimuths"][-1] - lsn["azimuths"][0]) - 13 * 9.3) < 1e-9,
          f"{lsn['azimuths'][0]:.2f}..{lsn['azimuths'][-1]:.2f}")
    bad = nwbmod.load_lsn_template(make_nwb(lsn_uniform=False))
    check("non-uniform blocks reported, not raised",
          bad["blocks_uniform"] is False and "error" in bad and not bad["reduced"])
    check("failure keeps the native template for inspection", bad["images"].shape[1:] == (16, 28))
