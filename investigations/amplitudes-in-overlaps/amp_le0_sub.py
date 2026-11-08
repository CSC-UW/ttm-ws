"""For the sampled spikes of amp_le0_snippets.csv: the amplitude SpikeInterface reads (main channel, spike frame plus
peak shift) once the other sorted spikes of the tetrode near it are subtracted. The snippet (3 ms before to 5 ms
after the spike, the tetrode's channels) is fitted jointly by least squares with the analyzer's mean template of
every sorted spike of the tetrode whose template reaches it, the spike itself included; the neighbours' fitted
contribution at the amplitude sample is taken off the raw value."""
import sys

import numpy as np
import pandas as pd
import spikeinterface as si
from spikeinterface.core.template_tools import get_template_peak_shift_on_main_channel

OUT = sys.argv[1]
TD = "/nvme/neuropixels/tetrode_data"
S = {"TTM-003-48H": "2026-10-20_09-02-12", "TTM-003-48H-2": "2026-11-03_09-01-38"}
B, A = 90, 150
N = pd.read_csv(f"{OUT}/amp_le0_snippets.csv")
rows = []
for exp, ses in S.items():
    a = si.load_sorting_analyzer(f"{TD}/{ses}/sortings_carry/main/analyzer.zarr")
    s, rec = a.sorting, a.recording
    chans = list(rec.channel_ids)
    ids = list(s.unit_ids)
    mua = np.asarray(s.get_property("is_mua"))
    grp = np.asarray(s.get_property("group")).astype(int)
    rgrp = np.asarray(rec.get_property("group")).astype(int)
    main_id = np.asarray(s.get_property("main_channel_id"))
    text = a.get_extension("templates")
    T = text.get_templates(operator="average").astype(np.float64)
    nb = text.nbefore
    L = T.shape[1]
    shifts = get_template_peak_shift_on_main_channel(a, with_dict=True)
    sv = s.to_spike_vector()
    tet_spikes = {}
    for g in np.unique(grp):
        ks = np.flatnonzero((grp == g) & ~mua)
        m = np.isin(sv["unit_index"], ks)
        tet_spikes[g] = (sv["sample_index"][m], sv["unit_index"][m])
    for _, q in N[N.exp == exp].iterrows():
        k = ids.index(q.unit)
        g = grp[k]
        tet = np.flatnonzero(rgrp == g)
        f = int(q.frame)
        tr = rec.get_traces(start_frame=f - B, end_frame=f + A, channel_ids=[chans[t] for t in tet],
                            return_in_uV=True).astype(np.float64)
        fr, un = tet_spikes[g]
        lo, hi = np.searchsorted(fr, f - B - (L - nb)), np.searchsorted(fr, f + A + nb)
        cols, labels, is_self = [], [], []
        for fj, kj in zip(fr[lo:hi], un[lo:hi]):
            col = np.zeros_like(tr)
            st = fj - nb - (f - B)
            a0, a1 = max(st, 0), min(st + L, tr.shape[0])
            if a1 <= a0:
                continue
            col[a0:a1] = T[kj][a0 - st:a1 - st][:, tet]
            cols.append(col.ravel())
            labels.append((int(fj), int(kj)))
            is_self.append(fj == f and kj == k)
        X = np.stack(cols, axis=1)
        coef = np.linalg.lstsq(X, tr.ravel(), rcond=None)[0]
        mcl = list(tet).index(chans.index(main_id[k]))
        t0 = B + shifts[ids[k]]
        nbr = np.array([not x for x in is_self])
        contrib = (X[:, nbr] @ coef[nbr]).reshape(tr.shape)[t0, mcl] if nbr.any() else 0.0
        raw = tr[t0, mcl]
        sign = 1.0 if abs(q.amp - raw) < abs(q.amp + raw) else -1.0
        self_coef = coef[np.array(is_self)][0] if any(is_self) else np.nan
        rows.append({"exp": exp, "unit": q.unit, "kind": q.kind, "frame": f, "amp": q.amp, "raw_flipped": sign * raw, "nbr_flipped": sign * contrib,
                     "corrected": sign * (raw - contrib), "self_scale": self_coef, "n_nbr": int(nbr.sum()),
                     "n_nbr_other_unit": int(sum(1 for (fj, kj), sf in zip(labels, is_self) if not sf and kj != k))})
    print(exp, "done", flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/amp_le0_sub.csv", index=False)
