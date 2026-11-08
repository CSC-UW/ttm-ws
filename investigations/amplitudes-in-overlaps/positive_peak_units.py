"""Per sorted unit of both analyzers: the analyzer's peak sign, the unit's main channel, its peak offset, and the
main channel's mean template at the peak offset, at its minimum and at its maximum (µV)."""
import sys

import numpy as np
import pandas as pd
import spikeinterface as si
from spikeinterface.core.template_tools import get_template_peak_shift_on_main_channel

OUT = sys.argv[1]
TD = "/nvme/neuropixels/tetrode_data"
S = {"TTM-003-48H": "2026-10-20_09-02-12", "TTM-003-48H-2": "2026-11-03_09-01-38"}
rows = []
for exp, ses in S.items():
    a = si.load_sorting_analyzer(f"{TD}/{ses}/sortings_carry/main/analyzer.zarr")
    s = a.sorting
    mua = np.asarray(s.get_property("is_mua"))
    grp = np.asarray(s.get_property("group")).astype(int)
    main_id = np.asarray(s.get_property("main_channel_id"))
    chans = list(a.channel_ids)
    text = a.get_extension("templates")
    T = text.get_templates(operator="average")
    nb = text.nbefore
    shifts = get_template_peak_shift_on_main_channel(a, with_dict=True)
    for k in np.flatnonzero(~mua):
        u = s.unit_ids[k]
        t = T[k][:, chans.index(main_id[k])]
        rows.append({"exp": exp, "unit": int(u), "tetrode": int(grp[k]) + 1, "peak_sign": a.peak_sign,
                     "main_ch": main_id[k], "peak_offset": int(shifts[u]), "at_offset_uv": float(t[nb + shifts[u]]),
                     "min_uv": float(t.min()), "max_uv": float(t.max())})
d = pd.DataFrame(rows)
d["reads_positive_peak"] = d.at_offset_uv > 0
d.to_csv(f"{OUT}/positive_peak_units.csv", index=False)
print(d[d.reads_positive_peak].to_string(index=False))
