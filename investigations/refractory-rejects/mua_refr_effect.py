"""Multi-unit counts per tetrode with and without the carry's refractory rejects, from the saved outputs of
the runs argv[2:]; one table per run goes to argv[1]."""
import sys
import time

import numpy as np
import pandas as pd
import spikeinterface as si

from ttm import analyzer as an

FS = 30000.0
for run in sys.argv[2:]:
    t0 = time.time()
    rejects = an.refractory_rejects(f"{run}/carry")
    t1 = time.time()
    p = an.AnalyzerParams()
    old, _ = an.multi_unit_events(f"{run}/multi_unit", p)
    new, tab = an.multi_unit_events(f"{run}/multi_unit", p, rejects=rejects, fs=FS)
    s = si.load(f"{run}/sorting")
    g = np.asarray(s.get_property("group")).astype(int)
    c = s.count_num_spikes_per_unit(outputs="array")
    sorted_per = pd.Series(c).groupby(g).sum()
    tab["sorted"] = tab.group.map(sorted_per).fillna(0).astype(int)
    tab["refractory_share_of_sorted_%"] = 100 * tab.refractory / tab.sorted.replace(0, np.nan)
    tab.to_csv(f"{sys.argv[1]}/multi_unit_counts_{run.split('/')[4]}.csv", index=False)
    print(f"\n{run.split('/')[4]}: rejects read in {t1 - t0:.0f} s")
    print(tab[["tetrode", "sorted", "admitted", "refractory", "admitted_near_refractory", "events",
               "refractory_share_of_sorted_%"]].round(2).to_string(index=False))
    print("totals: sorted", int(tab.sorted.sum()), "detected admitted", int(tab.admitted.sum()), "rejects",
          int(tab.refractory.sum()), "detected giving way", int(tab.admitted_near_refractory.sum()),
          "multi-unit events", int(tab.events.sum()), "| old multi-unit", sum(v.size for v in old.values()))
