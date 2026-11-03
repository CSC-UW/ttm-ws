import pathlib, time, numpy as np, pandas as pd, spikeinterface as si
from ttm import bad_periods as bp
S="/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12"
R=pathlib.Path(__file__).parent/"results"
rec=si.load(f"{S}/sortings_seed42_pcafix/mp_aligned/binary")
t=time.time(); st=bp.bin_stats(rec,bin_s=0.1,n_jobs=48); print(f"scan {time.time()-t:.0f} s", flush=True)
f=bp.bin_features(st,bin_s=1.0); b=bp.baselines(f)
names={int(g):f"TT{int(g)+1}" for g in f["group_ids"]}
tab=bp.detect_bad_periods(f,b,group_names=names); tab.to_csv(R/"bad_periods_aligned_table.csv",index=False)
old=pd.read_csv(f"{S}/bad_periods/global_median/bad_periods.csv")
dur=float(f["t"][-1]+f["bin_s"])
va=bp.valid_seconds(tab,dur,f["group_ids"]); vo=bp.valid_seconds(old,dur,f["group_ids"])
rows=[]
for g in f["group_ids"]:
    g=int(g)
    def kinds(t):
        s=t[t.group==g]; return {k:round(float((s.end_s-s.start_s)[s.kind.str.split('+').str[0]==k].sum())) for k in ("loud","common","global")}
    rows.append({"tetrode":names[g],"global_ref_s":round(dur-vo[g]),"aligned_ref_s":round(dur-va[g]),**{f"g_{k}":v for k,v in kinds(old).items()},**{f"a_{k}":v for k,v in kinds(tab).items()}})
d=pd.DataFrame(rows); print(d.to_string(index=False)); d.to_csv(R/"bad_periods_aligned.csv",index=False)
