import time, numpy as np, spikeinterface as si, spikeinterface.preprocessing as spre
from ttm.reference import aligned_median_reference
S="/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12"
z=np.load(f"{S}/sortings_seed42_pcafix/mp_aligned/wire_timing.npz")
timing={k:z[k] for k in ("starts","shifts","gains","filled")}
raw=si.read_zarr(f"{S}/2026-10-20_09-02-12.blosc-zstd.zarr"); raw.reset_times()
pp=aligned_median_reference(spre.bandpass_filter(raw,freq_min=300.0,freq_max=6000.0,dtype="float32"),timing,dtype="float32")
crop=pp.frame_slice(0,4*1800*30000)
q=spre.astype(spre.clip(crop,a_min=-32768,a_max=32767),dtype="int16")
t=time.time(); q.save(format="binary",folder=f"{S}/sortings_carry/int16test/binary",dtype="int16",n_jobs=64,chunk_duration="1s",progress_bar=False,overwrite=True)
print(f"int16 binary in {(time.time()-t)/60:.1f} min",flush=True)
f=si.load(f"{S}/sortings_seed42_pcafix/mp_aligned/binary"); g=si.load(f"{S}/sortings_carry/int16test/binary")
for s in (1000*30000, 3000*30000, 6000*30000):
    a=f.get_traces(start_frame=s,end_frame=s+30000); b=g.get_traces(start_frame=s,end_frame=s+30000).astype(np.float32)
    print(f"at {s//30000} s: max |int16 - float32| {np.abs(a-b).max():.3f} counts, rms {np.sqrt(((a-b)**2).mean()):.3f}",flush=True)
