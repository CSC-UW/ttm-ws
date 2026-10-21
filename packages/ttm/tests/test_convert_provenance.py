import json

import numpy as np
import spikeinterface as si
import spikeinterface.core as sc

from ttm.spikeinterface import convert_recording


def test_store_provenance_has_no_time_vector(tmp_path):
    inner = sc.generate_recording(num_channels=8, sampling_frequency=30000.0, durations=[2.0], seed=0)
    times = 5.0 + np.arange(inner.get_num_frames()) / 30000.0
    inner.set_times(times, with_warning=False)
    rec = inner.select_channels(inner.channel_ids)
    rec.set_times(times + 1.0, with_warning=False)
    assert "times_kwargs" in json.dumps(sc.core_tools.check_json(rec.to_dict(recursive=True)))
    out = tmp_path / "store.zarr"

    saved = convert_recording(rec, out, compressor="blosc-zstd", n_jobs=1)

    attrs = json.loads((out / ".zattrs").read_text())
    assert attrs["provenance"]["class"].endswith("ChannelSliceRecording")
    assert "times_kwargs" not in json.dumps(attrs["provenance"])
    assert (out / ".zattrs").stat().st_size < 100_000
    assert rec._serializability["json"] is True
    reloaded = si.read_zarr(str(out))
    np.testing.assert_allclose(reloaded.get_times(), times + 1.0)
    np.testing.assert_allclose(saved.get_times(), times + 1.0)
    np.testing.assert_array_equal(reloaded.get_traces(), rec.get_traces())
