import json
import pathlib

import pytest
from tetrode_analyses import config
from tetrode_analyses import experiment as ep
from tetrode_analyses import subject as sj

WORKSPACE = pathlib.Path(__file__).resolve().parents[3]


@pytest.fixture(autouse=True)
def _config(monkeypatch, tmp_path):
    monkeypatch.setenv("TTM_CONFIG", str(tmp_path / "none.toml"))
    monkeypatch.setenv("TTM_DATA_ROOT", str(tmp_path / "data"))
    monkeypatch.setenv("TTM_REMOTE_ARCHIVE", str(tmp_path / "archive"))
    monkeypatch.delenv("TTM_SUBJECT", raising=False)
    monkeypatch.delenv("TTM_EXPERIMENT", raising=False)
    monkeypatch.setenv("TTM_SUBJECTS_DIR", str(WORKSPACE / "subjects"))
    config.get.cache_clear()
    yield
    config.get.cache_clear()


def _write_subject(root: pathlib.Path, sid="TTM-009", exp="EXP", **extra) -> pathlib.Path:
    d = root / sid.lower()
    (d / "experiments").mkdir(parents=True)
    (d / "subject.toml").write_text(f'id = "{sid}"\ntimezone = "America/Chicago"\n')
    body = 'session = "2026-01-01_00-00-00"\nlightsOn = ["2026-01-01T09:00:00"]\nlightsOff = ["2026-01-01T21:00:00"]\n'
    for k, v in extra.items():
        body += f"{k} = {json.dumps(v)}\n"
    (d / "experiments" / f"{exp}.toml").write_text(body)
    return d / "experiments" / f"{exp}.toml"


def test_committed_manifests_load() -> None:
    subjects = sj.list_subjects()
    assert {s.id for s in subjects} >= {"TTM-001", "TTM-002"}
    nod = sj.find_subject("ttm-001").experiment("TTM-NOD")
    assert nod.session == "2026-05-27_09-07-52"
    assert nod.acquisition == nod.session
    assert len(nod.lightsOn) == 3 and len(nod.lightsOff) == 2
    assert nod.experiment_dir.parts[-2:] == ("TTM-001", "TTM-NOD")
    assert nod.lfp_path.name == "2026-05-27_09-07-52.lfp.zarr"
    assert nod.store_path("wavpack", 2.25).name == "2026-05-27_09-07-52.wavpack-bps2.25.zarr"
    assert nod.acquisition_dir.name == "2026-05-27_09-07-52"
    sd = sj.find_subject("TTM-002").experiment("TTM-002-SD")
    assert sd.acquisition == "2026-09-24_09-03-22/2026-09-24_09-03-22"


def test_manifest_defaults_and_validation(tmp_path: pathlib.Path) -> None:
    path = _write_subject(tmp_path / "subjects")
    exp = sj.Experiment.from_toml(path)
    assert (exp.subject, exp.experiment, exp.timezone) == ("TTM-009", "EXP", "America/Chicago")
    assert exp.acquisition == exp.session and exp.oe_stream_index == 0
    assert exp.novel_objects_start is None and exp.badChannels == []

    bad = _write_subject(tmp_path / "s2", sid="TTM-010", bogus=1)
    with pytest.raises(sj.ManifestError, match="unknown keys"):
        sj.Experiment.from_toml(bad)
    wrong = _write_subject(tmp_path / "s3", sid="TTM-011", subject="TTM-999")
    with pytest.raises(sj.ManifestError, match="subject"):
        sj.Experiment.from_toml(wrong)


def test_resolution_order(monkeypatch, tmp_path: pathlib.Path) -> None:
    root = tmp_path / "subjects"
    path = _write_subject(root)
    monkeypatch.setenv("TTM_SUBJECTS_DIR", str(root))
    with pytest.raises(sj.ManifestError, match="No experiment selected"):
        sj.resolve_experiment()
    assert sj.resolve_experiment(manifest=path).subject == "TTM-009"
    assert sj.resolve_experiment(subject="ttm-009", experiment="EXP").manifest_path == path
    monkeypatch.setenv("TTM_SUBJECT", "TTM-009")
    monkeypatch.setenv("TTM_EXPERIMENT", "EXP")
    assert sj.resolve_experiment().experiment == "EXP"
    with pytest.raises(sj.ManifestError, match="no experiment"):
        sj.resolve_experiment(subject="TTM-009", experiment="NOPE")
    with pytest.raises(sj.ManifestError, match="0 matches"):
        sj.find_subject("TTM-404", root)


def test_experiment_params_round_trip(tmp_path: pathlib.Path) -> None:
    exp = sj.Experiment.from_toml(_write_subject(tmp_path / "subjects", novel_objects_start="2026-01-01T10:00:00",
                                                 novel_objects_end="2026-01-01T12:00:00"))
    params = ep.ExperimentParams.from_experiment(exp, t0_unix=1767279600.0)  # 2026-01-01 09:00 CST
    out = ep.save_experiment_params(params, exp.experiment_params_path)
    assert out == tmp_path / "data" / "TTM-009" / "EXP" / "experiment_params.json"
    raw = json.loads(out.read_text())
    assert set(raw) == set(sj.Experiment.manifest_fields()) | {"t0_unix", "t0_datetime"}
    assert not any(str(v).startswith("/") for v in raw.values() if isinstance(v, str))
    back = ep.load_experiment_params(out)
    assert back == params
    assert ep.dt2t(back, "2026-01-01T09:00:00") == pytest.approx(0.0)
    assert ep.get_deprivation_period(back) == pytest.approx((3600.0, 3 * 3600.0))
    intervals, labels = ep.get_light_dark_periods(back)
    assert labels == ["on"] and intervals[0] == pytest.approx((0.0, 12 * 3600.0))

    out.write_text(json.dumps({**raw, "lfp_zarr": "/old/abs/path"}))
    with pytest.raises(ValueError, match="unknown keys"):
        ep.load_experiment_params(out)


def test_read_t0_from_slice_table_and_sync_messages(tmp_path: pathlib.Path) -> None:
    (tmp_path / "slice_table.csv").write_text("oe_software_time_s,fs\n1779890872.188,30000\n1779991788.126,30000\n")
    assert ep.read_session_t0_unix(tmp_path / "slice_table.csv") == pytest.approx(1779890872.188)
    rec = tmp_path / "acq" / "Record Node 101" / "experiment1" / "recording1"
    rec.mkdir(parents=True)
    (tmp_path / "acq" / "Record Node 101" / "experiment1" / "recording10").mkdir()
    (rec / "sync_messages.txt").write_text(
        "Software Time (milliseconds since midnight Jan 1st 1970 UTC): 1779232191976\n"
    )
    assert ep.first_recording_dir(tmp_path / "acq") == rec
    assert ep.read_session_t0_unix(tmp_path / "acq") == pytest.approx(1779232191.976)
