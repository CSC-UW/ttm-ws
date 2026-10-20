import pathlib

import pytest
from ttm import config


@pytest.fixture(autouse=True)
def _isolate(monkeypatch, tmp_path):
    for k in list(__import__("os").environ):
        if k.startswith(config.ENV_PREFIX):
            monkeypatch.delenv(k)
    monkeypatch.setenv("TTM_CONFIG", str(tmp_path / "config.toml"))
    config.get.cache_clear()
    yield
    config.get.cache_clear()


def test_missing_config_raises(tmp_path: pathlib.Path) -> None:
    with pytest.raises(config.ConfigError, match="data_root"):
        config.load()


def test_file_values(tmp_path: pathlib.Path) -> None:
    (tmp_path / "config.toml").write_text(
        'data_root = "/data"\ncache_dir = "~/c"\n[remotes]\narchive = "/arch"\nserver = "box:/data"\n'
    )
    cfg = config.get()
    assert cfg.data_root == pathlib.Path("/data")
    assert cfg.cache_dir == pathlib.Path("~/c").expanduser()
    assert cfg.remote("archive") == "/arch"
    assert cfg.remote("server") == "box:/data"
    with pytest.raises(config.ConfigError, match="TTM_REMOTE_SCRATCH"):
        cfg.remote("scratch")


def test_env_overrides_file(monkeypatch, tmp_path: pathlib.Path) -> None:
    (tmp_path / "config.toml").write_text('data_root = "/data"\n[remotes]\narchive = "/arch"\n')
    monkeypatch.setenv("TTM_DATA_ROOT", "/other")
    monkeypatch.setenv("TTM_REMOTE_ARCHIVE", "/arch2")
    monkeypatch.setenv("TTM_REMOTE_SCRATCH", "/scratch")
    cfg = config.load()
    assert cfg.data_root == pathlib.Path("/other")
    assert cfg.remotes == {"archive": "/arch2", "scratch": "/scratch"}
    assert cfg.cache_dir is None


def test_env_only(monkeypatch) -> None:
    monkeypatch.setenv("TTM_DATA_ROOT", "/only")
    assert config.load().data_root == pathlib.Path("/only")
