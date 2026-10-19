"""Per-machine configuration: where data live on this host.

Read from ``$TTM_CONFIG`` or ``~/.config/ttm/config.toml``::

    data_root = "<local dir holding <session>/ and <subject>/<experiment>/>"
    cache_dir = "<optional cache dir>"

    [remotes]
    archive = "<local mount of the raw Open Ephys acquisitions>"
    server  = "<rclone remote>:<the server's data_root>"

Environment variables override the file: ``TTM_DATA_ROOT``, ``TTM_CACHE_DIR``,
``TTM_REMOTE_<NAME>`` (``TTM_REMOTE_ARCHIVE=...``). There is no default
``data_root``; :func:`get` raises :class:`ConfigError` until one is set.

A remote value is either a local path (a mount on this host) or an rclone
``remote:path``; the caller knows which it expects.
"""

from __future__ import annotations

import dataclasses
import functools
import os
import pathlib
import tomllib

ENV_PREFIX = "TTM_"
DEFAULT_PATH = pathlib.Path.home() / ".config" / "ttm" / "config.toml"

TEMPLATE = """\
data_root = "/path/to/tetrode_data"

[remotes]
archive = "/path/to/raw/acquisitions"
"""


class ConfigError(RuntimeError):
    """Raised when the per-machine configuration is missing or incomplete."""


@dataclasses.dataclass(frozen=True)
class Config:
    data_root: pathlib.Path
    remotes: dict[str, str] = dataclasses.field(default_factory=dict)
    cache_dir: pathlib.Path | None = None
    source: str = ""

    def remote(self, name: str) -> str:
        try:
            return self.remotes[name]
        except KeyError:
            raise ConfigError(
                f"Remote {name!r} is not configured ({self.source}). Add it under "
                f"[remotes] or set {ENV_PREFIX}REMOTE_{name.upper()}."
            ) from None


def config_path() -> pathlib.Path:
    env = os.environ.get(f"{ENV_PREFIX}CONFIG")
    return pathlib.Path(env).expanduser() if env else DEFAULT_PATH


def _env_remotes() -> dict[str, str]:
    prefix = f"{ENV_PREFIX}REMOTE_"
    return {
        k[len(prefix) :].lower(): v for k, v in os.environ.items() if k.startswith(prefix)
    }


def load(path: str | pathlib.Path | None = None) -> Config:
    """Read the config file (if any) and apply environment overrides."""
    path = pathlib.Path(path).expanduser() if path else config_path()
    raw: dict = {}
    sources: list[str] = []
    if path.exists():
        raw = tomllib.loads(path.read_text())
        sources.append(str(path))

    data_root = os.environ.get(f"{ENV_PREFIX}DATA_ROOT") or raw.get("data_root")
    cache_dir = os.environ.get(f"{ENV_PREFIX}CACHE_DIR") or raw.get("cache_dir")
    remotes = {**raw.get("remotes", {}), **_env_remotes()}
    if any(k in os.environ for k in (f"{ENV_PREFIX}DATA_ROOT", f"{ENV_PREFIX}CACHE_DIR")) or _env_remotes():
        sources.append("environment")

    if not data_root:
        raise ConfigError(
            f"No data_root configured. Write {path} with\n\n{TEMPLATE}\n"
            f"or set {ENV_PREFIX}DATA_ROOT."
        )
    return Config(
        data_root=pathlib.Path(data_root).expanduser(),
        remotes={k: str(v) for k, v in remotes.items()},
        cache_dir=pathlib.Path(cache_dir).expanduser() if cache_dir else None,
        source=" + ".join(sources) or "defaults",
    )


@functools.cache
def get() -> Config:
    """The process-wide config, loaded once. ``get.cache_clear()`` to reload."""
    return load()
