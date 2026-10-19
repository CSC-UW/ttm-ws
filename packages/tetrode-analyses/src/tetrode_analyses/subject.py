"""Subject and experiment manifests under ``subjects/<id>/``.

``subjects/<dir>/subject.toml`` holds subject-level facts; each
``subjects/<dir>/experiments/<EXPERIMENT>.toml`` holds one experiment's
acquisition identity, light/dark schedule and deprivation window. The manifest
keys are the fields of :class:`Experiment`. Data paths are never in a manifest:
they are derived from the session name and the per-machine
:mod:`tetrode_analyses.config`.

Resolution for command-line tools (:func:`resolve_experiment`):
``--manifest PATH`` > ``--subject/--experiment`` > ``$TTM_SUBJECT`` /
``$TTM_EXPERIMENT`` > error. There is no default subject.
"""

from __future__ import annotations

import dataclasses
import os
import pathlib
import tomllib

from tetrode_analyses import config

SUBJECT_FILENAME = "subject.toml"
EXPERIMENTS_DIRNAME = "experiments"
EXPERIMENT_PARAMS_FILENAME = "experiment_params.json"


class ManifestError(RuntimeError):
    """A manifest is missing, ambiguous, or malformed."""


@dataclasses.dataclass(frozen=True, kw_only=True)
class Subject:
    id: str
    dir: pathlib.Path
    timezone: str
    notes: str = ""

    @classmethod
    def from_toml(cls, path: str | pathlib.Path) -> Subject:
        path = pathlib.Path(path)
        raw = tomllib.loads(path.read_text())
        try:
            return cls(id=raw["id"], timezone=raw["timezone"], dir=path.parent,
                       notes=raw.get("notes", ""))
        except KeyError as e:
            raise ManifestError(f"{path}: missing key {e}") from None

    @property
    def experiments_dir(self) -> pathlib.Path:
        return self.dir / EXPERIMENTS_DIRNAME

    def experiment_names(self) -> list[str]:
        return sorted(p.stem for p in self.experiments_dir.glob("*.toml"))

    def experiment(self, name: str) -> Experiment:
        path = self.experiments_dir / f"{name}.toml"
        if not path.exists():
            raise ManifestError(
                f"Subject {self.id} has no experiment {name!r} "
                f"(available: {self.experiment_names()}) in {self.experiments_dir}"
            )
        return Experiment.from_toml(path, subject=self)


@dataclasses.dataclass(frozen=True, kw_only=True)
class Experiment:
    """One experiment of one subject, as declared in its manifest.

    ``session`` names the derived-products directory ``<data_root>/<session>/``.
    ``acquisition`` is the Open Ephys session directory relative to the remote
    named by ``archive_remote`` (default ``archive``); it defaults to
    ``session``. Datetime fields are naive ISO-8601 strings in ``timezone``,
    the encoding WNE uses.
    """

    subject: str
    experiment: str
    timezone: str
    session: str
    acquisition: str
    archive_remote: str = "archive"
    oe_stream_index: int = 0
    lightsOn: list[str] = dataclasses.field(default_factory=list)
    lightsOff: list[str] = dataclasses.field(default_factory=list)
    novel_objects_start: str | None = None
    novel_objects_end: str | None = None
    badChannels: list[str] = dataclasses.field(default_factory=list)
    badTetrodes: list[str] = dataclasses.field(default_factory=list)
    manifest_path: pathlib.Path | None = dataclasses.field(default=None, compare=False)

    @classmethod
    def manifest_fields(cls) -> tuple[str, ...]:
        return tuple(f.name for f in dataclasses.fields(Experiment) if f.name != "manifest_path")

    @classmethod
    def from_toml(cls, path: str | pathlib.Path, *, subject: Subject | None = None) -> Experiment:
        path = pathlib.Path(path)
        raw = tomllib.loads(path.read_text())
        if subject is None:
            subject = Subject.from_toml(path.parent.parent / SUBJECT_FILENAME)
        raw.setdefault("subject", subject.id)
        raw.setdefault("timezone", subject.timezone)
        raw.setdefault("experiment", path.stem)
        raw.setdefault("acquisition", raw.get("session"))
        unknown = set(raw) - set(cls.manifest_fields())
        if unknown:
            raise ManifestError(f"{path}: unknown keys {sorted(unknown)}")
        missing = {"session"} - set(raw)
        if missing:
            raise ManifestError(f"{path}: missing keys {sorted(missing)}")
        if raw["subject"] != subject.id:
            raise ManifestError(f"{path}: subject {raw['subject']!r} != {subject.id!r}")
        return cls(**raw, manifest_path=path)

    # --- paths on this machine (data_root from the per-machine config) ---------

    @property
    def data_root(self) -> pathlib.Path:
        return config.get().data_root

    @property
    def session_dir(self) -> pathlib.Path:
        return self.data_root / self.session

    @property
    def experiment_dir(self) -> pathlib.Path:
        return self.data_root / self.subject / self.experiment

    @property
    def acquisition_dir(self) -> pathlib.Path:
        """Raw Open Ephys session directory on the ``archive_remote`` (a local mount)."""
        return pathlib.Path(config.get().remote(self.archive_remote)) / self.acquisition

    @property
    def experiment_params_path(self) -> pathlib.Path:
        return self.experiment_dir / EXPERIMENT_PARAMS_FILENAME

    @property
    def slice_table_path(self) -> pathlib.Path:
        return self.session_dir / "slice_table.csv"

    def store_path(self, compressor: str = "blosc-zstd", bps: float | None = None) -> pathlib.Path:
        """Converted 30 kHz Zarr store, e.g. ``<session>.blosc-zstd.zarr``."""
        tag = compressor if bps is None else f"{compressor}-bps{bps}"
        return self.session_dir / f"{self.session}.{tag}.zarr"

    @property
    def lfp_path(self) -> pathlib.Path:
        return self.session_dir / f"{self.session}.lfp.zarr"

    def subsampled_lfp_path(self, rate: int = 125) -> pathlib.Path:
        return self.session_dir / f"{self.session}.lfp.{rate}hz.zarr"

    @property
    def figures_dir(self) -> pathlib.Path:
        return self.experiment_dir / "figures"

    @property
    def idelta_path(self) -> pathlib.Path:
        return self.experiment_dir / "delta.idelta.zarr"

    @property
    def stft_bandpowers_path(self) -> pathlib.Path:
        return self.experiment_dir / "stft_bandpowers.zarr"

    @property
    def emg_path(self) -> pathlib.Path:
        return self.experiment_dir / "synthetic_emg_methods.zarr"

    def sorting_dir(self, sorting: str, subdir: str = "sortings_seed42_pcafix") -> pathlib.Path:
        return self.session_dir / subdir / sorting

    @property
    def viz_bundle_dir(self) -> pathlib.Path:
        return self.session_dir / "viz_bundle"

    @property
    def curation_bundles_dir(self) -> pathlib.Path:
        return self.session_dir / "curation_bundles"


# --- discovery --------------------------------------------------------------

def _looks_like_workspace(path: pathlib.Path) -> bool:
    pyproject = path / "pyproject.toml"
    return (path / "subjects").is_dir() and pyproject.exists() and (
        "[tool.uv.workspace]" in pyproject.read_text()
    )


def subjects_dir() -> pathlib.Path:
    """The ``subjects/`` directory: ``$TTM_SUBJECTS_DIR``, else the enclosing workspace.

    Walks up from the working directory to a uv workspace root holding
    ``subjects/``; falls back to the checkout this module is installed from.
    """
    env = os.environ.get("TTM_SUBJECTS_DIR")
    if env:
        return pathlib.Path(env).expanduser()
    for candidate in [pathlib.Path.cwd(), *pathlib.Path.cwd().parents]:
        if _looks_like_workspace(candidate):
            return candidate / "subjects"
    here = pathlib.Path(__file__).resolve().parents[4]
    if _looks_like_workspace(here):
        return here / "subjects"
    raise ManifestError(
        "Cannot locate subjects/: run from inside the workspace or set TTM_SUBJECTS_DIR."
    )


def list_subjects(root: pathlib.Path | None = None) -> list[Subject]:
    root = root or subjects_dir()
    return [Subject.from_toml(p) for p in sorted(root.glob(f"*/{SUBJECT_FILENAME}"))]


def find_subject(subject_id: str, root: pathlib.Path | None = None) -> Subject:
    subjects = list_subjects(root)
    matches = [s for s in subjects if s.id.lower() == subject_id.lower()
               or s.dir.name.lower() == subject_id.lower()]
    if len(matches) != 1:
        raise ManifestError(
            f"Subject {subject_id!r}: {len(matches)} matches among "
            f"{[s.id for s in subjects]}"
        )
    return matches[0]


def resolve_experiment(
    *,
    manifest: str | pathlib.Path | None = None,
    subject: str | None = None,
    experiment: str | None = None,
) -> Experiment:
    """Pick the experiment a tool should operate on. See the module docstring."""
    if manifest:
        return Experiment.from_toml(manifest)
    subject = subject or os.environ.get("TTM_SUBJECT")
    experiment = experiment or os.environ.get("TTM_EXPERIMENT")
    if not (subject and experiment):
        raise ManifestError(
            "No experiment selected: pass --manifest PATH, or --subject ID "
            "--experiment NAME, or set TTM_SUBJECT and TTM_EXPERIMENT."
        )
    return find_subject(subject).experiment(experiment)
