---
name: setup-workspace
description: Set up this workspace on a new machine - check uv/git/Python, sync the environment, run tests, report GPU status, point at the data. Use when asked to set up, install, onboard, get running, or when a command fails because nothing is configured yet.
---

# Set up the ttm-ws workspace

Work through the steps in order. Report each result in one line before
moving on. Ask the user only where a step says to. Do not clone other
repositories or edit tracked files; the only file written is the per-machine
config in step 5. The external
packages are pinned in `pyproject.toml` and resolve on their own.

## 1. Prerequisites

```bash
git --version
uv --version          # need 0.9 or newer
uv python find 3.14   # any 3.14+ interpreter uv can see; uv can download one
```

- No `uv`: offer the official installer
  (`curl -LsSf https://astral.sh/uv/install.sh | sh`), then re-check.
- No Python 3.14+: `uv python install 3.14`.

## 2. Environment

```bash
uv sync --all-packages --group dev
```

Success: exits 0. The first run downloads torch (about 1 GB). If it fails
fetching a git source, show the failing URL verbatim and say it is a network
or GitHub problem; do not edit `pyproject.toml` sources to work around it.

## 3. Tests

```bash
uv run pytest -q
```

Report passed/failed/skipped counts. Two tests skip when production data
are not mounted; that is expected off the server.

## 4. GPU

```bash
uv run python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())"
```

Report the line. On `tononi-2` it must show a CUDA 12 build and `True`.
Elsewhere a CPU build and `False` is normal; say that Kilosort4 needs the
GPU and everything else does not.

## 5. Per-machine config

```bash
uv run python -c "from ttm import config; print(config.get())"
```

If this raises `ConfigError`, ask the user where the data are on this
machine and write `~/.config/ttm/config.toml` as `docs/SETUP.md` section 5
shows (server: `data_root` plus `remotes.archive`; laptop: `data_root` on the
external drive plus `remotes.server` as an rclone `remote:path`). Re-run the
check. On `tononi-2`, also report `df -h /nvme`.

## 6. Orientation

Point the user at, in this order: `README.md`, `CLAUDE.md`,
`docs/pipeline.md`, `packages/ttm/README.md`,
`investigations/README.md`. State the three
rules they will hit first: always `uv sync --all-packages`, every tool needs
an experiment (`--subject ID --experiment NAME`), and commit small changes
directly to `main`.
