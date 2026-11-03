---
title: Data products of a 48 h tetrode recording
updated: 2026-11-03
---

# Data products

Every product of one 48 h recording, what makes it, what it is made
from, and its size. Sizes are TTM-003-48H's (64 channels, 30 kHz, 48 h).
`<session>` is `<data_root>/<session>/`, `<experiment>` is
`<data_root>/<subject>/<experiment>/`, `<run>` is
`<session>/sortings_carry/<run>/`. Steps and parameters of the sort are
in [`sorting_pipeline.md`](sorting_pipeline.md); the commands in
[`pipeline.md`](pipeline.md).

Every box has four lines, always in this order:

1. what the product is;
2. its format and sample type;
3. its size for TTM-003-48H;
4. what makes it.

Solid boxes are kept. Dashed grey boxes can be regenerated from their
inputs and may be deleted once nothing in use reads them. Dotted boxes
do not exist yet. An arrow points from an input to what is made from it.

```mermaid
flowchart LR
    classDef src fill:#ecfdf3,stroke:#047857,color:#111
    classDef keep fill:#e8f0fe,stroke:#1a56db,color:#111
    classDef temp fill:#f4f4f4,stroke:#9a9a9a,stroke-dasharray:5 4,color:#666
    classDef plan fill:#fff,stroke:#b45309,stroke-dasharray:2 3,color:#7c4a03

    subgraph IN["inputs"]
        RAW["<b>Open Ephys session</b><br/>binary, int16<br/>697 GB<br/>acquisition"]:::src
        MAN["<b>manifest</b><br/>TOML<br/>a few KB<br/>by hand"]:::src
        PUL["<b>laser pulse table</b><br/>CSV<br/>under 1 MB<br/>subject script"]:::src
    end

    subgraph SES["session directory"]
        STORE["<b>Zarr store</b><br/>Zarr, int16, blosc-zstd<br/>355 GB<br/>ttm-convert"]:::keep
        ST["<b>slice table</b><br/>CSV<br/>under 1 KB<br/>ttm-convert"]:::keep
        LFP["<b>LFP, 625 Hz</b><br/>Zarr, float32, blosc-zstd<br/>21 GB<br/>ttm-make-lfp"]:::keep
        LFP125["<b>LFP, 125 Hz, 16 leads</b><br/>Zarr, float32, blosc-zstd<br/>1.1 GB<br/>ttm-make-lfp"]:::keep
        BP["<b>bad periods</b><br/>CSV, JSON, npz cache<br/>0.2 MB, 2 GB<br/>ttm-bad-periods"]:::keep
        GBIN["<b>global-median binary</b><br/>binary, float32<br/>1.3 TB<br/>mp_baseline.py"]:::temp
    end

    subgraph EXP["experiment directory"]
        PAR["<b>experiment params</b><br/>JSON<br/>2 KB<br/>ttm-experiment-params"]:::keep
        SWA["<b>delta, band powers</b><br/>Zarr, float32, blosc-lz4<br/>1.1 GB<br/>ttm-swa"]:::keep
        EMG["<b>synthetic EMG</b><br/>Zarr, float64, blosc-lz4<br/>42 MB<br/>ttm-emg"]:::keep
        FIG["<b>figures</b><br/>PNG<br/>160 MB<br/>ttm-swa, ttm-emg"]:::temp
    end

    subgraph RUN["sorting run directory"]
        BIN["<b>aligned-median binary</b><br/>binary, int16<br/>about 660 GB<br/>ttm-sort binary"]:::temp
        WT["<b>wire timing</b><br/>npz<br/>0.1 MB<br/>ttm-sort binary"]:::keep
        SEED["<b>seed sort, bank, whitening</b><br/>sorter folder, npz<br/>20 MB<br/>ttm-sort seed"]:::keep
        CW["<b>per-window spikes, banks</b><br/>npz<br/>51 GB<br/>ttm-sort carry"]:::temp
        RS["<b>re-seed sorts</b><br/>sorter folders<br/>130 MB<br/>ttm-sort carry"]:::keep
        SORT["<b>sorting</b><br/>SpikeInterface folder, npy<br/>3.3 GB<br/>ttm-sort export"]:::keep
    end

    subgraph NEXT["after sorting"]
        AN["<b>SortingAnalyzer</b><br/>Zarr<br/>about 0.5 GB<br/>not yet made"]:::plan
        HYP["<b>hypnogram</b><br/>table<br/>under 1 MB<br/>not yet made"]:::plan
        MET["<b>metrics per state</b><br/>table<br/>under 1 MB<br/>not yet made"]:::plan
        CUR["<b>curated sorting</b><br/>SpikeInterface folder<br/>about 3 GB<br/>not yet made"]:::plan
    end

    subgraph LAP["copies for a laptop"]
        VIZ["<b>viz bundle</b><br/>Zarr, Parquet, JSON<br/>2.0 GB<br/>ttm-viz-bundle"]:::temp
        CB["<b>curation bundle</b><br/>JSON list of files<br/>analyzer and store<br/>ttm-curation-manifest"]:::temp
    end

    RAW --> STORE
    RAW --> ST
    MAN --> PAR
    ST --> PAR
    STORE --> LFP --> LFP125
    LFP --> SWA
    STORE --> EMG
    SWA --> EMG
    SWA --> FIG
    EMG --> FIG
    STORE --> BP
    MAN --> BP
    STORE --> GBIN
    MAN --> BIN
    STORE --> BIN
    BIN --> WT
    BIN --> SEED
    BIN --> CW
    SEED --> CW
    BP --> CW
    PUL --> CW
    CW --> RS
    RS --> CW
    CW --> SORT
    SORT --> AN
    BIN --> AN
    AN --> MET
    HYP --> MET
    BP --> MET
    MET --> CUR
    AN --> CB
    CB --> CUR
    EMG --> VIZ
    LFP125 --> VIZ
    SORT --> VIZ
    VIZ -- "scored in ttm-loupe" --> HYP
```

## Products

| Product | Path | Format | Made by | From | Size | Kept |
|---|---|---|---|---|---|---|
| Open Ephys session | archive `SUBJECT/<session>/` | binary, int16 (0.195 uV per bit) | acquisition |  | 697 GB | yes, the original |
| Manifest | `subjects/<id>/experiments/<NAME>.toml` | TOML | by hand |  | a few KB | yes, in git |
| Laser pulse table | `subjects/<id>/results/` | CSV | subject script (TTM-003: `laser_pulses.py`) | session TTL events | under 1 MB | yes, in git |
| Zarr store | `<session>/<session>.blosc-zstd.zarr` | Zarr, int16, blosc-zstd (lossless) | `ttm-convert` | session | 355 GB | yes; every product below starts here |
| Slice table | `<session>/slice_table.csv` | CSV | `ttm-convert` | session | under 1 KB | yes |
| Experiment params | `<experiment>/experiment_params.json` | JSON | `ttm-experiment-params` | manifest, slice table | 2 KB | yes; rewritten when the manifest changes |
| LFP, 625 Hz | `<session>/<session>.lfp.zarr` | Zarr, float32 (uV), blosc-zstd | `ttm-make-lfp` | store | 21 GB | yes |
| LFP, 125 Hz, 16 leads | `<session>/<session>.lfp.125hz.zarr` | Zarr, float32 (uV), blosc-zstd | `ttm-make-lfp --subsample 125` | 625 Hz LFP | 1.1 GB | yes |
| Delta, band powers | `<experiment>/delta.idelta.zarr`, `stft_bandpowers.zarr` | Zarr, float32, blosc-lz4 | `ttm-swa all` | 625 Hz LFP | 1.1 GB, 29 MB | yes |
| Synthetic EMG | `<experiment>/synthetic_emg_methods.zarr` | Zarr, float64, blosc-lz4 | `ttm-emg --full` | store, slice table, delta | 42 MB | yes |
| Figures | `<experiment>/figures/` | PNG | `ttm-swa plot`, `ttm-emg` | the above | 160 MB | regenerable |
| Bad periods | `<session>/bad_periods/` | CSV table; JSON of its input and thresholds; npz cache of the scan | `ttm-bad-periods` | store (aligned median, lazily), manifest; or a binary | 0.2 MB table, 2 GB cache | the table yes; the cache only to rerun thresholds |
| Global-median binary | `<session>/sortings_seed42_pcafix/mp_baseline/binary/` | binary, float32 (ADC counts) | `mp_baseline.py` | store | 1.3 TB | only while sortings made on it are in use; the pipeline no longer reads it |
| Aligned-median binary | `<run>/binary/` | binary, int16 (ADC counts, rounded) | `ttm-sort binary` | store, manifest | about 660 GB | while its sorting is analysed or curated; 40 min to rewrite |
| Wire timing | `<run>/wire_timing.npz` | npz, float64 | `ttm-sort binary` | store | 0.1 MB | yes; the binary is reproducible from it |
| Seed sort, bank, whitening | `<run>/seed_sort/`, `bank0.npz`, `whitening.npz` | MountainSort5 folder; npz, float32 | `ttm-sort seed` | binary | 20 MB | yes |
| Per-window spikes and banks | `<run>/carry/` | npz | `ttm-sort carry` | binary, bank, bad periods, pulses | 51 GB | until the sorting is final; resume and re-export read it |
| Re-seed sorts | `<run>/reseed_sorts/` | MountainSort5 folders | `ttm-sort carry` | binary | 130 MB | yes |
| Sorting | `<run>/sorting/` | SpikeInterface folder (npy) | `ttm-sort export` | per-window spikes | 3.3 GB | yes |
| Viz bundle | `<session>/viz_bundle/` | Zarr, Parquet, JSON | `ttm-viz-bundle` | EMG, 125 Hz LFP, a sorting | 2.0 GB | regenerable |
| Investigation products | `<session>/sortings_seed42_pcafix/` (hybrids, arms, the first 6 h binaries, earlier sorts) | various | `investigations/` |  | 1.1 TB besides the two 48 h float32 binaries (global and aligned median, `mp_baseline/` and `mp_aligned/`, 1.3 TB each) | until the user decides |

Not yet produced: the SortingAnalyzer of the new sorting (built on the
binary it was sorted from; June production's was 0.5 GB), the hypnogram
(scored in `ttm-loupe`), quality metrics per vigilance state, and the
curated sorting.
