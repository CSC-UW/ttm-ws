"""Per-pulse raster galleries for one sorted TTM-003 session.

One panel per selected pulse (-50 ms to pulse end + 100 ms); by default
the first, middle and last pulse of every train (``--pulses all`` for every
pulse). Three galleries, as multi-page PDFs under
``<experiment>/figures/galleries/``:

* ``channels``: one row per channel in tetrode order; a unit's spikes are
  drawn on the row of its peak channel, colour-coded per unit within the
  tetrode; unassigned threshold crossings on that channel in grey.
* ``tetrodes``: one row per tetrode, same colouring, unassigned crossings
  of the tetrode's four wires in grey.
* ``units``: one row per sorted cluster (every cluster, whatever its
  quality), grouped by tetrode, plus one grey row per tetrode with its
  unassigned crossings.
* ``tt_traces``: the 300-6000 Hz traces of each tetrode in ``--trace-tetrode``
  (lab names, default TT6, or ``all``) around each pulse with the sorted
  spikes marked, one PDF per tetrode.

"Unassigned" means a negative crossing of ``--thresh-mad`` (default 4.5)
times the channel's MAD on the preprocessed (bandpass + CMR) recording with
no sorted spike of any cluster on the same tetrode within 0.5 ms; these
are the spikes the sorter dropped or never detected. Channels listed as
bad in the manifest contribute no crossings and are labelled "(bad)". Crossings are cached
in ``<session>/sortings_opto/<variant>/crossings_<thresh>mad.npz``.

Usage: python subjects/ttm-003/scripts/opto_pulse_gallery.py --subject TTM-003 --experiment TTM-003-OPTO1 [--pulses first,middle,last|all]
"""

from __future__ import annotations

import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import spikeinterface as si

from ttm.cli._common import add_experiment_args, experiment_from_args

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000
DEAD = 15          # samples (0.5 ms) between crossings on one channel
MATCH = 15         # samples (0.5 ms) tolerance to call a crossing "assigned"
PRE_MS, POST_MS = 50, 100
UNASSIGNED_COLOR = "0.65"


def detect_crossings(rec, thresh_mad: float, chunk_s: float = 60.0) -> list[np.ndarray]:
    """Negative threshold crossings per channel (frame indices), 0.5 ms dead time."""
    n = rec.get_num_frames(); nch = rec.get_num_channels()
    ref = rec.get_traces(start_frame=0, end_frame=min(n, FS * 120))
    mad = np.median(np.abs(ref), axis=0) / 0.6745
    thr = -thresh_mad * mad
    out = [[] for _ in range(nch)]
    step = int(chunk_s * FS)
    prev_last = np.zeros(nch, bool)
    for s in range(0, n, step):
        x = rec.get_traces(start_frame=s, end_frame=min(n, s + step))
        below = x < thr
        start = below & ~np.vstack([prev_last[None, :], below[:-1]])
        prev_last = below[-1]
        for c in range(nch):
            idx = np.flatnonzero(start[:, c]) + s
            out[c].append(idx)
    res = []
    for c in range(nch):
        idx = np.concatenate(out[c]) if out[c] else np.array([], int)
        if idx.size:
            idx = idx[np.r_[True, np.diff(idx) > DEAD]]
        res.append(idx)
    return res


def unassigned_mask(cross: np.ndarray, assigned: np.ndarray) -> np.ndarray:
    if assigned.size == 0:
        return np.ones(len(cross), bool)
    j = np.searchsorted(assigned, cross)
    lo = assigned[np.clip(j - 1, 0, len(assigned) - 1)]
    hi = assigned[np.clip(j, 0, len(assigned) - 1)]
    return np.minimum(np.abs(cross - lo), np.abs(cross - hi)) > MATCH


def window(frames: np.ndarray, t: np.ndarray, on: float, pre: float, post: float) -> np.ndarray:
    """Spike frames -> ms relative to onset, restricted to the panel window."""
    i0 = np.searchsorted(t, on - pre); i1 = np.searchsorted(t, on + post)
    a = np.searchsorted(frames, i0); b = np.searchsorted(frames, i1)
    return (t[frames[a:b]] - on) * 1000


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--variant", default="asis")
    p.add_argument("--thresh-mad", type=float, default=4.5)
    p.add_argument("--pulses", default="first,middle,last",
                   help="which pulses of each train: comma list of first/middle/last, or 'all'")
    p.add_argument("--trace-tetrode", nargs="+", default=["TT6"], help="lab tetrode names, or 'all'")
    p.add_argument("--galleries", nargs="+", default=["channels", "tetrodes", "units", "tt_traces"])
    a = p.parse_args()
    exp = experiment_from_args(a)
    out = exp.sorting_dir(a.variant, subdir="sortings_opto")
    gal_dir = exp.figures_dir / "galleries"; gal_dir.mkdir(parents=True, exist_ok=True)

    rec = si.load(str(out / "preprocessed"))
    sorting = si.load(str(out / "sorting"))
    analyzer = si.load_sorting_analyzer(str(out / "analyzer"))
    t = rec.get_times(); t = t - t[0]
    ids = np.asarray(rec.channel_ids); groups = np.asarray(rec.get_property("group")).astype(int)
    tt_name = {g: f"TT{g + 1}" for g in np.unique(groups)}
    ugroup = pd.Series(sorting.get_property("group"), index=sorting.unit_ids).astype(int)
    peak_ch = si.get_template_extremum_channel(analyzer, outputs="index")
    units = list(sorting.unit_ids)
    spikes = {u: sorting.get_unit_spike_train(u) for u in units}

    cache = out / f"crossings_{a.thresh_mad:g}mad.npz"
    if cache.exists():
        z = np.load(cache); cross = [z[f"ch{c}"] for c in range(len(ids))]; unas = [z[f"un{c}"] for c in range(len(ids))]
    else:
        cross = detect_crossings(rec, a.thresh_mad)
        unas = []
        for c in range(len(ids)):
            assigned = np.sort(np.concatenate([spikes[u] for u in units if ugroup[u] == groups[c]] or [np.array([], int)]))
            unas.append(cross[c][unassigned_mask(cross[c], assigned)])
        np.savez(cache, **{f"ch{c}": cross[c] for c in range(len(ids))}, **{f"un{c}": unas[c] for c in range(len(ids))})
    bad = set(exp.badChannels)
    for c in range(len(ids)):
        if ids[c] in bad:
            unas[c] = unas[c][:0]
    tet_unas = {}
    for g in np.unique(groups):
        idx = np.sort(np.concatenate([unas[c] for c in np.flatnonzero(groups == g)]))
        tet_unas[g] = idx[np.r_[True, np.diff(idx) > DEAD]] if idx.size else idx
    print("unassigned crossings per tetrode (/s):", {tt_name[g]: round(len(v) / t[-1], 1) for g, v in tet_unas.items()})

    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv"); pulses = pulses[~pulses.glitch]
    trains = pd.read_csv(RESULTS / "stimulation_trains.csv"); trains = trains[trains.session == exp.session].set_index("ttl_train")
    if a.pulses != "all":
        picks = []
        for _, g in pulses.groupby("train"):
            n = len(g)
            want = {"first": 0, "middle": n // 2, "last": n - 1}
            picks.append(g.iloc[sorted({want[k] for k in a.pulses.split(",")})])
        pulses = pd.concat(picks)
    pre, post = PRE_MS / 1000, POST_MS / 1000

    # colours: per tetrode, tab10 cycling over its units
    cmap = plt.get_cmap("tab10")
    ucolor = {}
    for g in np.unique(ugroup.values):
        for k, u in enumerate([u for u in units if ugroup[u] == g]):
            ucolor[u] = cmap(k % 10)

    def title(row):
        tr = trains.loc[row.train]
        return f"train {int(tr.n)} ({tr.ttl_pattern}, {tr.power_mW:g} mW, {tr.state}) pulse {int(row.pulse_in_train)}/{int(tr.n_pulses)}  t={row.onset_s:.1f}s"

    def panel_common(ax, dur_ms, ttl, nrows):
        ax.axvspan(0, dur_ms, color="c", alpha=0.15, lw=0)
        ax.set_xlim(-PRE_MS, dur_ms + POST_MS); ax.set_ylim(-0.6, nrows - 0.4)
        ax.invert_yaxis(); ax.set_title(ttl, fontsize=7); ax.tick_params(labelsize=6)

    def gallery(name, rows_fn, labels, per_page, figsize):
        pdf_path = gal_dir / f"{name}_{a.variant}.pdf"
        nrow = len(labels)
        with PdfPages(pdf_path) as pdf:
            rows = list(pulses.itertuples())
            for s in range(0, len(rows), per_page):
                page = rows[s:s + per_page]
                ncol = 2 if per_page > 1 else 1
                fig, axes = plt.subplots(int(np.ceil(per_page / ncol)), ncol, figsize=figsize, squeeze=False)
                axes = axes.ravel()
                for ax in axes[len(page):]:
                    ax.axis("off")
                for ax, row in zip(axes, page):
                    dur_ms = (row.offset_s - row.onset_s) * 1000
                    on = row.onset_s
                    events, colors, lineoffsets = rows_fn(on, dur_ms)
                    if events:
                        ax.eventplot(events, colors=colors, lineoffsets=lineoffsets, linelengths=0.8, linewidths=0.7)
                    panel_common(ax, dur_ms, title(row), nrow)
                    ax.set_yticks(range(nrow)); ax.set_yticklabels(labels, fontsize=4.5)
                    ax.set_rasterized(True)
                axes[-1].set_xlabel("ms from laser onset", fontsize=7)
                fig.suptitle(f"{exp.experiment} {name}: rows = {name}; grey = unassigned crossings (>{a.thresh_mad:g} MAD, no sorted spike within 0.5 ms)", fontsize=8)
                fig.tight_layout(rect=(0, 0, 1, 0.97)); pdf.savefig(fig, dpi=130); plt.close(fig)
        print("wrote", pdf_path)

    # ---- channels gallery -------------------------------------------------
    order = np.argsort(groups, kind="stable")
    ch_labels = [f"{ids[c]} {tt_name[groups[c]]}" + (" (bad)" if ids[c] in bad else "") for c in order]
    row_of_ch = {c: r for r, c in enumerate(order)}
    units_on_ch = {c: [u for u in units if peak_ch[u] == c] for c in range(len(ids))}

    def rows_channels(on, dur_ms):
        ev, col, off = [], [], []
        for c in order:
            w = window(unas[c], t, on, pre, dur_ms / 1000 + post)
            if w.size: ev.append(w); col.append(UNASSIGNED_COLOR); off.append(row_of_ch[c])
            for u in units_on_ch[c]:
                w = window(spikes[u], t, on, pre, dur_ms / 1000 + post)
                if w.size: ev.append(w); col.append(ucolor[u]); off.append(row_of_ch[c])
        return ev, col, off

    # ---- tetrodes gallery -------------------------------------------------
    tets = sorted(np.unique(groups))
    tet_labels = [tt_name[g] + (" (bad)" if tt_name[g] in exp.badTetrodes else "") for g in tets]

    def rows_tetrodes(on, dur_ms):
        ev, col, off = [], [], []
        for r, g in enumerate(tets):
            w = window(tet_unas[g], t, on, pre, dur_ms / 1000 + post)
            if w.size: ev.append(w); col.append(UNASSIGNED_COLOR); off.append(r)
            for u in [u for u in units if ugroup[u] == g]:
                w = window(spikes[u], t, on, pre, dur_ms / 1000 + post)
                if w.size: ev.append(w); col.append(ucolor[u]); off.append(r)
        return ev, col, off

    # ---- units gallery ----------------------------------------------------
    unit_rows = []   # (kind, key, label)
    for g in tets:
        for u in [u for u in units if ugroup[u] == g]:
            unit_rows.append(("unit", u, f"u{u} {tt_name[g]} {ids[peak_ch[u]]}"))
        unit_rows.append(("unas", g, f"{tt_name[g]} unassigned"))
    unit_labels = [r[2] for r in unit_rows]

    def rows_units(on, dur_ms):
        ev, col, off = [], [], []
        for r, (kind, key, _) in enumerate(unit_rows):
            src = spikes[key] if kind == "unit" else tet_unas[key]
            w = window(src, t, on, pre, dur_ms / 1000 + post)
            if w.size: ev.append(w); col.append(ucolor[key] if kind == "unit" else UNASSIGNED_COLOR); off.append(r)
        return ev, col, off

    if "channels" in a.galleries:
        gallery("channels", rows_channels, ch_labels, per_page=4, figsize=(16, 20))
    if "tetrodes" in a.galleries:
        gallery("tetrodes", rows_tetrodes, tet_labels, per_page=8, figsize=(16, 18))
    if "units" in a.galleries:
        gallery("units", rows_units, unit_labels, per_page=4, figsize=(16, 20))

    # ---- traces of one tetrode ---------------------------------------------
    trace_tts = [tt_name[g] for g in tets] if a.trace_tetrode == ["all"] else a.trace_tetrode
    for trace_tt in (trace_tts if "tt_traces" in a.galleries else []):
        g = int(trace_tt.removeprefix("TT")) - 1
        chans = np.flatnonzero(groups == g)
        tunits = [u for u in units if ugroup[u] == g]
        gain = float(rec.get_property("gain_to_uV")[chans[0]]) if "gain_to_uV" in rec.get_property_keys() else 1.0
        pdf_path = gal_dir / f"{trace_tt}_traces_{a.variant}.pdf"
        with PdfPages(pdf_path) as pdf:
            rows = list(pulses.itertuples())
            for s in range(0, len(rows), 4):
                fig, axes = plt.subplots(4, 1, figsize=(16, 20))
                for ax in axes[len(rows[s:s + 4]):]:
                    ax.axis("off")
                for ax, row in zip(axes, rows[s:s + 4]):
                    dur_ms = (row.offset_s - row.onset_s) * 1000
                    i0 = np.searchsorted(t, row.onset_s - pre); i1 = np.searchsorted(t, row.onset_s + dur_ms / 1000 + post)
                    x = rec.get_traces(start_frame=i0, end_frame=i1, channel_ids=ids[chans]) * gain
                    ts = (t[i0:i1] - row.onset_s) * 1000
                    step = max(np.percentile(np.abs(x), 99.9) * 2.2, 200)
                    for k in range(len(chans)):
                        ax.plot(ts, x[:, k] - k * step, lw=0.4, color="k")
                        ax.text(ts[0], -k * step, str(ids[chans[k]]), fontsize=6, ha="right", va="center")
                    ymin = -(len(chans) - 0.3) * step; ymax = step * 0.7
                    for j, u in enumerate(tunits):
                        w = window(spikes[u], t, row.onset_s, pre, dur_ms / 1000 + post)
                        ax.plot(w, np.full(len(w), ymax - (j + 1) * step * 0.12), "|", color=ucolor[u], ms=7, mew=1.2, label=f"u{u}")
                    w = window(tet_unas[g], t, row.onset_s, pre, dur_ms / 1000 + post)
                    ax.plot(w, np.full(len(w), ymax - (len(tunits) + 1) * step * 0.12), "|", color=UNASSIGNED_COLOR, ms=7, mew=1.2, label="unassigned")
                    ax.axvspan(0, dur_ms, color="c", alpha=0.15, lw=0); ax.set_xlim(-PRE_MS, dur_ms + POST_MS); ax.set_ylim(ymin, ymax + step * 0.2)
                    ax.set_yticks([]); ax.set_title(title(row), fontsize=7); ax.tick_params(labelsize=6)
                    ax.set_rasterized(True)
                axes[0].legend(fontsize=5, ncol=8, loc="upper right")
                axes[-1].set_xlabel("ms from laser onset", fontsize=7)
                fig.suptitle(f"{exp.experiment} {trace_tt}: 300-6000 Hz + CMR traces (uV) with sorted spikes (colour) and unassigned crossings (grey)", fontsize=8)
                fig.tight_layout(rect=(0, 0, 1, 0.97)); pdf.savefig(fig, dpi=130); plt.close(fig)
        print("wrote", pdf_path)


if __name__ == "__main__":
    main()
