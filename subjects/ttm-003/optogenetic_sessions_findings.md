---
title: TTM-003 optogenetic sessions of 2026-10-16, first analysis
updated: 2026-10-22
status: needs_verification
---

# TTM-003 optogenetic sessions: artifacts, sorting, tagging, population effects

Analysis of the two sessions described in `optogenetic_stimulation.md`,
run on 2026-10-21 with the scripts in `scripts/` (converted stores and
figures on the data volume; small tables in `results/`). Measurements are
stated as such; sentences marked *conjecture* are the analyst's reading
and have not been confirmed by the experimenter. Tetrodes are named as
the lab names them (`TT1`..`TT16`, the recording's `tetrode` property);
the `group` column in the CSVs is the 0-based index.

## 1. Light artifact

Method (`opto_onset_zoom.py`, `opto_explore.py`): 300 to 6000 Hz
pulse-triggered averages at light onset and offset for the lowest- and
highest-power trains; the raw (unfiltered) averages for every train. A
photoelectric transient is locked to the light edge (0 ms, again at
offset) and scales with power; evoked spikes start after about 1 ms and
jitter.

| Session, train | noise floor (uV) | at onset edge, median / max (uV) | at offset edge, median / max | evoked 1 to 6 ms, max |
|---|---|---|---|---|
| OPTO1 train 1, 0.5 mW | 2.1 | 4.5 / 7.9 | 5.6 / 9.7 | 28.8 |
| OPTO1 train 6, 9.8 mW | 2.1 | 3.1 / 8.6 | 4.3 / 8.0 | 42.5 |
| OPTO2 train 1, 9.8 mW | 2.0 | 2.5 / 7.2 | 2.7 / 4.8 | 60.4 |
| OPTO2 train 3, 4.2 mW | 2.0 | 2.0 / 6.9 | 1.5 / 4.4 | 47.0 |

Nothing stands above the noise floor at either light edge at any power,
on any channel, so no spike-band artifact was detected and no blanking was
applied before sorting. The only edge-locked spike-band signal is a
synchronous volley at 1 to 4 ms after onset: on TT4 it peaks at 2.7 ms
(about -40 uV in the 200-pulse average at 9.8 mW, -28 uV at 0.5 mW); on
TT6 (4 ms) and TT10 (4.5 ms) it appears only at the high powers. Single
trials on TT4 show individual spikes at 2 to 4 ms.

The raw averages carry a slow potential of several hundred microvolts
that lasts the whole pulse and reverses after offset: median sustained
shift across channels 109 uV (max 272) at 0.5 mW and 193 uV (max 374) at
9.8 mW, negative on TT1 to TT3 and TT13 to TT16, positive on TT4 to
TT12. *Conjecture:* its weak dependence on a twenty-fold power change,
its sign reversal across tetrodes and its coincidence with the firing
suppression below point to a physiological evoked potential rather than
a photovoltaic artifact, which would scale with intensity; an
opsin-negative or fibre-disconnected control recording would settle it.

## 2. Bad channels

Raw per-channel standard deviation over 60 s (`opto_explore` run):
TT14 (CH10, CH12, CH14, CH16) about 1.4 mV on all four wires with
large slow waves and no spikes; CH2 on TT13 1.3 mV with slow drift.
Every other channel is 165 to 470 uV. Recorded in both manifests as
`badChannels` and `badTetrodes = ["TT14"]`.

## 3. Sorting

`opto_sort.py`: bandpass 300 to 6000 Hz, global common median reference,
MountainSort5 scheme 2 per tetrode, detect threshold 5.5, whitening seed
42; analyzer with waveforms, templates, amplitudes, correlograms and
quality metrics under `<session>/sortings_opto/asis/`. 56 units in each
session (per-tetrode counts in `results/units_*.csv`; metrics in
`results/unit_metrics_*.csv`). The global median changes by less than
15 % in standard deviation during pulses, so the reference does not
inject pulse-locked activity.

## 4. Population effect of the pulses

`opto_tag.py` classes each unit per pulse pattern (set `B`: 64 ms pulses
at 9.8 mW; `A`: 207 ms pulses at 4.2 mW; `A_low`: 207 ms at 0.5 to
0.55 mW) and `opto_population.py` pools the suppressed units.

| | OPTO1 (56 units) | OPTO2 (56 units) |
|---|---|---|
| suppressed during pattern B | 40 | 44 |
| activated or weakly activated | 2 | 1 |
| unaffected | 14 | 11 |
| median in-pulse rate of suppressed units, fraction of baseline | 0.027 | 0.029 |
| median rate 0 to 50 ms after offset | 0.031 | 0.037 |

Pooled over the suppressed units, firing falls below half its baseline
in the first 5 ms bin after onset in every train of both sessions, stays
at 1 to 5 % of baseline for the rest of the pulse, and recovers to half
baseline 45 to 105 ms after light offset (pattern A) or 70 to 105 ms
(pattern B); a rebound of 1.2 to 1.9 times baseline follows. The two
0.5 mW trains produce the same picture (in-pulse 3 to 5 %, recovery 50
and 100 ms after offset) as the 4.2 and 9.8 mW trains. Awake and drowsy
trains do not differ in these measures.

The suppression is spatial: every unit on TT5 to TT8 and TT10 to TT12 is
suppressed, two of three on TT4, none on TT1 to TT3, TT9, TT13, TT15 or
TT16 (those tetrodes hold few, low-rate units). Multi-unit threshold
crossings per tetrode (`opto_explore`, 5 MAD on each wire) show the same
split, with the pulse-period rate of TT5 to TT12 at 0.5 to 0.8 of
baseline and TT1 to TT3, TT13, TT15, TT16 at about 1.0; the crossing
counts sit on a floor of noise events (TT14, which has no spikes, still
gives 15 crossings/s), which is why they understate the suppression the
sorted units show. *Conjecture:* the suppressed
tetrodes are the ones within reach of the fibre, and the complete,
saturating silence even at 0.5 mW means the lowest power used was
already far above threshold.

The hand log says the 0.5 and 0.55 mW trains gave no or doubtful neural
response. By these measures they gave the full response.

## 5. Tagging candidates

Criteria (`opto_tag.py`, 64 ms pulses primarily): SALT first-spike
latency test with twenty 10 ms baseline windows (p floor 0.005), early
(1 to 10 ms) rate over 3 times baseline, reliability over 0.3, in-pulse
vs out-of-train waveform correlation over 0.85. Two cells pass in each
session, the same two by tetrode, waveform and rate:

| | TT6 cell (OPTO1 unit 21, OPTO2 unit 22) | TT10 cluster (OPTO1 unit 34, OPTO2 unit 35) |
|---|---|---|
| rate, SNR, trough amplitude | 13 Hz, 9.5, -110 uV on CH49 | 8 to 10 Hz, 7.2, -80 uV on CH19 |
| peak-to-trough, trough half-width | 0.63 ms, 0.23 ms (population median 0.64, 0.19) | 0.55 to 0.64 ms, 0.15 ms |
| ISI violation ratio, rp contamination | 0.11 to 0.17, 0.05 to 0.11 | 0.60 to 0.65, 0.46 to 0.49 |
| first-spike latency (median), jitter | 4.4 to 4.8 ms (modes at 3.5 and 8 ms), 3 to 4 ms | 6.1 to 7.4 ms, 4.5 ms |
| reliability (a spike in 1 to 10 ms) | 0.40 to 0.48 | 0.53 to 0.72 |
| during the rest of the pulse | silent (0.06 to 0.17 of baseline) | keeps firing, decaying over ~50 ms (1.2 to 1.5 of baseline) |
| 0 to 50 ms after offset | silent | 0.25 to 0.35 of baseline, then rebound |
| at 0.5 mW | same onset response, same silence | same |
| in-pulse vs outside waveform r | 0.96 to 0.98 | 0.83 to 0.91 |

The TT6 cell fires one or two spikes 3 to 9 ms after every second pulse,
then nothing for the remainder of a 64 ms or a 207 ms pulse and for
about 100 ms after offset, like its neighbours. Its cross-correlograms
with the other 55 units outside the trains (`optotag_asis_unit21_ccgs.png`)
show no trough or peak at 1 to 5 ms except positive peaks with two other
TT6 clusters (units 14 and 20, z about 6.5), which may be split
fragments of the same cell. Spike-triggered LFP on its own tetrode shows
only the spike's own leak into the 1 kHz signal; phase locking outside
the trains is weak in every band (mean resultant length 0.03 to 0.09,
Rayleigh p small only because n is 18,000). Its autocorrelogram has a
wide (about 10 ms) relative refractory dip and a broad ISI distribution
with a mode near 10 to 20 ms. *Conjecture:* onset spike followed by
block is what Herman et al. 2014 describe for regular-spiking SST
interneurons under long pulses (see the literature survey), but a cell
that receives a brief disynaptic kick before the general suppression
would look the same; the correlograms do not show it inhibiting anyone,
so these sessions do not identify it.

The TT10 cluster is the only activity that follows the light for tens of
milliseconds while everything around it is silent, and it is depressed
after offset rather than rebounding, which is the direct-activation
signature in the survey. It is not a single unit: ISI violations and
contamination are high, and its correlograms with the other two TT10
clusters have broad central peaks. *Conjecture:* TT10 (and TT4, whose
2.7 ms volley the sorter assigned to no unit) carry the directly driven
cells, as more than one unit the sorter has merged; a sort with a lower
detection threshold, or Kilosort4, on TT4 and TT10 alone, followed by
manual curation, is the next step for isolating them.

### 5a. Does the TT6 cell keep firing through the pulse? (added later on 2026-10-21)

Graham Findlay suspected the TT6 cluster keeps only the first one or two
spikes of a train that continues through the pulse. `opto_tt6_inpulse.py`
detects every crossing above 3.5 MAD on the four TT6 wires around the
9.8 mW pulses, labels each as claimed by a cluster or not, and compares
its four-wire amplitude vector with the cell's template (cosine; a
fingerprint independent of overall amplitude). `results/tt6_event_rates_asis.csv`
and `results/tt_inpulse_<experiment>_TT6.json` hold the numbers.

| TT6, 9.8 mW pulses | baseline | onset 0 to 10 ms | in pulse (10 ms to 64 ms) | 0 to 100 ms after |
|---|---|---|---|---|
| all crossings (/s), OPTO1 / OPTO2 | 126 / 171 | 266 / 319 | 156 / 200 | 26 / 37 |
| claimed by any TT6 cluster (/s) | 38 / 48 | 54 / 61 | 2.5 / 2.0 | 9.5 / 13 |
| unclaimed, cosine to the cell's fingerprint > 0.95 (/s) | 20 / 26 | 126 / 139 | 133 / 160 | 3 / 4 |
| amplitude of those, 10th / 50th / 90th pct (uV) | -142 / -68 / -52 | -92 / -69 / -55 | -71 / -58 / -50 | -148 / -82 / -52 |
| the cell's own spikes (uV) | -156 / -111 / -80 | | | |

So the sorter is not dropping full-size spikes: during the pulse there
are only 2 to 2.5 claimed events per second on the whole tetrode and
nothing unclaimed above about -75 uV. What the sorter drops is a train of
small events, 130 to 160 per second, with the cell's wire dominance (CH49
largest), its spike width (trough-to-peak 0.60 vs 0.63 ms) and a dominant-
wire shape correlation of 0.99 with its template, at half its amplitude.
The train is regular (ISI median 4.5 ms, 5th percentile 3.2 ms, CV 0.6 to
0.7, so a refractory floor exists) and decays from about 210 per second in
the first 10 ms of the pulse to about 70 per second in its last 15 ms,
with the amplitude settling at -55 uV. The wire ratio is not identical:
in-pulse events are more CH49-dominant (CH49:CH53 of 3.3) than the
template (2.4) or than baseline look-alikes (2.3). *Conjecture:* the
simplest reading is that this is the same cell firing at 200 Hz with
spikes shrunk by sustained depolarisation, which changes the spatial
pattern as well as the amplitude; the alternative is a second, smaller
cell near CH49 that the light drives. A sort of TT6 with a lower threshold
that lets the small events form their own cluster, then a comparison of
that cluster's waveform against the cell outside the trains, is the test;
MountainSort5 at 5.5 MAD can never recover them.

The per-pulse galleries (first, middle and last pulse of each train) that
show this are `figures/galleries/{channels,tetrodes,units}_asis.pdf` and
`figures/galleries/TT6_traces_asis.pdf` (`opto_pulse_gallery.py`); grey
ticks are unclaimed crossings at 4.5 MAD.

### 5b. TT6 re-sorted at a lower threshold (2026-10-21, later)

`opto_tt_resort.py` sorts TT6 alone from the same preprocessed recording
with MountainSort5 scheme 2 at thresholds 4.0 and 3.5 (the default is
5.5; MountainSort5 applies it to the whitened signal). It yields 8 to 10
clusters per session. One cluster per sort reproduces 97 to 99 % of the
TT6 cell's spikes and has a fingerprint cosine of 0.998 to 0.999 to its
template; the sorter does not split the small in-pulse events off, it
folds them into that cluster, which grows from 13 Hz to 36 Hz (OPTO1) and
46 Hz (OPTO2) at threshold 3.5 and now fires 141 and 126 per second
through the 9.8 mW pulse (340 per second in the first bins, decaying to
about 70), then nothing for 100 ms after offset. The merged cluster's
trough amplitude outside the trains is a single continuous distribution
from -50 to -160 uV; inside the pulse it is -46 to -72 uV.
`results/resort_<experiment>_TT6.csv`, `figures/opto_TT6_resort_thr*.png`.

`opto_tt_amplitude_split.py` then asks whether the large (<= -85 uV) and
small (>= -70 uV) spikes of that cluster are one cell, using their
cross-correlogram outside the trains:

| | OPTO1 | OPTO2 |
|---|---|---|
| large / small spikes outside trains | 18,229 / 12,450 | 10,903 / 10,621 |
| small spikes within 0.75 ms of a large one, observed / expected if independent | 35 / 316 | 31 / 303 |
| small spikes 2 to 3 ms after a large one, ratio to independence | 1.6 | 1.5 |
| ISI floor of in-pulse small spikes (5th percentile) | 3.3 ms | 3.5 ms |

The two populations avoid each other at zero lag (one tenth of the
independent expectation) and the small spikes follow the large ones at
2 to 3 ms. That is one cell whose spike amplitude falls within a burst,
and the in-pulse train is the same cell firing at 100 to 300 per second
with every spike shrunk. The experimenter's suspicion holds: the main sort
kept the first one or two full-size spikes of each pulse and dropped the
rest. A second CH49-dominant cell (re-sort cluster 9, -180 uV, 3.4 Hz)
shares the wire but is silenced by the light and never activated.
`results/ampsplit_*.json`, `figures/opto_TT6_ampsplit_thr3.5_c*.png`.

For tagging on these tetrodes, the practical consequence is that the
sort must run at a lower threshold on the tetrodes of interest, or the
in-pulse spikes must be assigned afterwards by shape to the cell's
template with the amplitude left free. The TTM-001 production sort
(`investigations/spike-sorting/97_production_pipeline.py`) does not do
this: its matcher accepts a spike only at 0.8 of the template amplitude
or more, and its residual pass only looks at unclaimed events of 10 MAD
or more, so it would drop these half-amplitude spikes as well.

### 5c. TT10 re-sorted at a lower threshold (2026-10-22)

Graham Findlay saw the same problem in `galleries/TT10_traces_asis.pdf`:
small spikes during the pulse left grey (unassigned), and in other panels
in-pulse spikes with no tick at all. The second case is a display
threshold: the TT10 wires have a noise MAD of about 12 uV, the in-pulse
spikes are mostly 40 to 60 uV, and the galleries only ticked crossings
above 4.5 MAD (54 uV). `opto_tt6_inpulse.py --tetrode TT10` (crossings
above 3.5 MAD, fingerprint of unit 34/35) gives, per second, over the
9.8 mW pulses:

| | baseline | onset 0 to 10 ms | rest of pulse | 100 ms after |
|---|---|---|---|---|
| OPTO1: crossings / claimed by unit 34 / fingerprint look-alikes left unsorted | 189 / 7.7 / 33 | 373 / 51 / 119 | 207 / 13 / 115 | 130 / 5.6 / 23 |
| OPTO2: same for unit 35 | 213 / 8.4 / 41 | 390 / 45 / 114 | 250 / 12 / 139 | 111 / 4.0 / 22 |

So the light raises the rate of CH19-dominant events three to four times
and the 5.5 sort claims about a tenth of them; their median trough is
-52 and -49 uV (about 4.3 MAD), just under the sorter's threshold. In the
first pulse of OPTO2's first train, CH19 carries a train of spikes every
3.3 ms for the first 30 ms, each at -50 uV, of which the main sort keeps
three. Unlike TT6, the in-pulse events are not the reference unit's
waveform shrunk: their mean waveform is broader (trough to peak 0.83 vs
0.63 ms in OPTO1, 0.73 vs 0.60 in OPTO2; dominant-wire shape correlation
to the template 0.92 to 0.95, vs 0.99 for the baseline look-alikes).

`opto_tt_resort.py --tetrode TT10` at thresholds 4.0, 3.5 and 3.0 gives
three to five clusters per session, one per wire, every one with 0.2 to
0.5 refractory contamination: TT10 is a multi-unit tetrode and the sorter
is dividing its hash by dominant wire. The CH19 cluster absorbs the
in-pulse events (OPTO1 at 3.5: 43 per second baseline, 264 at onset, 178
through the pulse, 29 after; OPTO2: 44, 232, 170, 20), so the galleries
now show the pulse filled with its ticks (`TT10_traces_asis_resort3.5.pdf`,
and `TT10_traces_asis_resort3.5_3.5mad.pdf` with grey ticks down to 3.5
MAD). In OPTO2 only, threshold 4.0 splits the CH19 activity in two:
cluster 4 (-59 uV, 21 Hz, ISI violation ratio 0.12, refractory
contamination 0.05) fires 13 per second at baseline, 174 at onset, 114
through the pulse and 4 in the 100 ms after, with a broad trough and a
large positive afterwave; cluster 3 (-68 uV, 17 Hz, contamination 1.0)
is silenced by the light like the rest of the tetrode. Outside the trains
their cross-correlogram is flat at 0.5 to 1 ms (83 and 107 pairs against
83 expected) while cluster 4's autocorrelogram dips to half of
expectation out to 1.5 ms, so cluster 4 is a credible single unit, and
*conjecture:* the directly driven cell on TT10. The split does not
survive thresholds 3.5 or 3.0, does not appear in OPTO1 at any
threshold, and disappears again when the snippet window is lengthened to
1.5 ms after the trough (`--ms5 snippet_T2=45 --tag t2_45`, meant to give
the afterwave more weight). `results/resort_<experiment>_TT10*.csv`,
`figures/opto_TT10_resort_thr*.png`, `TT10_traces_asis_resort4.pdf`
(OPTO2).

What the cell needs is a sort that uses waveform shape with the amplitude
left free, at a detection threshold near 3.5 MAD, with the refractory
period of the candidate used to reject the hash; MountainSort5's
threshold and PCA clustering only reach it by chance.

## 6. Figures

On the data volume under `<data_root>/TTM-003/<experiment>/figures/`:
`opto_train*.png` (per-train raw and LFP averages, spike-band average,
MUA PSTH), `opto_onset_zoom_*.png` and `opto_onset_single_trials_*.png`
(artifact check), `opto_session_mua.png`, `optotag_asis_population_{A,B}.png`
(z-scored PSTH of every unit), `optotag_asis_unit*_{summary,ccgs,lfp}.png`
(candidates), `opto_population_asis.png` (pooled suppression per train),
`opto_TT6_inpulse_asis.png` (section 5a), `opto_TT10_inpulse_asis.png`
and `opto_TT10_resort_thr*.png` (section 5c), `galleries/*.pdf` (per-pulse
rasters; `_resort<thr>` suffix = that re-sort's clusters drawn in place of
the main sort's on TT6 and TT10).

## 7. What would help next

Short pulses (1 to 5 ms) at several powers including well below 0.5 mW,
for latency-based tagging and a dose-response curve; a control recording
with the fibre disconnected or an opsin-negative animal for the slow
evoked potential; a shape-based, amplitude-free assignment pass on TT4,
TT6 and TT10 (section 5c); and
sleep scoring for the drowsy trains, which the `state` column of the hand
log only approximates.
