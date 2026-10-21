---
title: Light artifacts and optotagging on tetrodes, literature survey
updated: 2026-10-21
status: needs_verification
---

# Light artifacts and optotagging on tetrodes: literature survey

Compiled by a web-research agent on 2026-10-21 for the TTM-003 optogenetic
sessions (`subjects/ttm-003/optogenetic_stimulation.md`): 16 tetrodes, 64 ch,
30 kHz Open Ephys, implanted fibre with a blue laser at 0.5 to 9.8 mW,
207 ms pulses at 1 Hz and 64 ms pulses at 2.86 Hz, SOM-ChR2 mouse. Numbers
are quoted from the cited sources as the agent read them; citations have not
been re-checked by a person, hence `needs_verification`. Gaps are flagged at
the end.

## 1. Nature of the light artifact on metal electrodes

- Mechanism and shape. Kozai and Vazquez 2015 (J Mater Chem B 3:4965,
  doi:10.1039/C5TB00108K, PMC4494685) distinguish photovoltaic (Becquerel),
  photoelectrochemical, photoconductive, photothermal and fluorescence modes.
  For the Becquerel effect at a metal/electrolyte interface the charge
  transfer is transient and non-Faradaic at light on and off; under steady
  illumination no net charge transfer is observed. On metal wires, expect
  onset and offset transients rather than a sustained DC step
  (semiconductor substrates behave differently). Their Fig. 6 gives the
  polarity: the leading artifact is negative and a positive deflection
  follows the end of the pulse. Current is proportional to light intensity
  and exposed surface area; distance from the source attenuates it; blue
  (473 nm) light produces artifacts despite many electrode materials having
  a threshold well below 473 nm. The review does not rank Pt/Ir vs NiCr vs
  gold.
- Time course on metal, in vivo. Cardin et al. 2010 (Nat Protoc 5:247,
  PMC3655719): artifact onset is exactly aligned with light onset, whereas
  evoked spikes lag by several ms; with metal electrodes artifacts lasted at
  least 8 ms after a pulse ended; for a 1 ms pulse the offset decay is about
  2 to 10 ms; size is proportional to total laser power and exposed metal.
  Twisted 12.5 um stereotrode wires gave artifacts 84.6 +/- 16.3 % smaller
  than 75 um tungsten tips; angling electrodes so light does not strike
  them nearly eliminates it. Safe irradiance up to ~75 mW/mm2 for 0.5 to
  50 ms pulses; >100 mW/mm2 sustained (>500 ms) caused damage; 2 mW out of
  a 200 um fibre is ~63.7 mW/mm2.
- Spike-band contamination. Mikulovic et al. 2016 (Neurophotonics 3:015002,
  doi:10.1117/1.NPh.3.1.015002, PMC4717288): square pulses produced sharp
  deflections detected as spikes after 500 to 5000 Hz band-pass filtering;
  artifact power grows about exponentially with intensity (b ~ 1.9 to 2.3
  per mW at 1 to 5 mW); artifacts appear at the stimulation frequency and
  its harmonics; sinusoidal drive produces no spike-like transients.
  Recommends opsin-negative controls.
- Fibre-tetrode geometry. Anikeeva et al. 2012 optetrode (Nat Neurosci
  15:163, PMC4164695): 12 um NiCr gold-plated to 150 to 350 kOhm, tetrodes
  300 to 1000 um below a 200 um fibre; 60 to 160 mW/mm2 at the tip gives
  ~5 to 12 mW/mm2 at 0.5 mm. They report increased multi-unit background
  and changed spike amplitudes during ChR2 trains that degraded isolation.
  Kvitsiani et al. 2013 (Nature 498:363, PMC4349584): tetrodes parallel to
  the fibre, within 500 um of the tip; when an artifact appeared, light
  intensity was lowered; 2 to 6 mW at the tip.
- Low-power local light abolishes it. Stark, Koos and Buzsaki 2012
  (J Neurophysiol 108:349, PMC3434617): diode-probe light under 40 uW gave
  no artifacts, which are often seen at 1 to 10 mW. Roux, Stark, Sjulson
  and Buzsaki (arXiv:2311.10410): photoelectric artifacts distort spike
  waveforms, especially with short pulses; high intensity causes
  synchronous firing and superimposed waveforms that cannot be sorted; with
  interneuron ChR2, post-stimulus rebound spiking may occur. Neuropixels
  Opto (Nat Methods 2026, PMC11838571): surface illumination with sharp
  onsets gives a photoelectric artifact of 1 mV or more; on-shank emitters
  gave ~30 uV, removed by standard preprocessing before Kilosort.

Expectation for 0.5 to 10 mW through an implanted fibre next to tetrodes:
mV-scale, spike-band, onset- and offset-locked transients lasting ~2 to
10 ms, proportional to power and exposed tip area, varying wire to wire.

## 2. Preprocessing and sorting around pulses

- Blanking windows. SpikeInterface `remove_artifacts` defaults
  `ms_before=0.5, ms_after=3.0, mode="zeros"`; modes `zeros` (only for
  zero-centred traces, i.e. after a high-pass), `linear`, `cubic`
  (`fit_sample_spacing`; check the fit visually), `average`/`median`
  (template subtraction, needs `list_labels`; `scale_amplitude`,
  `time_jitter`). `silence_periods` zeros, noise-fills or apodizes whole
  periods. The cibrrig pipeline blanks both onsets and offsets with
  `ms_before=0.125, ms_after=0.25` before KS4 (Neuropixels); tetrodes near a
  fibre need wider windows (Cardin's 2 to 10 ms decay).
- Kilosort. KS4 `artifact_threshold` (default inf) zeroes the entire batch
  containing a value above it; users report it removing real spikes
  (MouseLand/Kilosort#1009). Kilosort#179: a spike-count dip ~5 to 10 ms
  before and after LED onset, scaling with intensity, present in
  opsin-negative controls; the maintainer's advice is to remove the artifact
  from the raw data before Kilosort because filtering spreads it in both
  directions. KS2 crashed on per-pulse laser artifacts not removed by
  whitening (#154). For <= 64 channels: `nblocks=0`, larger `batch_size`,
  `nearest_templates` <= channel count, `dminx` = median lateral spacing.
  Kilosort paper: CAR after high-pass filtering.
- Referencing. Kozai and Vazquez: common-average referencing of identical,
  identically illuminated electrodes can dramatically attenuate the
  Becquerel effect but not remove it; also PCA component removal or
  regression on an out-of-brain reference. Illumination differs across
  tetrodes, so a global 64-channel median leaves residuals on the
  near-fibre tetrodes; a per-tetrode median is closer to identically
  illuminated but subtracts common spike energy.
- Joint vs separate sorting. Every pipeline found sorts stim and non-stim
  data together after artifact suppression, then checks waveforms: the
  microLED optoelectrode study (Nat Commun 2020, PMC7188816) uses
  spontaneous and light-induced waveforms without discrimination; Neuropixels
  Opto likewise. Anikeeva re-assessed cluster quality during illumination
  (46 of 50 PV units kept isolation). CellExplorer sorts everything but
  excludes `*.manipulation.mat` intervals from waveform, ACG, monosynaptic
  and PSTH metrics by default; no built-in statistical tagging test.

## 3. Optotagging criteria

| Source | Pulse | Test | Thresholds |
|---|---|---|---|
| Kvitsiani 2013 (tetrodes, PV/SOM, ACC) | 1 to 2 ms at 4 to 100 Hz, 2 to 6 mW | SALT: first-spike latency distribution in a 10 ms post-light window vs equal-length baseline epochs (Jensen-Shannon divergence; Supplement) | P<0.01; waveform r>0.85 spontaneous vs evoked; latency <4 ms, low jitter; isolation distance >=18, L-ratio <=0.2 |
| Cohen 2012 (Nature 482:85, PMC3271183) | 20 and 50 Hz trains | waveform correlation vs response energy | r>0.9; within a few ms; follows 50 Hz |
| Pi 2013 (PMC4017628) | 1 ms, 1.5 to 3 mW | SDF half-peak epoch, Mann-Whitney | P<0.05; waveform r in the top 5 % of a bootstrap; ID>=20, L-ratio<=0.1; latency 2.5 +/- 0.2 ms, jitter 1.4 +/- 0.2 ms, reliability 0.62 |
| Stark 2012 (diode probes) | 2 to 50 ms, 3 to 180 uW | Mann-Whitney stim vs baseline bins; Poisson/Bonferroni for latency | P<0.05; median latency 5.5 ms, jitter 2 ms, spike probability 0.7; waveform consistency median 0.98 with a label-permutation test |
| Allen SDK optotagging notebook | 10 ms single, 2.5 ms at 10 Hz, 1 s raised cosine | evoked 1 to 9 ms vs baseline -10 to -2 ms | evoked/(baseline+1) > 2; ignore spikes at 0 to 0.5 and 9.5 to 10 ms (likely artifactual); at least 1 ms to generate a real light-evoked spike; Sst units fire at a range of latencies with mixed FS/RS waveforms |
| Neuropixels Opto (Allen) | 10 ms at 20 Hz, 100 uW | rate increase on >= 4 of 5 pulses, Holm-Sidak | p<0.05; latency <8 ms; reliability >0.3 |
| Lima 2009 PINP (PLoS ONE 4:e6099) | 3 to 10 ms, 17 mW | follows a 5 Hz train; NBQX control | direct units 2.8 +/- 1.3 ms; synaptic units fail to follow trains; latency alone did not separate direct from synaptic |

Direct vs synaptic: Roux/Buzsaki note indirect activation is minor for
GABAergic ChR2 (no excitatory output) but rebound spiking occurs; Lima's
train-following test and collision/NBQX logic are the only hard
discriminators found. Network effects of SOM activation: Kvitsiani's
suppression criterion is firing <50 % of baseline within 100 ms,
Mann-Whitney P<0.05, with SOM inhibition longer and more variable than PV.
Funk et al. 2017 (J Neurosci 37:9132, PMC5607463): 20 and 50 ms pulses,
10 mW, in M2 via optrode; OFF period = >=100 ms of pooled-MUA silence from
pulse onset; evoked NREM OFF 190 +/- 32 ms vs spontaneous 174 +/- 10 ms;
all 72 single units silenced, outlasting the pulse; PV gave more transient
silencing; median first post-pulse spike latency SOM+ 128 +/- 48 ms vs PV+
59 +/- 19 ms (tagged-unit latencies after the suppression, not
direct-activation latencies). Niethard et al. 2018 (PNAS) is two-photon
imaging: SOM activity rises before the down state and falls during up
states.

## 4. Long (200 ms) pulses: desensitization, depolarization block, rebound

- Herman et al. 2014 (eLife, PMC3904216): at ~40 mW/mm2, 20 Hz trains,
  regular-spiking SST interneurons fired an initial spike or short burst
  then entered depolarization block at 49 ms pulses, while fast-spiking SST
  interneurons kept firing; in vivo, RS-SST silenced by long pulses let
  presumptive pyramidal firing return to baseline; median latencies 1.8 ms
  (RS-SST) and 6 ms (FS-SST); block occurred regardless of intensity; they
  recommend short pulses for interneurons.
- Mattis et al. 2012 (Nat Methods 9:159, PMC4165888): ChR2 steady-state to
  peak ratio 0.30 +/- 0.01 at 1 s of light; tau_off 11.6 +/- 0.4 ms; H134R
  recovers from desensitization more slowly. In a 200 ms pulse the
  photocurrent falls to ~30 % within tens of ms: expect an onset burst then
  adaptation or silence in RS-SOM cells.
- Rebound: Li et al. 2019 (eLife 48622): rebound after 500 ms
  photoinhibition ~10 % of baseline (1.5 to 7 mW), >30 % for 4 s pulses,
  decaying over hundreds of ms; they ramp intensity down over the last 100
  to 200 ms. A 2025 V1 study (PMC12369306): with 500 to 1000 ms pulses,
  PV-driven pyramidal rebound ~2x SOM-driven; directly activated SOM/PV
  cells show post-stimulus suppression 10 to 100 ms after offset, opposite
  in sign to the pyramidal rebound. A tagged SOM cell should be depressed,
  not rebounding, right after offset.

## 5. Open Ephys TTLs and alignment

Open Ephys synchronization docs: align events by sample number (subtract
the first continuous `sample_numbers` value); unsynchronized streams carry
timestamps of -1. Acquisition Board digital inputs share the FPGA clock, so
a laser TTL into the board's digital port needs no sync line.
open-ephys-python-tools exposes `sample_number`, `timestamp`, `line`,
`state` (1 rising, 0 falling). Only one event is recorded if two land on
the same sample. Measure the laser-vs-TTL delay (shutter/AOM) rather than
assuming zero; per Cardin the artifact itself is the best onset fiducial.

## Recommended pipeline (tetrodes, SpikeInterface, KS4 or MS5)

1. Characterize first: average raw 30 kHz traces on each wire around TTL
   rising and falling edges, lowest-power train as the near-control. Record
   per-wire artifact amplitude, polarity and decay to the noise floor.
2. Events: `sample_number` from `events/.../TTL`, state 1 onset, 0 offset;
   verify pulse widths and the artifact latency relative to the TTL.
3. Preprocess: high-pass (300 Hz), `remove_artifacts(mode="cubic" or
   "linear")` on onsets and offsets with windows from step 1 (start ~0.5 ms
   before, extend to the measured decay), then common median reference.
   Keep KS4 `artifact_threshold=inf`.
4. Sort everything jointly; KS4 with `nblocks=0`, `dminx`/`dmin` from the
   tetrode geometry, larger `batch_size`; MS5 needs no drift settings.
   Cross-check waveform r between stim-epoch and non-stim spikes per unit
   (>= 0.85 to 0.9) and isolation metrics within stim epochs.
5. Tag with the short-pulse trains primarily: SALT-style first-spike
   latency test (10 ms window, P<0.01) plus latency <4 to 8 ms, jitter
   ~1 to 2 ms, reliability >0.3 to 0.6, excluding spikes in the blanked
   windows. Treat 200 ms pulses as confirmatory: onset burst then
   adaptation, sustained firing only in FS-SOM; after offset a tagged SOM
   cell should be depressed while untagged units may rebound.
6. Network check: MUA or pyramidal suppression on the same drive
   (Kvitsiani <50 % baseline within 100 ms; Funk >=100 ms OFF periods);
   units that follow the suppression with long latency are candidates for
   synaptic or rebound responders, not tagged cells.

## Gaps

No source quantified the artifact on NiCr vs Pt-Ir vs gold specifically;
the SALT Jensen-Shannon details are in the Kvitsiani 2013 Supplement (not
retrieved); the "Zwitter test" returned nothing and may be a misremembered
name; CellExplorer has no built-in tagging statistic.
