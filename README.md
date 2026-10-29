<div align="center">

# An Experimental Study of BMOCZ over Underwater Acoustic Channels

**A MATLAB modem for Binary Modulation On Conjugate-reciprocal Zeros (BMOCZ), simulated and then tested in a real underwater acoustic channel**

![MATLAB](https://img.shields.io/badge/MATLAB-R2023a%2B-orange?logo=mathworks)
![Toolboxes](https://img.shields.io/badge/Toolboxes-Comms%20%7C%20Signal%20%7C%20DAQ%20%7C%20Image-blue)
![Hardware](https://img.shields.io/badge/Hardware-NI%20DAQ-lightgrey)
![Status](https://img.shields.io/badge/status-completed%20project-success)

Communication Lab, Faculty of Electrical & Computer Engineering, Technion – Israel Institute of Technology<br>
Presented at the TAU Underwater Acoustics Symposium 2026

</div>

---

## Table of contents

1. [About the project](#about-the-project)
2. [Key findings](#key-findings)
3. [Modem architecture](#modem-architecture)
4. [Repository structure](#repository-structure)
5. [Requirements](#requirements)
6. [Quick start: simulation](#quick-start-simulation)
7. [Running a pool experiment](#running-a-pool-experiment)
8. [Function reference](#function-reference)
9. [Parameters](#parameters)
10. [Results](#results)
11. [Known issues](#known-issues)
12. [Team](#team)
13. [References](#references)
14. [Citation](#citation)


---

## About the project

Underwater acoustic channels are among the hardest channels in communications. They are slow, narrow-band, and full of multipath echoes, so coherent receivers need long pilots and channel estimation before they can decode anything. For **short packets** (sensor readings, commands, acknowledgments) that overhead can cost more than the payload itself.

**MOCZ** (Modulation On Conjugate-reciprocal Zeros) avoids channel estimation altogether. Information is carried in the **zeros of a polynomial** whose coefficients are the transmitted samples. A linear time-invariant channel multiplies the signal's z-transform by its own, which adds new zeros but **cannot move the message zeros**. The receiver can therefore decode **blindly**, with no pilots and no channel state information.

Previous work studied MOCZ in theory, in simulation and over the air with SDRs. To the best of our knowledge, **this project is the first to test BMOCZ in a real underwater acoustic channel.** We:

- built a complete BMOCZ modem in MATLAB, covering pulse shaping, Zadoff–Chu synchronization, and three decoders (**DiZeT**, **Greedy**, **ML**);
- validated it in simulation over multipath channels;
- built two MATLAB App Designer GUIs that drive NI DAQ hardware: one to transmit and record, one to inspect received packets;
- ran pool experiments that swept SNR, payload size `K`, and guard length `L`, and sounded the channel with chirps;
- found and explained an **asymmetric bias in DiZeT and ML decoders** that appears when the guard interval is longer than the channel needs.

## Key findings

| # | Finding |
|---|---------|
| 1 | **BMOCZ works blindly underwater.** With a well-chosen guard length, packets were decoded in the pool at low BER without any channel estimation. |
| 2 | **The guard length `L` is a design parameter, not a safety margin.** It must cover the channel's delay spread, `L ≥ L_ch − 1 ≈ ⌈T_d / T⌉`, or the convolution tail is cut off and the zeros move. Making it *longer* than needed hurts too. |
| 3 | **Excess guards cause an asymmetric bias.** Every extra noise-only sample adds a spurious "noise zero" near the unit circle. Bits sent as `0` (zero inside the unit circle) are increasingly decoded as `1`, while bits sent as `1` stay almost error-free. Trimming the decoding window back to `K + L_ch` removes the bias. |
| 4 | **ML beats DiZeT in simulation, but DiZeT is close in practice.** Over a harsh 5-tap simulated channel, ML is about 6 dB better. In the pool, both decoders follow the same trend versus `L`. |
| 5 | **PAPR stays bounded as `K` grows.** The 1% PAPR stays below about 9 dB for `K = 7…12`, because the zero radius `R` approaches 1. |
| 6 | **Received gain does not matter.** Attenuation and spreading loss scale every coefficient equally, and zeros depend only on coefficient ratios. Both DiZeT and ML are therefore invariant to a constant channel gain. |


## Modem architecture

![BMOCZ modem scheme](Simulation%20Modem/BMOCZ%20Communication%20Modem%20Scheme.png)

**Frame format.** One packet carries `B` BMOCZ blocks. Each sample lasts `Tsym = 1 ms`.

```
| Zadoff–Chu (N=63, u=2) | L zeros | x_1 (K+1) | L zeros | x_2 (K+1) | L zeros | ... | x_B (K+1) | L zeros |
```

**Transmitter.**
1. **Bits → zeros → coefficients** (`generateFunc.m`).
2. **Raised-cosine pulse shaping** (`pulseSHP.m`, `rcosdesign(β=1, span=6, sps=200)`).
3. **Up-conversion** to the carrier.
4. **Scaling** to the DAQ voltage range (±5 V).

**Receiver.**
1. **Down-conversion.**
2. **Matched filter.**
3. **Zadoff–Chu cross-correlation** for exact timing. The receiver warns if the peak-to-mean ratio is below 15.
4. **Sampling** once per symbol.
5. **Reshaping** into `B` blocks of `K + 1 + L` samples.
6. **Optional truncation** to `n_trun` samples.
7. **DiZeT, ML or Greedy decoding.**

## Repository structure

```
MOCZ_project/
├── README.md
├── Simulation Modem/
│   ├── BMOCZ Communication Modem Scheme.png
│   ├── Sim Code/
│   │   ├── berSim.m              ← entry point: BER vs Eb/N0 for DiZeT, Greedy and ML
│   │   ├── PAPRstatistics.m      ← entry point: PAPR CCDF and PAPR vs K
│   │   ├── makeParams.m          ← all modem parameters, plus precomputed ML/Greedy matrices
│   │   ├── generateFunc.m        ← random packet → passband signal
│   │   ├── pulseSHP.m            ← raised-cosine modulate/demodulate
│   │   ├── MOCZRealDecoding.m    ← simulated channel + AWGN + full receiver chain
│   │   ├── DiZeT.m · Greedy.m · ML.m   ← decoders
│   └── Sim Results/              ← figures (see Results)
└── Experiments/
    ├── Experiments GUIs/
    │   ├── Guis Code/
    │   │   ├── Apps/
    │   │   │   ├── Pool_transmit_GUI.mlapp   ← transmit + record with NI DAQ
    │   │   │   └── recieveGui.mlapp          ← inspect received packets
    │   │   ├── makeParams.m                  ← experiment defaults (K=5, L=25, fc=25 kHz, ±5 V)
    │   │   ├── generateFunc.m                ← packet generator (random / all-0 / all-1 / image)
    │   │   ├── generate_TX_signal_pool_5Vm.m ← builds the full set of TX signals for an experiment
    │   │   ├── generate_TX_signal_pool_10Vm.m
    │   │   ├── decode_mocz.m                 ← hardware receiver: sync, sample, truncate, decode, BER
    │   │   ├── runMLandGparams.m             ← rebuilds ML/Greedy matrices for a given K, L
    │   │   ├── chirpFnc.m                    ← chirp channel sounding (generate / decode / simulate)
    │   │   ├── image_tx_rx.m                 ← send a 256×256 binary image over BMOCZ
    │   │   ├── pulseSHP.m · DiZeT.m · Greedy.m · ML.m
    │   └── Guis Images/              ← annotated screenshots of both GUIs
    └── Experiments Results/
        ├── Results Analyzing Code/
        │   ├── AnalyzeScript_L_changes_Channel_Profilling.m  ← BER vs L, split by bit value
        │   └── AnalyzeScript_SNR_changes.m                   ← BER vs SNR at several effective L
        └── Results Graphs/           ← figures (see Results)
```

> **Note:** the simulation and experiment folders each contain their own copies of `makeParams.m`, `generateFunc.m` and the decoders, and the copies differ. Keep only one of the two folders on the MATLAB path at a time.

## Requirements

| Item | Needed for |
|------|-----------|
| MATLAB **R2023a or newer**. The GUIs were built in R2023a and R2023b. | Everything |
| Communications Toolbox (`zadoffChuSeq`, `awgn`) | Everything |
| Signal Processing Toolbox (`rcosdesign`, `xcorr`, `findpeaks`, `pspectrum`, `chirp`) | Everything |
| Statistics and Machine Learning Toolbox (`ff2n`) | Packet generation, PAPR |
| Data Acquisition Toolbox + NI-DAQmx support package | Pool GUI (hardware only) |
| Image Processing Toolbox (`imbinarize`, `imresize`) | Image transmission only |
| NI DAQ device(s) with `ao0` output and `ai0`/`ai1` inputs, an underwater transducer and hydrophone(s), and a power amplifier as needed | Pool experiments |

The simulation needs no hardware.

## Quick start: simulation

```matlab
cd 'Simulation Modem/Sim Code'
berSim            % BER vs Eb/N0 (0–18 dB) for DiZeT, Greedy and ML
PAPRstatistics    % PAPR CCDF for K = 5 and 12, and PAPR at 1% vs K
```

**Things you can change:**

- **Payload, guards, packet size.** Edit `K`, `L` and `B` in `makeParams.m`. ML runs only when `K ≤ 8`, because it searches all `2^K` messages.
- **Channel.** Edit `h_channel` in `MOCZRealDecoding.m`. The default `[0.8 −0.45 0.2 −0.1 0.003]` is the "harsh" channel; `[0.8 −0.1 0.003]` is the "soft" one. Keep `L ≥ length(h_channel) − 1`.
- **Constellation plots.** Set `figFlag = 1` in `makeParams.m` to plot the received zeros against the transmitted ones for every block. Use a small `B` and a single SNR when you do this.
- **Run time.** It is controlled by `attemptNum` in `berSim.m` (200 + 30·i packets at SNR step i).

**Using the modem in your own script:**

```matlab
p = makeParams();                                        % default modem
[P, x_pb] = generateFunc(p);                             % P: K×B bits, x_pb: passband signal
[Pd, Pg, Pm] = MOCZRealDecoding(p, 10, x_pb, P);         % channel + 10 dB AWGN + decode
ber_DiZeT = mean(P(:) ~= Pd(:));
```

## Running a pool experiment

All of the following runs from `Experiments/Experiments GUIs/Guis Code` (add the `Apps` subfolder to the path).

### 1. Generate the transmit pool

```matlab
generate_TX_signal_pool_5Vm
```

The script creates `TX_signals/fiveVmax/MOCZ_packets/` with one folder per sweep. Each folder holds 5 runs per setting, and each file contains `signal_pb_total`, `P`, `x_bb` and `simParams`.

| Folder | Sweep | Fixed |
|--------|-------|-------|
| `K_changes` | `K = 5…20` | `L = 25` |
| `L_changes` | `L = 1…50` | `K = 5` |
| `SNR_changes` | TX power × 0.5…1.0 (21 levels) | `K = 5`, `L = 25` |
| `(N_over_K)Ratio_changes` | `L = 5, 10, …, 50` | `K = 5` |

You can also create a single signal from the transmit GUI with **Generate Signal**: choose `R`, `K`, `L` and the packet type (all-0, all-1 or random).

To create the chirp used for channel sounding (50 ms, logarithmic sweep from 20 to 45 kHz), run `chirpFnc(1)`. It saves `Chirp_Tx.mat`.

### 2. Transmit and record: `Pool_transmit_GUI`

![Transmit GUI](Experiments/Experiments%20GUIs/Guis%20Images/Overview%20of%20the%20GUI%20used%20for%20Transmitting.png)

1. **DAQ setup (1).** Choose *Single* (one device for TX and RX) or *Dual* (separate TX and RX devices), and set their indices from `daqlist`.
2. **Queue (6).** Load a single `.mat` or a whole folder from step 1.
3. **Chirp (5).** Optionally tick *Transmit Chirp for test* and point it to `Chirp_Tx.mat`.
4. **Output (7).** Choose the recordings folder.
5. **Session details (8).** Enter Tx–Rx distance, location (*middle of pool*, *near wall*, *near wall and surface*), pause between packets, and extra listening time.
6. **Run (2).** Press **START TRANSMIT**. Progress appears in field 3.

Each packet is saved as `<signal name>_recorded_result.mat`, holding `recordedData` (a timetable with `ai0` and `ai1`) and `currentParams` (all modem parameters, the transmitted bits `P`, location and distance). A session log goes to `transmission_session.txt`.

The GUI has two extra panels:

- **Analyze Recording (9)** decodes a recording on the spot, with optional truncation, and shows the BER.
- **Image transmission (4)** sends a 256×256 black-and-white image with `K = 8` and shows the decoded image.

### 3. Inspect received packets: `recieveGui`

![Receive GUI](Experiments/Experiments%20GUIs/Guis%20Images/Overview%20of%20the%20GUI%20used%20for%20Inspecting%20received%20packets.png)

1. **Load (1).** Load a `*_recorded_result.mat`.
2. **Message data (2).** Shows `L`, `K`, `B`, location, and the BER of DiZeT and ML.
3. **Selection (3).** Choose the RX channel and which of the `B` packets to examine.
4. **Truncation (4).** Set how many samples the decoder uses, from `K + 1` to `K + 1 + L`. The GUI suggests the optimal truncation, defined as the shortest window with the lowest DiZeT BER.
5. **Constellation and weights (5).** Received zeros vs. transmitted zeros in the z-plane, with the ML weight matrix $(V^HV)^{-1/2}$ for the examined packet or averaged over the packet.
6. **Signal views (6).** The passband signal, its spectrum and its spectrogram.

### 4. Decode from the command line

```matlab
S = load('signal_Vm5_K5_L11_Power100_run1_recorded_result.mat');
p = runMLandGparams(S.currentParams, S.currentParams.K, S.currentParams.L);
x = S.recordedData.Variables;  x = x(:,1);               % first hydrophone (ai0)
[bits, BER] = decode_mocz(x, p, 1, p.K + 1 + p.L);        % 1 = DiZeT, 2 = ML
```

### 5. Batch analysis

1. Set `folderPath` at the top of the relevant script to your recordings folder. It currently points to the lab's network share.
2. Run one of the scripts:
   - `AnalyzeScript_L_changes_Channel_Profilling.m` gives the aggregate BER vs `L` and the BER split by the transmitted bit value. It averages both receivers and all 5 runs.
   - `AnalyzeScript_SNR_changes.m` re-decodes every recording from one physical `L = 25` at several *effective* guard lengths (`n_trun = K + 1 + L_eff`, default `L_eff = 1, 4, 15`) and plots BER vs relative SNR.
3. To profile the channel, decode a recorded chirp with `chirpFnc(2, 'chirp_recorded_result.mat')`. It plots the cross-correlation (the channel's delay profile) and the zeros of the estimated channel.

## Function reference

| Function | Signature | Purpose |
|----------|-----------|---------|
| `makeParams` | `p = makeParams()` | Builds the parameter struct and precomputes `V_bank` (Greedy) and `V_all`/`Cinv_all` (ML) |
| `runMLandGparams` | `p = runMLandGparams(p, K, L)` | Rebuilds the decoder matrices after changing `K` or `L` |
| `generateFunc` | `[P, x_pb, x_bb, gd, p] = generateFunc(p)` | Random unique messages → ZC + guarded BMOCZ blocks → passband |
| `pulseSHP` | `y = pulseSHP(x, p, 'modulate' or 'demodulate')` | Raised-cosine shaping plus carrier mixing |
| `MOCZRealDecoding` | `[Pd, Pg, Pm] = MOCZRealDecoding(p, SNR, x_pb, P)` | *(Simulation)* channel + AWGN + full receiver, all three decoders |
| `decode_mocz` | `[P̂, BER, Y, W] = decode_mocz(x_pb, p, mode, n_trun)` | *(Hardware)* receiver; `mode` 1 = DiZeT, 2 = ML; truncates each block to `n_trun` samples |
| `DiZeT` | `m = DiZeT(R, theta, y, K)` | Zero-testing decoder, O(K·N) |
| `Greedy` | `m = Greedy(y, p)` | DiZeT-initialized greedy subspace decoder |
| `ML` | `[m, W] = ML(y, p)` | Exhaustive whitened-projection ML decoder, 2^K hypotheses |
| `chirpFnc` | `chirpFnc(mode, path)` | 1 = generate chirp, 2 = decode recorded chirp, 3 = simulate a 5-tap channel |
| `image_tx_rx` | `[out, p] = image_tx_rx('tx' or 'rx', data, p)` | Image ↔ BMOCZ signal (`K = 8`, 256×256 binary) |

## Parameters

| Symbol | Field | Simulation | Experiment | Meaning |
|--------|-------|-----------|------------|---------|
| `K` | `K` | 7 | 5 (swept 5–20) | Bits per block; the polynomial has degree `K` |
| `L` | `L` | 5 | 25 (swept 1–50) | Zero guards after each block |
| `B` | `B` | 20 | 20 | Blocks per packet |
| `λ` | `lambda` | 1 | 1 | Sets the zero radius `R = sqrt(1 + 2λ sin(π/K))` |
| `N_ZC`, `u` | `zadoffChuPair` | 63, 2 | 63, 2 | Sync sequence length and root |
| `T` | `Tsym` | 1 ms | 1 ms | Symbol duration |
| `f_c` | `fc` | 15 kHz | 25 kHz | Carrier frequency |
| `F_s` | `Fs` | 200 kHz | 200 kHz | Sampling rate (DAQ rate) |
| `β` | `betha` | 1 | 1 | Raised-cosine roll-off (span 6, 200 samples per symbol) |
| — | `a`, `b` | ±10 (hard-coded) | ±5 V | Output voltage range |
| — | `packet_type` | random | 0, 1 or 2 | All-0, all-1 or random bits |

**Choosing `L`:** take the channel's delay spread `T_d` from a chirp sounding and use `L ≈ ⌈T_d / T⌉`. Don't add a large safety margin; see the results below.

## Results

### Simulation

| BER vs Eb/N0, harsh 5-tap channel | Zeros: soft 3-tap channel | Zeros: harsh 5-tap channel |
|:---:|:---:|:---:|
| ![BER sim](Simulation%20Modem/Sim%20Results/BER%20vs%20EbN0%20of%203%20decoders%20with%20harsh%205%20tap%20channel.png) | ![soft](Simulation%20Modem/Sim%20Results/Constellation%20exmp%20soft%203%20tap%20channel%20and%20SNR%20-3.59%20dB.png) | ![harsh](Simulation%20Modem/Sim%20Results/Constellation%20exmp%20harsh%205%20tap%20channel%20and%20SNR%20-3.59%20dB.png) |

- **Decoder ranking.** Over the harsh channel (`K = 7`, `L = 5`), ML reaches a BER of about 10⁻⁵ near 9 dB. Greedy and DiZeT need about 15–16 dB. Greedy improves on DiZeT mostly at low SNR.
- **Constellations.** At −3.6 dB SNR, the message zeros still sit at radius `R` or `1/R`. The channel zeros appear elsewhere in the plane.

| PAPR CCDF, K = 5 vs K = 12 | PAPR at 1% vs K | PAPR at 15% vs K |
|:---:|:---:|:---:|
| ![ccdf](Simulation%20Modem/Sim%20Results/CCDF%20of%20PAPR%20for%20two%20polynomial%20orders%20%28payload%20sizes%29.png) | ![papr1](Simulation%20Modem/Sim%20Results/PAPR%20at%200.01%20CCDF%20probability%20as%20a%20function%20of%20block%20size.png) | ![papr15](Simulation%20Modem/Sim%20Results/PAPR%20at%200.15%20CCDF%20probability%20as%20a%20function%20of%20block%20size.png) |

PAPR does not grow with the block size: the 1% PAPR stays between about 7.9 and 9 dB for `K = 7…12`.

### Pool experiments

**Channel sounding.** Near the pool walls, the chirp cross-correlation shows a direct path followed by strong echoes about 1.5–3.5 ms later. This delay spread sets the guard length the channel needs.

<p align="center"><img src="Experiments/Experiments%20Results/Results%20Graphs/Chirp%20cross%20correlation%20in%20harsh%20channel%20-%20near%20pool%20walls.png" width="520"></p>

**BER vs guard length: the asymmetric bias.**

| Aggregate BER vs L | BER vs L, split by transmitted bit |
|:---:|:---:|
| ![agg](Experiments/Experiments%20Results/Results%20Graphs/BER_aggregate_vs_L.png) | ![split](Experiments/Experiments%20Results/Results%20Graphs/BER_vs_L_bitSplit.png) |

- **Short guards are best.** The BER is lowest for the shortest guards and rises once `L` exceeds what the channel needs.
- **The errors are one-sided.** Splitting by the transmitted bit shows that almost all errors are `0 → 1`. Bits sent outside the unit circle are almost never wrong, while bits sent inside it approach and exceed 50% error.
- **Why it happens.** Each noise-only guard sample adds a spurious zero on a ring of radius $\rho \approx \mathrm{SNR}^{1/(2L)}$, just outside the unit circle. As `L` grows the ring approaches `|z| = 1`, and the test point $R^{-1}e^{j\theta_k}$ stops looking special to DiZeT and ML. The trend is the same for both decoders.

**BER vs SNR at different effective guard lengths.** These recordings come from a difficult channel. The same recordings (physical `L = 25`) were truncated at the receiver:

<p align="center"><img src="Experiments/Experiments%20Results/Results%20Graphs/BER_vs_SNR_by_effectiveL_DiZeT.png" width="620"></p>

With `L_eff = 1` or `4` the BER stays at about 0–2% over the whole 3 dB range. With `L_eff = 15` it is stuck at 6–9% no matter the SNR. **Trimming the window at the receiver removes most of the bias without retransmitting.**

**What the decoder sees.**

| Suitable guard length | Many excess guards | ML weight matrix $(VV^H)^{-1}$ |
|:---:|:---:|:---:|
| ![good](Experiments/Experiments%20Results/Results%20Graphs/Symbols%20in%20Z-plane%20under%20suitable%20conditions.png) | ![bad](Experiments/Experiments%20Results/Results%20Graphs/Symbols%20in%20Z-plane%20with%20many%20excess%20guard%20units.png) | ![W](Experiments/Experiments%20Results/Results%20Graphs/Heatmap%20Example%20of%20a%20Weight%20Matrix%20%28VV%5EH%20%29%5E%28-1%29.png) |

- **Suitable guards.** The received zeros (blue) land on the transmitted ones (green).
- **Excess guards.** A ring of noise zeros forms around the unit circle, crowding the inner candidates.
- **Weight matrix.** The ML weight matrix is close to diagonal. This explains why DiZeT, which ignores the off-diagonal terms, performs close to ML.

## Known issues

- **`chirpFnc.m`, modes 2 and 3.** The peak-extraction section uses `corr_positive`, which is never defined. Add `corr_positive = corr_out(lags >= 0);` (and the matching `lags`) before the threshold line.
- **Hard-coded paths.** The analysis scripts point to a lab network share. Edit `folderPath` before running them.
- **Expected input layout.** The receive GUI and the analysis scripts assume the recording timetable's first data column is the hydrophone of interest (`Dev4_ai0`).
- **ML cost.** ML is limited to `K ≤ 8` because its cost grows as `2^K`.

## Team

**Students**
- O. Simchy
- S. Marom

**Supervisors**
- D. Dikarov
- A. Sinai
- M. Zayde
- A. Amar

The project was carried out at the Communication Lab, ECE Faculty, Technion.

**Acknowledgments.** We thank the Communication Lab staff for the hardware and pool access, and the organizers of the **Tel Aviv University Underwater Acoustics Symposium 2026**, where this work was presented.

## References

[1]  P. Walk, P. Jung, and B. Hassibi, “MOCZ for Blind Short-Packet Communication: Basic Principles,” IEEE Transactions on Wireless Communications, vol. 18, no. 11, pp. 5080–5097, Nov. 2019.
[2]  D. A. Huffman, “The Generation of Impulse-Equivalent Pulse Trains,” IRE Transactions on Information Theory, vol. 8, no. 5, pp. 10–16, Sept. 1962.
[3]  M. H. Ackroyd, “The Design of Huffman Sequences,” IEEE Transactions on Aerospace and Electronic Systems, no. 6, pp. 790–796, 1970.
[4] B. Conrey, A. Granville, B. Poonen, and K. Soundararajan, “Zeros of Fekete polynomials,” Ann. Inst. Fourier (Grenoble), vol. 50, no. 3, pp. 865–889, 2000.
[5] “Zadoff–Chu sequence,” Wikipedia, The Free Encyclopedia. [Online]. Available: https://en.wikipedia.org/wiki/Zadoff%E2%80%93Chu_sequence.
[6] O. Schaul, A. Amram, M. Zeyde, and A. Amar, “Greedy Blind Equalization for Binary Modulation on Conjugate-Reciprocal Zeros,” in EUSIPCO 2026, pp. 2066–2070, 2026.
[7] P. Huggins and A. Şahin, “Jutted BMOCZ for Non-Coherent OFDM,” IEEE Transactions on Wireless Communications, vol. 25, pp. 17864–17878, May 2026, doi: 10.1109/TWC.2026.3692740.


## Citation

If you use this code or these results, please cite:

```bibtex
@misc{simchy_marom_2026_bmocz_underwater,
  author       = {O. Simchy and S. Marom},
  title        = {An Experimental Study of {BMOCZ} over Underwater Acoustic Channels},
  howpublished = {Project report, Communication Lab, Technion -- Israel Institute of Technology. Presented at the TAU Underwater Acoustics Symposium},
  year         = {2026},
  note         = {Supervisors: D. Dikarov, A. Sinai, M. Zayde, A. Amar. \url{https://github.com/orsimchy7/MOCZ_project}}
}
```

