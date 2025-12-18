# Drone Communication Anti-Jamming Simulation

## Overview
This MATLAB project simulates a robust uplink communication scenario between a drone and a ground receiver equipped with a **Uniform Linear Array (ULA)**. The simulation demonstrates advanced **Anti-Jamming (AJ)** capabilities using beamforming techniques, focusing on identifying and nulling strong interference while maintaining the integrity of the desired signal under realistic hardware constraints.

The system features a complete processing chain: from **OFDM signal generation** and **channel modeling**, through **DOA estimation** and **Target Identification**, to **Spatial Smoothing Beamforming** and **Hardware Quantization modeling**.

## Key Features

### 📡 Signal & Interference Modeling
* **OFDM Waveform:** Generates QPSK-modulated OFDM signals with a cyclic preamble structure for robust synchronization and correlation.
* **Packet Structure:** Supports frame-based processing with configurable guard intervals and packet lengths.
* **Versatile Jammer Generator:**
    * `Spot`: Narrowband burst interference (time-gated).
    * `Barrage`: Wideband noise jamming.
    * `Sweep`: Frequency-hopping interference.
    * `MultiTone`: Multiple simultaneous CW tones.
    * `CW`: Continuous Wave interference.

### 🎛️ Array Processing & Beamforming
* **DOA Estimation:** Utilizes **MVDR** (Minimum Variance Distortionless Response) for initial angle-of-arrival detection.
* **Spatial Smoothing:** Implements `pc_beamformer_ss` to handle coherent signals and improve covariance matrix rank, ensuring effective nulling even with correlated interference.
* **Principal Component (PC) Beamformer:** Projects the signal onto a subspace orthogonal to the interference.

### 🧠 Smart Target Identification
* **Preamble Correlation:** Distinguishes the desired user from jammers (even when $SIR < 0$ dB) by correlating beamformed outputs with a known preamble sequence.
* **Raw Data Buffering:** Uses a sliding buffer mechanism to detect preambles that cross frame boundaries.

### ⚙️ Hardware Realism
* **Quantization Modeling:** Simulates the effects of limited-resolution phase shifters and attenuators (Digital Step Attenuators).
    * Configurable bits for **Phase** (e.g., 6-8 bits).
    * Configurable bits for **Gain** (e.g., 5-6 bits).

## File Structure

### Root
* `main.m`: **Entry point.** Configures global parameters (`gl_params`), executes the frame-based simulation loop, and visualizes performance (SNR, Spectrum, Time-domain).

### 📁 /signals
* `generate_ofdm_signal_multi.m`: Generates the baseband OFDM signal (Packets + Preamble + Guard). Uses a cyclic QPSK preamble vector for low PAPR.

### 📁 /jammers
* `generate_jammer.m`: Factory function for creating interference.
* `generate_spot_jammer.m`, `generate_barrage_jammer.m`, etc.: Specific jammer implementations.

### 📁 /utilities
* `pc_beamformer_ss.m`: **Core Algorithm.** Principal Component Beamformer with Spatial Smoothing.
* `select_desired_doa_by_preamble.m`: Logic for selecting the correct DOA using matched filtering on buffered data.
* `quantize_weight_vector.m`: Simulates hardware quantization effects on the beamforming weights.
* `fix_angle_indexing_5.m`: Tracker to maintain consistent DOA indexing across frames.
* `steering_vec_ula.m`: Generates the steering vector for the ULA.

### 📁 /TESTS
Standalone scripts for performance sensitivity analysis:
* `Test_Frame_Duration.m`: Analyzes output SNR vs. processing frame duration.
* `Test_Preamble_Duration.m`: Analyzes performance vs. preamble length.
* `Test_Quantization_vs_Osnr.m`: Evaluates the degradation caused by limited phase/gain bits.
* `Test_Sir_vs_Osnr.m`: Tests the system's limit by sweeping input SIR (from -30dB to +10dB).

## Getting Started

1.  **Prerequisites:**
    * MATLAB (R2021b+ recommended).
    * Phased Array System Toolbox.

2.  **Configuration (`gl_params` in `main.m`):**
    Adjust the simulation parameters at the top of the script:
    ```matlab
    gl_params.SNR_in_dB = 20;          % Thermal Noise floor
    gl_params.SIR_in_dB = -20;         % Strong Jammer scenario
    gl_params.jammerType = 'Spot';     % Interference type
    gl_params.use_quantization = true; % Enable hardware simulation
    gl_params.bits_phase = 8;          % Phase shifter resolution
    gl_params.bits_gain = 6;           % Attenuator resolution
    ```

3.  **Run:**
    Execute `main.m`. The script will output:
    * **Figure 1:** Per-frame Output SNR.
    * **Figure 1 (Subplots):** Time-domain comparison of input vs. output signals (Desired, Jammer, Noise).

## Simulation Insights

* **Quantization Effect:** When `use_quantization` is enabled, the null depth is limited by the bit resolution. A 6-bit gain quantization typically limits the maximum achievable SIR improvement to ~35-40 dB.
* **Spot Jammer:** The Spot Jammer operates with a duty cycle. Note that `main.m` normalizes power over the *entire* duration, meaning the jammer's peak power during the "On" time is significantly higher than the average SIR indicates, providing a stress test for the beamformer.
* **Spatial Smoothing:** This is crucial when the jammer is coherent or when the covariance matrix is ill-conditioned. It effectively creates a "sub-array" averaging effect.

## Author
[Yuval Ginzberg - Technion]