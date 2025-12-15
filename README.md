# Drone Communication Anti-Jamming Simulation with Spatial Smoothing

## Overview
This MATLAB project simulates an uplink communication scenario between a drone (or ground station) and a receiver equipped with a Uniform Linear Array (ULA). The primary goal of the simulation is to demonstrate robust **Anti-Jamming (AJ)** capabilities using beamforming techniques.

The system is designed to distinguish a desired OFDM signal from various types of interference (jammers) using **DOA (Direction of Arrival) estimation**, **Preamble Correlation**, and **Spatial Smoothing**.

## Key Features

* **OFDM Signal Generation:** Generates QPSK-modulated OFDM signals with cyclic preambles for synchronization and identification.
* **Interference Modeling:** Includes a versatile jammer generator supporting multiple modes:
    * `CW` (Continuous Wave)
    * `Barrage` (Wideband Noise)
    * `Spot` (Narrowband Burst)
    * `Sweep` (Chirped/Hopping)
    * `MultiTone`
* **Array Processing:**
    * **ULA Configuration:** Simulation of a 4-element array (configurable).
    * **Spatial Smoothing:** Implemented to decorrelate coherent signals or improve covariance matrix rank estimation.
    * **MVDR Estimation:** Uses Minimum Variance Distortionless Response for initial DOA detection.
    * **PC Beamforming:** Principal Component Beamforming with spatial smoothing (`pc_beamformer_ss`) to null nullify interference.
* **Smart Source Selection:**
    * Algorithm to automatically identify the desired user vs. the jammer.
    * Uses **Preamble Correlation** (Matched Filter) on raw buffered data to associate the correct DOA with the desired signal, even when the jammer is stronger (negative SIR).
    * Includes a **Frame-Based Buffer** mechanism to handle preambles crossing frame boundaries.

## File Structure

### Root Directory
* `main.m`: The entry point of the simulation. Configures parameters (`gl_params`), runs the frame loop, performs beamforming, and plots SNR/Signal results.

### /signals
* `generate_ofdm_signal_multi.m`: Generates the desired baseband OFDM waveform with packet structure and guard intervals.

### /jammers
* `generate_jammer.m`: Wrapper function to generate interference.
* `generate_barrage_jammer.m`, `generate_spot_jammer.m`, etc.: Specific implementations for different jamming profiles.

### /utilities
* `select_desired_doa_by_preamble.m`: Core logic for distinguishing the desired signal from the jammer using correlation metrics and historical raw data buffering.
* `pc_beamformer_ss.m`: Implementation of the PC Beamformer with Spatial Smoothing.
* `steering_vec_ula.m`: Computes steering vectors for the ULA.
* `fix_angle_indexing_5.m`: Tracks and stabilizes DOA estimates across frames to prevent index swapping.

## Installation & Usage

1.  **Prerequisites:**
    * MATLAB (R2021b or later recommended).
    * Phased Array System Toolbox (required for `phased.ULA` and `phased.MVDREstimator`).

2.  **Running the Simulation:**
    * Open `main.m`.
    * Adjust the global parameters (`gl_params`) section to test different scenarios:
        ```matlab
        gl_params.SNR_in_dB = 15;        % Signal-to-Noise Ratio
        gl_params.SIR_in_dB = -20;       % Signal-to-Interference Ratio (Negative = Strong Jammer)
        gl_params.theta_desired_deg = -30;
        gl_params.theta_jammer_deg  = 19;
        gl_params.jammerType = 'spot';   % Choose: 'cw', 'barrage', 'spot', etc.
        ```
    * Run the script.

3.  **Outputs:**
    * **SNR Plot:** Shows the output SINR per processed frame.
    * **Beamformer Output:** Time-domain visualization of the recovered signal vs. residual interference and noise.

## Algorithm Logic
1.  **Frame Processing:** The continuous signal is sliced into short processing frames.
2.  **Covariance Estimation:** The sample covariance matrix is computed for each frame.
3.  **DOA Estimation:** MVDR estimates candidate angles of arrival.
4.  **Target Identification:** The `select_desired_doa_by_preamble` function steers a beam towards each candidate angle using historical raw data (Pre-Beamforming Buffer) and correlates the output with the known preamble sequence. The angle yielding the highest correlation peak is identified as the **Desired Signal**.
5.  **Beamforming:** A Principal Component beamformer (with Spatial Smoothing) is calculated to maximize gain towards the desired angle while suppressing the jammer.

