%% Test: Frame Duration vs. Output SNR
clear; clc; close all;
addpath(genpath('../')); % Ensure access to parent folders

%% 1. Parameters Setup
jammerTypeToTest = 'barrage'; 

% Define range of Frame durations to test [seconds]
% e.g., 10us to 200us. (At 30MHz, 5us = 150 samples)
frameDurs = [5e-6, 10e-6, 20e-6, 50e-6, 100e-6]; 
results_oSNR = zeros(size(frameDurs));

fprintf('Starting Frame Duration Test with %s Jammer...\n', jammerTypeToTest);

%% 2. Loop
for i = 1:length(frameDurs)
    % Initialize params exactly as in main.m
    gl_params = get_main_params(jammerTypeToTest);
    
    % --- Override Frame Duration ---
    gl_params.frameDur = frameDurs(i);
    
    % --- Recalculate dependent frame parameters ---
    gl_params.frameLen  = round(gl_params.frameDur * gl_params.fs);
    gl_params.numFrames = ceil(gl_params.N / gl_params.frameLen);
    gl_params.Nproc     = gl_params.numFrames * gl_params.frameLen;
    
    % Run Simulation
    avgSNR = run_simulation_core(gl_params);
    results_oSNR(i) = avgSNR;
    
    fprintf('Run %d/%d: FrameDur=%.1f us -> oSNR=%.2f dB\n', ...
        i, length(frameDurs), frameDurs(i)*1e6, avgSNR);
end

%% 3. Plot
figure;
semilogx(frameDurs*1e6, results_oSNR, '-o', 'LineWidth', 2);
xlabel('Frame Duration [\mus]');
ylabel('Average Output SNR [dB]');
title(['Frame Duration vs. Output SNR (' jammerTypeToTest ')']);
grid on;

%% --- Helper Functions ---
function gl_params = get_main_params(jType)
    gl_params = struct();
    
    % --- 1) DEFINE GLOBAL PARAMETERS ---
    % Physical constants
    gl_params.c        = 3e8;
    gl_params.f0       = 1.5e9;
    gl_params.lambda   = gl_params.c / gl_params.f0;
    
    % Array configuration
    gl_params.numElements = 4;
    gl_params.d           = 0.5 * gl_params.lambda;
    
    % Sampling and timing
    gl_params.fs      = 30e6;                 
    gl_params.Tsim    = 1e-3; % Kept short for testing speed
    gl_params.t       = 0 : 1/gl_params.fs : gl_params.Tsim - 1/gl_params.fs;
    gl_params.N       = numel(gl_params.t);
    
    % --- Quantization Parameters ---
    gl_params.use_quantization = true; 
    gl_params.bits_phase       = 8;    
    gl_params.bits_gain        = 6;
    gl_params.sqnr_phase = 6.02*gl_params.bits_phase + 1.76;
    gl_params.sqnr_gain = 6.02*gl_params.bits_gain + 1.76;
    
    % --- FRAME-BASED PROCESSING (Default values, overwritten in loop) ---
    gl_params.frameDur  = 5e-6;                                
    gl_params.frameLen  = round(gl_params.frameDur * gl_params.fs);
    gl_params.numFrames = ceil(gl_params.N / gl_params.frameLen);
    gl_params.Nproc     = gl_params.numFrames * gl_params.frameLen;
    
    % Covariance memory factor
    gl_params.lambda_mem = 0.0;
    
    % Signal bandwidth
    gl_params.bw_sig = 20e6;
    gl_params.bw_tx  = gl_params.bw_sig;
    
    % SNR / SIR
    gl_params.SNR_in_dB = 20;
    gl_params.SIR_in_dB = -20;
    
    % DOAs
    gl_params.theta_desired_deg = 30;
    gl_params.theta_jammer_deg  = -10;
    gl_params.num_signals = 2;
    gl_params.scanAngles = -90:0.5:90;
    
    % Jammer type selection
    gl_params.jammerType = jType;
    
    % Jammer parameters
    switch lower(gl_params.jammerType)
        case 'cw'
            gl_params.jamParams = struct('fOffsetHz', 5e6);
        case 'barrage'
            gl_params.jamParams = struct('bwHz', 100e6);
        case 'spot'
            gl_params.jamParams = struct('bwHz', 20e6, 'onTimeSec', 600e-6, 'startTimeSec', 250e-6);
        case 'sweep'
            gl_params.jamParams = struct('numTones', 20, 'dwellTimeSec', 100e-6, 'fSpanHz', 40e6);
        case 'multitone'
            gl_params.jamParams = struct('numTones', 5, 'fSpanHz', 40e6);
    end
    
    % --- 3) DESIRED SIGNAL GENERATION PARAMETERS ---
    gl_params.numPackets = 20;
    gl_params.dataSymbolsPerPacket = 10;
    gl_params.guardIntervalSec = 0e-6;
    gl_params.guardIntervalN = round(gl_params.guardIntervalSec*gl_params.fs);
    gl_params.preambleDurSec = 0.3e-3;
    gl_params.randSeed = rand; 
    gl_params.NFFT = 128; 
    gl_params.preambleShiftPerSym = 10;
end

function avgSNR = run_simulation_core(gl_params)
    % Initialize objects (Section 2)
    ula = phased.ULA('NumElements', gl_params.numElements, 'ElementSpacing', gl_params.d);
    mvdrEstimator = phased.MVDREstimator('SensorArray', ula, 'OperatingFrequency', gl_params.f0, ...
        'NumSignals', gl_params.num_signals, 'DOAOutputPort', true, 'ForwardBackwardAveraging', true);
    
    % Steering vectors
    theta_des_rad = deg2rad(gl_params.theta_desired_deg);
    theta_jam_rad = deg2rad(gl_params.theta_jammer_deg);
    elemIdx = (0:gl_params.numElements-1).';
    a_desired = exp(1j * 2*pi * gl_params.d / gl_params.lambda * elemIdx * sin(theta_des_rad));
    a_jammer  = exp(1j * 2*pi * gl_params.d / gl_params.lambda * elemIdx * sin(theta_jam_rad));
    
    % Generate Signals (Section 3)
    [s_bb, ofdm_params, preamble_td] = generate_ofdm_signal_multi(gl_params);
    s_bb = s_bb / sqrt(mean(abs(s_bb).^2)); % Normalize unit RMS
    P_sig = 1; 
    s_bb = s_bb * sqrt(P_sig);
    
    % Generate Jammer (Section 4)
    jammer_unit = generate_jammer(gl_params.jammerType, gl_params.t, gl_params.fs, gl_params.f0, gl_params.bw_tx, gl_params.jamParams);
    P_jam_unit = mean(abs(jammer_unit).^2);
    
    % SIR Scaling (Section 5)
    sir_lin = 10^(gl_params.SIR_in_dB/10);
    P_jam_target = P_sig / sir_lin;
    alpha_jam = sqrt(P_jam_target / (P_jam_unit + 1e-12));
    jammer_scaled = alpha_jam * jammer_unit;
    
    % Array Signal Construction (Section 6)
    S_array = a_desired * s_bb;
    J_array = a_jammer * jammer_scaled;
    noisePower = P_sig / (10^(gl_params.SNR_in_dB/10));
    noise = sqrt(noisePower/2)*(randn(gl_params.numElements, gl_params.N) + 1j*randn(gl_params.numElements, gl_params.N));
    X_total = S_array + J_array + noise;
    
    % Frame Loop Setup (Section 8)
    gl_params.preamble = preamble_td(:).';
    gl_params.Lp = numel(gl_params.preamble);
    gl_params.L = 2;
    
    % Buffer calculation
    bufferLenSamples = 1*ofdm_params.packetSamples + gl_params.guardIntervalN;
    corrBuf = zeros(gl_params.numElements, bufferLenSamples);
    
    hist1 = NaN(1,5); hist2 = NaN(1,5);
    snr_vec = [];
    
    R = eye(4);
    
    for k = 1:gl_params.numFrames
        idxStart = (k-1)*gl_params.frameLen + 1;
        idxEnd = k*gl_params.frameLen;
        
        % Check bounds
        if idxEnd > size(X_total, 2)
            break; 
        end
        
        Xk = X_total(:, idxStart:idxEnd);
        Xk_sig = S_array(:, idxStart:idxEnd);
        Xk_jam = J_array(:, idxStart:idxEnd);
        Xk_noise = noise(:, idxStart:idxEnd);
        Lk = size(Xk, 2);
        
        R_inst = (Xk*Xk')/Lk;
        % Apply memory factor (simplification used in user script)
        R = R*gl_params.lambda_mem + (1-gl_params.lambda_mem)*R_inst; 
        
        % MVDR
        [~, EstAngs] = mvdrEstimator(Xk.');
        [theta_trk, hist1, hist2] = fix_angle_indexing_5(EstAngs, hist1, hist2, 15);
        
        % Preamble Selection
        [estTheta_des, ~, corrBuf] = select_desired_doa_by_preamble(Xk, theta_trk, gl_params, corrBuf);
        
        npc = mdltest(R_inst);
        if npc < 2, estTheta_des = max(theta_trk); end
        if ~npc, npc=1; end
        
        % Beamforming
        w = pc_beamformer_ss(R_inst, npc, gl_params.numElements, estTheta_des, gl_params.L);
        w = w/norm(w);
        
        % Quantization
        if gl_params.use_quantization
            w = quantize_weight_vector(w, gl_params.bits_phase, gl_params.bits_gain);
        end
        
        % Apply (Spatial Smoothing dimension adjustment)
        M_final = gl_params.numElements - gl_params.L + 1;
        yk_sig = w' * Xk_sig(1:M_final, :);
        yk_interf = w' * (Xk_jam(1:M_final, :) + Xk_noise(1:M_final, :));
        
        P_s = mean(abs(yk_sig).^2);
        P_i = mean(abs(yk_interf).^2);
        
        % Avoid Log of zero
        if P_i < 1e-15, P_i = 1e-15; end
        snr_vec(end+1) = 10*log10(P_s/P_i);
    end
    
    % Average SNR (linear mean of dB values or mean of power ratios)
    % The user's script used pow2db(mean(db2pow(...))).
    avgSNR = pow2db(mean(db2pow(snr_vec(1:end)))); 
end