%% DRONE COMMUNICATION SIMULATION – MAIN SCRIPT (UPDATED WITH gl_params)

clear; clc; close all;
addpath ('C:\Users\ginyu\OneDrive - Technion\Desktop\טכניון\תואר שני\Thesis\PMO_project\Drone_comm\*')
%% ========================= 1) DEFINE GLOBAL PARAMETERS ===================

gl_params = struct();

% Physical constants
gl_params.c        = 3e8;
gl_params.f0       = 1.5e9;                          % Carrier [Hz]
gl_params.lambda   = gl_params.c / gl_params.f0;

% Array configuration
gl_params.numElements = 4;
gl_params.d           = 0.5 * gl_params.lambda;      % Element spacing

% Sampling and timing
gl_params.fs      = 1e9;                           % Dev sampling rate
gl_params.Tsim    = 10e-3;                            % Total sim time [s]
gl_params.t       = (0:1/gl_params.fs:gl_params.Tsim - 1/gl_params.fs);
gl_params.N       = numel(gl_params.t);

% Frame parameters
gl_params.frameDur   = 100e-6;                       % Frame duration [s]
gl_params.frameLen   = round(gl_params.frameDur * gl_params.fs);
gl_params.numFrames  = floor(gl_params.N / gl_params.frameLen);
gl_params.lambda_mem = 0.95;                         % Memory factor

% Signal band
gl_params.bw_sig  = 20e6;                            % Desired bandwidth
gl_params.bw_tx   = gl_params.bw_sig;                % Jammer reference span

% SNR and SIR
gl_params.SNR_in_dB = 10;                            % Desired vs noise
gl_params.SIR_in_dB = -20;                           % Desired vs jammer

% DOAs (uniform random in [-60,60])
gl_params.theta_desired_deg = 30;
gl_params.theta_jammer_deg  = -20;
gl_params.num_signals = 2;

% MUSIC scanning grid
gl_params.scanAngles = -90:0.5:90;

% Jammer type selection
gl_params.jammerType = 'Barrage';                   % 'CW','Barrage','Spot','Sweep','MultiTone'

% Jammer parameters
switch lower(gl_params.jammerType)
    case 'cw'
        gl_params.jamParams = struct('fOffsetHz', 5e6);

    case 'barrage'
        gl_params.jamParams = struct('bwHz', 100e6);

    case 'spot'
        gl_params.jamParams = struct( ...
            'bwHz', 20e6, ...
            'onTimeSec', 500e-6, ...
            'startTimeSec', 300e-6 );

    case 'sweep'
        gl_params.jamParams = struct( ...
            'numTones', 20, ...
            'dwellTimeSec', 50e-6, ...
            'fSpanHz', 40e6 );

    case 'multitone'
        gl_params.jamParams = struct( ...
            'numTones', 5, ...
            'fSpanHz', 40e6 );

    otherwise
        error('Unsupported jammer type.');
end

%% ========================= 2) ULA + STEERING VECTORS ======================

ula = phased.ULA('NumElements', gl_params.numElements, ...
                 'ElementSpacing', gl_params.d);

theta_des_rad = deg2rad(gl_params.theta_desired_deg);
theta_jam_rad = deg2rad(gl_params.theta_jammer_deg);

elemIdx = (0:gl_params.numElements-1).';     % [4×1]

a_desired = exp(1j * 2*pi * gl_params.d / gl_params.lambda * elemIdx * sin(theta_des_rad));
a_jammer  = exp(1j * 2*pi * gl_params.d / gl_params.lambda * elemIdx * sin(theta_jam_rad));

%% ========================= 3) DESIRED SIGNAL GENERATION ===================

fs     = gl_params.fs;
t      = gl_params.t;
bw_sig = gl_params.bw_sig;

% White Gaussian complex baseband
[s_bb, ofdm_params, preamble_td] = generate_ofdm_signal(gl_params);

% Band-limit using bandfilter()
s0_filt = bandfilter(s_bb.', 1, fs, -bw_sig/2, bw_sig/2);
s_bb = s0_filt.';                    % row vector

% Normalize desired signal to unit RMS
s_bb = s_bb / sqrt(mean(abs(s_bb).^2));

P_sig = 1;                           % Defined unit power
s_bb = s_bb * sqrt(P_sig);

%% ========================= 4) GENERATE UNIT-POWER JAMMER ==================

jammer_unit = generate_jammer(gl_params.jammerType, ...
                              gl_params.t, ...
                              gl_params.fs, ...
                              gl_params.f0, ...
                              gl_params.bw_tx, ...
                              gl_params.jamParams);

P_jam_unit = mean(abs(jammer_unit).^2);

%% ========================= 5) SIR-BASED JAMMER SCALING ====================

sir_lin = 10^(gl_params.SIR_in_dB/10);
P_jam_target = P_sig / sir_lin;

alpha_jam = sqrt(P_jam_target / P_jam_unit);
jammer_scaled = alpha_jam * jammer_unit;   % Correct SIR at input

%% ========================= 6) ARRAY SIGNAL CONSTRUCTION ===================

% Desired mapped to array
S_array = a_desired * s_bb;         % [M x N]

% Jammer mapped to array
J_array = a_jammer * jammer_scaled; % [M x N]

% Noise power from SNR
noisePower = P_sig / (10^(gl_params.SNR_in_dB/10));
noise = sqrt(noisePower/2) * (randn(gl_params.numElements, gl_params.N) + ...
                              1j*randn(gl_params.numElements, gl_params.N));

% Total signal at array
X_total = S_array + J_array + noise;

%% ========================= 7) INITIALIZE COVARIANCE =======================

R = 1e-6 * eye(gl_params.numElements);

%% ========================= 8) FRAME LOOP ==================================

y_out      = zeros(1, gl_params.N);
SNR_out_dB = zeros(1, gl_params.numFrames);

for k = 1:gl_params.numFrames
    idxStart = (k-1)*gl_params.frameLen + 1;
    idxEnd   = k * gl_params.frameLen;
    idx      = idxStart:idxEnd;

    Xk_total = X_total(:, idx);          % [M x L]
    Xk_sig   = S_array(:, idx);
    Xk_jam   = J_array(:, idx);
    Xk_noise = noise(:, idx);
    Lk       = size(Xk_total, 2);

    % ---------- Instantaneous covariance ----------
    R_inst = (Xk_total * Xk_total') / Lk;

    % ---------- Exponential memory ----------
    R = gl_params.lambda_mem * R + (1 - gl_params.lambda_mem) * R_inst;

    % ---------- DOA estimation (placeholder for MUSIC) ----------
    estTheta_des = gl_params.theta_desired_deg;
    estTheta_jam = gl_params.theta_jammer_deg;

    % ---------- Beamforming weights (placeholder for pc_bf) ----------
    % Replace with your actual function:
    w = pc_beamformer(R, gl_params.num_signals, gl_params.numElements, gl_params.theta_desired_deg);
    % w = a_desired / norm(a_desired);    % Temporary MRC

    w_q = w;   % If quantizing later: apply fi() here

    % ---------- Apply beamformer ----------
    yk_total = w_q' * Xk_total;
    yk_sig   = w_q' * Xk_sig;
    yk_jam   = w_q' * Xk_jam;
    yk_noise = w_q' * Xk_noise;

    % Store output
    y_out(idx) = yk_total;

    % ---------- Compute output SNR ----------
    P_sig_out   = mean(abs(yk_sig).^2);
    P_noise_out = mean(abs(yk_jam + yk_noise).^2);
    SNR_out_dB(k) = 10*log10(P_sig_out / P_noise_out);
end

%% ========================= 9) PLOTS ======================================

figure;
subplot(2,1,1);
plot(1:gl_params.numFrames, SNR_out_dB, '-o');
xlabel('Frame index');
ylabel('Output SNR [dB]');
title('Per-frame output SNR');
grid on;

subplot(2,1,2);
plot(gl_params.t*1e3, real(y_out));
xlabel('Time [ms]');
ylabel('Real\{y(t)\}');
title('Beamformer Output – Real Part');
grid on;
