%% DRONE COMMUNICATION SIMULATION – MAIN SCRIPT (UPDATED WITH gl_params)

clear; clc; close all;
addpath(genpath('C:\Users\ginyu\OneDrive - Technion\Desktop\טכניון\תואר שני\Thesis\PMO_project\Drone_comm\'))
%% ========================= 1) DEFINE GLOBAL PARAMETERS ===================

gl_params = struct();

% Physical constants
gl_params.c        = 3e8;
gl_params.f0       = 1.5e9;
gl_params.lambda   = gl_params.c / gl_params.f0;

% Array configuration
gl_params.numElements = 4;
gl_params.d           = 0.5 * gl_params.lambda;

% Sampling and timing
gl_params.fs      = 60e6;                 % Sampling frequency [Hz]
gl_params.Tsim    = 1e-3;                % Total simulation time [s]
gl_params.t       = 0 : 1/gl_params.fs : gl_params.Tsim - 1/gl_params.fs;
gl_params.N       = numel(gl_params.t);

% ---------------- FRAME-BASED PROCESSING (FIXED) ----------------
gl_params.frameDur  = 10e-6;                                % Frame duration [s]
gl_params.frameLen  = round(gl_params.frameDur * gl_params.fs);

% Number of frames needed to cover ALL samples
gl_params.numFrames = ceil(gl_params.N / gl_params.frameLen);

% Safety: number of processed samples (always >= N)
gl_params.Nproc     = gl_params.numFrames * gl_params.frameLen;
% ---------------------------------------------------------------

% Covariance memory factor
gl_params.lambda_mem = 0.0;

% Signal bandwidth
gl_params.bw_sig = 20e6;
gl_params.bw_tx  = gl_params.bw_sig;

% SNR / SIR
gl_params.SNR_in_dB = 15;
gl_params.SIR_in_dB = -20;

% DOAs
gl_params.theta_desired_deg = -30;
gl_params.theta_jammer_deg  = 19;
gl_params.num_signals = 2;

% Scan grid
gl_params.scanAngles = -90:0.5:90;


% Jammer type selection
gl_params.jammerType = 'Spot';                   % 'CW','Barrage','Spot','Sweep','MultiTone'

% Jammer parameters
switch lower(gl_params.jammerType)
    case 'cw'
        gl_params.jamParams = struct('fOffsetHz', 5e6);

    case 'barrage'
        gl_params.jamParams = struct('bwHz', 100e6);

    case 'spot'
        gl_params.jamParams = struct( ...
            'bwHz', 20e6, ...
            'onTimeSec', 600e-6, ...
            'startTimeSec', 100e-6 );

    case 'sweep'
        gl_params.jamParams = struct( ...
            'numTones', 20, ...
            'dwellTimeSec', 100e-6, ...
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

mvdrEstimator = phased.MVDREstimator('SensorArray', ula, ...
                                      'OperatingFrequency', gl_params.f0, ...
                                      'NumSignals', gl_params.num_signals, ...
                                      'DOAOutputPort',true, 'ForwardBackwardAveraging',true);

theta_des_rad = deg2rad(gl_params.theta_desired_deg);
theta_jam_rad = deg2rad(gl_params.theta_jammer_deg);

elemIdx = (0:gl_params.numElements-1).';     % [4×1]

a_desired = exp(1j * 2*pi * gl_params.d / gl_params.lambda * elemIdx * sin(theta_des_rad));
a_jammer  = exp(1j * 2*pi * gl_params.d / gl_params.lambda * elemIdx * sin(theta_jam_rad));

%% ========================= 3) DESIRED SIGNAL GENERATION ===================

fs     = gl_params.fs;
t      = gl_params.t;
bw_sig = gl_params.bw_sig;
gl_params.numPackets = 20;
gl_params.dataSymbolsPerPacket = 10;
gl_params.guardIntervalSec = 100e-6;
gl_params.guardIntervalN = round(100e-6*gl_params.fs);
gl_params.preambleDurSec = 0.3e-3;
gl_params.randSeed = rand; 
gl_params.NFFT = 128; 
gl_params.preambleShiftPerSym = 10;


% White Gaussian complex baseband
[s_bb, ofdm_params, preamble_td] = generate_ofdm_signal_multi(gl_params);

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

R = eye(4, 4);

%% ========================= 8) FRAME LOOP ==================================

y_out       = zeros(1, gl_params.N);
y_sig_out   = y_out;
y_jam_out   = y_out;
y_noise_out = y_out;
SNR_out_dB  = zeros(1, gl_params.numFrames);
FirstTimeNPC2 = 0;
FirstTimeNPC1 = 0;
% Correlation configuration
gl_params.preamble = preamble_td(:).';           % row vector
gl_params.Lp = numel(gl_params.preamble);
gl_params.L = 2;
% Buffer to handle preamble across frame boundaries
bufferLenSamples = 2*ofdm_params.packetSamples + gl_params.guardIntervalN;
corrBuf = zeros(gl_params.numElements, bufferLenSamples);         

hist1 = NaN(1,5);
hist2 = NaN(1,5);

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
    R = R_inst;

    
    % ---------- DOA estimation using phased.MVDREstimator ----------


    [Y, Estimated_angs] = mvdrEstimator(Xk_total.');   % Estimated_angs contains 2 angles (deg)
    [theta_trk, hist1, hist2] = fix_angle_indexing_5(Estimated_angs, hist1, hist2, 15);
    Estimated_angs = theta_trk;   % now index is stable
    % Estimated_angs
    % Select which estimated DOA corresponds to the desired signal using preamble correlation
    [estTheta_des, corrMetric, corrBuf] = select_desired_doa_by_preamble( ...
        (Xk_total), Estimated_angs, gl_params, corrBuf);
     % plot(abs(corrBuf))
     % drawnow
    estTheta_des_v(k) = estTheta_des;
    % estTheta_des = 30;
    % The other DOA is treated as the interferer (if two distinct angles exist)
    estTheta_jam = pick_other_angle(Estimated_angs, estTheta_des);

    % ---------- Exponential memory ----------
    npc = mdltest(R_inst);
    if npc < 2
        estTheta_des = max(Estimated_angs);
    end

    if npc==2 && FirstTimeNPC2==0
        R_2 = R_inst;
        FirstTimeNPC2 = 1;
    end
    if npc==1 && FirstTimeNPC1==0
        R_1 = R_inst;
        FirstTimeNPC1 = 1;
    end

  
    if npc == 2
        R_2 = gl_params.lambda_mem * R_2 + (1 - gl_params.lambda_mem) * R_inst;
        R = R_2;
    elseif npc == 1
        R_1 = gl_params.lambda_mem * R_1 + (1 - gl_params.lambda_mem) * R_inst;
        R = R_1;
    else
        R = R_inst;
    end

    if ~npc
        npc = 1;
    end
    % ---------- Beamforming weights (placeholder for pc_bf) ----------
    
    % w = pc_beamformer(R, npc, gl_params.numElements, estTheta_des);
    w = pc_beamformer_ss(R, npc, gl_params.numElements, estTheta_des, gl_params.L);
    w = (w)/norm(w);
    

    % ula_pat = phased.ULA('NumElements', gl_params.numElements-gl_params.L+1, ...
    %              'ElementSpacing', gl_params.d);
    % pattern(ula_pat,gl_params.f0 ,gl_params.scanAngles,0,...
    %     'Weights',w,'CoordinateSystem','rectangular',...
    %     'Type','directivity');
    % hold on
    % xline(gl_params.theta_desired_deg, 'Color','green')
    % xline(gl_params.theta_jammer_deg, 'Color','red')
    % ylim([-50, 10]);
    % xlim([-90, 90]);
    % drawnow
    % hold off
    w_q = w;   % If quantizing later: apply fi() here

    % ---------- Apply beamformer ----------
    yk_total = w_q' * Xk_total(1:gl_params.numElements-gl_params.L+1, :);
    yk_sig   = w_q' * Xk_sig(1:gl_params.numElements-gl_params.L+1, :);
    yk_jam   = w_q' * Xk_jam(1:gl_params.numElements-gl_params.L+1, :);
    yk_noise = w_q' * Xk_noise(1:gl_params.numElements-gl_params.L+1, :);

    % Store output
    y_out(idx) = yk_total;
    y_sig_out(idx) = yk_sig;
    y_jam_out(idx) = yk_jam;
    y_noise_out(idx) = yk_noise;


    % ---------- Compute output SNR ----------
    P_sig_out   = mean(abs(yk_sig).^2);
    P_noise_out = mean(abs(yk_jam + yk_noise).^2);
    SNR_out_dB(k) = 10*log10(P_sig_out / P_noise_out);
end

%% ========================= 9) PLOTS ======================================

figure;

% ===================== Subplot 1: Output SNR =====================
subplot(3,1,1);
plot(1:gl_params.numFrames, SNR_out_dB, '-o');
xlabel('Frame index');
ylabel('Output SNR [dB]');
title('Per-frame output SNR');
ylim([-10, 26]);
grid on;
legend('Beamformer output SNR','Location','best');

% ===================== Subplot 2: Beamformer Output Components =====================
subplot(3,1,2);
hold on;
plot(gl_params.t*1e3, db(medfilt1(abs(y_sig_out)/norm(y_sig_out), 200)));
plot(gl_params.t*1e3, db(medfilt1(abs(y_jam_out)/norm(y_sig_out), 200)));
plot(gl_params.t*1e3, db(medfilt1(abs(y_noise_out)/norm(y_sig_out), 200)));
hold off;

xlabel('Time [ms]');
ylabel('Normalized magnitude');
title('Beamformer output components');
grid on;
legend({'Desired signal','Interference','Noise'}, 'Location','best');

% ===================== Subplot 3: Beamformer Input Components =====================
subplot(3,1,3);
hold on;
plot(gl_params.t*1e3, db(medfilt1(abs(s_bb)/norm(s_bb), 200)));
plot(gl_params.t*1e3, db(medfilt1(abs(jammer_scaled)/norm(s_bb), 200)));
plot(gl_params.t*1e3, db(medfilt1(abs(noise(1,:))/norm(s_bb), 200)));
hold off;

xlabel('Time [ms]');
ylabel('Normalized magnitude');
title('Beamformer input components (element 1)');
grid on;
legend({'Desired signal','Interference','Noise'}, 'Location','best');


