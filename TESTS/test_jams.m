%% Test script for jammer generation and spectrum verification

clear; clc; close all;

%% Global simulation parameters

fs    = 1.5e9;          % Sampling frequency [Hz] (reduced for testing)
Tsim  = 3e-3;           % Total simulation time [s]
t     = (0:1/fs:Tsim-1/fs);   % Time vector [1 x N]
N     = numel(t);
f0    = 1.5e9;          % Carrier frequency [Hz] (not directly used here)
bw_tx = 20e6;           % Nominal TX bandwidth [Hz]

%% Helper function for FFT plotting (nested function)
plot_dB_spectrum = @(sig, fs, titleStr) ...
    local_plot_fft(sig, fs, titleStr);

%% 1) CW jammer test

params_cw = struct();
params_cw.fOffsetHz = 5e6;   % 5 MHz offset from baseband

jammer_cw = generate_jammer('CW', t, fs, f0, bw_tx, params_cw);

figure;
plot_dB_spectrum(jammer_cw, fs, 'CW jammer (5 MHz offset)');

%% 2) Barrage jammer test

params_barrage = struct();
params_barrage.bwHz = bw_tx;  % 50 MHz barrage bandwidth

jammer_barrage = generate_jammer('Barrage', t, fs, f0, bw_tx, params_barrage);

figure;
plot_dB_spectrum(jammer_barrage, fs, 'Barrage jammer (50 MHz BW)');

%% 3) Spot jammer test

params_spot = struct();
params_spot.bwHz        = 20e6;   % 20 MHz spot bandwidth
params_spot.onTimeSec   = 200e-6; % 200 us burst
params_spot.startTimeSec = 300e-6; % start at 300 us

jammer_spot = generate_jammer('Spot', t, fs, f0, bw_tx, params_spot);

figure;
subplot(2,1,1);
plot(t*1e3, real(jammer_spot));
xlabel('Time [ms]');
ylabel('Real part');
title('Spot jammer (time domain)');
grid on;

subplot(2,1,2);
plot_dB_spectrum(jammer_spot, fs, 'Spot jammer (20 MHz BW, burst)');

%% 4) Sweep jammer test

params_sweep = struct();
params_sweep.numTones     = 10;      % 10 frequencies
params_sweep.dwellTimeSec = 1e-3;   % 1ms per tone
params_sweep.fSpanHz      = 1e9;    % sweep over +/-500 MHz

jammer_sweep = generate_jammer('Sweep', t, fs, f0, bw_tx, params_sweep);

figure;
plot_dB_spectrum(jammer_sweep, fs, 'Sweep jammer (20 tones)');

%% 5) MultiTone jammer test

params_mt = struct();
params_mt.numTones     = 5;        % 5 tones
params_mt.fSpanHz      = 20e6;     % tones across +/-20 MHz

jammer_multitone = generate_jammer('MultiTone', t, fs, f0, bw_tx, params_mt);

figure;
plot_dB_spectrum(jammer_multitone, fs, 'MultiTone jammer (5 tones)');

%% Local helper function for FFT plotting
function local_plot_fft(sig, fs, titleStr)
    % Compute magnitude spectrum in dB and plot vs frequency in MHz

    sig = sig(:).';                    % ensure row
    N   = numel(sig);

    SigF = fftshift(fft(sig));
    mag  = abs(SigF);

    % Avoid log of zero
    mag_dB = 20*log10(mag + 1e-12);
    mag_dB = mag_dB - max(mag_dB);     % normalize to 0 dB peak

    f_axis = (-N/2:N/2-1)/N * fs;      % frequency axis [Hz]

    plot(f_axis/1e6, mag_dB, 'LineWidth', 1.5);
    xlabel('Frequency [MHz]');
    ylabel('Magnitude [dB]');
    title(titleStr);
    grid on;
    xlim([min(f_axis) max(f_axis)]/1e6);
end
