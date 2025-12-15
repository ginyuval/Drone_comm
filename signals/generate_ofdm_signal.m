function [s_bb, ofdm_params, preamble_td] = generate_ofdm_signal(gl_params)
%GENERATE_OFDM_SIGNAL Generate baseband OFDM waveform with multiple packets.
%
% Each packet structure:
%   [preamble (Np symbols)] + [data (Nd symbols)] + [guard zeros (Ng samples)]
%
% Required fields (gl_params):
%   fs, t, N, bw_sig
%
% Optional fields (gl_params):
%   preambleDurSec        : preamble duration [s], default 0.3e-3
%   numPackets            : number of packets, default 3
%   dataSymbolsPerPacket  : number of OFDM data symbols per packet, default 10
%   guardIntervalSec      : zero-guard between packets [s], default 50e-6
%   Nfft                  : IFFT size, default 1024
%   cpFrac                : CP fraction, default 1/8
%   randSeed              : set RNG seed for repeatability (optional)
%
% Outputs:
%   s_bb        : [1 x N] complex waveform, normalized to unit RMS over ACTIVE part
%   ofdm_params : struct (packet layout + OFDM config)
%   preamble_td : [1 x Lp] time-domain preamble (one packet), normalized consistently

% ---------------------- Basic parameters ----------------------------------
fs     = gl_params.fs;
N      = gl_params.N;
bw_sig = gl_params.bw_sig;

% Optional RNG seed
if isfield(gl_params, 'randSeed')
    rng(gl_params.randSeed);
end

% Packetization parameters
if isfield(gl_params, 'numPackets'),           numPackets = gl_params.numPackets; else, numPackets = 3; end
if isfield(gl_params, 'dataSymbolsPerPacket'), Nd = gl_params.dataSymbolsPerPacket; else, Nd = 10; end
if isfield(gl_params, 'guardIntervalSec'),     guardSec = gl_params.guardIntervalSec; else, guardSec = 50e-6; end
if isfield(gl_params, 'preambleDurSec'),       preambleSec = gl_params.preambleDurSec; ...
elseif isfield(gl_params, 'preambleDur'),      preambleSec = gl_params.preambleDur; ...
else,                                          preambleSec = 0.3e-3;
end

% OFDM parameters
if isfield(gl_params, 'Nfft'),   Nfft = gl_params.Nfft; else, Nfft = 1024; end
if isfield(gl_params, 'cpFrac'), cpFrac = gl_params.cpFrac; else, cpFrac = 1/8; end
cpLen = round(cpFrac * Nfft);

deltaF = fs / Nfft;

% Number of used subcarriers
Nused = floor(bw_sig / deltaF);
Nused = min(Nused, Nfft - 2);
if mod(Nused,2)==1, Nused = Nused - 1; end
usedIdx_dc = (-Nused/2 : Nused/2-1);

NsymSamples = Nfft + cpLen;
Tsym        = NsymSamples / fs;

% Preamble symbols count
Np = max(1, round(preambleSec / Tsym));   % preamble symbols per packet

% Guard samples
Ng = round(guardSec * fs);

% ---------------------- Build deterministic preamble (KNOWN) --------------
% Use a fixed QPSK pattern for preamble so correlation works across runs.
% (If you want a specific standard later, replace this.)
bitsI_p = randi([0 1], Nused, Np);
bitsQ_p = randi([0 1], Nused, Np);
preambleSym = ((2*bitsI_p - 1) + 1j*(2*bitsQ_p - 1)) / sqrt(2); % [Nused x Np]

% ---------------------- Allocate full waveform buffer ----------------------
packetSamples = (Np + Nd) * NsymSamples;      % samples in preamble+data (no guard)
totalSamples  = numPackets * packetSamples + (numPackets-1) * Ng;

s_full = zeros(1, totalSamples);

% ---------------------- Helper to build one OFDM symbol --------------------
    function x_cp = ofdm_symbol_from_subcarriers(symVecUsed)
        % symVecUsed: [Nused x 1] values for active subcarriers in DC-centered indexing
        X_dc = zeros(Nfft, 1);
        centerIdx = Nfft/2 + 1;
        X_dc(centerIdx + usedIdx_dc) = symVecUsed;

        X_ifft = ifftshift(X_dc);
        x_no_cp = ifft(X_ifft);
        x_cp = [x_no_cp(end-cpLen+1:end); x_no_cp];  % [NsymSamples x 1]
    end

% ---------------------- Construct packets ----------------------------------
writePos = 1;
for p = 1:numPackets

    % --- Data symbols for this packet (random QPSK) ---
    bitsI_d = randi([0 1], Nused, Nd);
    bitsQ_d = randi([0 1], Nused, Nd);
    dataSym = ((2*bitsI_d - 1) + 1j*(2*bitsQ_d - 1)) / sqrt(2); % [Nused x Nd]

    % --- Concatenate [preamble | data] in frequency domain ---
    SymFD = [preambleSym, dataSym];  % [Nused x (Np+Nd)]

    % --- Convert to time domain symbol-by-symbol ---
    pkt_td = zeros(1, packetSamples);
    for m = 1:(Np+Nd)
        x_cp = ofdm_symbol_from_subcarriers(SymFD(:,m));
        idxSym = (m-1)*NsymSamples + 1 : m*NsymSamples;
        pkt_td(idxSym) = x_cp.';
    end

    % --- Write packet into global buffer ---
    s_full(writePos : writePos + packetSamples - 1) = pkt_td;
    writePos = writePos + packetSamples;

    % --- Insert guard zeros between packets (not after last) ---
    if p < numPackets
        writePos = writePos + Ng; % zeros already in s_full
    end
end

% ---------------------- Extract time-domain preamble (one packet) ----------
preambleLenSamples = Np * NsymSamples;

% Build preamble time-domain explicitly from preambleSym (more robust)
preamble_td_full = zeros(1, preambleLenSamples);
for m = 1:Np
    x_cp = ofdm_symbol_from_subcarriers(preambleSym(:,m));
    idxSym = (m-1)*NsymSamples + 1 : m*NsymSamples;
    preamble_td_full(idxSym) = x_cp.';
end

% ---------------------- Normalize -----------------------------------------
% Normalize using ACTIVE samples only (exclude guard zeros) to keep SIR meaningful.
activeMask = abs(s_full) > 0; % guard zeros excluded
P_est = mean(abs(s_full(activeMask)).^2);

if ~isempty(P_est) && P_est > 0
    scale = 1 / sqrt(P_est);
    s_full = s_full * scale;
    preamble_td_full = preamble_td_full * scale;
end

% ---------------------- Trim/pad to requested N ----------------------------
if totalSamples >= N
    s_bb = s_full(1:N);
else
    s_bb = [s_full, zeros(1, N-totalSamples)];
end

% Final preamble (normalized)
preamble_td = preamble_td_full;

% ---------------------- Fill ofdm_params ----------------------------------
ofdm_params = struct();
ofdm_params.Nfft = Nfft;
ofdm_params.cpLen = cpLen;
ofdm_params.cpFrac = cpFrac;
ofdm_params.deltaF = deltaF;
ofdm_params.Nused = Nused;
ofdm_params.usedIdx_dc = usedIdx_dc;
ofdm_params.NsymSamples = NsymSamples;
ofdm_params.Tsym = Tsym;

ofdm_params.numPackets = numPackets;
ofdm_params.numPreambleSymbols = Np;
ofdm_params.dataSymbolsPerPacket = Nd;
ofdm_params.guardIntervalSec = guardSec;
ofdm_params.guardSamples = Ng;

ofdm_params.packetSamples = packetSamples;
ofdm_params.totalSamplesGenerated = totalSamples;
ofdm_params.preambleLenSamples = preambleLenSamples;
ofdm_params.preambleDurSec = Np * Tsym;

end
