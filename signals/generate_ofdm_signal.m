function [s_bb, ofdm_params, preamble_td] = generate_ofdm_signal(gl_params)
%GENERATE_OFDM_SIGNAL Generate baseband OFDM signal with leading preamble.
%
%   [s_bb, ofdm_params, preamble_td] = GENERATE_OFDM_SIGNAL(gl_params)
%
%   Inputs (from gl_params):
%       gl_params.fs       : sampling frequency [Hz]
%       gl_params.t        : time vector [1 x N]
%       gl_params.N        : number of samples
%       gl_params.bw_sig   : desired OFDM occupied bandwidth [Hz]
%       gl_params.preambleDur (optional): preamble duration [s],
%                                        default = 0.3e-3 (0.3 ms)
%
%   Outputs:
%       s_bb        : complex baseband OFDM signal [1 x N], unit RMS power
%       ofdm_params : struct with OFDM configuration
%       preamble_td : time-domain preamble samples [1 x Lp], where
%                    Lp = ofdm_params.preambleLenSamples (after normalization)
%
%   Notes:
%       1) OFDM uses QPSK symbols on active subcarriers.
%       2) Preamble is implemented as the first Ns_pream OFDM symbols.
%       3) The first preambleLenSamples of s_bb are exactly preamble_td.
%

% Extract basic parameters
fs     = gl_params.fs;
t      = gl_params.t;
N      = gl_params.N;
bw_sig = gl_params.bw_sig;

% Preamble duration
if isfield(gl_params, 'preambleDur')
    preambleDur = gl_params.preambleDur;
else
    preambleDur = 0.3e-3;  % default: 0.3 ms
end

% ---------------------- OFDM parameter definition -------------------------

% IFFT size
Nfft   = 1024;
cpFrac = 1/8;                        % CP length as fraction of Nfft
cpLen  = round(cpFrac * Nfft);       % CP in samples

% Subcarrier spacing
deltaF = fs / Nfft;                  % [Hz]

% Number of used subcarriers to approximate desired bandwidth
Nused  = floor(bw_sig / deltaF);
Nused  = min(Nused, Nfft - 2);       % margin away from edges
if mod(Nused, 2) == 1
    Nused = Nused - 1;               % enforce even number
end

% DC-centered indices of active subcarriers: [-Nused/2 ... Nused/2-1]
usedIdx_dc = (-Nused/2 : Nused/2-1);     % length = Nused

% Samples per OFDM symbol including CP
NsymSamples = Nfft + cpLen;

% Symbol duration
Tsym = NsymSamples / fs;             % [s]

% Number of preamble symbols to approximate preambleDur
numPreambleSymbols = max(1, round(preambleDur / Tsym));

% Number of symbols needed to cover total duration
numSymbols = ceil(N / NsymSamples);
numPreambleSymbols = min(numPreambleSymbols, numSymbols);  % safety

% ---------------------- Generate QPSK data --------------------------------
% Generate QPSK for all active subcarriers and OFDM symbols
bitsI   = randi([0 1], Nused, numSymbols);
bitsQ   = randi([0 1], Nused, numSymbols);
dataSym = ((2*bitsI - 1) + 1j*(2*bitsQ - 1)) / sqrt(2);  % QPSK, unit power

% For now, the preamble is simply the first numPreambleSymbols OFDM symbols
% using the first columns of dataSym. If later you want a fixed known
% pattern, you can overwrite dataSym(:,1:numPreambleSymbols) with a
% deterministic sequence.
preambleSym = dataSym(:, 1:numPreambleSymbols);

% ---------------------- Time-domain OFDM construction ---------------------

s_bb_full = zeros(1, numSymbols * NsymSamples);

for m = 1:numSymbols
    % Frequency-domain vector in DC-centered indexing
    X_dc = zeros(Nfft, 1);
    
    % DC is at index Nfft/2 + 1 in DC-centered indexing
    centerIdx = Nfft/2 + 1;
    X_dc(centerIdx + usedIdx_dc) = dataSym(:, m);
    
    % Convert from DC-centered order to standard IFFT order
    X_ifft = ifftshift(X_dc);
    
    % IFFT to obtain OFDM symbol without CP
    x_no_cp = ifft(X_ifft);                     % [Nfft x 1]
    
    % Add cyclic prefix
    x_cp = [x_no_cp(end - cpLen + 1:end); x_no_cp];  % [NsymSamples x 1]
    
    % Insert into global buffer
    idxSym = (m-1)*NsymSamples + 1 : m*NsymSamples;
    s_bb_full(idxSym) = x_cp.';                 % row insertion
end

% ---------------------- Extract preamble in time domain -------------------

preambleLenSamples = numPreambleSymbols * NsymSamples;
preamble_td_full   = s_bb_full(1:preambleLenSamples);   % before normalization

% ---------------------- Trim and normalize --------------------------------

% Trim to exact simulation length
s_bb = s_bb_full(1:N);                          % [1 x N]

% Normalize to unit RMS power
P_sig_est = mean(abs(s_bb).^2);
if P_sig_est > 0
    scale = 1 / sqrt(P_sig_est);
    s_bb = s_bb * scale;
    preamble_td_full = preamble_td_full * scale;
else
    scale = 1;
end

% Final preamble_td (normalized); trim in case total N is shorter
Lp = min(preambleLenSamples, N);
preamble_td = preamble_td_full(1:Lp);

% ---------------------- Fill ofdm_params struct ---------------------------

ofdm_params = struct();
ofdm_params.Nfft                = Nfft;
ofdm_params.cpLen               = cpLen;
ofdm_params.Nused               = Nused;
ofdm_params.usedIdx_dc          = usedIdx_dc;
ofdm_params.deltaF              = deltaF;
ofdm_params.NsymSamples         = NsymSamples;
ofdm_params.numSymbols          = numSymbols;
ofdm_params.numPreambleSymbols  = numPreambleSymbols;
ofdm_params.preambleLenSamples  = preambleLenSamples;
ofdm_params.preambleDur         = numPreambleSymbols * Tsym;

end
