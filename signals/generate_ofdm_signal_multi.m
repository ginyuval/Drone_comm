function [s_bb, ofdm_params, preamble_td] = generate_ofdm_signal_multi(gl_params)
%GENERATE_OFDM_SIGNAL Generate baseband OFDM waveform with multiple packets and a cyclic QPSK preamble.
%
% Packet structure:
%   [preamble (Np OFDM symbols)] + [data (Nd OFDM symbols)] + [guard zeros (Ng samples)]
%
% Preamble design (NEW):
%   - Create a short cyclic sequence of QPSK symbols defined by 2-bit words
%     (random per generation unless user provides gl_params.preamblePatternBits)
%   - Repeat this sequence until filling all preamble resource elements
%     (Nused_eff * Np), then reshape to [Nused_eff x Np]
%   - DC subcarrier is explicitly NOT used
%
% Required gl_params fields:
%   gl_params.fs, gl_params.t, gl_params.N, gl_params.bw_sig
%
% Optional gl_params fields:
%   gl_params.numPackets            (default 3)
%   gl_params.dataSymbolsPerPacket  (default 10)
%   gl_params.guardIntervalSec      (default 50e-6)
%   gl_params.preambleDurSec        (default 0.3e-3) OR gl_params.preambleDur
%   gl_params.Nfft                  (default 1024)
%   gl_params.cpFrac                (default 1/8)
%   gl_params.randSeed              (optional; affects DATA + preamble pattern randomness)
%   gl_params.preamblePatternLen    (default 5)
%   gl_params.preamblePatternBits   (optional explicit [L x 2] bits, rows are [bI bQ])
%
% Outputs:
%   s_bb        : [1 x N] complex waveform, unit RMS over ACTIVE (non-zero) samples
%   ofdm_params : struct with OFDM and packetization parameters
%   preamble_td : [1 x Lp] time-domain preamble of ONE packet (for correlation)

% ---------------------- Basic parameters ----------------------------------
fs     = gl_params.fs;
N      = gl_params.N;
bw_sig = gl_params.bw_sig;

% Optional RNG seed (affects DATA + preamble pattern randomness)
if isfield(gl_params, 'randSeed')
    rng(gl_params.randSeed);
end

% ---------------------- Packetization parameters --------------------------
if isfield(gl_params, 'numPackets'),           numPackets = gl_params.numPackets; else, numPackets = 3; end
if isfield(gl_params, 'dataSymbolsPerPacket'), Nd = gl_params.dataSymbolsPerPacket; else, Nd = 10; end
if isfield(gl_params, 'guardIntervalSec'),     guardSec = gl_params.guardIntervalSec; else, guardSec = 50e-6; end

if isfield(gl_params, 'preambleDurSec')
    preambleSec = gl_params.preambleDurSec;
elseif isfield(gl_params, 'preambleDur')
    preambleSec = gl_params.preambleDur;
else
    preambleSec = 0.3e-3;
end

% ---------------------- OFDM parameters -----------------------------------
if isfield(gl_params, 'Nfft'),   Nfft = gl_params.Nfft; else, Nfft = 1024; end
if isfield(gl_params, 'cpFrac'), cpFrac = gl_params.cpFrac; else, cpFrac = 1/8; end
cpLen = round(cpFrac * Nfft);

deltaF = fs / Nfft;

% Requested number of used carriers to approximate bw_sig
Nused_req = floor(bw_sig / deltaF);
Nused_req = min(Nused_req, Nfft - 2);
if mod(Nused_req,2)==1, Nused_req = Nused_req - 1; end

% Active carrier indices in DC-centered indexing, EXCLUDING DC bin
usedIdx_dc = (-Nused_req/2 : Nused_req/2-1);
usedIdx_dc(usedIdx_dc == 0) = [];         % ---- remove DC explicitly ----
Nused_eff = numel(usedIdx_dc);            % effective number of active subcarriers

NsymSamples = Nfft + cpLen;
Tsym        = NsymSamples / fs;

% Preamble symbols count
Np = max(1, round(preambleSec / Tsym));

% Guard samples
Ng = round(guardSec * fs);

% ---------------------- Build CYCLIC QPSK preamble (LOWER PAPR) ------------
% Key change: do NOT repeat a short pattern across subcarriers.
% Instead create one QPSK symbol per active subcarrier (length Nused_eff),
% then cyclic-shift this full vector across preamble OFDM symbols.

if isfield(gl_params, 'preambleShiftPerSym')
    shiftPerSym = gl_params.preambleShiftPerSym;
else
    shiftPerSym = 1;
end

% One QPSK value per active carrier (random per generation)
bitsI0 = randi([0 1], Nused_eff, 1);
bitsQ0 = randi([0 1], Nused_eff, 1);
baseVec = ((2*bitsI0 - 1) + 1j*(2*bitsQ0 - 1)) / sqrt(2);   % [Nused_eff x 1]

% Store base preamble bits/symbols (one per active subcarrier)
preambleBitsI0 = bitsI0;
preambleBitsQ0 = bitsQ0;
preambleBaseVec = baseVec;


% Cyclic across symbols (frequency-domain circular shift)
preambleSym = zeros(Nused_eff, Np);
for m = 1:Np
    sh = mod(shiftPerSym*(m-1), Nused_eff);
    preambleSym(:,m) = circshift(baseVec, sh);
end


% ---------------------- Packet/sample sizes --------------------------------
packetSamples = (Np + Nd) * NsymSamples;             % preamble+data (no guard)
totalSamples  = numPackets * packetSamples + (numPackets-1) * Ng;

s_full = zeros(1, totalSamples);

% ---------------------- Helper: build one OFDM symbol ----------------------
    function x_cp = ofdm_symbol_from_subcarriers(symVecUsed)
        % symVecUsed: [Nused_eff x 1] values for active subcarriers in DC-centered indexing (no DC)
        X_dc = zeros(Nfft, 1);
        centerIdx = Nfft/2 + 1;
        X_dc(centerIdx + usedIdx_dc) = symVecUsed;

        X_ifft  = ifftshift(X_dc);
        x_no_cp = ifft(X_ifft);
        x_cp    = [x_no_cp(end-cpLen+1:end); x_no_cp];  % [NsymSamples x 1]
    end

% ---------------------- Construct packets ----------------------------------
writePos = 1;
for p = 1:numPackets

    % Data symbols for this packet (random QPSK on active carriers)
    bitsI_d = randi([0 1], Nused_eff, Nd);
    bitsQ_d = randi([0 1], Nused_eff, Nd);
    dataSym = ((2*bitsI_d - 1) + 1j*(2*bitsQ_d - 1)) / sqrt(2);

    % Concatenate [preamble | data] in frequency domain
    SymFD = [preambleSym, dataSym];  % [Nused_eff x (Np+Nd)]

    % Convert packet to time domain
    pkt_td = zeros(1, packetSamples);
    for m = 1:(Np+Nd)
        x_cp = ofdm_symbol_from_subcarriers(SymFD(:,m));
        idxSym = (m-1)*NsymSamples + 1 : m*NsymSamples;
        pkt_td(idxSym) = x_cp.';
    end

    % Write packet into waveform buffer
    s_full(writePos : writePos + packetSamples - 1) = pkt_td;
    writePos = writePos + packetSamples;

    % Insert guard zeros between packets (not after last)
    if p < numPackets
        writePos = writePos + Ng; % zeros already in s_full
    end
end

% ---------------------- Build time-domain preamble explicitly --------------
preambleLenSamples = Np * NsymSamples;
preamble_td_full = zeros(1, preambleLenSamples);
for m = 1:Np
    x_cp = ofdm_symbol_from_subcarriers(preambleSym(:,m));
    idxSym = (m-1)*NsymSamples + 1 : m*NsymSamples;
    preamble_td_full(idxSym) = x_cp.';
end

% Ensure zero time-domain mean (removes any residual DC offset)
preamble_td_full = preamble_td_full - mean(preamble_td_full);

% ---------------------- Normalize (active samples only) --------------------
activeMask = abs(s_full) > 0;                 % excludes guard zeros
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
    s_bb = [s_full, zeros(1, N - totalSamples)];
end

preamble_td = preamble_td_full;

% ---------------------- Fill ofdm_params ----------------------------------
ofdm_params = struct();
ofdm_params.fs                    = fs;
ofdm_params.Nfft                  = Nfft;
ofdm_params.cpLen                 = cpLen;
ofdm_params.cpFrac                = cpFrac;
ofdm_params.deltaF                = deltaF;

ofdm_params.Nused_requested        = Nused_req;
ofdm_params.Nused_effective        = Nused_eff;
ofdm_params.usedIdx_dc             = usedIdx_dc;        % DC removed
ofdm_params.NsymSamples            = NsymSamples;
ofdm_params.Tsym                   = Tsym;

ofdm_params.numPackets             = numPackets;
ofdm_params.numPreambleSymbols     = Np;
ofdm_params.dataSymbolsPerPacket   = Nd;
ofdm_params.guardIntervalSec       = guardSec;
ofdm_params.guardSamples           = Ng;

ofdm_params.packetSamples          = packetSamples;
ofdm_params.totalSamplesGenerated  = totalSamples;

ofdm_params.preambleLenSamples     = preambleLenSamples;
ofdm_params.preambleDurSec         = Np * Tsym;

ofdm_params.preambleType      = 'CyclicFullVectorQPSK';
ofdm_params.preambleShiftPerSym = shiftPerSym;

ofdm_params.preambleBitsI0    = preambleBitsI0;    % [Nused_eff x 1]
ofdm_params.preambleBitsQ0    = preambleBitsQ0;    % [Nused_eff x 1]
ofdm_params.preambleBaseVec   = preambleBaseVec;   % [Nused_eff x 1]


end
