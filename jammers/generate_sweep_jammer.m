function jammer = generate_sweep_jammer(t, fs, bw_tx, params)
% params fields:
%   .numTones     : number of frequencies in sweep, e.g. 20
%   .dwellTimeSec : dwell time per frequency [s], e.g. 1e-3
%   .fSpanHz      : span to sweep over (default: bw_tx)
%
% Output jammer has unit RMS over the non-zero part.

if ~isfield(params, 'numTones'), params.numTones = 20; end
if ~isfield(params, 'dwellTimeSec'), params.dwellTimeSec = 1e-3; end
if ~isfield(params, 'fSpanHz'), params.fSpanHz = bw_tx; end

jammer = zeros(size(t));                        % row vector

fVec = linspace(-params.fSpanHz/2, params.fSpanHz/2, params.numTones);

for k = 1:params.numTones
    tStart = (k-1)*params.dwellTimeSec;
    tEnd   = k*params.dwellTimeSec;
    idx = (t >= tStart) & (t < tEnd);

    if ~any(idx)
        continue;
    end

    fk = fVec(k);
    jammer(idx) = exp(1j*2*pi*fk*t(idx));
end

% Normalize RMS over non-zero samples
idxNZ = abs(jammer) > 0;
if any(idxNZ)
    rms_val = sqrt(mean(abs(jammer(idxNZ)).^2));
    if rms_val > 0
        jammer(idxNZ) = jammer(idxNZ) / rms_val;
    end
end

end
