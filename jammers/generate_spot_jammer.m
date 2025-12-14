function jammer = generate_spot_jammer(t, fs, params)
% params fields:
%   .bwHz        : spot bandwidth [Hz], e.g. 10e6 or 20e6
%   .onTimeSec   : total on-time [s], e.g. 500e-6 to 5e-3
%   .startTimeSec: start time of burst [s] within t-vector
%
% Output jammer has unit RMS over its active interval (where it is non-zero).

if ~isfield(params, 'bwHz'), params.bwHz = 10e6; end
if ~isfield(params, 'onTimeSec'), params.onTimeSec = 1e-3; end
if ~isfield(params, 'startTimeSec'), params.startTimeSec = 0; end

N = numel(t);

% 1) Generate complex white Gaussian noise
w0 = (randn(N,1) + 1j*randn(N,1)) / sqrt(2);   % [N x 1]

% 2) Band-limit to the specified spot bandwidth
bw = params.bwHz;
lo = -bw/2;
hi =  bw/2;
w_filt = bandfilter(w0, 1, fs, lo, hi);        % [N x 1]

% 3) Time-gate according to on-time window
tStart = params.startTimeSec;
tEnd   = tStart + params.onTimeSec;
idxOn  = (t.' >= tStart) & (t.' < tEnd);        % logical [N x 1]

jammer_vec = zeros(N,1);
jammer_vec(idxOn) = w_filt(idxOn);

% 4) Normalize RMS only over the active samples
if any(idxOn)
    rms_val = sqrt(mean(abs(jammer_vec(idxOn)).^2));
    if rms_val > 0
        jammer_vec(idxOn) = jammer_vec(idxOn) / rms_val;
    end
end

jammer = jammer_vec.';                          % row vector

end
