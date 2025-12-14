function jammer = generate_barrage_jammer(t, fs, bw_tx, params)
% params fields:
%   .bwHz : barrage bandwidth [Hz], e.g. 100e6 or 200e6
%
% Output jammer has unit RMS in-band.

if ~isfield(params, 'bwHz'), params.bwHz = bw_tx; end

N = numel(t);

% 1) Generate complex white Gaussian noise CN(0,1)
w0 = (randn(N,1) + 1j*randn(N,1)) / sqrt(2);   % [N x 1]

% 2) Band-limit to the desired bandwidth around baseband
bw = params.bwHz;
lo = -bw/2;
hi =  bw/2;
w_filt = bandfilter(w0, 1, fs, lo, hi);        % [N x 1]

% 3) Normalize to unit RMS
rms_val = sqrt(mean(abs(w_filt).^2));
if rms_val > 0
    w_norm = w_filt / rms_val;
else
    w_norm = w_filt;
end

jammer = w_norm.';                              % row vector

end
