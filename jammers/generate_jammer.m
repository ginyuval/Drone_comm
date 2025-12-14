function jammer = generate_jammer(jammerType, t, fs, f0, bw_tx, params)
%GENERATE_JAMMER Generate unit-power jammer signal for drone comm simulation.
%
%   jammer = GENERATE_JAMMER(jammerType, t, fs, f0, bw_tx, params)
%
%   Inputs:
%       jammerType : string, type of jammer:
%                    'CW', 'Barrage', 'Spot', 'Sweep', 'MultiTone'
%       t          : time vector [1 x N] or [N x 1] in seconds
%       fs         : sampling frequency [Hz]
%       f0         : carrier (central) frequency [Hz] of comm system
%                    (kept for completeness; not directly used here)
%       bw_tx      : TX signal bandwidth [Hz] (used as default span)
%       params     : struct with fields depending on jammerType.
%
%   Output:
%       jammer     : complex baseband jammer signal [1 x N]
%
%   Notes:
%       1) The signal is generated in complex baseband (relative to f0).
%       2) The returned jammer is normalized so that its RMS power in the
%          active interval (where it is non-zero) is approximately 1.
%       3) SIR-based power scaling is done outside this function in the
%          main simulation code.
%

t = t(:).';              % ensure row vector
jammer = zeros(size(t)); % complex baseband

switch lower(jammerType)
    case 'cw'
        jammer = generate_cw_jammer(t, params);

    case 'barrage'
        jammer = generate_barrage_jammer(t, fs, bw_tx, params);

    case 'spot'
        jammer = generate_spot_jammer(t, fs, params);

    case 'sweep'
        jammer = generate_sweep_jammer(t, fs, bw_tx, params);

    case 'multitone'
        jammer = generate_multitone_jammer(t, fs, bw_tx, params);

    otherwise
        error('Unknown jammer type: %s', jammerType);
end

end
