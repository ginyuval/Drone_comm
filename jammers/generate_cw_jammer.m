function jammer = generate_cw_jammer(t, params)
% params fields:
%   .fOffsetHz : frequency offset from carrier [Hz]
%
% Output jammer has approximately unit RMS.

if ~isfield(params, 'fOffsetHz'), params.fOffsetHz = 0; end

fOff = params.fOffsetHz;

j = exp(1j * 2*pi * fOff * t);

rms_val = sqrt(mean(abs(j).^2));
if rms_val > 0
    jammer = j / rms_val;
else
    jammer = j;
end

end
