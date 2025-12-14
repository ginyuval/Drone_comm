function jammer = generate_multitone_jammer(t, fs, bw_tx, params)
    % params fields:
    %   .numTones : number of tones (default 5)
    %   .fSpanHz  : frequency span (default: bw_tx)
    %
    % Multi-tone jammer = sum of several simultaneous (not sequential) tones.
    % Output is normalized to have unit RMS.
    
    if ~isfield(params, 'numTones'), params.numTones = 5; end
    if ~isfield(params, 'fSpanHz'), params.fSpanHz = bw_tx; end
    
    Ntones = params.numTones;
    span   = params.fSpanHz;
    
    % Tone frequencies in baseband
    fVec = linspace(-span/2, span/2, Ntones);  % [Hz]
    
    % Generate multi-tone
    j = zeros(size(t));
    for k = 1:Ntones
        j = j + exp(1j * 2*pi * fVec(k) * t);
    end
    
    % Normalize to unit RMS
    rms_val = sqrt(mean(abs(j).^2));
    if rms_val > 0
        jammer = j / rms_val;
    else
        jammer = j;
    end

end
