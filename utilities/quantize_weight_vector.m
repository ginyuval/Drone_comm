function w_q = quantize_weight_vector(w, bits_phase, bits_gain)
% QUANTIZE_WEIGHT_VECTOR Simulates hardware quantization for beamforming weights.
%
%   w_q = quantize_weight_vector(w, bits_phase, bits_gain)
%
%   Inputs:
%       w          : Complex weight vector [M x 1] (Double precision)
%       bits_phase : Number of bits for the Phase Shifter (e.g., 6 for 64 states)
%       bits_gain  : Number of bits for the Gain Shifter/Attenuator (e.g., 5)
%
%   Output:
%       w_q        : Quantized complex weight vector (Double precision approximation)
%
%   The function uses MATLAB's Fixed-Point Designer 'fi' object to simulate
%   limited resolution in both magnitude and phase.

    % 1. Extract magnitude and phase
    mag_w = abs(w);
    phase_w = angle(w); % Returns values between -pi and pi

    % ----------------------- PHASE QUANTIZATION -----------------------
    % Phase Shifters typically operate in a 0 to 360 degree range.
    % We normalize the phase to the [0, 1) range representing a full circle.
    
    norm_phase = (phase_w + pi) / (2*pi); 
    
    % Create Fixed-Point object (Unsigned)
    % FractionLength = WordLength ensures the value is a fraction in [0, 1).
    % This simulates dividing the circle into 2^bits_phase equal sectors.
    F_phase = fi(norm_phase, 0, bits_phase, bits_phase);
    
    % Convert back to Double and rescale to radians [-pi, pi)
    quant_phase = double(F_phase) * 2*pi - pi;

    % ----------------------- GAIN QUANTIZATION ------------------------
    % Simulating Digital Attenuators.
    % Assume the maximum possible gain is the maximum value in the current vector.
    max_val = max(mag_w);
    % max_val = 70;
    if max_val > 0
        % Normalize values to [0, 1) relative to the maximum
        norm_mag = mag_w / max_val;
        
        % Create Fixed-Point object for gain (Unsigned)
        F_gain = fi(norm_mag, 0, bits_gain, bits_gain);
        
        % Reconstruct the quantized magnitude
        quant_mag = double(F_gain) * max_val;
    else
        quant_mag = mag_w; % If the vector is all zeros
    end

    % ----------------------- RECONSTRUCTION ---------------------------
    w_q = quant_mag .* exp(1j * quant_phase);

end