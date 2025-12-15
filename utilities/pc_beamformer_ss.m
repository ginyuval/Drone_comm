function [w, R_r, R_ss] = pc_beamformer_ss(R, npc, numelements, angle, ss_params)
%% PC Beamformer with Spatial Smoothing
%
% Inputs:
%   R          : [M x M] sample covariance (ULA)
%   npc        : number of principal components kept (signal+interferers subspace)
%   numelements: M (number of array elements)
%   angle      : look direction (deg)
%   ss_params  : struct with fields (optional)
%       .enable        (default true)
%       .L             subarray length (<= M). Typical: L = M-1 or M-2
%       .mode          'forward' (default) or 'fb' (forward-backward)
%
% Outputs:
%   w    : [M x 1] beamforming weights
%   R_r  : [M x M] reduced covariance (projected onto signal+interferer subspace)
%   R_ss : [M x M] spatially-smoothed covariance actually used

% ---------------- Defaults ----------------
if nargin < 5 || isempty(ss_params)
    ss_params = struct();
end
if ~isfield(ss_params, 'enable'), ss_params.enable = true; end
if ~isfield(ss_params, 'mode'),   ss_params.mode   = 'forward'; end
if ~isfield(ss_params, 'L')
    % Default subarray length: M-1 (common choice for 4 elements -> L=3)
    ss_params.L = max(2, numelements - 1);
end

M = numelements;
L = ss_params.L;

% ---------------- Spatial smoothing ----------------
if ~ss_params.enable
    R_ss = R;
else
    if L > M
        error('Spatial smoothing: L must satisfy L <= numelements.');
    end
    K = M - L + 1;  % number of overlapping subarrays

    % Forward smoothing on covariance: average principal submatrices
    Rf = zeros(M, M);
    for k = 1:K
        idx = k:(k+L-1);
        Rk  = R(idx, idx);     % [L x L]
        Rf(idx, idx) = Rf(idx, idx) + Rk;
    end
    Rf = Rf / K;

    % Forward-backward option
    if strcmpi(ss_params.mode, 'fb')
        J = fliplr(eye(M));
        Rb = J * conj(Rf) * J;
        R_ss = 0.5 * (Rf + Rb);
    else
        R_ss = Rf;
    end

    % Enforce Hermitian numerically
    R_ss = 0.5 * (R_ss + R_ss');
end

% ---------------- PC reduction (same as your original) ----------------
v_m = exp(1j*pi*((0:M-1)')*sin(deg2rad(angle)));

[V, D] = eig(R_ss, 'vector');
[D, ind] = sort(D, "descend");
V = V(:, ind);

% Keep npc components
npc = min(npc, M);
SSI = diag(D(1:npc));          % [npc x npc]
U_SI = V(:, 1:npc);

% Reduced covariance (signal+interferer subspace)
R_r = U_SI * SSI * U_SI';

% ---------------- PC beamformer weights ----------------
% Equivalent stable form: w = A*v / (v'*A*v) with A = U*inv(SSI)*U'
A = U_SI * (SSI \ eye(npc)) * U_SI';          % U * inv(SSI) * U'
den = (v_m' * A * v_m);
w = (A * v_m) / (den + 1e-12);               % [M x 1]

end
