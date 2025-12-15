function [w, R_r, R_ss] = pc_beamformer_ss(R, npc, numelements, angle, L)

% ---------------- PC reduction (same as your original) ----------------
M = numelements-L+1;
v_m = exp(1j*pi*((0:M-1)')*sin(deg2rad(angle)));
R_ss = spsmooth(R, L);
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
