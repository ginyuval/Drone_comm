function [theta_trk, hist1, hist2] = fix_angle_indexing_5(Estimated_angs, hist1, hist2, gate_deg)
% Keeps consistent indexing of 2 DOA estimates using last 5 values.
% - hist1, hist2 are 1x5 buffers (NaN allowed)
% - gate_deg: max jump allowed before declaring "missing" (e.g., 10-20 deg)

if nargin < 4 || isempty(gate_deg), gate_deg = 15; end

angs = Estimated_angs(:).';                     % 1x2
angs(~isfinite(angs)) = NaN;                    % remove inf/nan

% Previous representative angles (median over last 5)
p1 = median(hist1, 'omitnan');
p2 = median(hist2, 'omitnan');

% If first call: initialize tracks directly
if ~isfinite(p1) && ~isfinite(p2)
    hist1 = push5(hist1, angs(1));
    hist2 = push5(hist2, angs(2));
    theta_trk = [median(hist1,'omitnan'), median(hist2,'omitnan')];
    return;
end

% Compute distances to previous tracks
d11 = abs(angs(1) - p1);  d12 = abs(angs(1) - p2);
d21 = abs(angs(2) - p1);  d22 = abs(angs(2) - p2);

% Assign by minimal total distance (2x2 assignment)
if (d11 + d22) <= (d21 + d12)
    a1 = angs(1); a2 = angs(2);
else
    a1 = angs(2); a2 = angs(1);
end

% Gate: if update is too far (or NaN), treat as missing -> hold (median)
if ~isfinite(a1) || abs(a1 - p1) > gate_deg, a1 = p1; end
if ~isfinite(a2) || abs(a2 - p2) > gate_deg, a2 = p2; end

% Update histories
hist1 = push5(hist1, a1);
hist2 = push5(hist2, a2);

theta_trk = [median(hist1,'omitnan'), median(hist2,'omitnan')];

end

function h = push5(h, x)
h = [h(2:end), x];
end
