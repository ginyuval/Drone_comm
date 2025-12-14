function theta_other = pick_other_angle(candAnglesDeg, theta_chosen)
    %PICK_OTHER_ANGLE Return the candidate angle that is not theta_chosen.
    candAnglesDeg = candAnglesDeg(:).';
    candAnglesDeg = candAnglesDeg(~isnan(candAnglesDeg));
    
    if numel(candAnglesDeg) < 2
        theta_other = NaN;
        return;
    end
    
    % Choose the angle farthest from theta_chosen
    [~, idx] = max(abs(candAnglesDeg - theta_chosen));
    theta_other = candAnglesDeg(idx);
    
    % If that accidentally picks the same (rare), fall back
    if abs(theta_other - theta_chosen) < 1e-6
        theta_other = candAnglesDeg(3-idx);
    end
end
