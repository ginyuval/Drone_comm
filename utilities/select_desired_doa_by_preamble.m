function [theta_des_hat, metric, corrBufOut] = select_desired_doa_by_preamble(Xk, candAnglesDeg, gl_params, corrBufIn)
    % SELECT_DESIRED_DOA_BY_PREAMBLE
    % Now re-processes historical Raw Data to enable a fair comparison between angles.
    candAnglesDeg = candAnglesDeg(:).';
    candAnglesDeg = candAnglesDeg(~isnan(candAnglesDeg));
    
    % If there are no candidates, return the buffer as is (updated with new Xk at the end)
    if isempty(candAnglesDeg)
        theta_des_hat = NaN;
        metric = -Inf;
        % Update FIFO buffer even when no decision is made
        bufLen = size(corrBufIn, 2);
        X_combined = [corrBufIn, Xk];
        if size(X_combined, 2) > bufLen
            corrBufOut = X_combined(:, end-bufLen+1:end);
        else
            corrBufOut = X_combined;
        end
        return;
    end
    
    p = gl_params.preamble; 
    % Create Matched Filter
    h = conj(fliplr(p));
    
    metrics = -Inf(1, numel(candAnglesDeg));
    
    % Concatenate raw data: history + current frame
    % X_proc will be of size [M x (BufferLen + CurrentFrameLen)]
    X_proc = [corrBufIn, Xk];
    
    for i = 1:numel(candAnglesDeg)
        th = candAnglesDeg(i);
    
        % Create steering vector for the candidate
        a = steering_vec_ula(th, gl_params);
        w = a / norm(a);
    
        % Re-apply Beamforming to all history and present data
        % Now, the Preamble located in the Buffer "feels" the new weights w
        y_proc = w' * X_proc; 
    
        % Perform correlation on the processed signal
        r = conv(y_proc, h, 'full');

        % if i == 1
        %     plot(abs(r), '-o', 'Color', 'red')
        %     hold on;
        % else
        %     plot(abs(r), '-*', 'Color', 'blue')
        % end
        

        % The metric is the peak magnitude
        metrics(i) = max(abs(r));
    end
    % hold off
    % drawnow
    % Select the angle with the highest correlation
    [metric, idxBest] = max(metrics);
    theta_des_hat = candAnglesDeg(idxBest);
    
    % Update buffer for the next frame - storing RAW DATA
    bufLen = size(corrBufIn, 2);
    
    % Take the last samples from the combined matrix
    if size(X_proc, 2) >= bufLen
        corrBufOut = X_proc(:, end-bufLen+1:end);
    else
        % Pad with zeros on the left if the buffer is not yet full (at the beginning)
        padLen = bufLen - size(X_proc, 2);
        corrBufOut = [zeros(gl_params.numElements, padLen), X_proc];
    end
end