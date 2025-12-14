function [signal] = bandfilter(signal,num_channels, fs, lo, hi)
    N = length(signal(:,1));
    nbands = length(lo);
    freq_axis = linspace(-fs/2, fs/2-fs/N, N);
    rect_filter = zeros(size(freq_axis));
    for i=1:nbands
        rect_filter = rect_filter + ((freq_axis) >= (lo(i)) & (freq_axis) <= (hi(i)));
    end
    for ch=1:num_channels
        TX_FFT = fftshift(fft(signal(:,ch)));
        TX_FFT_filtered = TX_FFT .* (rect_filter.');
        % Convert back to time domain
        signal(:, ch) = ifft(ifftshift(TX_FFT_filtered));
    end 
end