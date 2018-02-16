function [ H_EST ] = fFD_channel_est(fft_tx,fft_rx)
    H_EST_ML1 = ( conj(fft_tx(2:27)).*fft_rx(2:27) ) ./ ( conj(fft_tx(2:27)).*fft_tx(2:27) );
    H_EST_ML2 = ( conj(fft_tx(39:end)).*fft_rx(39:end) ) ./ ( conj(fft_tx(39:end)).*fft_tx(39:end) );
    H_EST = [0;H_EST_ML1;zeros(11,1);H_EST_ML2];

%     h_est = ifft(H_EST,64);
%     h_est = h_est(1:5);
%     H_EST = fft(h_est,64);
end