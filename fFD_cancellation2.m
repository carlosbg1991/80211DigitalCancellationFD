function [ cancellation] = fFD_cancellation2( tx_block,rx_block,H_EST,blocks )
    cancellation = zeros(80,size(tx_block,2));
    for nob = 1:blocks
        rx_pred = fft(tx_block(end-63:end,nob),64).* H_EST;
        rx_pred = ifft(rx_pred,64);
        rx_pred = [rx_pred(end-15:end);rx_pred];
        cancellation(:,nob) = rx_block(:,nob) - rx_pred;
    end
end

