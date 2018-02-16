function [ t_ofdm_blocks ] = fFD_blocks_extraction2(data,sampXblock,blocks)
t_ofdm_blocks = zeros(sampXblock,blocks);
    for nob = 1:blocks
        t_ofdm_blocks(:,nob) = data(1 + sampXblock*(nob-1) : sampXblock*nob);
    end
end

