function [ output ] = fFD_channel_equal(input_symbols,H_EST)
    output1 = input_symbols(2:27,nob,np)./H_EST(2:27);
    output2 = input_symbols(39:end,nob,np)./H_EST(39:end);
    output = [output2;output1];
end