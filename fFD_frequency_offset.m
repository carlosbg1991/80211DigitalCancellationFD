function [ output ] = fFD_frequency_offset(fs,fo,input)
    CFO = fo/fs;
    phase_freq_offset = exp(1i*2*pi*(0:length(input)-1)*CFO).';
    output = input .* phase_freq_offset;
end