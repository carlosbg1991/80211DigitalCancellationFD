function [ f_off ] = fFD_frequency_offset_estimation(fs,input,type)
switch(type)
    case 'short_preamble'
        sh_dat = input(1:96);         %3 cycles
        for m=0:32
            theta(m+1) = angle(sh_dat(1+m:32+m)'*sh_dat(32+1+m:64+m))/(2*pi);
        end
        theta = mean(theta);
        f_off = theta*fs/32;        
        
    case 'long_preamble'
        lo_dat = input(33:160);
        theta = angle(lo_dat(1:64)'*lo_dat(64+1:128))/(2*pi);
        f_off = theta*fs/64;

    otherwise
        display('Wrong case');
end