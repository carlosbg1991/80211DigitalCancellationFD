function [decoded_signal,r,d_len,tail_bits] = fFD_decode_signal( eq_signal )
    inter=[0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,...
        1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,2,5,8,...
        11,14,17,20,23,26,29,32,35,38,41,44,47];

    bits=-real(eq_signal);

    bits=bits(inter+1); % deinterleave

    trellis = poly2trellis(7,[133 171]);

    hVitDec = comm.ViterbiDecoder(trellis, 'InputFormat', 'Unquantized', ...
         'TracebackDepth', 7, 'TerminationMethod', 'Terminated');

    decoded_signal = step(hVitDec,bits);
    parity=mod(sum(decoded_signal(1:17)),2);
    r=0;
    d_len=0;
    for in=1:17
        if (in<5 && decoded_signal(in))
            r=bitor(r,bitshift(1,(in-1)));
        end
        if (decoded_signal(in) && in>5)
            d_len=bitor(d_len,bitshift(1,(in-6)));
        end
    end

    tail_bits=sum(decoded_signal(19:24));

%     if parity ~= decoded_signal(18)
%         display('parity wrong');
%     end

%     switch r
%         case 11
%             display('3 Mbits/sec');
%         case 15
%             display('4.5 Mbits/sec');
%         case 10
%             display('6 Mbits/sec');
%         case 14
%             display('9 Mbits/sec');
%         case 9
%             display('12 Mbits/sec');
%         case 13
%             display('18 Mbits/sec');
%         case 8
%             display('24 Mbits/sec');
%         case 12
%             display('27 Mbits/sec');
%         otherwise
%             display('Unknown encoding, or wrong packet');
%     end    

    % [r d_len tail_bits]
    % decoded'
end