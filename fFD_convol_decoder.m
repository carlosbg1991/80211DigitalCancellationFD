function [ decoded_bits ] = fFD_convol_decoder( coded_bits )
    trellis = poly2trellis(7,[133 171]);

    hVitDec = comm.ViterbiDecoder(trellis, 'InputFormat', 'Unquantized', ...
        'TracebackDepth', 7, 'TerminationMethod', 'Terminated');

    decoded_bits = step(hVitDec,coded_bits);
end