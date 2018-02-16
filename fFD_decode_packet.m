function [descram_bits,r,d_len,message,crc_pass] = fFD_decode_packet( qam_packet , nsamplesXblock, num_packet,...
                             nblocks, n_bpsc, n_cbps, n_dbps, psdu_size)
    n_sym           = ceil((16 + 8*psdu_size + 6)/n_dbps);
    n_data          = n_sym * n_dbps;
    n_encoded_bits  = n_sym*n_cbps;

    descram_bits = zeros(n_data,num_packet);
    for packet = 1:num_packet
        % BLOCK EXTRACTION
        [ sym ] = fFD_blocks_extraction2(qam_packet(:,packet),nsamplesXblock,nblocks);

        % EQUALIZATION
        [ eq_sym ] = fFD_equalize_symbols(sym,nsamplesXblock,nblocks);

        % RESHAPE
        eq_symb_serie = reshape(eq_sym, 1, 48*nblocks);

        % DECODE SIGNAL AND DATA
        eq_signal = eq_symb_serie(1:n_cbps);
        [ ~ , r , d_len , ~] = fFD_decode_signal( eq_signal.' );

        % DEMODULATE DATA BPSK r=1/2
        eq_data = eq_symb_serie(49:end);
        bits_inter=-real(eq_data);

        % DATA DEINTERLEAVE
        [ cod_bits ] = fFD_deinterleave( bits_inter, n_cbps , n_bpsc, n_sym, n_encoded_bits);

        % DATA CONVOLUCIONAL DECODE
        [ scrambled_bits ] = fFD_convol_decoder( cod_bits );

        % DATA DESCRAMBLING
        [descram_bytes, descram_bits(:,packet)] = fFD_descramble( scrambled_bits );
        
        % MESSAGE DECODING
        message{packet} = char(descram_bytes(3:41));
        coded=logical(descram_bits(:,packet)');
        CRCgen = comm.CRCGenerator([32 26 23 22 16 12 11 10 8 7 5 4 2 1 0],'InitialConditions',1,'DirectMethod',1);
        [data]=step(CRCgen,coded(17:end)');
        crc_pass(packet)=sum(data(end-31:end)'.*(2.^(31:-1:0)))==3338984827;
    end
end