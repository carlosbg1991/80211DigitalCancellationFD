function [ ofdm_sym_eq ] = fFD_equalize_symbols( packet , symbol_len, num_ofdm_sym )
    ofdm_sym_eq = zeros(48,num_ofdm_sym);
    for nob = 1:num_ofdm_sym
        sig=packet(1:symbol_len,nob);
        sig=sig(17:80);
        sig=fftshift(fft(sig));
        mag=(abs(sig(12))+abs(sig(26))+abs(sig(40))+abs(sig(54)))/4;    

        Polarity = [
            1, 1, 1, 1,-1,-1,-1, 1,-1,-1,-1,-1, 1, 1,-1, 1,...
            -1,-1, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1,-1, 1,...
            1, 1,-1, 1, 1,-1,-1, 1, 1, 1,-1, 1,-1,-1,-1, 1,...
            -1, 1,-1,-1, 1,-1,-1, 1, 1, 1, 1, 1,-1,-1, 1, 1,...
            -1,-1, 1,-1, 1,-1, 1, 1,-1,-1,-1, 1, 1,-1,-1,-1,...
            -1, 1,-1,-1, 1,-1, 1, 1, 1, 1,-1, 1,-1, 1,-1, 1,...
            -1,-1,-1,-1,-1, 1,-1, 1, 1,-1, 1,-1, 1, 1, 1,-1,...
            -1, 1,-1,-1,-1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1 ];

        p=Polarity(nob);

        p1 = angle( p * sig(12));
        p2 = angle( p * sig(26) * conj(p * sig(12))) + p1;
        p3 = angle( p * sig(40) * conj(p * sig(26))) + p2;
        p4 = angle(-p * sig(54) * conj(p * sig(40))) + p3;

        my=(p1+p2+p3+p4)/4;

        mx=(11+25+39+53)/4;

        var = (((11.0*11.0)+(25.0*25.0)+(39.0*39.0)+(53.0*53.0))/4)-(mx*mx);

        cov =  (((p1*11) + (p2*25) + (p3*39) + (p4*53))/4)-(mx*my);

        beta = cov / var;
        alpha = my - beta * mx;

        sfd=[sig(7:11) ;sig(13:25); sig(27:32); sig(34:39); sig(41:53); sig(55:59)];
        inds=[6:10 12:24 26:31 33:38 40:52 54:58]';
        ofdm_sym_eq(:,nob)=(sfd.*exp(complex(0,-inds*beta-alpha)))/mag;
    end
end

