function [ deinter ] = fFD_deinterleave( bits_inter, n_cbps , n_bpsc, n_sym, n_encoded_bits)
    first=zeros(n_cbps,1);
	second=zeros(n_cbps,1);
	s = max(n_bpsc / 2, 1);
    deinter=zeros(n_encoded_bits,1);
	for j=1:n_cbps 
		first(j) = (s * ((j-1) / s) + mod((j + (floor(16.0 * (j-1) / n_cbps))) , s))+1;
    end

	for i=1:n_cbps
		second(i) = (16 * (i-1) - (n_cbps - 1) * (floor(16.0 * (i-1) / n_cbps)))+1;
    end

	for i = 1:n_sym 
		for k=1:n_cbps
			deinter((i-1) * n_cbps + second(first(k))) = bits_inter((i-1) * n_cbps + k);
        end
    end
end

