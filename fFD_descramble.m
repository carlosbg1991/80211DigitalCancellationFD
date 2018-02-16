function [out_bytes, out_bits] = fFD_descramble( decoded_bits )
    state = 0;
    dsize=length(decoded_bits);
	for i = 1:8
		if(decoded_bits(i)) 
			state = bitor(state , 2^(8 - i));
        end
    end

	feedback=0;

	for i = 9:dsize
		feedback = bitxor((~~bitand(state , 64)) , (~~bitand(state , 8)));
		out_bits(i) = bitxor(feedback , decoded_bits(i));
		state = bitor(bitand(bitshift(state, 1) , 126) , feedback);
    end

	for i = 1:dsize
		bit = mod(i-1, 8);
		byte = floor((i-1) / 8)+1;
		if(bit == 0) 
			out_bytes(byte) = 0;
        end

		if (out_bits(i))
			out_bytes(byte) = bitor(out_bytes(byte) , bitshift(1,bit));
        end
    end
end

