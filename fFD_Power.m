function pow = fFD_Power( input )
    pow = sum(input.*conj(input))./length(input);
end