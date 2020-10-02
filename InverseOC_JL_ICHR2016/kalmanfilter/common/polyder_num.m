function diff_coeff = polyder_num(coeff)
    crossmulti = (length(coeff)-1):-1:1;
    diff_coeff = crossmulti .* coeff(1:end-1);
    
    % validate = polyder(coeff)
end