function outMtx = skewSym3(inMtx)
    % generate a skew symmetric matrix
    outMtx = [0 -inMtx(3)  inMtx(2);
       inMtx(3)         0 -inMtx(1);
      -inMtx(2)  inMtx(1)         0];
end