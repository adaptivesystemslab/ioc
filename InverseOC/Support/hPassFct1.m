function [b, have, want] = hPassFct1(Hhat, dimWeights)
    [rows, cols] = size(Hhat);
    have = rows;
    want = dimWeights * cols;
    b = rows >= dimWeights * cols;
end