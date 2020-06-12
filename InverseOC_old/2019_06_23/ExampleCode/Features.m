function [Jx,Ju]=Features(x,u)

Jx=zeros(2,length(x));

Ju=[u(1) 0;
    0  u(2)];

