function stop=Judgement(H,gamma,delta,n)

stop=0;
[r,c]=size(H);
s=svd(H);
s=sort(s);
if (s(1)<delta && abs(s(2)/s(1))>gamma && r>=n*c)
    stop=1;
end


