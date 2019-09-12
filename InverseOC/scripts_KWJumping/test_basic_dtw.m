t = 0:0.01:1;
b = 5;

x = sin(b*t);
y = sin(b*t) + 0.1;

figure(1); clf;
subplot(1,3,1); hold on; grid on;
plot(x,'r'); 
plot(y,'b');

[dist,ix,iy] = dtw(x,y);
subplot(1,3,2); hold on; grid on;
plot(x,'r'); 
plot(y,'b');
plot(x(ix),'r--');
plot(y(iy),'b--');

[acor,lag] = xcorr(x,y);
[~,I] = max(abs(acor));
lagDiff = lag(I);

subplot(1,3,2); hold on; grid on;
if(lagDiff > 0)
    x_lag = x(1+lagDiff:end);
    y_lag = y(1:end-lagDiff);
else
    x_lag = x(1:end-lagDiff);
    y_lag = y(1+lagDiff:end);
end

plot(x_lag,'r'); 
plot(y_lag,'b');

