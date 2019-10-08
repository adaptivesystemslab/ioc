function shift = shiftMatchJAs(x,y,maxShift)
% tests all amplitude matches of signals x and y up to +/- maxsamp
% x and y must be vectors of the same length

debug = 0;

x = x./max(abs(y)); %normalize signals to approx -1 < x < 1
y = y./max(abs(y));

if(maxShift >= length(x)-1)
    maxShift = length(x)-2;
end

shiftErr = zeros(maxShift*2+1,1);

%first shift y back relative to x
for i = maxShift:-1:1
    tmp_x = x(1:end-i);
    tmp_y = y(1+i:end);
    % NOTE: both the equations below were tested, first one aligns better
    %       (1st equation "ignores" variable trajectory sections better)
    shiftErr(maxShift-i + 1) = sum(abs(tmp_x - tmp_y))/length(tmp_x);
%     shiftErr(maxShift-i + 1) = sum((tmp_x - tmp_y).^2)/length(tmp_x);
    
    if(debug)
        figure(1); clf; hold on;
        plot(tmp_x,'k');
        plot(tmp_y,'r');
    end
    
end

%error for no shifting
shiftErr(maxShift + 1) = sum(abs(x - y))/length(x);

if(debug)
        figure(1); clf; hold on;
        plot(tmp_x,'k');
        plot(tmp_y,'r');
end
    
%finally, shift y forward relative to x
for i = 1:maxShift
    tmp_x = x(1+i:end);
    tmp_y = y(1:end-i);
    shiftErr(i + maxShift + 1) = sum(abs(tmp_x - tmp_y))/length(tmp_x);
    
    if(debug)
        figure(1); clf; hold on;
        plot(tmp_x,'k');
        plot(tmp_y,'r');
    end
    
end

[~,shift] = min(shiftErr);
shift = shift - maxShift - 1;

