function [x1, y1] = mapUnitCircle(x0, y0)
    % take the angles of x, y and map it to the unit circle by extracting
    % the joint angles
    
    r = 1;
    
    th = atan2(y0, x0);
    y1 = r*sin(th);
    x1 = r*cos(th); 