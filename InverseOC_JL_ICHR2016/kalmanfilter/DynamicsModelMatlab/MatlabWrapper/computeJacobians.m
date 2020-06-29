function [J, dJ] = computeJacobians(mdl,ee_frame_name)

ee = mdl.getFrameByName(ee_frame_name);

Jv = zeros(3,numel(mdl.joints));
Jw = zeros(3,numel(mdl.joints));
Jv_dot = zeros(3,numel(mdl.joints));
Jw_dot = zeros(3,numel(mdl.joints));
Oei_mat = zeros(3,numel(mdl.joints));

T0e = ee.t;
for i=1:numel(mdl.joints)
    
    T0i = mdl.joints(i).frame_in.t;
    Oei = T0e(1:3,4) - T0i(1:3,4);
    Oei_mat(:,i) = Oei;
    Zi = T0i(1:3,3);
    if(strcmp(mdl.joints(i).type,'revolute'))
        Jv(:,i) = cross(Zi,Oei);
        Jw(:,i) = Zi;
    else
        Jv(:,i) = Zi;
        Jw(:,i) = zeros(3,1);
    end
end

for i=1:numel(mdl.joints)
    
    T0i = mdl.joints(i).frame_in.t;
    Oei = T0e(1:3,4) - T0i(1:3,4);
    
    Zi = T0i(1:3,3);
    W0im1 = Jw(1:3,1:i-1)*mdl.velocity(1:i-1);
    W0i = Jw(1:3,1:i)*mdl.velocity(1:i);
    Zi_dot = cross(W0im1,Zi);
    Oei_dot = cross(W0i,Oei) + Jv(1:3,i+1:end)*mdl.velocity(i+1:end);
    
    if(strcmp(mdl.joints(i).type,'revolute'))
        Jv_dot(:,i) = cross(Zi_dot,Oei) + cross(Zi,Oei_dot);
        Jw_dot(:,i) = Zi_dot;
    else
        Jv_dot(:,i) = Zi_dot;
        Jw_dot(:,i) = zeros(3,1);
    end
end

J = [Jv;Jw];
dJ = [Jv_dot; Jw_dot];

end