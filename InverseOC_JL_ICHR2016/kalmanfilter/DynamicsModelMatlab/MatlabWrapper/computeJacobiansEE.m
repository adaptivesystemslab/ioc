function [J, dJ, Zi_dot_mat,Oei_dot_mat] = computeJacobiansEE(mdl,ee_frame_name)

ee = mdl.getFrameByName(ee_frame_name);

Jv = zeros(3,numel(mdl.joints));
Jw = zeros(3,numel(mdl.joints));
Jv_dot = zeros(3,numel(mdl.joints));
Jw_dot = zeros(3,numel(mdl.joints));
Oei_mat = zeros(3,numel(mdl.joints));

T0e = ee.t;
Te0 = T0e';
Te0(4,1:3) = 0;
Te0(1:3,4) = -Te0(1:3,1:3)*T0e(1:3,4);

for i=1:numel(mdl.joints)
    T0i = mdl.joints(i).frame_in.t;
    Tei = Te0*T0i;
    Oei = -Tei(1:3,4);
    Oei_mat(:,i) = Oei;
    Zi = Te0(1:3,1:3)*T0i(1:3,3);
    
    if(strcmp(mdl.joints(i).type,'revolute'))
        Jv(:,i) = cross(Zi,Oei);
        Jw(:,i) = Zi;
    else
        Jv(:,i) = Zi;
        Jw(:,i) = zeros(3,1);
    end
end

Oei_dot_mat =zeros(3,numel(mdl.joints));
Zi_dot_mat =zeros(3,numel(mdl.joints));
for i=1:numel(mdl.joints)
    
    T0i = mdl.joints(i).frame_in.t;
    Tei = Te0*T0i;
    Oei = -Tei(1:3,4);
    Zi = Te0(1:3,1:3)*T0i(1:3,3);
    
    W0im1 = Jw(1:3,1:i-1)*mdl.velocity(1:i-1);
    W0i = Jw(1:3,1:i)*mdl.velocity(1:i);
    Zi_dot = cross(W0im1,Zi);
    Zi_dot_mat(:,i) = Zi_dot;
    Oei_dot = cross(W0i,Oei) + Jv(1:3,i+1:end)*mdl.velocity(i+1:end);
    Oei_dot = ee.v(4:6) - Tei(1:3,1:3)*mdl.joints(i).frame_in.v(4:6);
    Oei_dot_mat(:,i) = Oei_dot;
    
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