function [S2sc_optim_L, S2sc_optim_R, normUpArm_L, normUpArm_R] = getS2scDrop(trc,offset3D,window) %,S2sc_estim)
% finds offset between shoulder marker on acromion and shoulder joint
% center that minimizes CHANGE in norm distance between shoulder and elbow 
% joint centers. Units in METERS.

% if "3Doffset" = 0, assume shoulder center is directly below SL/SR marker
% if "3Doffset" = 0, find offset

% fmincon Variables
offsetZ_lb = -0.15;
offsetZ_ub = 0; % cannot be above shoulder marker
offset3D_lb = [-0.1,-0.1,-0.15];
offset3D_ub = [0.1,0.1,0];
S2sc_estim = -0.060; %average global Z offset from SL/SR marker to shoulder center

options = optimoptions('fmincon','Display','off');

sf = window(1); %start frame
ef = window(2); %end frame

% Left Shoulder
S_marker = trc.data.SL(sf:ef,:)/1000; %shoulder marker pos
E_center = ((trc.data.ELLat(sf:ef,:) + trc.data.ELMed(sf:ef,:))/2)/1000; %elbow center pos
% fun = @(x)sum(abs(norm( (S_marker + [x(1),x(2),x(3)]) - E_center) - mean(norm( (S_marker + [x(1),x(2),x(3)]) - E_center))));
if(~offset3D)
    fun = @(x)sum(abs(sqrt(sum(( (S_marker + [0,0,x]) - E_center).^2, 2)) - mean(sqrt(sum(( (S_marker + [0,0,x]) - E_center).^2, 2)))));
    S2sc_optim_L = fmincon(fun,S2sc_estim,[],[],[],[],offsetZ_lb,offsetZ_ub,[],options);
    normUpArm_L = sqrt(sum(( (S_marker + [0,0,S2sc_optim_L]) - E_center).^2, 2));
else
    fun = @(x)sum(abs(sqrt(sum(( (S_marker + [x(1),x(2),x(3)]) - E_center).^2, 2)) - mean(sqrt(sum(( (S_marker + [x(1),x(2),x(3)]) - E_center).^2, 2)))));
    S2sc_optim_L = fmincon(fun,[0,0,S2sc_estim],[],[],[],[],offset3D_lb,[100,100,0],[],options);
    normUpArm_L = sqrt(sum(( (S_marker + S2sc_optim_L) - E_center).^2, 2));
end

% Right Shoulder
S_marker = trc.data.SR(sf:ef,:)/1000; %shoulder marker pos
E_center = ((trc.data.ERLat(sf:ef,:) + trc.data.ERMed(sf:ef,:))/2)/1000; %elbow center pos
if(~offset3D)
    fun = @(x)sum(abs(sqrt(sum(( (S_marker + [0,0,x]) - E_center).^2, 2)) - mean(sqrt(sum(( (S_marker + [0,0,x]) - E_center).^2, 2)))));
    S2sc_optim_R = fmincon(fun,S2sc_estim,[],[],[],[],offsetZ_lb,offsetZ_ub,[],options);
    normUpArm_R = sqrt(sum(( (S_marker + [0,0,S2sc_optim_R]) - E_center).^2, 2));
else
    fun = @(x)sum(abs(sqrt(sum(( (S_marker + [x(1),x(2),x(3)]) - E_center).^2, 2)) - mean(sqrt(sum(( (S_marker + [x(1),x(2),x(3)]) - E_center).^2, 2)))));
    S2sc_optim_R = fmincon(fun,[0,0,S2sc_estim],[],[],[],[],offset3D_lb,offset3D_ub,[],options);
    normUpArm_R = sqrt(sum(( (S_marker + S2sc_optim_R) - E_center).^2, 2));
end

