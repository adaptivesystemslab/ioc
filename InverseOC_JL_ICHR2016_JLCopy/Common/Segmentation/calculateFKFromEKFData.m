function YrecoverEF = calculateFKFromEKFData(q)

	link1Marker =  [0, 0,    -0.428];
    link2Marker =  [0, 0,    -0.428];

    ulLenRot = rotz(pi)*roty(pi/2); % proper
    llLenRot = rotz(-pi/2)*rotx(pi);
    ulMarker = ulLenRot'*link1Marker'; % joint 1 to joint 2 length
    llMarker = llLenRot'*link2Marker'; % joint 1 to shimmer length
    
    YrecoverEF = zeros(size(q, 1), 6);
    
    for ind = 1:size(q, 1)
        currQ = q(ind, 1:5);
        [rotMtxPost1, rotMtxPost2] = rotmtxCombined(currQ);
        YrecoverEF(ind, 1:3) = rotMtxPost1*ulMarker;
        YrecoverEF(ind, 4:6) = rotMtxPost1*ulMarker + rotMtxPost1*rotMtxPost2*llMarker;
    end

function [R03, R35] = rotmtxCombined(q)
% calculate FK for EKF leg data
R01 = hfun_R01(q(1));
R12 = hfun_R12(q(2));
R23 = hfun_R23(q(3));
R34 = hfun_R34(q(4));
R45 = hfun_R45(q(5));

R03 = R01*R12*R23;
R35 = R34*R45;

function R = hfun_R01(q)
R = [-sin(q),  0, cos(q);
    cos(q),  0, sin(q);
    0,  1,      0];

function R = hfun_R12(q)
R = [ -sin(q),  0, -cos(q);
    cos(q),  0, -sin(q);
    0, -1,       0];

function R = hfun_R23(q)
R = [  sin(q),  0, cos(q);
    -cos(q),  0, sin(q);
    0, -1,      0];

function R = hfun_R34(q)
R = [ -sin(q), 0, cos(q);
    cos(q), 0, sin(q);
    0, 1,      0];

function R = hfun_R45(q)
R = [ cos(q), -sin(q), 0;
    sin(q),  cos(q), 0;
    0,       0, 1];