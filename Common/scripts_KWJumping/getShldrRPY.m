function shldrRPY = getShldrRPY(model,shldrJointName,side)

jointNames = {model.joints.name};
trNum = find(ismember(jointNames,shldrJointName)==1);

R_0_preShldr = model.joints(trNum).frame_in.t(1:3,1:3);
R_0_postShldr = model.joints(trNum+2).frame_out.t(1:3,1:3);

R_preShldr_postShldr = (R_0_preShldr')*R_0_postShldr;

% Debugging: choose correct frame decomposition order
% shldrRPY = zeros(12,3);
% rotOrder = {'ZYX', 'ZYZ', 'ZXY', 'ZXZ', 'YXZ', 'YXY', 'YZX', 'YZY', 'XYZ', 'XYX', 'XZY', 'XZX'};
% for i = 1:numel(rotOrder)
%     [a,b,c] = dcm2angle(R_preShldr_postShldr,rotOrder{i});
%     shldrRPY(i,:) = [a,b,c];
% end


% [a,b,c] = dcm2angle(R_preShldr_postShldr,'XYZ');
% if(side == 'L') % left shoulder
%     shldrRPY = [-wrapToPi(c+pi/2),wrapToPi(a+pi/2),b];
% elseif(side == 'R') % right shoulder
%     shldrRPY = [-wrapToPi(c+pi/2),-wrapToPi(a-pi/2),-b];
% end

% [a,b,c] = dcm2angle(R_preShldr_postShldr,'ZYZ');
% if(side == 'L') % left shoulder
%     shldrRPY = wrapToPi([-(c-pi),-(b-pi/2),-(a-pi/2)]);
% elseif(side == 'R') % right shoulder
%     shldrRPY = wrapToPi([-c,-(b-pi/2),-(a+pi/2)]);
% end

[a,b,c] = dcm2angle(R_preShldr_postShldr,'XYZ');
if(side == 'L') % left shoulder
    shldrRPY = [-wrapToPi(c+pi/2),a+pi/2,b];
elseif(side == 'R') % right shoulder
    shldrRPY = [-wrapToPi(c+pi/2),-(a-pi/2),-b];
end


% Debugging: compare calculated and actual shoulder rotations
% if(side == 'R')
%     transformNames = {model.transforms.name};
%     trNum = find(ismember(transformNames,shldrJointName)==1);
%     R_shldr_decompose = rotz(shldrRPY(1))...
%                         * model.transforms(trNum+1).t(1:3,1:3)...
%                         * rotz(shldrRPY(2))...
%                         * model.transforms(trNum+3).t(1:3,1:3)...
%                         * rotz(shldrRPY(3));
% 
%     figure(3); clf; hold on; grid on;
%     trplot(eye(3),'thick',1,'color','k','length',2,'view',[-20,20]);
%     trplot(R_preShldr_postShldr,'thick',2,'color','r','length',2,'view',[-20,20]);
%     trplot(R_shldr_decompose,'thick',2,'color','b','length',1.5,'view',[-20,20]);
%     axis equal
%     pause(0.05);
% end
