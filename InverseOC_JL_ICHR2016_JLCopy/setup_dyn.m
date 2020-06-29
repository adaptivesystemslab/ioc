function param = setup_dyn(param)
% o=ones(1,param.knot_count+1);
o = ones(1, param.NbSample);
switch param.jointChainStart 
    case 'ankle'
        switch param.datasetTag
            case 'sim'
                Parameters_sim;
                
            otherwise
                Parameters_AnkleToHip;
        end
        
    case 'elbow'
        switch param.datasetTag
            otherwise
                Parameters_HipToAnkle;
        end
end

param.link_dh_d =[-param.L1;param.L3;param.L4;param.L5;param.L6;param.L7;-param.L8;0]';%d
param.link_dh_r =[ param.L2;0;0;0;0;0;0;0]';%R
param.link_mass = [param.M1;param.M2;param.M3;param.M4;param.M5;param.M6;param.M7;param.M8];
param.link_com = {param.COM1;param.COM2;param.COM3;param.COM4;param.COM5;param.COM6;param.COM7};
param.link_i = {IMat_huygens1;IMat_huygens2;IMat_huygens3;IMat_huygens4;IMat_huygens5;IMat_huygens6;IMat_huygens7;IMat_huygens8;IMat_huygens9;IMat_huygens10};

% matrix_conversion;
% BASE0;
% base0.gravity = gravityVal;
% base0.G2=gravityVal*o;
% z=zeros(1,1);
% for i=1:7
%     eval(['ext_wrenches.FX(' num2str(i) ',:)=z;'])
%     eval(['ext_wrenches.FY(' num2str(i) ',:)=z;'])
%     eval(['ext_wrenches.FZ(' num2str(i) ',:)=z;'])
%     eval(['ext_wrenches.CX(' num2str(i) ',:)=z;'])
%     eval(['ext_wrenches.CY(' num2str(i) ',:)=z;'])
%     eval(['ext_wrenches.CZ(' num2str(i) ',:)=z;'])
% end
% 
% param.base0 = base0;
% param.ext_wrenches = ext_wrenches;