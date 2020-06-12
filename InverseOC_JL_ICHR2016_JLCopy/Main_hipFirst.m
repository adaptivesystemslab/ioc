clc;clear all;close all;
rmse=@(a,b) sqrt(sum((a-b).^2)/length(a));
nrmse=@(a,b) 100*(sqrt(sum((a-b).^2)/length(a)))/sqrt(sum((a).^2)/length(a));
CC=@(a,b) corr2(a,b);

addpath('../Model/')
addpath('../Model/animate_robot_functions/')
addpath('../Model/calc_HJC/')
% addpath('../Model/data/Subject01/Session 1/')
% addpath('../Human_data/Subject06/Session 1/')

Model_name='SQUAT_7DOF';

%% Import data
Subject=1;
Trial=2;

load(['C:\Documents\aslab\data\7DOF_Planar\Human_data\Elaborated_data\Subject0',num2str(Subject),'\Kinematics_meas_Sb',num2str(Subject),'_Tr',num2str(Trial)]);
load(['C:\Documents\aslab\data\7DOF_Planar\Human_data\Elaborated_data\Subject0',num2str(Subject),'\Dynamics_meas_Sb',num2str(Subject),'_Tr',num2str(Trial)]);

%% Calcualte/import kinematics
NbSample=floor(length(Markers.RTOE));
o=ones(1,NbSample);
TimeSpan=1:NbSample;
q=q(:,TimeSpan);
dq=dq(:,TimeSpan);
ddq=ddq(:,TimeSpan);
param.Te=0.01;


figure;hold on;
plot(q(1,:),'r')
plot(q(2,:),'g')
plot(q(3,:),'b')
plot(q(4,:),'y')
plot(q(5,:),'m')
plot(q(6,:),'k')
plot(q(7,:),'c')
ylabel('Joint angle [rad]')

figure;subplot(211);hold on;
plot(dq(1,:),'r')
plot(dq(2,:),'g')
plot(dq(3,:),'b')
plot(dq(4,:),'y')
plot(dq(5,:),'m')
plot(dq(6,:),'k')
plot(dq(7,:),'c')
ylabel('Joint velocity [rad.s-1]')
subplot(212);hold on;
plot(ddq(1,:),'r')
plot(ddq(2,:),'g')
plot(ddq(3,:),'b')
plot(ddq(4,:),'y')
plot(ddq(5,:),'m')
plot(ddq(6,:),'k')
plot(ddq(7,:),'c')
ylabel('Joint acc [rad.s-2]')

% %% Calculate q
% % calculate HJC
% Markers.Mid_ASIS=(Markers.LASI+Markers.RASI)/2;
% Markers.Mid_PSIS=(Markers.LPSI+Markers.RPSI)/2;
% Markers.Mid_Pelvis=(Markers.Mid_ASIS+Markers.Mid_PSIS)/2;
% Markers.Mid_LTorso=(Markers.T10+Markers.STRN)/2;
% Markers.Mid_HTorso=(Markers.C7+Markers.CLAV)/2;
% 
% for i=1:NbSample
%     PW(i,:)=norm(Markers.Mid_ASIS(i,:)-Markers.Mid_PSIS(i,:));% Pelvis width
% end
% PW=mean(PW(1:100));% Only take the 100 first samples during no motion since marker data need to be improved
% 
% % Hip joint center position in pelvis (HJCpv) reference frame [Bell, 1990]
% % right hip joint center
% HJCpv(1)=-0.14*PW;
% HJCpv(2)=-0.30*PW;
% HJCpv(3)= 0.24*PW;
% 
% %%Determine orientation of the pelvis in the global reference frame
% [gOl gRl] = ARF(Markers, 'pelvis');
% 
% % transform HJC from pelvic to global reference frame
% Markers.RHJC = transform(gOl, gRl, repmat(HJCpv,size(gRl,1),1),2);%%%%%%%
% 
% 
% Markers.RANKproj=Markers.RANK;
% Markers.RANKproj(:,1)=Markers.RANKproj(:,1)+0.1*ones(size(Markers.RANKproj,1),1);
% Markers.base0=zeros(size(Markers.RANK,1),size(Markers.RANK,2));
% Markers.baseX0=Markers.base0;
% Markers.baseX0(:,1)=-ones(size(Markers.RANK,1),1);
% for ii_mk=1:NbSample
%     
%     q(1,ii_mk)=0*(angle2d([  Markers.base0(ii_mk,1:2) Markers.baseX0(ii_mk,1:2) Markers.RTOE(ii_mk,1:2) Markers.RHEE(ii_mk,1:2)]));
%     q(2,ii_mk)=(angle2d([Markers.RANKproj(ii_mk,1:2) Markers.RANK(ii_mk,1:2) Markers.RKNE(ii_mk,1:2) Markers.RANK(ii_mk,1:2)]));
%     q(3,ii_mk)=(angle2d([Markers.RANK(ii_mk,1:2) Markers.RKNE(ii_mk,1:2)   Markers.RKNE(ii_mk,1:2) Markers.RHJC(ii_mk,1:2)  ]));
%     q(4,ii_mk)=(angle2d([Markers.RHJC(ii_mk,1:2) Markers.RKNE(ii_mk,1:2) Markers.Mid_Pelvis(ii_mk,1:2) Markers.RHJC(ii_mk,1:2) ]));
%     q(5,ii_mk)=(angle2d([Markers.Mid_Pelvis(ii_mk,1:2) Markers.RHJC(ii_mk,1:2) Markers.Mid_LTorso(ii_mk,1:2) Markers.Mid_Pelvis(ii_mk,1:2) ]));
%     q(6,ii_mk)=(angle2d([Markers.Mid_HTorso(ii_mk,1:2) Markers.Mid_LTorso(ii_mk,1:2)  Markers.Mid_LTorso(ii_mk,1:2) Markers.Mid_Pelvis(ii_mk,1:2)]));
%     q(7,ii_mk)=(angle2d([Markers.Mid_LTorso(ii_mk,1:2) Markers.Mid_HTorso(ii_mk,1:2) Markers.RELB(ii_mk,1:2) Markers.Mid_HTorso(ii_mk,1:2)  ]));
% end
% 
% 
% param.Te=0.01;%Sample time 100Hz
% TimeSpan=1:NbSample;
% cut_off=5;
% %Filter joint angles
% q=lowpass_filter(q,1/param.Te,cut_off,5);
% % Compute dq and ddq
% deriv=First_order_diff(q,param.Te,2,1);
% dq=deriv(1:7,:);
% ddq=deriv(8:end,:);


%% Geometrical parameters and base0
param.L1=abs(Markers.RANK(1,1));
param.L2=abs(Markers.RANK(1,2));
param.L3=norm(Markers.RANK(1,1:2)-Markers.RKNE(1,1:2));
param.L4=norm(Markers.RHJC(1,1:2)-Markers.RKNE(1,1:2));
param.L5=norm(Markers.RHJC(1,1:2)-Markers.Mid_Pelvis(1,1:2));
param.L6=norm(Markers.Mid_LTorso(1,1:2)-Markers.Mid_Pelvis(1,1:2));
param.L7=norm(Markers.Mid_LTorso(1,1:2)-Markers.Mid_HTorso(1,1:2));
param.L8=norm(Markers.RELB(1,1:2)-Markers.Mid_HTorso(1,1:2));
param.body_weight=body_weight;
param.body_height=body_height;

Parameters;
matrix_conversion;
BASE0
base0.G2=-9.81*o;




%% Regressor
R=eval([Model_name '_dim2(q,dq,ddq,base0,param)']);
%% GAMMA and EN
z=zeros(1,1);
for i=1:7
    eval(['ext_wrenches.FX(' num2str(i) ',:)=z;'])
    eval(['ext_wrenches.FY(' num2str(i) ',:)=z;'])
    eval(['ext_wrenches.FZ(' num2str(i) ',:)=z;'])
    eval(['ext_wrenches.CX(' num2str(i) ',:)=z;'])
    eval(['ext_wrenches.CY(' num2str(i) ',:)=z;'])
    eval(['ext_wrenches.CZ(' num2str(i) ',:)=z;'])
end

[GAMMA EN]=eval([Model_name '_dyn3(q,dq,ddq,base0,ext_wrenches,param)']);

Est_IDM=[EN(1,:) EN(2,:) EN(6,:)];


figure;hold on;
plot(GAMMA(1,:),'r')
plot(GAMMA(2,:),'g')
plot(GAMMA(3,:),'b')
plot(GAMMA(4,:),'y')
plot(GAMMA(5,:),'m')
plot(GAMMA(6,:),'k')
ylabel('Joint Torque [N.m]')

%%%%%%% Check in Forward_kin if the right model_name for DH
%% COM
for ii=1:NbSample
    Tr=forward_kin([0;0;0],q(:,ii),param);
    CoM=eval(['COM_' Model_name '(Tr,param);']);
    CoM(:,8)*param.Mtot*9.81;
end
%% COP
COP=CoP_calculation(EN);
%% plot FX,Fy,Fz,Mx,My,Mz
Est_IDM=Est_IDM';
figure;hold on;
plot(Est_IDM, 'g')
Fp_plotting=[-FP.Force(TimeSpan,1)' -FP.Force(TimeSpan,2)' -(FP.Moment(TimeSpan,3))'];
plot(Fp_plotting, 'b')
legend('Model dyn','Plate dyn');


%% Animation

base00=zeros(3,NbSample);
animate_robot(base00,[q],param,param,Markers,FP,TimeSpan,NbSample)
