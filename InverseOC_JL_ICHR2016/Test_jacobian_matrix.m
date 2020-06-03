clc;close all;clear all;


param.body_height=1.70;
param.body_weight=70;
param.L1=param.body_height*(1/1.77)*0.219;
param.L2=param.body_height*(1/1.77)*0.06;
param.L3=param.body_height*(1/1.77)*0.433;
param.L4=param.body_height*(1/1.77)*0.432;
param.L5=param.body_height*(1/1.77)*0.094;
param.L6=param.body_height*(1/1.77)*0.151;
param.L7=param.body_height*(1/1.77)*0.334;
param.L8=param.body_height*(1/1.77)*0.271;
 Parameters
param.L=[param.L1 param.L2 param.L3 param.L4 param.L5...
 param.L6 param.L7 param.L8 param.L9 param.L10 param.L11];

Te=0.01;
time=0:Te:2;
NbSample=length(time);
for ii=1:7
q(ii,:)=0.1*cos(2*pi*time);
dq(ii,:)=-0.1*2*pi*sin(2*pi*time);
ddq(ii,:)=-0.1*(2*pi)^2*cos(2*pi*time);

end

J=SQUAT_7DOF_Ext_wrenches_Jacobian(q,param);
 
Tr = FKM_SQUAT_7DOF(repmat(eye(4),[1 1 NbSample]),q',param.L);
 

dJ = dJ_SQUAT_7DOF_Ext_wrenches_Jacobian(q,dq,param.L);


%% calculation of JpQp
 dXnum=zeros(3,NbSample);

for ii=2:NbSample
    dXnum(:,ii)= (Tr(1:3,4,end,ii) -Tr(1:3,4,end,ii-1) )/Te;
    ddXnum(:,ii)= (dXnum(:,ii) -dXnum(:,ii-1) )/Te;
    dXJ(:,ii)=J([1 2 6],:,ii)*dq(:,ii);
    ddX(:,ii)=J([1 2 6],:,ii)*ddq(:,ii)+dJ([1 2 6],:,ii)*dq(:,ii); 
    
end


figure;
subplot(211);hold on;
plot(dXJ(1,:),'r')
plot(dXJ(2,:),'g')
%plot(dXJ(3,:),'b')

plot(dXnum(1,:),'--r')
plot(dXnum(2,:),'--g')
subplot(212);hold on;
plot(ddX(1,:),'r')
plot(ddX(2,:),'g')
%plot(dXJ(3,:),'b')

plot(ddXnum(1,:),'--r')
plot(ddXnum(2,:),'--g')
%  








