clear
clc

%implementation of cost weights tracking

%loading trajectory
data=load('traj1.mat');
traj_u=data.control;
traj_x=data.state;
traj_t=data.time;
clear data

%sampling the trajectory with white noise
sig=1e-4;
traj_st=linspace(0,1,2001);
traj_su=interp1(traj_t,traj_u,traj_st,'spline')+sig*randn(length(traj_st),size(traj_u,2));
traj_sx=interp1(traj_t,traj_x,traj_st,'spline')+sig*randn(length(traj_st),size(traj_x,2));
%sampling rate
dt=traj_st(2)-traj_st(1);
%samples count
maxlength=length(traj_st);

%initialization
t=2; %start from sample 2
traj_costs=[]; %store the trajectory cost
traj_size=[]; %store the window size trajectory
traj_ct=[];

%feature count
fcount=2;


%threshold
gamma=1300;
delta=1e-6;
n=20;
maxcapture=2000;

%stop flag
stop=0;
while stop==0
    %initialize the standing point and windowsize
    x=traj_sx(t,:)';
    u=traj_su(t-1,:)';
    x0=traj_sx(t-1,:)';
    u1=traj_su(t,:)';
    s=0; %window size
    %initilize window-information matrix H
    [px,pu]=Features(x,u);
    px=px';
    pu=pu';
    [xnext,fx,fu]=Arm_Dyn(x,u1,dt); %#ok<ASGLU>
    df_dx=fx';
    [xnext,fx,fu]=Arm_Dyn(x0,u,dt);
    df_du=fu';        
    H2=df_du*df_dx;
    H1=df_du*px+pu;
    H=[H1,-H2];
    Hhat=H/norm(H,'fro'); %normalize
    
    while Judgement(Hhat,gamma,delta,n)==0 && stop==0 && s<maxcapture
        s=s+1;
        tt=s+t;
        if tt> maxlength
            stop=1;
        else
            %capture the next sample and update H
            x=traj_sx(tt,:)';
            u=traj_su(tt-1,:)';
            x0=traj_sx(tt-1,:)';
            u1=traj_su(tt,:)';
            [px,pu]=Features(x,u);
            px=px';
            pu=pu';
            [xnext,fx,fu]=Arm_Dyn(x,u1,dt);%#ok<ASGLU>
            df_dx=fx';
            [xnext,fx,fu]=Arm_Dyn(x0,u,dt);
            df_du=fu';
            H1=[H1+H2*px;
                df_du*px+pu];
            H2=[H2*df_dx;
                df_du*df_dx];
            H=[H1 -H2]; 
            Hhat=H/norm(H,'fro'); %normalize H
        end
        
    end
    
    if (s==maxcapture)
        traj_st(t)
        v=traj_costs(end,:)';
    else
        % nomalize the costs
        Hhat=H/norm(H,'fro'); %normalize H
        [U,S,V]=svd(Hhat);
        v=V(:,end);
        v=v(1:fcount);
        v=v/sum(v);
    end

    %store solutions
    traj_costs(end+1,:)=v';
    traj_size(end+1)=s;
    traj_ct(end+1)=t-1;
    %move the window
    t=t+1
end

groundtruth=[0.6,0.4];

% endi=1450;
% traj_ct=traj_ct(1:endi);
% traj_costs=traj_costs(1:endi,:);
% traj_size=traj_size(1:endi);


% figure(1)
% subplot(2,1,1)
% plot(traj_ct,traj_costs(:,1))
% hold on
% plot(traj_ct,traj_costs(:,2))
% legend('$\hat{\omega}_1$','$\hat{\omega}_2$','Interpreter','latex')
% xlabel('Observation Starting Time')
% ylabel('Recovered Weights')
% axis([min(traj_ct),max(traj_ct),0,0.8]);
% 
% subplot(2,1,2)
% plot(traj_ct,traj_size)
% xlabel('Observation Starting Time')
% ylabel('Minimal Observation Length')
% axis([min(traj_ct),max(traj_ct),0,max(traj_size)*1.2]);

%compute error
error=Compute_Error(traj_costs,groundtruth);
final_error=mean(error)
final_length=mean(traj_size)