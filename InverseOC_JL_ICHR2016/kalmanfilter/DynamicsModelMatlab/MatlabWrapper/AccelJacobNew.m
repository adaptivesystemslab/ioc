% Build and test Accelerometer Jacobian

%These are the values we test with
pos_init = [0.6787    0.7577    0.7431    0.3922    0.6555]'*10;
%pos_init = [0 0 0 0 0]';
vel_init = [0.6787    0.7577    0.7431    0.3922    0.6555]'*10;
accel_init = [0.9787    1.7577    2.7431    0.3322    0.7555]'*10;

mdl = rlCModel('Models/arm_fixed.xml');
sens = SensorCore('IMU');
sens.addDecorator('accelerometer');
mdl.forwardPosition();
mdl.addSensor(sens,'post:rradius',eye(4));
mdl.forwardPosition();
mdl.g = [0 0 0]';

vis = rlVisualizer('vis',640,480);
vis.addModel(mdl);
vis.update();

%% Lets compute the Jacobian Numerically, Assuming joint Velocity is 0
dxdd_dqdd = zeros(numel(sens.measurement),numel(mdl.joints));
mdl.position = pos_init;
mdl.velocity(:) = 0;
mdl.acceleration = accel_init;
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mes_init = sens.measurement;
eps = 1e-7;
for i=1:numel(mdl.joints)
    accel = accel_init;
    accel(i) = accel(i) + eps;
    mdl.acceleration = accel;
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    dxdd_dqdd(:,i) = (sens.measurement-mes_init)./eps;
end
disp('Derivative of acceleration wrt joint acceleration Jacobian Differences');
mdl.calculateSensorJacobians();
disp(norm(sens.obsJacobian(:,2*end/3+1:end)-dxdd_dqdd));

%% Now lets look at derivative WRT q assuming gravity is zero
mdl.velocity(:) = 0;
mdl.position = pos_init;
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();

%Now lets try to do it analytically
S = skew([0 0 -1]);
T_0ee = sens.transform;
T_ee0 = T_0ee;
T_ee0(1:3,1:3) = T_ee0(1:3,1:3)';
T_ee0(1:3,4) = -T_ee0(1:3,1:3)*T_ee0(1:3,4);

%Compute dri_dthetac and dei_dthetac

%Derivative of ri in end effector frame wrt joint positions
dr_dtheta = zeros(3,numel(mdl.joints),numel(mdl.joints));
dr_dtheta_new = dr_dtheta;
%Derivative of ei in end effector frame wrt joint positions
de_dtheta = zeros(3,numel(mdl.joints),numel(mdl.joints));
%all of the ri in end effector frame
ris_ee = zeros(3,numel(mdl.joints));
ris_b = zeros(3,numel(mdl.joints));
%all of the ei in end effector frame
eis_ee = zeros(3,numel(mdl.joints));
eis_b = zeros(3,numel(mdl.joints));

for i=1:numel(mdl.joints)
    
    T_0i = mdl.joints(i).frame_in.t;
    T_i0 = T_0i';
    T_i0(4,1:3) = 0;
    T_i0(1:3,4) = -T_i0(1:3,1:3)*T_0i(1:3,4);
    
    T_eei = T_ee0*T_0i;
    T_iee = T_i0*T_0ee;
    
    ri = -T_eei(1:3,4);
    ris_ee(:,i) = ri;
    
    ei = T_eei(1:3,3);
    eis_ee(:,i) = ei;
end

for i=1:numel(mdl.joints)
    T_0i = mdl.joints(i).frame_in.t;
    T_i0 = T_0i';
    T_i0(4,1:3) = 0;
    T_i0(1:3,4) = -T_i0(1:3,1:3)*T_0i(1:3,4);
    
    T_eei = T_ee0*T_0i;
    T_iee = T_i0*T_0ee;
    %The Jacobian Jv is based on ri and ei. Since accelerometer is in end
    %effector only the joints after i change ri and ei when moved
    for c=i:numel(mdl.joints)
        T_0c = mdl.joints(c).frame_in.t;
        T_c0 = T_0c';
        T_c0(4,1:3) = 0;
        T_c0(1:3,4) = -T_c0(1:3,1:3)*T_0c(1:3,4);
        T_eec = T_ee0*T_0c;
        T_cee = T_c0*T_0ee;
        R_eec = T_eec(1:3,1:3);
        R_c0 = T_c0(1:3,1:3);
        
        T_ci = T_c0*T_0i;
        T_ic = T_i0*T_0c;
        R_ci = T_ci(1:3,1:3);
        %From joint i to joint c in frame i
        if(strcmp(mdl.joints(c).type,'revolute'))
            
            %joint i to joint c in end effector frame
            r_ci_ee = ris_ee(:,i) - ris_ee(:,c);
            r_ci_i = T_iee(1:3,1:3)*r_ci_ee;
            r_ci_i = T_ic(1:3,4);
            
            
            dri_dc = R_eec*skew([0 0 1])*R_ci*r_ci_i;
            dr_dtheta(:,i,c) = dri_dc;
            
            dei_dc = R_eec*skew([0 0 1])*R_ci*[0 0 1]';
            de_dtheta(:,i,c) = dei_dc;
        else
            dr_dtheta(:,i,c) = eis_ee(:,c);
        end
        
    end
end

%Numerically verify dri_dtheta
dr_dtheta_num = zeros(size(dr_dtheta));
de_dtheta_num = zeros(size(dr_dtheta));
for i=1:numel(mdl.joints)
    
    T_0ee = sens.transform;
    T_ee0 = T_0ee;
    T_ee0(1:3,1:3) = T_ee0(1:3,1:3)';
    T_ee0(1:3,4) = -T_ee0(1:3,1:3)*T_ee0(1:3,4);
    
    T_0i = mdl.joints(i).frame_in.t;
    T_i0 = T_0i';
    T_i0(4,1:3) = 0;
    T_i0(1:3,4) = -T_i0(1:3,1:3)*T_0i(1:3,4);
    
    T_eei = T_ee0*T_0i;
    %From joint i to end effector in end effector frame
    r_i_ee_init = T_eei(1:3,4);
    e_i_ee_init = T_eei(1:3,3);
    for c=i:numel(mdl.joints)
        pos = pos_init;
        pos(c) = pos(c) + eps;
        mdl.position = pos;
        mdl.forwardPosition;
        
        T_0ee = sens.transform;
        T_ee0 = T_0ee';
        T_ee0(4,1:3) = 0;
        T_ee0(1:3,4) = -T_ee0(1:3,1:3)*T_0ee(1:3,4);
        
        T_0i = mdl.joints(i).frame_in.t;
        T_i0 = T_0i;
        T_i0(4,1:3) = 0;
        T_i0(1:3,4) = -T_i0(1:3,1:3)*T_0i(1:3,4);
        r_i_ee = T_ee0*T_0i*[0 0 0 1]';
        T_eei = T_ee0*T_0i;
        e_i_ee = T_eei(1:3,1:3)*[0 0 1]';
        
        dr_dtheta_num(:,i,c) = (r_i_ee(1:3)-r_i_ee_init(1:3))./eps;
        de_dtheta_num(:,i,c) = (e_i_ee(1:3)-e_i_ee_init(1:3))./eps;
    end
    %Reset
    mdl.position = pos_init;
    mdl.forwardPosition();
end


%% Compare the full thing

%We set the acceleration to something random
mdl.acceleration(:) = 5;
dxdd_dq_num = zeros(numel(sens.measurement),numel(mdl.joints));
mdl.g = [0 0 0]';
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mes_init = sens.measurement;
pos_init = mdl.position;
eps = 1e-5;
for i=1:numel(mdl.joints)
    pos = pos_init;
    pos(i) = pos(i) + eps;
    mdl.position = pos;
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    dxdd_dq_num(:,i) = (sens.measurement-mes_init)./eps;
end

mdl.position = pos_init;
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();

%Now analytical
dxdd_dq = zeros(numel(sens.measurement),numel(mdl.joints));
for i=1:numel(mdl.joints)
    H_tmp = zeros(3,numel(mdl.joints));
    for c=1:numel(mdl.joints)
        if(strcmp(mdl.joints(c).type,'revolute'))
            H_tmp(:,c) = (cross(de_dtheta_num(:,c,i),ris_ee(:,c)) + cross(eis_ee(:,c),-dr_dtheta_num(:,c,i)));
            disp(['Htmp: ' num2str(i-1) ' col: ' num2str(c-1)]);
            disp('First Cross: ');
            disp(cross(de_dtheta(:,c,i),ris_ee(:,c)));
            disp('Second Cross: ');
            disp(cross(eis_ee(:,c),dr_dtheta(:,c,i)));
        else
            H_tmp(:,c) = de_dtheta_num(:,c,i);
        end
    end
    dxdd_dq(:,i) = H_tmp*mdl.acceleration;
end

disp('Derivative of acceleration wrt joint position Jacobian Differences');
mdl.calculateSensorJacobians();
disp(norm(dxdd_dq_num-dxdd_dq));

%% Lets take a look at J_dot, how do we get it numerically ??
mdl.position = pos_init;
mdl.acceleration(:) = 0;
mdl.velocity(1) = 5;
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();

T_0ee_init = sens.transform;
T_ee0_init = T_0ee_init;
T_ee0_init(1:3,1:3) = T_ee0_init(1:3,1:3)';
T_ee0_init(1:3,4) = -T_ee0_init(1:3,1:3)*T_ee0_init(1:3,4);

%Figure out ris_b and eis s
[eis_b_init, ris_b_init] = computeEiRi(mdl,'post:rradius');
eis_ee_init = sens.transform(1:3,1:3)'*eis_b_init;
ris_ee_init = sens.transform(1:3,1:3)'*ris_b_init;

[J_init_ee, ~] =  computeJacobiansEE(mdl,'post:rradius');
[J_init_b, ~] =  computeJacobians(mdl,'post:rradius');

dJ_dq_ee_num = zeros(size(J_init_ee,1),size(J_init_ee,2),size(J_init_ee,2));
dJ_dq_b_num = zeros(size(J_init_ee,1),size(J_init_ee,2),size(J_init_ee,2));
dris_b_dq_num = zeros(size(ris_b,1),size(ris_b,2),size(ris_b,2));
deis_b_dq_num = zeros(size(eis_b,1),size(eis_b,2),size(eis_b,2));
dris_ee_dq_num = zeros(size(ris_b,1),size(ris_b,2),size(ris_b,2));
deis_ee_dq_num = zeros(size(eis_b,1),size(eis_b,2),size(eis_b,2));

eps = 1e-8;
for i=1:numel(mdl.joints)
    pos = pos_init;
    pos(i) = pos(i) + eps;
    mdl.position = pos;
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration
    
    T_0ee = sens.transform;
    T_ee0 = T_0ee;
    T_ee0(1:3,1:3) = T_ee0(1:3,1:3)';
    T_ee0(1:3,4) = -T_ee0(1:3,1:3)*T_ee0(1:3,4);

    
    [J_eps_ee, ~] = computeJacobiansEE(mdl,'post:rradius');
    [J_eps_b, ~] = computeJacobians(mdl,'post:rradius');
    dJ_dq_ee_num(1:3,:,i) = T_0ee_init(1:3,1:3)'*(J_eps_b(1:3,:)-J_init_b(1:3,:))./eps;
    dJ_dq_ee_num(4:6,:,i) = (J_eps_ee(4:6,:)-J_init_ee(4:6,:))./eps;
    dJ_dq_b_num(:,:,i) = (J_eps_b-J_init_b)./eps;
    [eis_b, ris_b] = computeEiRi(mdl,'post:rradius');
    eis_ee = T_0ee_init(1:3,1:3)'*eis_b;
    ris_ee = T_0ee_init(1:3,1:3)'*ris_b;
    dris_b_dq_num(:,:,i) = (ris_b-ris_b_init)./eps;
    deis_b_dq_num(:,:,i) = (eis_b-eis_b_init)/eps;
    dris_ee_dq_num(:,:,i) = (ris_ee-ris_ee_init)./eps;
    deis_ee_dq_num(:,:,i) = (eis_ee-eis_ee_init)/eps;
end

%Notice dris_b_dq_num is computed as d[r1 r2 r3 r4]/dq1 d[r1 r2 r3 r4]/dq2
%We want it in dr1/d[q1 q2 q3 q4]... to get r1_dot so permute
dris_b_dq_num = permute(dris_b_dq_num,[1,3,2]);
deis_b_dq_num = permute(deis_b_dq_num,[1,3,2]);
dris_ee_dq_num = permute(dris_ee_dq_num,[1,3,2]);
deis_ee_dq_num = permute(deis_ee_dq_num,[1,3,2]);

mdl.position = pos_init;
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();

J_dot_ee_num = zeros(size(J_init_ee));
ris_dot_ee_num = ris_ee;
eis_dot_ee_num = eis_ee;
ris_dot_b_num = ris_ee;
eis_dot_b_num = eis_ee;
eis_dot_ee = eis_ee;
ris_dot_ee = ris_ee;
for i=1:numel(mdl.joints)
    J_dot_ee_num(:,i) = dJ_dq_ee_num(:,:,i)*mdl.velocity;
    ris_dot_b_num(:,i) = dris_b_dq_num(:,:,i)*mdl.velocity;
    eis_dot_b_num(:,i) = deis_b_dq_num(:,:,i)*mdl.velocity;
    ris_dot_ee_num(:,i) = dris_ee_dq_num(:,:,i)*mdl.velocity;
    eis_dot_ee_num(:,i) = deis_ee_dq_num(:,:,i)*mdl.velocity;
    eis_dot_ee(:,i) = de_dtheta(:,:,i)*mdl.velocity;
    ris_dot_ee(:,i) = dr_dtheta(:,:,i)*mdl.velocity;
end

[J, J_dot] = computeJacobians(mdl,'post:rradius');
[J_ee, J_dot_ee, ei_dots, ri_dots] = computeJacobiansEE(mdl,'post:rradius');


%Lets take a look if dri_dot_dtheta can be expressed as dJ_dtheta*q_dot -
%dReej_dtheta*joint.v

%Here is the numeric dri_dot_dtheta
dris_dot_dtheta_num = zeros(size(ri_dots,1),size(ri_dots,2),size(ri_dots,2));
deis_dot_dtheta_num = zeros(size(ri_dots,1),size(ri_dots,2),size(ri_dots,2));

djv_dot_dtheta_num = zeros(size(ri_dots,1),size(ri_dots,2),size(ri_dots,2));

T_0ee = sens.transform;
T_ee0 = T_0ee;
T_ee0(1:3,1:3) = T_ee0(1:3,1:3)';
T_ee0(1:3,4) = -T_ee0(1:3,1:3)*T_ee0(1:3,4);

%Joint frame in velocity in end effector frame
joint_frames_in_v_init = zeros(3,numel(mdl.joints));
for j=1:numel(mdl.joints)
    jv_ee =  T_ee0(1:3,1:3)'*mdl.joints(j).frame_in.t(1:3,1:3)*mdl.joints(j).frame_in.v(4:6);
    joint_frames_in_v_init(:,j) = jv_ee;
end

for i=1:numel(mdl.joints)
    pos = pos_init;
    pos(i) = pos(i) + eps;
    mdl.position = pos;
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration
    [J_ee, J_dot_ee, ei_dots_eps, ri_dots_eps] = computeJacobiansEE(mdl,'post:rradius');
    dris_dot_dtheta_num(:,:,i) = (ri_dots_eps-ri_dots)/eps;
    deis_dot_dtheta_num(:,:,i) = (ei_dots_eps-ei_dots)/eps;
    joint_frames_in_v = zeros(size(joint_frames_in_v_init));
    T_0ee = sens.transform;
    T_ee0 = T_0ee;
    T_ee0(1:3,1:3) = T_ee0(1:3,1:3)';
    T_ee0(1:3,4) = -T_ee0(1:3,1:3)*T_ee0(1:3,4);
    for j=1:numel(mdl.joints)
       jv_ee =  T_ee0(1:3,1:3)'*mdl.joints(j).frame_in.t(1:3,1:3)*mdl.joints(j).frame_in.v(4:6);
       joint_frames_in_v(:,j) = jv_ee;
    end
    djv_dot_dtheta_num(:,:,i) = (joint_frames_in_v- joint_frames_in_v_init)/eps;
end


dris_dot_dtheta_num = permute(dris_dot_dtheta_num,[1,3,2]);
deis_dot_dtheta_num = permute(deis_dot_dtheta_num,[1,3,2]);
%Here is the analytic one
dris_dot_dtheta = zeros(size(dris_dot_dtheta_num));
deis_dot_dtheta = zeros(size(deis_dot_dtheta_num));


T_0ee = sens.transform;
T_ee0 = T_0ee;
T_ee0(1:3,1:3) = T_ee0(1:3,1:3)';
T_ee0(1:3,4) = -T_ee0(1:3,1:3)*T_ee0(1:3,4);

for i=1:numel(mdl.joints)
    dJ_dthetai_mul_q_dot = dJ_dq_ee_num(1:3,:,i)*mdl.velocity;
    
    for c=i+1:numel(mdl.joints)
        T_0c = mdl.joints(c).frame_in.t;
        T_c0 = T_0c';
        T_c0(4,1:3) = 0;
        T_c0(1:3,4) = -T_c0(1:3,1:3)*T_0c(1:3,4);
        T_eec = T_ee0*T_0c;
        T_cee = T_c0*T_0ee;
        R_eec = T_eec(1:3,1:3);
        R_c0 = T_c0(1:3,1:3);
        
        T_ci = T_c0*T_0i;
        T_ic = T_i0*T_0c;
        R_ci = T_ci(1:3,1:3);
        dris_dot_dtheta(:,c,i) = dJ_dthetai_mul_q_dot - R_eec*skew([0 0 1])*R_ci*mdl.joints(i).frame_in.v(4:6);
    end
end

%So now if we have qd = 0 and qdd != 0 sens.mes = J*qdd + J_dot*qd = J*qdd
%Thus d_sens.mes_dtheta = dJ_dq_ee_num*qdd

dxdd_dq = zeros(numel(sens.measurement),numel(mdl.joints));
mdl.position = pos_init;
mdl.velocity(:) = 5;
mdl.acceleration(:) = 0;
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mes_init = sens.measurement;
eps = 1e-8;
for i=1:numel(mdl.joints)
    pos = pos_init;
    pos(i) = pos(i) + eps;
    mdl.position = pos;
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    dxdd_dq(:,i) = (sens.measurement-mes_init)./eps;
end

%% Now for derivative of acceleration wrt q when velocity isnt zero but acceleration is
mdl.position = pos_init;
mdl.acceleration(:) = 0;
mdl.velocity = vel_init;
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mes_init = sens.measurement;

[~, Jv_dot_init, eis_dot_init] = computeJacobiansEE(mdl,'post:rradius');
Jv_dot_init = Jv_dot_init(1:end/2,:);
dJvdot_dq_num = zeros(size(Jv_dot_init,1),size(Jv_dot_init,2),size(Jv_dot_init,2));
deis_dot_dq_num = zeros(size(eis_dot_init,1),size(eis_dot_init,2),size(eis_dot_init,2));
eps = 1e-8;
for i=1:numel(mdl.joints)
    pos = pos_init;
    pos(i) = pos(i) + eps;
    mdl.position = pos;
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    dxdd_dq_num(:,i) = (sens.measurement-mes_init)./eps;
    [~, Jv_dot, eis_dot, ~] = computeJacobiansEE(mdl,'post:rradius');
    Jv_dot = Jv_dot(1:end/2,:);
    dJvdot_dq_num(:,:,i) = (Jv_dot-Jv_dot_init)./eps;
    deis_dot_dq_num(:,:,i) = (eis_dot- eis_dot_init)./eps;
end

mdl.position = pos_init;
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();

dxdd_dq_num_semi = dxdd_dq_num;
for i=1:numel(mdl.joints)
    dxdd_dq_num_semi(:,i) = dJvdot_dq_num(:,:,i)*mdl.velocity;
end





%% Lets look at derivative of ei dot wrt to q in end effector frame

%First lets try to get ri_dot and ei_dot in a numerical way
mdl.position = pos_init;
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();

%Compute eis and ris
T_0ee = sens.transform;
T_ee0 = T_0ee';
T_ee0(4,1:3) = 0;
T_ee0(1:3,4) = -T_ee0(1:3,1:3)*T_0ee(1:3,4);
eis_ee_init = eis_ee;
ris_ee_init = ris_ee;
for i=1:numel(mdl.joints)
    
    T_0i = mdl.joints(i).frame_in.t;
    T_i0 = T_0i';
    T_i0(4,1:3) = 0;
    T_i0(1:3,4) = -T_i0(1:3,1:3)*T_0i(1:3,4);
    
    T_eei = T_ee0*T_0i;
    T_iee = T_i0*T_0ee;
    
    ri = -T_eei(1:3,4);
    ris_ee_init(:,i) = ri;
    
    ei = T_eei(1:3,3);
    eis_ee_init(:,i) = ei;
end

%Use dt to update position from velocity
dt = 0.001;
mdl.position = mdl.position + mdl.velocity*dt;
mdl.forwardPosition();

T_0ee = sens.transform;
T_ee0 = T_0ee';
T_ee0(4,1:3) = 0;
T_ee0(1:3,4) = -T_ee0(1:3,1:3)*T_0ee(1:3,4);

%Recompute eis and ris
for i=1:numel(mdl.joints)
    
    T_0ee = sens.transform;
    T_ee0 = T_0ee';
    T_ee0(4,1:3) = 0;
    T_ee0(1:3,4) = -T_ee0(1:3,1:3)*T_0ee(1:3,4);
    
    
    T_0i = mdl.joints(i).frame_in.t;
    T_i0 = T_0i';
    T_i0(4,1:3) = 0;
    T_i0(1:3,4) = -T_i0(1:3,1:3)*T_0i(1:3,4);
    
    T_eei = T_ee0*T_0i;
    T_iee = T_i0*T_0ee;
    
    ri = -T_eei(1:3,4);
    ris_ee(:,i) = ri;
    
    ei = T_eei(1:3,3);
    eis_ee(:,i) = ei;
end

ris_dot_ee_num = (ris_ee-ris_ee_init)/dt;
eis_dot_ee_num = (eis_ee-eis_ee_init)/dt;

%Reset Model
mdl.position = pos_init;
mdl.forwardPosition();

%Numerical
d_ei_dot_d_theta_num = zeros(3,numel(mdl.joints),numel(mdl.joints));
[J_init, dJ_init] = computeJacobiansEE(mdl,'post:rradius');

eis_dot = zeros(3,numel(mdl.joints));
for i=1:numel(mdl.joints)
    %We pull out ei from J
    ei_dot_init = cross(J_init(4:6,1:i)*mdl.velocity(1:i),J_init(4:6,i));
    eis_dot(:,i) = ei_dot_init;
    for c=1:numel(mdl.joints)
        pos = pos_init;
        pos(c) = pos(c) + eps;
        mdl.position = pos;
        mdl.forwardPosition();
        [J, dJ] = computeJacobiansEE(mdl,'post:rradius');
        ei_dot = cross(J(4:6,1:i)*mdl.velocity(1:i),J(4:6,i));
        d_ei_dot_d_theta_num(:,i,c) = (ei_dot - ei_dot_init)./eps;
    end
    %Reset
    mdl.position = pos_init;
    mdl.forwardPosition();
end


%Analytical
d_ei_dot_d_theta = d_ei_dot_d_theta_num;
for i=2:numel(mdl.joints)
    %Build d_ei_dot_d_theta_c
    for c=1:numel(mdl.joints)
        d_ei_dot_d_theta_c = cross(de_dtheta(:,1:i-1,c)*mdl.velocity(1:i-1),eis_ee(:,i)) + ...
            cross(eis_ee(:,1:i-1)*mdl.velocity(1:i-1),de_dtheta(:,i,c));
        d_ei_dot_d_theta(:,i,c) = d_ei_dot_d_theta_c;
    end
end

%% Lets look at derivative of ri dot wrt to q in end effector frame

%Numerical
d_ri_dot_d_theta_num = zeros(3,numel(mdl.joints),numel(mdl.joints));
ris_dot = eis_dot;
[J_init, dJ_init] = computeJacobiansEE(mdl,'post:rradius');
for i=1:numel(mdl.joints)
    %We pull out ei from J
    ri_dot_init = cross(J_init(4:6,1:i)*mdl.velocity(1:i),ris_ee(:,i));
    for j = i+1:numel(mdl.joints)
        ri_dot_init = ri_dot_init + cross(J_init(4:6,j)*mdl.velocity(j),ris_ee(:,j));
    end
    ris_dot(:,i) = ri_dot_init;
    for c=1:numel(mdl.joints)
        pos = pos_init;
        pos(c) = pos(c) + eps;
        mdl.position = pos;
        mdl.forwardPosition();
        
        T_0ee = sens.transform;
        T_ee0 = T_0ee;
        T_ee0(1:3,1:3) = T_ee0(1:3,1:3)';
        T_ee0(1:3,4) = -T_ee0(1:3,1:3)*T_ee0(1:3,4);
        T_0i = mdl.joints(i).frame_in.t;
        T_i0 = T_0i;
        T_i0(1:3,1:3) = T_0i(1:3,1:3)';
        T_i0(1:3,4) = -T_0i(1:3,1:3)*T_0i(1:3,4);
        T_eei = T_ee0*T_0i;
        ri = T_eei(1:3,4);
        [J, dJ] = computeJacobiansEE(mdl,'post:rradius');
        ri_dot = cross(J(4:6,1:i)*mdl.velocity(1:i),ri);
        for j = i+1:numel(mdl.joints)
            T_0j = mdl.joints(j).frame_in.t;
            T_j0 = T_0i;
            T_j0(1:3,1:3) = T_0j(1:3,1:3)';
            T_j0(1:3,4) = -T_0j(1:3,1:3)*T_0j(1:3,4);
            T_eej = T_ee0*T_0j;
            rj = T_eej(1:3,4);
            ri_dot = ri_dot + cross(J(4:6,j)*mdl.velocity(j),rj);
        end
        d_ri_dot_d_theta_num(:,i,c) = (ri_dot - ri_dot_init)./eps;
    end
    %Reset
    mdl.position = pos_init;
    mdl.forwardPosition();
end

%Analytical
d_ri_dot_d_theta = d_ri_dot_d_theta_num;
[J, dJ] = computeJacobiansEE(mdl,'post:rradius');
for i=1:numel(mdl.joints)
    for l=1:numel(mdl.joints)
        %First the (d_ej_d_q * q_dot) x ri part
        %We need the sum 1:i of d_ej_d_q * q_dot then cross it with ri
        d_ri_dot_dl = cross(de_dtheta(:,1:i,l)*mdl.velocity(1:i),ris_ee(:,i));
        %Next the w cross derivative of ri wrt theta
        d_ri_dot_dl = d_ri_dot_dl + cross(J(4:6,1:i)*mdl.velocity(1:i),dr_dtheta(:,i,l));
        
        for j=i+1:numel(mdl.joints)
            d_ri_dot_dl = d_ri_dot_dl + cross(de_dtheta(:,j,l)*mdl.velocity(j),ris_ee(:,j));
            d_ri_dot_dl = d_ri_dot_dl + cross(eis_ee(:,j)*mdl.velocity(j),dr_dtheta(:,j,l));
        end
        d_ri_dot_d_theta(:,i,l) = d_ri_dot_dl;
    end
end

%% Now put everything together into one nice dJ_dot_dq
dJvdot_dq = dJvdot_dq_num;
for c = 1:numel(mdl.joints)
    for i=1:numel(mdl.joints)
        dJvdot_dq(:,i,c) = cross(d_ei_dot_d_theta(:,i,c),ris_ee(:,i)) + ...
            cross(eis_dot(:,i),dr_dtheta(:,i,c)) + ...
            cross(de_dtheta(:,i,c),ris_dot(:,i)) + ...
            cross(eis_ee(:,i),d_ri_dot_d_theta(:,i,c));
    end
end


%% Lets try the full comparisson with velocity and acceleration

mdl.position = rand(numel(mdl.joints),1);
mdl.velocity = rand(numel(mdl.joints),1);
mdl.acceleration = rand(numel(mdl.joints),1);
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();

%Compute Numerical Jacobian of accelerometer WRT q
mes_init = sens.measurement;
pos_init = mdl.position;
eps = 1e-4;
for i=1:numel(mdl.joints)
    pos = pos_init;
    pos(i) = pos(i) + eps;
    mdl.position = pos;
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    dxdd_dq_num(:,i) = (sens.measurement-mes_init)./eps;
end
%Reset Model
mdl.position = pos_init;
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();

%Now our analytical version
T_0ee = sens.transform;
T_ee0 = T_0ee;
T_ee0(1:3,1:3) = T_ee0(1:3,1:3)';
T_ee0(1:3,4) = -T_ee0(1:3,1:3)*T_ee0(1:3,4);

%Derivative of ri in end effector frame wrt joint positions
dr_dtheta = zeros(3,numel(mdl.joints),numel(mdl.joints));
%Derivative of ei in end effector frame wrt joint positions
de_dtheta = zeros(3,numel(mdl.joints),numel(mdl.joints));
%all of the ri in end effector frame
ris_ee = zeros(3,numel(mdl.joints));
ris_dot = zeros(3,numel(mdl.joints));
%all of the ei in end effector frame
eis_ee = zeros(3,numel(mdl.joints));
eis_dot = zeros(3,numel(mdl.joints));
for i=1:numel(mdl.joints)
    
    T_0i = mdl.joints(i).frame_in.t;
    T_i0 = T_0i;
    T_i0(1:3,1:3) = T_0i(1:3,1:3)';
    T_i0(1:3,4) = -T_0i(1:3,1:3)*T_0i(1:3,4);
    
    T_eei = T_ee0*T_0i;
    
    ri = T_eei(1:3,4);
    ris_ee(:,i) = ri;
    
    ei = T_eei(1:3,3);
    eis_ee(:,i) = ei;
    
    %The Jacobian Jv is based on ri and ei. Since accelerometer is in end
    %effector only the joints after i change ri and ei when moved
    for c=i+1:numel(mdl.joints)
        T_0c = mdl.joints(c).frame_in.t;
        T_c0 = T_0c;
        T_c0(1:3,1:3) = T_0c(1:3,1:3)';
        T_c0(1:3,4) = -T_c0(1:3,1:3)*T_c0(1:3,4);
        T_eec = T_ee0*T_0c;
        R_eec = T_eec(1:3,1:3);
        
        T_ci = T_c0*T_0i;
        T_ic = T_i0*T_0c;
        R_ci = T_ci(1:3,1:3);
        %From joint i to joint c in frame i
        r_ci_i = T_ic(1:3,4);
        dri_dc = R_eec*skew([0 0 1])*R_ci*r_ci_i;
        dr_dtheta(:,i,c) = dri_dc;
        
        dei_dc = R_eec*skew([0 0 -1])*R_ci*[0 0 1]';
        de_dtheta(:,i,c) = dei_dc;
    end
end
%Now compute derivative of Jacobian wrt q
d_Jqdd_d_theta = zeros(3,numel(mdl.joints));
for i=1:numel(mdl.joints)
    H_tmp = zeros(3,numel(mdl.joints));
    for c=1:numel(mdl.joints)
        H_tmp(:,c) = (cross(de_dtheta(:,c,i),ris_ee(:,c)) + cross(eis_ee(:,c),dr_dtheta(:,c,i)));
    end
    d_Jqdd_d_theta(:,i) = H_tmp*mdl.acceleration;
end

[J, dJ] = computeJacobiansEE(mdl,'post:rradius');
for i=1:numel(mdl.joints)
    %We pull out ei from J
    ei_dot = cross(J(4:6,1:i)*mdl.velocity(1:i),J(4:6,i));
    ri_dot = cross(J(4:6,1:i)*mdl.velocity(1:i),ris_ee(:,i));
    for j = i+1:numel(mdl.joints)
        ri_dot = ri_dot + cross(J(4:6,j)*mdl.velocity(j),ris_ee(:,j));
    end
    ris_dot(:,i) = ri_dot;
    eis_dot(:,i) = ei_dot;
end

%We need derivative of ei_dot wrt theta
d_ei_dot_d_theta = d_ei_dot_d_theta_num;
for i=1:numel(mdl.joints)
    %Build d_ei_dot_d_theta_c
    for c=1:numel(mdl.joints)
        d_ei_dot_d_theta_c = cross(de_dtheta(:,1:i-1,c)*mdl.velocity(1:i-1),eis_ee(:,i)) + ...
            cross(eis_ee(:,1:i-1)*mdl.velocity(1:i-1),de_dtheta(:,i,c));
        d_ei_dot_d_theta(:,i,c) = d_ei_dot_d_theta_c;
    end
end

%We need derivative of ri_dot wrt theta
d_ri_dot_d_theta = d_ri_dot_d_theta_num;
for i=1:numel(mdl.joints)
    for l=1:numel(mdl.joints)
        %First the (d_ej_d_q * q_dot) x ri part
        %We need the sum 1:i of d_ej_d_q * q_dot then cross it with ri
        d_ri_dot_dl = cross(de_dtheta(:,1:i,l)*mdl.velocity(1:i),ris_ee(:,i));
        %Next the w cross derivative of ri wrt theta
        d_ri_dot_dl = d_ri_dot_dl + cross(J(4:6,1:i)*mdl.velocity(1:i),dr_dtheta(:,i,l));
        
        for j=i+1:numel(mdl.joints)
            d_ri_dot_dl = d_ri_dot_dl + cross(de_dtheta(:,j,l)*mdl.velocity(j),ris_ee(:,j));
            d_ri_dot_dl = d_ri_dot_dl + cross(eis_ee(:,j)*mdl.velocity(j),dr_dtheta(:,j,l));
        end
        d_ri_dot_d_theta(:,i,l) = d_ri_dot_dl;
    end
end

%Now for the derivative of Jv_dot wrt theta multiplied by q dot
dJvdotqd_dtheta = d_Jqdd_d_theta;
dJvdot_dq = dJvdot_dq_num;
for c = 1:numel(mdl.joints)
    for i=1:numel(mdl.joints)
        
        dJvdot_dq_col = cross(d_ei_dot_d_theta(:,i,c),ris_ee(:,i)) + ...
            cross(eis_dot(:,i),dr_dtheta(:,i,c)) + ...
            cross(de_dtheta(:,i,c),ris_dot(:,i)) + ...
            cross(eis_ee(:,i),d_ri_dot_d_theta(:,i,c));
        
        dJvdot_dq(:,i,c) = dJvdot_dq_col;
    end
    dJvdotqd_dtheta(:,c) = dJvdot_dq(:,:,c)*mdl.velocity;
end
%Finally the whole Jacobian
dxdd_dq = - dJvdotqd_dtheta - d_Jqdd_d_theta;

disp('Derivative of acceleration wrt joint position Jacobian Differences');
disp(norm(dxdd_dq-dxdd_dq_num));

%% Lets see if we got everything right

%We set the acceleration to something random
mdl.position = pos_init;
mdl.velocity = vel_init;
%mdl.acceleration = rand(numel(mdl.joints),1);
mdl.acceleration = accel_init;
J_num = zeros(numel(sens.measurement),numel(mdl.joints)*3);
mdl.g = [0 0 9.81]';
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mes_init = sens.measurement;
eps = 1e-5;
for i=1:numel(mdl.joints)
    pos = pos_init;
    pos(i) = pos(i) + eps;
    mdl.position = pos;
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    J_num(:,i) = (sens.measurement-mes_init)./eps;
end
mdl.position = pos_init;
eps = 1e-5;
for i=1:numel(mdl.joints)
    vel = vel_init;
    vel(i) = vel(i) + eps;
    mdl.velocity = vel;
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    J_num(:,i+numel(mdl.joints)) = (sens.measurement-mes_init)./eps;
end

mdl.velocity = vel_init;
eps = 1e-5;
for i=1:numel(mdl.joints)
    acc = accel_init;
    acc(i) = acc(i) + eps;
    mdl.acceleration = acc;
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    J_num(:,i+numel(mdl.joints)*2) = (sens.measurement-mes_init)./eps;
end
mdl.acceleration = accel_init;
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mdl.calculateSensorJacobians();
disp('FINAL C++ TO NUMERIC DIFF');
disp(norm(J_num - sens.obsJacobian));



