function [xnext,df_dx,df_du]=Arm_Dyn(x,u,dt)


%dynamic paramter setting
m1=1; %mass of the 1st link
m2=1; %mass of the 2nd link
l1=1; %length of the 1st link
l2=1; %length of the 2nd link
lc1=0.5; %CoM of the 1st link
lc2=0.5; %CoM of the 2nd link
I1=0.5; %moment of initia of the 1st link
I2=0.5; %moment of initia of the 2nd link
g=9.8; %gravity costant

x1=x(1);
x2=x(2);
x3=x(3);
x4=x(4);
u1=u(1);
u2=u(2);

x1next=x1+dt*x2;
x3next=x3+dt*x4;

M11=m1*lc1^2+m2*(l1*l1+lc2*lc2+2*l1*lc2*cos(x3))+I1+I2;
M12=m2*(lc2*lc2+l1*lc2*cos(x3))+I2;
M21=M12;
M22=(m2*lc2*lc2+I2);

C1=-m2*l1*lc2*sin(x3)*x4*x4-2*m2*l1*lc2*sin(x3)*x2*x4;
C2=m2*l1*lc2*sin(x3)*x2*x2;

G1=(m1*lc1+m2*l1)*g*cos(x1)+m2*g*lc2*cos(x1+x3);
G2=m2*g*lc2*cos(x1+x3);

x2next=x2+dt*(M22*(u1-C1-G1)-M12*(u2-C2-G2))/(M11*M22-M12*M21);
x4next=x4+dt*(-M21*(u1-C1-G1)+M11*(u2-C2-G2))/(M11*M22-M12*M21);

xnext=[x1next,x2next,x3next,x4next]';

%compute df/dx
df1_dx=[1, dt, 0, 0];
df2_dx=[ (dt*((m2*lc2^2 + I2)*(g*sin(x1)*(l1*m2 + lc1*m1) + g*lc2*m2*sin(x1 + x3)) - g*lc2*m2*sin(x1 + x3)*(I2 + m2*(lc2^2 + l1*cos(x3)*lc2))))/((m2*lc2^2 + I2)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2), (dt*(2*l1*lc2*m2*x4*sin(x3)*(m2*lc2^2 + I2) + 2*l1*lc2*m2*x2*sin(x3)*(I2 + m2*(lc2^2 + l1*cos(x3)*lc2))))/((m2*lc2^2 + I2)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2) + 1, - (dt*((- l1*lc2*m2*cos(x3)*x2^2 + g*lc2*m2*sin(x1 + x3))*(I2 + m2*(lc2^2 + l1*cos(x3)*lc2)) - (m2*lc2^2 + I2)*(l1*lc2*m2*cos(x3)*x4^2 + 2*l1*lc2*m2*x2*cos(x3)*x4 + g*lc2*m2*sin(x1 + x3)) + l1*lc2*m2*sin(x3)*(l1*lc2*m2*sin(x3)*x2^2 - u2 + g*lc2*m2*cos(x1 + x3))))/((m2*lc2^2 + I2)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2) - (dt*(2*l1*lc2*m2*sin(x3)*(I2 + m2*(lc2^2 + l1*cos(x3)*lc2)) - 2*l1*lc2*m2*sin(x3)*(m2*lc2^2 + I2))*((I2 + m2*(lc2^2 + l1*cos(x3)*lc2))*(l1*lc2*m2*sin(x3)*x2^2 - u2 + g*lc2*m2*cos(x1 + x3)) + (m2*lc2^2 + I2)*(l1*lc2*m2*sin(x3)*x4^2 + 2*l1*lc2*m2*x2*sin(x3)*x4 + u1 - g*cos(x1)*(l1*m2 + lc1*m1) - g*lc2*m2*cos(x1 + x3))))/((m2*lc2^2 + I2)*(m1*lc1^2 + I1 + I2 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2)^2, (dt*(m2*lc2^2 + I2)*(2*l1*lc2*m2*x2*sin(x3) + 2*l1*lc2*m2*x4*sin(x3)))/((m2*lc2^2 + I2)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2)];
df3_dx=[ 0, 0, 1, dt];
df4_dx=[ -(dt*((g*sin(x1)*(l1*m2 + lc1*m1) + g*lc2*m2*sin(x1 + x3))*(I2 + m2*(lc2^2 + l1*cos(x3)*lc2)) - g*lc2*m2*sin(x1 + x3)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2))))/((m2*lc2^2 + I2)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2), -(dt*(2*l1*lc2*m2*x2*sin(x3)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) + 2*l1*lc2*m2*x4*sin(x3)*(I2 + m2*(lc2^2 + l1*cos(x3)*lc2))))/((m2*lc2^2 + I2)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2), (dt*((- l1*lc2*m2*cos(x3)*x2^2 + g*lc2*m2*sin(x1 + x3))*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))*(l1*lc2*m2*cos(x3)*x4^2 + 2*l1*lc2*m2*x2*cos(x3)*x4 + g*lc2*m2*sin(x1 + x3)) + l1*lc2*m2*sin(x3)*(l1*lc2*m2*sin(x3)*x4^2 + 2*l1*lc2*m2*x2*sin(x3)*x4 + u1 - g*cos(x1)*(l1*m2 + lc1*m1) - g*lc2*m2*cos(x1 + x3)) + 2*l1*lc2*m2*sin(x3)*(l1*lc2*m2*sin(x3)*x2^2 - u2 + g*lc2*m2*cos(x1 + x3))))/((m2*lc2^2 + I2)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2) + (dt*(2*l1*lc2*m2*sin(x3)*(I2 + m2*(lc2^2 + l1*cos(x3)*lc2)) - 2*l1*lc2*m2*sin(x3)*(m2*lc2^2 + I2))*((I2 + m2*(lc2^2 + l1*cos(x3)*lc2))*(l1*lc2*m2*sin(x3)*x4^2 + 2*l1*lc2*m2*x2*sin(x3)*x4 + u1 - g*cos(x1)*(l1*m2 + lc1*m1) - g*lc2*m2*cos(x1 + x3)) + (l1*lc2*m2*sin(x3)*x2^2 - u2 + g*lc2*m2*cos(x1 + x3))*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2))))/((m2*lc2^2 + I2)*(m1*lc1^2 + I1 + I2 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2)^2, 1 - (dt*(2*l1*lc2*m2*x2*sin(x3) + 2*l1*lc2*m2*x4*sin(x3))*(I2 + m2*(lc2^2 + l1*cos(x3)*lc2)))/((m2*lc2^2 + I2)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2)];
df_dx=[df1_dx;df2_dx;df3_dx;df4_dx];
%compute df/dx
df1_du=[ 0, 0];
df2_du=[ (dt*(m2*lc2^2 + I2))/((m2*lc2^2 + I2)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2),...
          -(dt*(I2 + m2*(lc2^2 + l1*cos(x3)*lc2)))/((m2*lc2^2 + I2)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2)];      
df3_du=[ 0, 0];
df4_du=[ -(dt*(I2 + m2*(lc2^2 + l1*cos(x3)*lc2)))/((m2*lc2^2 + I2)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2),...
        (dt*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)))/((m2*lc2^2 + I2)*(I1 + I2 + lc1^2*m1 + m2*(l1^2 + 2*cos(x3)*l1*lc2 + lc2^2)) - (I2 + m2*(lc2^2 + l1*cos(x3)*lc2))^2)];
    
df_du=[df1_du;df2_du;df3_du;df4_du];







