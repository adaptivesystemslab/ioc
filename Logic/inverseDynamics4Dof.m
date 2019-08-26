% (********************************************)
% (** SYMORO+ : SYmbolic MOdelling of RObots **)
% (**========================================**)
% (**      IRCCyN-ECN - 1, rue de la Noe     **)
% (**      B.P.92101                         **)
% (**      44321 Nantes cedex 3, FRANCE      **)
% (**      www.irccyn.ec-nantes.fr           **)
% (********************************************)


%    Name of file : C:\Users\dell.dell-PC\Desktop\arm_ioc\v3\arm_4dofv3.dyn




%      Geometric parameters   

% j     ant   mu    sigma gamma b     alpha d     theta r
% 1     0     1     0     0     0     a1    0     t1    0
% 2     1     1     0     0     0     a2    0     t2    0
% 3     2     1     0     0     0     a3    0     t3    0
% 4     3     1     0     0     0     a4    d4    t4    r4
% 5     4     0     2     0     0     a5    d5    t5    r5



%              Inertial parameters

% j     XX    XY    XZ    YY    YZ    ZZ    MX    MY    MZ    M     Ia

% 1     0     0     0     0     0     0     0     0     0     0     0

% 2     0     0     0     0     0     0     0     0     0     0     0

% 3     0     0     0     0     0     0     0     0     0     0     0

% 4     XX4   XY4   XZ4   YY4   YZ4   ZZ4   MX4   MY4   MZ4   M4    0

% 5     XX5   XY5   XZ5   YY5   YZ5   ZZ5   MX5   MY5   MZ5   M5    0



%  External forces,friction parameters, joint velocities and accelerations

% j      FX     FY     FZ     CX     CY     CZ     FS     FV     QP     QDP

% 1      0      0      0      0      0      0      0      0      QP1    QDP1

% 2      0      0      0      0      0      0      0      0      QP2    QDP2

% 3      0      0      0      0      0      0      0      0      QP3    QDP3

% 4      0      0      0      0      0      0      0      0      QP4    QDP4

% 5      0      0      0      0      0      0      0      0      0      0

% Base velocity, base accelerations, and gravity

% j     W0    WP0   V0    VP0   G

% 1     0     0     0     0     0

% 2     0     0     0     0     0

% 3     0     0     0     0     G3

%  Dynamic model: Newton Euler method
% Equations:

% Declaration of the function
function torques = inverseDynamics4Dof(armModel, accelerations)

    % Arm state
     state = armModel.getState();

     % Arm current state
     [t1, t2, t3, t4] = feval(@(x) x{:}, num2cell(state(1,:)));
     [QP1, QP2, QP3, QP4] = feval(@(x) x{:}, num2cell(state(2,:)));
     [QDP1, QDP2, QDP3, QDP4] = feval(@(x) x{:}, num2cell(accelerations));

     joint1 = armModel.joints(1);
     joint2 = armModel.joints(2);
     joint3 = armModel.joints(3);
     joint4 = armModel.joints(4);
     joint5 = armModel.joints(5);

     % Arm MDH parameters
     a1 = joint1.twist;
     a2 = joint2.twist;
     a3 = joint3.twist;
     a4 = joint4.twist;
     a5 = joint5.twist;
     r3 = joint3.offset;
     r5 = joint5.offset;

     % First moment of inertia
     [MX3, MY3, MZ3] = feval(@(x) x{:}, num2cell((joint3.com * joint3.mass)));
     [MX5, MY5 ,MZ5] = feval(@(x) x{:}, num2cell((joint5.com * joint5.mass)));
     M3 = joint3.mass;
     M5 = joint5.mass;

     % Second moment of inertia, i.e., inertia tensor
     tensor3 = joint3.inertiaTensor;
     tensor5 = joint5.inertiaTensor;

     XX3 = tensor3(1,1);
     XY3 = tensor3(1,2);
     XZ3 = tensor3(1,3);
     YY3 = tensor3(2,2);
     YZ3 = tensor3(2,3);
     ZZ3 = tensor3(3,3);

     XX5 = tensor5(1,1);
     XY5 = tensor5(1,2);
     XZ5 = tensor5(1,3);
     YY5 = tensor5(2,2);
     YZ5 = tensor5(2,3);
     ZZ5 = tensor5(3,3);

     % Gravity
     [~, ~, G3] = feval(@(x) x{:}, num2cell(armModel.gravity));

% Function description:

	S1=sin(t1);
	C1=cos(t1);
	Sa1=sin(a1);
	Ca1=cos(a1);
	A311=S1.*Sa1;
	A321=C1.*Sa1;
	S2=sin(t2);
	C2=cos(t2);
	Sa2=sin(a2);
	Ca2=cos(a2);
	A212=Ca2.*S2;
	A222=C2.*Ca2;
	A312=S2.*Sa2;
	A322=C2.*Sa2;
	S3=sin(t3);
	C3=cos(t3);
	Sa3=sin(a3);
	Ca3=cos(a3);
	A213=Ca3.*S3;
	A223=C3.*Ca3;
	A313=S3.*Sa3;
	A323=C3.*Sa3;
	LOO23=-(r3.*Sa3);
	LOO33=Ca3.*r3;
	S4=sin(t4);
	C4=cos(t4);
	Sa4=sin(a4);
	Ca4=cos(a4);
	A214=Ca4.*S4;
	A224=C4.*Ca4;
	A314=S4.*Sa4;
	A324=C4.*Sa4;
	Sa5=sin(a5);
	Ca5=cos(a5);
	LOO25=-(r5.*Sa5);
	LOO35=Ca5.*r5;
	VP11=-(A311.*G3);
	VP21=-(A321.*G3);
	VP31=-(Ca1.*G3);
	WI12=A312.*QP1;
	WI22=A322.*QP1;
	WI32=Ca2.*QP1;
	W32=QP2 + WI32;
	WP12=A312.*QDP1 + QP2.*WI22;
	WP22=A322.*QDP1 - QP2.*WI12;
	WP32=Ca2.*QDP1 + QDP2;
	DV112=-WI12.^2;
	DV222=-WI22.^2;
	DV332=-W32.^2;
	DV122=WI12.*WI22;
	DV132=W32.*WI12;
	DV232=W32.*WI22;
	U122=DV122 - WP32;
	U132=DV132 + WP22;
	U222=DV112 + DV332;
	U232=DV232 - WP12;
	U322=DV232 + WP12;
	U332=DV112 + DV222;
	VP12=C2.*VP11 + A212.*VP21 + A312.*VP31;
	VP22=-(S2.*VP11) + A222.*VP21 + A322.*VP31;
	VP32=-(Sa2.*VP21) + Ca2.*VP31;
	WI13=A313.*W32 + C3.*WI12 + A213.*WI22;
	WI23=A323.*W32 - S3.*WI12 + A223.*WI22;
	WI33=Ca3.*W32 - Sa3.*WI22;
	W33=QP3 + WI33;
	WP13=QP3.*WI23 + C3.*WP12 + A213.*WP22 + A313.*WP32;
	WP23=-(QP3.*WI13) - S3.*WP12 + A223.*WP22 + A323.*WP32;
	WP33=QDP3 - Sa3.*WP22 + Ca3.*WP32;
	DV113=-WI13.^2;
	DV223=-WI23.^2;
	DV333=-W33.^2;
	DV123=WI13.*WI23;
	DV133=W33.*WI13;
	DV233=W33.*WI23;
	U113=DV223 + DV333;
	U123=DV123 - WP33;
	U133=DV133 + WP23;
	U213=DV123 + WP33;
	U223=DV113 + DV333;
	U233=DV233 - WP13;
	U313=DV133 - WP23;
	U323=DV233 + WP13;
	U333=DV113 + DV223;
	VSP13=LOO23.*U122 + LOO33.*U132 + VP12;
	VSP23=LOO23.*U222 + LOO33.*U232 + VP22;
	VSP33=LOO23.*U322 + LOO33.*U332 + VP32;
	VP13=C3.*VSP13 + A213.*VSP23 + A313.*VSP33;
	VP23=-(S3.*VSP13) + A223.*VSP23 + A323.*VSP33;
	VP33=-(Sa3.*VSP23) + Ca3.*VSP33;
	F13=MX3.*U113 + MY3.*U123 + MZ3.*U133 + M3.*VP13;
	F23=MX3.*U213 + MY3.*U223 + MZ3.*U233 + M3.*VP23;
	F33=MX3.*U313 + MY3.*U323 + MZ3.*U333 + M3.*VP33;
	PIS13=-YY3 + ZZ3;
	PIS23=XX3 - ZZ3;
	PIS33=-XX3 + YY3;
	No13=DV233.*PIS13 + WP13.*XX3 - U313.*XY3 + U213.*XZ3 + (-DV223 + DV333).*YZ3;
	No23=DV133.*PIS23 + U323.*XY3 + (DV113 - DV333).*XZ3 + WP23.*YY3 - U123.*YZ3;
	No33=DV123.*PIS33 + (-DV113 + DV223).*XY3 - U233.*XZ3 + U133.*YZ3 + WP33.*ZZ3;
	WI14=A314.*W33 + C4.*WI13 + A214.*WI23;
	WI24=A324.*W33 - S4.*WI13 + A224.*WI23;
	WI34=Ca4.*W33 - Sa4.*WI23;
	W34=QP4 + WI34;
	WP14=QP4.*WI24 + C4.*WP13 + A214.*WP23 + A314.*WP33;
	WP24=-(QP4.*WI14) - S4.*WP13 + A224.*WP23 + A324.*WP33;
	WP34=QDP4 - Sa4.*WP23 + Ca4.*WP33;
	DV114=-WI14.^2;
	DV224=-WI24.^2;
	DV334=-W34.^2;
	DV124=WI14.*WI24;
	DV134=W34.*WI14;
	DV234=W34.*WI24;
	U124=DV124 - WP34;
	U134=DV134 + WP24;
	U224=DV114 + DV334;
	U234=DV234 - WP14;
	U324=DV234 + WP14;
	U334=DV114 + DV224;
	VP14=C4.*VP13 + A214.*VP23 + A314.*VP33;
	VP24=-(S4.*VP13) + A224.*VP23 + A324.*VP33;
	VP34=-(Sa4.*VP23) + Ca4.*VP33;
	WI25=Sa5.*W34 + Ca5.*WI24;
	WI35=Ca5.*W34 - Sa5.*WI24;
	WP25=Ca5.*WP24 + Sa5.*WP34;
	WP35=-(Sa5.*WP24) + Ca5.*WP34;
	DV115=-WI14.^2;
	DV225=-WI25.^2;
	DV335=-WI35.^2;
	DV125=WI14.*WI25;
	DV135=WI14.*WI35;
	DV235=WI25.*WI35;
	U115=DV225 + DV335;
	U125=DV125 - WP35;
	U135=DV135 + WP25;
	U215=DV125 + WP35;
	U225=DV115 + DV335;
	U235=DV235 - WP14;
	U315=DV135 - WP25;
	U325=DV235 + WP14;
	U335=DV115 + DV225;
	VSP15=LOO25.*U124 + LOO35.*U134 + VP14;
	VSP25=LOO25.*U224 + LOO35.*U234 + VP24;
	VSP35=LOO25.*U324 + LOO35.*U334 + VP34;
	VP25=Ca5.*VSP25 + Sa5.*VSP35;
	VP35=-(Sa5.*VSP25) + Ca5.*VSP35;
	F15=MX5.*U115 + MY5.*U125 + MZ5.*U135 + M5.*VSP15;
	F25=MX5.*U215 + MY5.*U225 + MZ5.*U235 + M5.*VP25;
	F35=MX5.*U315 + MY5.*U325 + MZ5.*U335 + M5.*VP35;
	PIS15=-YY5 + ZZ5;
	PIS25=XX5 - ZZ5;
	PIS35=-XX5 + YY5;
	No15=DV235.*PIS15 + WP14.*XX5 - U315.*XY5 + U215.*XZ5 + (-DV225 + DV335).*YZ5;
	No25=DV135.*PIS25 + U325.*XY5 + (DV115 - DV335).*XZ5 + WP25.*YY5 - U125.*YZ5;
	No35=DV125.*PIS35 + (-DV115 + DV225).*XY5 - U235.*XZ5 + U135.*YZ5 + WP35.*ZZ5;
	N15=No15 - MZ5.*VP25 + MY5.*VP35;
	N25=No25 - MX5.*VP35 + MZ5.*VSP15;
	N35=No35 + MX5.*VP25 - MY5.*VSP15;
	FDI25=Ca5.*F25 - F35.*Sa5;
	FDI35=Ca5.*F35 + F25.*Sa5;
	N14=FDI35.*LOO25 - FDI25.*LOO35 + N15;
	N24=F15.*LOO35 + Ca5.*N25 - N35.*Sa5;
	N34=-(F15.*LOO25) + Ca5.*N35 + N25.*Sa5;
	FDI14=C4.*F15 - FDI25.*S4;
	FDI24=A214.*F15 + A224.*FDI25 - FDI35.*Sa4;
	FDI34=A314.*F15 + A324.*FDI25 + Ca4.*FDI35;
	E13=F13 + FDI14;
	E23=F23 + FDI24;
	E33=F33 + FDI34;
	N13=C4.*N14 + No13 - N24.*S4 - MZ3.*VP23 + MY3.*VP33;
	N23=A214.*N14 + A224.*N24 + No23 - N34.*Sa4 + MZ3.*VP13 - MX3.*VP33;
	N33=A314.*N14 + A324.*N24 + Ca4.*N34 + No33 - MY3.*VP13 + MX3.*VP23;
	FDI13=C3.*E13 - E23.*S3;
	FDI23=A213.*E13 + A223.*E23 - E33.*Sa3;
	FDI33=A313.*E13 + A323.*E23 + Ca3.*E33;
	N12=FDI33.*LOO23 - FDI23.*LOO33 + C3.*N13 - N23.*S3;
	N22=FDI13.*LOO33 + A213.*N13 + A223.*N23 - N33.*Sa3;
	N32=-(FDI13.*LOO23) + A313.*N13 + A323.*N23 + Ca3.*N33;
	N31=A312.*N12 + A322.*N22 + Ca2.*N32;
	GAM1=N31;
	GAM2=N32;
	GAM3=N33;
	GAM4=N34;

    torques(1) = GAM1; torques(2) = GAM2; torques(3) = GAM3; torques(4) = GAM4;
    
    
% *=*
% Number of operations : 218 '+' or '-', 274 '*' or '/'
