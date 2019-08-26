% (********************************************)
% (** SYMORO+ : SYmbolic MOdelling of RObots **)
% (**========================================**)
% (**      IRCCyN-ECN - 1, rue de la Noe     **)
% (**      B.P.92101                         **)
% (**      44321 Nantes cedex 3, FRANCE      **)
% (**      www.irccyn.ec-nantes.fr           **)
% (********************************************)

%    Name of file : E:\simpleArm\simpleArm.dyn




%      Geometric parameters   

% j     ant   mu    sigma gamma b     alpha d     theta r
% 1     0     1     0     0     0     a1    0     t1    0
% 2     1     1     0     0     0     a2    0     t2    0
% 3     2     1     0     0     0     a3    d3    t3    0
% 4     3     0     2     0     0     a4    d4    t4    0



%              Inertial parameters

% j     XX    XY    XZ    YY    YZ    ZZ    MX    MY    MZ    M     Ia

% 1     0     0     0     0     0     0     0     0     0     0     0

% 2     XX2   XY2   XZ2   YY2   YZ2   ZZ2   MX2   MY2   MZ2   M2    0

% 3     XX3   XY3   XZ3   YY3   YZ3   ZZ3   MX3   MY3   MZ3   M3    0

% 4     0     0     0     0     0     0     0     0     0     0     0



%  External forces,friction parameters, joint velocities and accelerations

% j      FX     FY     FZ     CX     CY     CZ     FS     FV     QP     QDP

% 1      0      0      0      0      0      0      0      0      QP1    QDP1

% 2      0      0      0      0      0      0      0      0      QP2    QDP2

% 3      0      0      0      0      0      0      0      0      QP3    QDP3

% 4      0      0      0      0      0      0      0      0      0      0

% Base velocity, base accelerations, and gravity

% j     W0    WP0   V0    VP0   G

% 1     0     0     0     0     G1

% 2     0     0     0     0     G2

% 3     0     0     0     0     G3

%  Dynamic model: Newton Euler method
% Equations:

% Declaration of the function
function torques = inverseDynamics3Dof(armModel, accelerations)

    % Arm state
    state = armModel.getState();
   
    % Arm current state
    [t1, t2, t3] = feval(@(x) x{:}, num2cell(state(1,:)));
    [QP1, QP2, QP3] = feval(@(x) x{:}, num2cell(state(2,:)));
    [QDP1, QDP2, QDP3] = feval(@(x) x{:}, num2cell(accelerations));
    
    joint1 = armModel.joints(1);
    joint2 = armModel.joints(2);
    joint3 = armModel.joints(3);
    
    % Arm MDH parameters
    a1 = joint1.twist;
    a2 = joint2.twist;
    a3 = joint3.twist;
    d3 = joint3.length;
        
    % First moment of inertia
    [MX2, MY2, MZ2] = feval(@(x) x{:}, num2cell((joint2.com * joint2.mass)));
    [MX3, MY3 ,MZ3] = feval(@(x) x{:}, num2cell((joint3.com * joint3.mass)));
    M3 = joint3.mass;
    
    % Second moment of inertia, i.e., inertia tensor
    tensor2 = joint2.inertiaTensor;
    tensor3 = joint3.inertiaTensor;
    
    XX2 = tensor2(1,1);
    XY2 = tensor2(1,2);
    XZ2 = tensor2(1,3);
    YY2 = tensor2(2,2);
    YZ2 = tensor2(2,3);
    ZZ2 = tensor2(3,3);
    
    XX3 = tensor3(1,1);
    XY3 = tensor3(1,2);
    XZ3 = tensor3(1,3);
    YY3 = tensor3(2,2);
    YZ3 = tensor3(2,3);
    ZZ3 = tensor3(3,3);
     
    % Gravity
    [G1, G2, G3] = feval(@(x) x{:}, num2cell(armModel.gravity));

% Function description:
    S1=sin(t1);
	C1=cos(t1);
	Sa1=sin(a1);
	Ca1=cos(a1);
	A211=Ca1.*S1;
	A221=C1.*Ca1;
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
	VP11=-(C1.*G1) - A211.*G2 - A311.*G3;
	VP21=-(A221.*G2) - A321.*G3 + G1.*S1;
	VP31=-(Ca1.*G3) + G2.*Sa1;
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
	U112=DV222 + DV332;
	U122=DV122 - WP32;
	U132=DV132 + WP22;
	U212=DV122 + WP32;
	U232=DV232 - WP12;
	U312=DV132 - WP22;
	U322=DV232 + WP12;
	VP12=C2.*VP11 + A212.*VP21 + A312.*VP31;
	VP22=-(S2.*VP11) + A222.*VP21 + A322.*VP31;
	VP32=-(Sa2.*VP21) + Ca2.*VP31;
	PIS12=-YY2 + ZZ2;
	PIS22=XX2 - ZZ2;
	PIS32=-XX2 + YY2;
	No12=DV232.*PIS12 + WP12.*XX2 - U312.*XY2 + U212.*XZ2 + (-DV222 + DV332).*YZ2;
	No22=DV132.*PIS22 + U322.*XY2 + (DV112 - DV332).*XZ2 + WP22.*YY2 - U122.*YZ2;
	No32=DV122.*PIS32 + (-DV112 + DV222).*XY2 - U232.*XZ2 + U132.*YZ2 + WP32.*ZZ2;
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
	VSP13=d3.*U112 + VP12;
	VSP23=d3.*U212 + VP22;
	VSP33=d3.*U312 + VP32;
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
	N13=No13 - MZ3.*VP23 + MY3.*VP33;
	N23=No23 + MZ3.*VP13 - MX3.*VP33;
	N33=No33 - MY3.*VP13 + MX3.*VP23;
	FDI23=A213.*F13 + A223.*F23 - F33.*Sa3;
	FDI33=A313.*F13 + A323.*F23 + Ca3.*F33;
	N12=C3.*N13 + No12 - N23.*S3 - MZ2.*VP22 + MY2.*VP32;
	N22=-(d3.*FDI33) + A213.*N13 + A223.*N23 + No22 - N33.*Sa3 + MZ2.*VP12 - MX2.*VP32;
	N32=d3.*FDI23 + A313.*N13 + A323.*N23 + Ca3.*N33 + No32 - MY2.*VP12 + MX2.*VP22;
	N31=A312.*N12 + A322.*N22 + Ca2.*N32;
	GAM1=N31;
	GAM2=N32;
	GAM3=N33;

    torques(1) = GAM1; torques(2) = GAM2; torques(3) = GAM3;
    
    
% *=*
% Number of operations : 218 '+' or '-', 274 '*' or '/'
