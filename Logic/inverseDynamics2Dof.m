% (********************************************)
% (** SYMORO+ : SYmbolic MOdelling of RObots **)
% (**========================================**)
% (**      IRCCyN-ECN - 1, rue de la Noe     **)
% (**      B.P.92101                         **)
% (**      44321 Nantes cedex 3, FRANCE      **)
% (**      www.irccyn.ec-nantes.fr           **)
% (********************************************)


%    Name of file : E:\2DofArm_Symoro\Arm2Dof.dyn




%      Geometric parameters   

% j     ant   mu    sigma gamma b     alpha d     theta r
% 1     0     1     0     0     0     a1    0     t1    0
% 2     1     1     0     0     0     a2    d2    t2    0
% 3     2     0     2     0     0     a3    d3    t3    0



%              Inertial parameters

% j     XX    XY    XZ    YY    YZ    ZZ    MX    MY    MZ    M     Ia

% 1     XX1   XY1   XZ1   YY1   YZ1   ZZ1   MX1   MY1   MZ1   M1    0

% 2     XX2   XY2   XZ2   YY2   YZ2   ZZ2   MX2   MY2   MZ2   M2    0

% 3     0     0     0     0     0     0     0     0     0     0     0



%  External forces,friction parameters, joint velocities and accelerations

% j      FX     FY     FZ     CX     CY     CZ     FS     FV     QP     QDP

% 1      0      0      0      0      0      0      0      0      QP1    QDP1

% 2      0      0      0      0      0      0      0      0      QP2    QDP2

% 3      0      0      0      0      0      0      0      0      0      0

% Base velocity, base accelerations, and gravity

% j     W0    WP0   V0    VP0   G

% 1     0     0     0     0     G1

% 2     0     0     0     0     G2

% 3     0     0     0     0     G3

function torques = inverseDynamics2Dof(armModel, accelerations)

    QDP1 = accelerations(1);
    QDP2 = accelerations(2);
    
     % Arm current state
    state = armModel.getState();
    [t1, t2] = feval(@(x) x{:}, num2cell(state(1,:)));
    [QP1, QP2] = feval(@(x) x{:}, num2cell(state(2,:)));
    
    % Arm MDH parameters
    joint1 = armModel.joints(1);
    joint2 = armModel.joints(2);
    a1 = joint1.twist;
    a2 = joint2.twist;
    d2 = joint2.length;
    
    % Definition of first moments of inertia
    [MX1, MY1, ~] = feval(@(x) x{:}, num2cell((joint1.com * joint1.mass)));
    [MX2, MY2 ,MZ2] = feval(@(x) x{:}, num2cell((joint2.com * joint2.mass)));
    M2 = joint2.mass;
    
    % Second moment of inertia, i.e., inertia tensor
    tensor1 = joint1.inertiaTensor;
    tensor2 = joint2.inertiaTensor;
    
    ZZ1 = tensor1(3,3);
    
    XX2 = tensor2(1,1);
    XY2 = tensor2(1,2);
    XZ2 = tensor2(1,3);
    YY2 = tensor2(2,2);
    YZ2 = tensor2(2,3);
    ZZ2 = tensor2(3,3);
    
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
	DV331=-QP1.^2;
	VP11=-(C1.*G1) - A211.*G2 - A311.*G3;
	VP21=-(A221.*G2) - A321.*G3 + G1.*S1;
	VP31=-(Ca1.*G3) + G2.*Sa1;
	No31=QDP1.*ZZ1;
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
	U222=DV112 + DV332;
	U232=DV232 - WP12;
	U312=DV132 - WP22;
	U322=DV232 + WP12;
	U332=DV112 + DV222;
	VSP12=d2.*DV331 + VP11;
	VSP22=d2.*QDP1 + VP21;
	VP12=A312.*VP31 + C2.*VSP12 + A212.*VSP22;
	VP22=A322.*VP31 - S2.*VSP12 + A222.*VSP22;
	VP32=Ca2.*VP31 - Sa2.*VSP22;
	F12=MX2.*U112 + MY2.*U122 + MZ2.*U132 + M2.*VP12;
	F22=MX2.*U212 + MY2.*U222 + MZ2.*U232 + M2.*VP22;
	F32=MX2.*U312 + MY2.*U322 + MZ2.*U332 + M2.*VP32;
	PIS12=-YY2 + ZZ2;
	PIS22=XX2 - ZZ2;
	PIS32=-XX2 + YY2;
	No12=DV232.*PIS12 + WP12.*XX2 - U312.*XY2 + U212.*XZ2 + (-DV222 + DV332).*YZ2;
	No22=DV132.*PIS22 + U322.*XY2 + (DV112 - DV332).*XZ2 + WP22.*YY2 - U122.*YZ2;
	No32=DV122.*PIS32 + (-DV112 + DV222).*XY2 - U232.*XZ2 + U132.*YZ2 + WP32.*ZZ2;
	N12=No12 - MZ2.*VP22 + MY2.*VP32;
	N22=No22 + MZ2.*VP12 - MX2.*VP32;
	N32=No32 - MY2.*VP12 + MX2.*VP22;
	FDI22=A212.*F12 + A222.*F22 - F32.*Sa2;
	N31=d2.*FDI22 + A312.*N12 + A322.*N22 + Ca2.*N32 + No31 - MY1.*VP11 + MX1.*VP21;
	GAM1=N31;
	GAM2=N32;

    torques = [GAM1 GAM2];
    
% *=*
% Number of operations : 66 '+' or '-', 84 '*' or '/'
end

