% (********************************************)
% (** SYMORO+ : SYmbolic MOdelling of RObots **)
% (**========================================**)
% (**      IRCCyN-ECN - 1, rue de la Noe     **)
% (**      B.P.92101                         **)
% (**      44321 Nantes cedex 3, FRANCE      **)
% (**      www.irccyn.ec-nantes.fr           **)
% (********************************************)


%    Name of file : E:\2DofArm_Symoro\Arm2Dof.ddm




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

function accelerations = forwardDynamics2Dof(armModel, tau)

    % Torques to apply
    GAM1 = tau(1);
    GAM2 = tau(2);
    
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
	S2=sin(t2);
	C2=cos(t2);
	Sa2=sin(a2);
	Ca2=cos(a2);
	A211=Ca1.*S1;
	A221=C1.*Ca1;
	A311=S1.*Sa1;
	A321=C1.*Sa1;
	A212=Ca2.*S2;
	A222=C2.*Ca2;
	A312=S2.*Sa2;
	A322=C2.*Sa2;
	WI12=A312.*QP1;
	WI22=A322.*QP1;
	WI32=Ca2.*QP1;
	W32=QP2 + WI32;
	JPR132=A212.*d2;
	JPR232=A222.*d2;
	JPR332=-(d2.*Sa2);
	JW12=WI12.*XX2 + WI22.*XY2 + W32.*XZ2;
	JW22=WI12.*XY2 + WI22.*YY2 + W32.*YZ2;
	JW32=WI12.*XZ2 + WI22.*YZ2 + W32.*ZZ2;
	KW12=-(JW22.*W32) + JW32.*WI22;
	KW22=JW12.*W32 - JW32.*WI12;
	KW32=JW22.*WI12 - JW12.*WI22;
	SW12=-(W32.*(MX2.*W32 - MZ2.*WI12)) + WI22.*(MY2.*WI12 - MX2.*WI22);
	SW22=-(WI12.*(MY2.*WI12 - MX2.*WI22)) + W32.*(-(MY2.*W32) + MZ2.*WI22);
	SW32=WI12.*(MX2.*W32 - MZ2.*WI12) - WI22.*(-(MY2.*W32) + MZ2.*WI22);
	WQ12=QP2.*WI22;
	WQ22=-(QP2.*WI12);
	LW12=-(C2.*d2.*QP1.^2);
	LW22=d2.*QP1.^2.*S2;
	JD2=1./ZZ2;
	JU12=JD2.*XZ2;
	JU22=JD2.*YZ2;
	JU32=JD2.*ZZ2;
	JU42=-(JD2.*MY2);
	JU52=JD2.*MX2;
	GW2=GAM2 - KW32;
	GK112=XX2 - JU12.*XZ2;
	GK122=XY2 - JU12.*YZ2;
	GK132=XZ2 - JU12.*ZZ2;
	GK142=JU12.*MY2;
	GK152=-(JU12.*MX2) - MZ2;
	GK212=XY2 - JU22.*XZ2;
	GK222=YY2 - JU22.*YZ2;
	GK232=YZ2 - JU22.*ZZ2;
	GK242=JU22.*MY2 + MZ2;
	GK252=-(JU22.*MX2);
	GK312=XZ2 - JU32.*XZ2;
	GK322=YZ2 - JU32.*YZ2;
	GK332=ZZ2 - JU32.*ZZ2;
	GK342=-MY2 + JU32.*MY2;
	GK352=MX2 - JU32.*MX2;
	GK412=-(JU42.*XZ2);
	GK422=MZ2 - JU42.*YZ2;
	GK432=-MY2 - JU42.*ZZ2;
	GK442=M2 + JU42.*MY2;
	GK452=-(JU42.*MX2);
	GK512=-MZ2 - JU52.*XZ2;
	GK522=-(JU52.*YZ2);
	GK532=MX2 - JU52.*ZZ2;
	GK542=JU52.*MY2;
	GK552=M2 - JU52.*MX2;
	NG12=GK142.*LW12 + GK152.*LW22 + GK112.*WQ12 + GK122.*WQ22;
	NG22=GK242.*LW12 + GK252.*LW22 + GK212.*WQ12 + GK222.*WQ22;
	NG32=GK342.*LW12 + GK352.*LW22 + GK312.*WQ12 + GK322.*WQ22;
	NG42=GK442.*LW12 + GK452.*LW22 + GK412.*WQ12 + GK422.*WQ22;
	NG52=GK542.*LW12 + GK552.*LW22 + GK512.*WQ12 + GK522.*WQ22;
	NG62=MY2.*WQ12 - MX2.*WQ22;
	VS12=GW2.*JU12 + NG12;
	VS22=GW2.*JU22 + NG22;
	VS32=GW2.*JU32 + NG32;
	VS42=GW2.*JU42 + NG42;
	VS52=GW2.*JU52 + NG52;
	AP12=KW12 + VS12;
	AP22=KW22 + VS22;
	AP32=KW32 + VS32;
	AP42=SW12 + VS42;
	AP52=SW22 + VS52;
	AP62=NG62 + SW32;
	GX312=A312.*GK112 + A322.*GK212 + Ca2.*GK312 + GK412.*JPR132 + GK512.*JPR232 + JPR332.*MY2;
	GX322=A312.*GK122 + A322.*GK222 + Ca2.*GK322 + GK422.*JPR132 + GK522.*JPR232 - JPR332.*MX2;
	GX332=A312.*GK132 + A322.*GK232 + Ca2.*GK332 + GK432.*JPR132 + GK532.*JPR232;
	GX342=A312.*GK142 + A322.*GK242 + Ca2.*GK342 + GK442.*JPR132 + GK542.*JPR232;
	GX352=A312.*GK152 + A322.*GK252 + Ca2.*GK352 + GK452.*JPR132 + GK552.*JPR232;
	GX362=JPR332.*M2 - A322.*MX2 + A312.*MY2;
	GX412=C2.*GK412 - GK512.*S2;
	GX422=C2.*GK422 - GK522.*S2;
	GX432=C2.*GK432 - GK532.*S2;
	GX442=C2.*GK442 - GK542.*S2;
	GX452=C2.*GK452 - GK552.*S2;
	GX512=A212.*GK412 + A222.*GK512 - MY2.*Sa2;
	GX522=A212.*GK422 + A222.*GK522 + MX2.*Sa2;
	GX532=A212.*GK432 + A222.*GK532;
	GX542=A212.*GK442 + A222.*GK542;
	GX552=A212.*GK452 + A222.*GK552;
	GX562=-(M2.*Sa2);
	GX612=A312.*GK412 + A322.*GK512 + Ca2.*MY2;
	GX622=A312.*GK422 + A322.*GK522 - Ca2.*MX2;
	GX632=A312.*GK432 + A322.*GK532;
	GX642=A312.*GK442 + A322.*GK542;
	GX652=A312.*GK452 + A322.*GK552;
	GX662=Ca2.*M2;
	TKT332=A312.*GX312 + A322.*GX322 + Ca2.*GX332 + GX342.*JPR132 + GX352.*JPR232 + GX362.*JPR332;
	TKT432=A312.*GX412 + A322.*GX422 + Ca2.*GX432 + GX442.*JPR132 + GX452.*JPR232;
	TKT532=A312.*GX512 + A322.*GX522 + Ca2.*GX532 + GX542.*JPR132 + GX552.*JPR232 + GX562.*JPR332;
	TKT632=A312.*GX612 + A322.*GX622 + Ca2.*GX632 + GX642.*JPR132 + GX652.*JPR232 + GX662.*JPR332;
	MJE331=TKT332 + ZZ1;
	MJE431=-MY1 + TKT432;
	MJE531=MX1 + TKT532;
	VBE31=-(A312.*AP12) - A322.*AP22 - AP32.*Ca2 - AP42.*JPR132 - AP52.*JPR232 - AP62.*JPR332;
	JD1=1./MJE331;
	JU41=JD1.*MJE431;
	JU51=JD1.*MJE531;
	JU61=JD1.*TKT632;
	GW1=GAM1 + VBE31;
	VR41=-(C1.*G1) - A211.*G2 - A311.*G3;
	VR51=-(A221.*G2) - A321.*G3 + G1.*S1;
	VR61=-(Ca1.*G3) + G2.*Sa1;
	GU1=JU41.*VR41 + JU51.*VR51 + JU61.*VR61;
	QDP1=-GU1 + GW1.*JD1;
	VR12=A312.*QDP1 + WQ12;
	VR22=A322.*QDP1 + WQ22;
	VR32=Ca2.*QDP1;
	VR42=LW12 + JPR132.*QDP1 + C2.*VR41 + A212.*VR51 + A312.*VR61;
	VR52=LW22 + JPR232.*QDP1 - S2.*VR41 + A222.*VR51 + A322.*VR61;
	GU2=JU12.*VR12 + JU22.*VR22 + JU32.*VR32 + JU42.*VR42 + JU52.*VR52;
	QDP2=-GU2 + GW2.*JD2;
	QDP3= 0;

    accelerations = [QDP1 ,QDP2];

% *=*
% Number of operations : 160 '+' or '-', 241 '*' or '/'

end

