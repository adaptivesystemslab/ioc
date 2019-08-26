% (********************************************)
% (** SYMORO+ : SYmbolic MOdelling of RObots **)
% (**========================================**)
% (**      IRCCyN-ECN - 1, rue de la Noe     **)
% (**      B.P.92101                         **)
% (**      44321 Nantes cedex 3, FRANCE      **)
% (**      www.irccyn.ec-nantes.fr           **)
% (********************************************)

%    Name of file : E:\simpleArm\simpleArm.ddm




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


% Equations:

% Declaration of the function
function accelerations = fowardDynamics3Dof(armModel, tau)

    % Forces to apply
    GAM1 = tau(1);
    GAM2 = tau(2);
    GAM3 = tau(3);

    % Arm state
    state = armModel.getState();
    
    % Arm current state
    [t1, t2, t3] = feval(@(x) x{:}, num2cell(state(1,:)));
    [QP1, QP2, QP3] = feval(@(x) x{:}, num2cell(state(2,:)));
    
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
	S2=sin(t2);
	C2=cos(t2);
	Sa2=sin(a2);
	Ca2=cos(a2);
	S3=sin(t3);
	C3=cos(t3);
	Sa3=sin(a3);
	Ca3=cos(a3);
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
	A213=Ca3.*S3;
	A223=C3.*Ca3;
	A313=S3.*Sa3;
	A323=C3.*Sa3;
	WI13=A313.*W32 + C3.*WI12 + A213.*WI22;
	WI23=A323.*W32 - S3.*WI12 + A223.*WI22;
	WI33=Ca3.*W32 - Sa3.*WI22;
	W33=QP3 + WI33;
	JPR123=-(A313.*d3);
	JPR133=A213.*d3;
	JPR223=-(A323.*d3);
	JPR233=A223.*d3;
	JPR323=-(Ca3.*d3);
	JPR333=-(d3.*Sa3);
	JW12=WI12.*XX2 + WI22.*XY2 + W32.*XZ2;
	JW22=WI12.*XY2 + WI22.*YY2 + W32.*YZ2;
	JW32=WI12.*XZ2 + WI22.*YZ2 + W32.*ZZ2;
	KW12=-(JW22.*W32) + JW32.*WI22;
	KW22=JW12.*W32 - JW32.*WI12;
	KW32=JW22.*WI12 - JW12.*WI22;
	WQ12=QP2.*WI22;
	WQ22=-(QP2.*WI12);
	JW13=WI13.*XX3 + WI23.*XY3 + W33.*XZ3;
	JW23=WI13.*XY3 + WI23.*YY3 + W33.*YZ3;
	JW33=WI13.*XZ3 + WI23.*YZ3 + W33.*ZZ3;
	KW13=-(JW23.*W33) + JW33.*WI23;
	KW23=JW13.*W33 - JW33.*WI13;
	KW33=JW23.*WI13 - JW13.*WI23;
	SW13=-(W33.*(MX3.*W33 - MZ3.*WI13)) + WI23.*(MY3.*WI13 - MX3.*WI23);
	SW23=-(WI13.*(MY3.*WI13 - MX3.*WI23)) + W33.*(-(MY3.*W33) + MZ3.*WI23);
	SW33=WI13.*(MX3.*W33 - MZ3.*WI13) - WI23.*(-(MY3.*W33) + MZ3.*WI23);
	WQ13=QP3.*WI23;
	WQ23=-(QP3.*WI13);
	LW13=A313.*d3.*W32.*WI12 + A213.*d3.*WI12.*WI22 + C3.*(-(d3.*W32.^2) - d3.*WI22.^2);
	LW23=A323.*d3.*W32.*WI12 + A223.*d3.*WI12.*WI22 - S3.*(-(d3.*W32.^2) - d3.*WI22.^2);
	LW33=Ca3.*d3.*W32.*WI12 - d3.*Sa3.*WI12.*WI22;
	JD3=1./ZZ3;
	JU13=JD3.*XZ3;
	JU23=JD3.*YZ3;
	JU33=JD3.*ZZ3;
	JU43=-(JD3.*MY3);
	JU53=JD3.*MX3;
	GW3=GAM3 - KW33;
	GK113=XX3 - JU13.*XZ3;
	GK123=XY3 - JU13.*YZ3;
	GK133=XZ3 - JU13.*ZZ3;
	GK143=JU13.*MY3;
	GK153=-(JU13.*MX3) - MZ3;
	GK213=XY3 - JU23.*XZ3;
	GK223=YY3 - JU23.*YZ3;
	GK233=YZ3 - JU23.*ZZ3;
	GK243=JU23.*MY3 + MZ3;
	GK253=-(JU23.*MX3);
	GK313=XZ3 - JU33.*XZ3;
	GK323=YZ3 - JU33.*YZ3;
	GK333=ZZ3 - JU33.*ZZ3;
	GK343=-MY3 + JU33.*MY3;
	GK353=MX3 - JU33.*MX3;
	GK413=-(JU43.*XZ3);
	GK423=MZ3 - JU43.*YZ3;
	GK433=-MY3 - JU43.*ZZ3;
	GK443=M3 + JU43.*MY3;
	GK453=-(JU43.*MX3);
	GK513=-MZ3 - JU53.*XZ3;
	GK523=-(JU53.*YZ3);
	GK533=MX3 - JU53.*ZZ3;
	GK543=JU53.*MY3;
	GK553=M3 - JU53.*MX3;
	NG13=GK143.*LW13 + GK153.*LW23 + LW33.*MY3 + GK113.*WQ13 + GK123.*WQ23;
	NG23=GK243.*LW13 + GK253.*LW23 - LW33.*MX3 + GK213.*WQ13 + GK223.*WQ23;
	NG33=GK343.*LW13 + GK353.*LW23 + GK313.*WQ13 + GK323.*WQ23;
	NG43=GK443.*LW13 + GK453.*LW23 + GK413.*WQ13 + GK423.*WQ23;
	NG53=GK543.*LW13 + GK553.*LW23 + GK513.*WQ13 + GK523.*WQ23;
	NG63=LW33.*M3 + MY3.*WQ13 - MX3.*WQ23;
	VS13=GW3.*JU13 + NG13;
	VS23=GW3.*JU23 + NG23;
	VS33=GW3.*JU33 + NG33;
	VS43=GW3.*JU43 + NG43;
	VS53=GW3.*JU53 + NG53;
	AP13=KW13 + VS13;
	AP23=KW23 + VS23;
	AP33=KW33 + VS33;
	AP43=SW13 + VS43;
	AP53=SW23 + VS53;
	AP63=NG63 + SW33;
	GX113=C3.*GK113 - GK213.*S3;
	GX123=C3.*GK123 - GK223.*S3;
	GX213=A213.*GK113 + A223.*GK213 + GK413.*JPR123 + GK513.*JPR223 + JPR323.*MY3 - GK313.*Sa3;
	GX223=A213.*GK123 + A223.*GK223 + GK423.*JPR123 + GK523.*JPR223 - JPR323.*MX3 - GK323.*Sa3;
	GX233=A213.*GK133 + A223.*GK233 + GK433.*JPR123 + GK533.*JPR223 - GK333.*Sa3;
	GX243=A213.*GK143 + A223.*GK243 + GK443.*JPR123 + GK543.*JPR223 - GK343.*Sa3;
	GX253=A213.*GK153 + A223.*GK253 + GK453.*JPR123 + GK553.*JPR223 - GK353.*Sa3;
	GX263=JPR323.*M3 - A223.*MX3 + A213.*MY3;
	GX313=A313.*GK113 + A323.*GK213 + Ca3.*GK313 + GK413.*JPR133 + GK513.*JPR233 + JPR333.*MY3;
	GX323=A313.*GK123 + A323.*GK223 + Ca3.*GK323 + GK423.*JPR133 + GK523.*JPR233 - JPR333.*MX3;
	GX333=A313.*GK133 + A323.*GK233 + Ca3.*GK333 + GK433.*JPR133 + GK533.*JPR233;
	GX343=A313.*GK143 + A323.*GK243 + Ca3.*GK343 + GK443.*JPR133 + GK543.*JPR233;
	GX353=A313.*GK153 + A323.*GK253 + Ca3.*GK353 + GK453.*JPR133 + GK553.*JPR233;
	GX363=JPR333.*M3 - A323.*MX3 + A313.*MY3;
	GX413=C3.*GK413 - GK513.*S3;
	GX423=C3.*GK423 - GK523.*S3;
	GX433=C3.*GK433 - GK533.*S3;
	GX443=C3.*GK443 - GK543.*S3;
	GX453=C3.*GK453 - GK553.*S3;
	GX513=A213.*GK413 + A223.*GK513 - MY3.*Sa3;
	GX523=A213.*GK423 + A223.*GK523 + MX3.*Sa3;
	GX533=A213.*GK433 + A223.*GK533;
	GX543=A213.*GK443 + A223.*GK543;
	GX553=A213.*GK453 + A223.*GK553;
	GX563=-(M3.*Sa3);
	GX613=A313.*GK413 + A323.*GK513 + Ca3.*MY3;
	GX623=A313.*GK423 + A323.*GK523 - Ca3.*MX3;
	GX633=A313.*GK433 + A323.*GK533;
	GX643=A313.*GK443 + A323.*GK543;
	GX653=A313.*GK453 + A323.*GK553;
	GX663=Ca3.*M3;
	TKT113=C3.*GX113 - GX123.*S3;
	TKT213=C3.*GX213 - GX223.*S3;
	TKT313=C3.*GX313 - GX323.*S3;
	TKT413=C3.*GX413 - GX423.*S3;
	TKT513=C3.*GX513 - GX523.*S3;
	TKT613=C3.*GX613 - GX623.*S3;
	TKT223=A213.*GX213 + A223.*GX223 + GX243.*JPR123 + GX253.*JPR223 + GX263.*JPR323 - GX233.*Sa3;
	TKT323=A213.*GX313 + A223.*GX323 + GX343.*JPR123 + GX353.*JPR223 + GX363.*JPR323 - GX333.*Sa3;
	TKT423=A213.*GX413 + A223.*GX423 + GX443.*JPR123 + GX453.*JPR223 - GX433.*Sa3;
	TKT523=A213.*GX513 + A223.*GX523 + GX543.*JPR123 + GX553.*JPR223 + GX563.*JPR323 - GX533.*Sa3;
	TKT623=A213.*GX613 + A223.*GX623 + GX643.*JPR123 + GX653.*JPR223 + GX663.*JPR323 - GX633.*Sa3;
	TKT333=A313.*GX313 + A323.*GX323 + Ca3.*GX333 + GX343.*JPR133 + GX353.*JPR233 + GX363.*JPR333;
	TKT433=A313.*GX413 + A323.*GX423 + Ca3.*GX433 + GX443.*JPR133 + GX453.*JPR233;
	TKT533=A313.*GX513 + A323.*GX523 + Ca3.*GX533 + GX543.*JPR133 + GX553.*JPR233 + GX563.*JPR333;
	TKT633=A313.*GX613 + A323.*GX623 + Ca3.*GX633 + GX643.*JPR133 + GX653.*JPR233 + GX663.*JPR333;
	MJE112=TKT113 + XX2;
	MJE212=TKT213 + XY2;
	MJE312=TKT313 + XZ2;
	MJE512=-MZ2 + TKT513;
	MJE612=MY2 + TKT613;
	MJE222=TKT223 + YY2;
	MJE322=TKT323 + YZ2;
	MJE422=MZ2 + TKT423;
	MJE622=-MX2 + TKT623;
	MJE332=TKT333 + ZZ2;
	MJE432=-MY2 + TKT433;
	MJE532=MX2 + TKT533;
	VBE12=-(AP13.*C3) - KW12 + AP23.*S3;
	VBE22=-(A213.*AP13) - A223.*AP23 - AP43.*JPR123 - AP53.*JPR223 - AP63.*JPR323 - KW22 + AP33.*Sa3;
	VBE32=-(A313.*AP13) - A323.*AP23 - AP33.*Ca3 - AP43.*JPR133 - AP53.*JPR233 - AP63.*JPR333 - KW32;
	JD2=1./MJE332;
	JU12=JD2.*MJE312;
	JU22=JD2.*MJE322;
	JU32=JD2.*MJE332;
	JU42=JD2.*MJE432;
	JU52=JD2.*MJE532;
	JU62=JD2.*TKT633;
	GW2=GAM2 + VBE32;
	GK112=MJE112 - JU12.*MJE312;
	GK122=MJE212 - JU12.*MJE322;
	GK132=MJE312 - JU12.*MJE332;
	GK212=MJE212 - JU22.*MJE312;
	GK222=MJE222 - JU22.*MJE322;
	GK232=MJE322 - JU22.*MJE332;
	GK312=MJE312 - JU32.*MJE312;
	GK322=MJE322 - JU32.*MJE322;
	GK332=MJE332 - JU32.*MJE332;
	GK412=-(JU42.*MJE312) + TKT413;
	GK422=-(JU42.*MJE322) + MJE422;
	GK432=-(JU42.*MJE332) + MJE432;
	GK512=-(JU52.*MJE312) + MJE512;
	GK522=-(JU52.*MJE322) + TKT523;
	GK532=-(JU52.*MJE332) + MJE532;
	GK612=-(JU62.*MJE312) + MJE612;
	GK622=-(JU62.*MJE322) + MJE622;
	GK632=-(JU62.*MJE332) + TKT633;
	NG12=GK112.*WQ12 + GK122.*WQ22;
	NG22=GK212.*WQ12 + GK222.*WQ22;
	NG32=GK312.*WQ12 + GK322.*WQ22;
	VS12=GW2.*JU12 + NG12;
	VS22=GW2.*JU22 + NG22;
	VS32=GW2.*JU32 + NG32;
	AP12=-VBE12 + VS12;
	AP22=-VBE22 + VS22;
	AP32=-VBE32 + VS32;
	GX312=A312.*GK112 + A322.*GK212 + Ca2.*GK312;
	GX322=A312.*GK122 + A322.*GK222 + Ca2.*GK322;
	GX332=A312.*GK132 + A322.*GK232 + Ca2.*GK332;
	GX412=C2.*GK412 - GK512.*S2;
	GX422=C2.*GK422 - GK522.*S2;
	GX432=C2.*GK432 - GK532.*S2;
	GX512=A212.*GK412 + A222.*GK512 - GK612.*Sa2;
	GX522=A212.*GK422 + A222.*GK522 - GK622.*Sa2;
	GX532=A212.*GK432 + A222.*GK532 - GK632.*Sa2;
	GX612=A312.*GK412 + A322.*GK512 + Ca2.*GK612;
	GX622=A312.*GK422 + A322.*GK522 + Ca2.*GK622;
	GX632=A312.*GK432 + A322.*GK532 + Ca2.*GK632;
	TKT332=A312.*GX312 + A322.*GX322 + Ca2.*GX332;
	TKT432=A312.*GX412 + A322.*GX422 + Ca2.*GX432;
	TKT532=A312.*GX512 + A322.*GX522 + Ca2.*GX532;
	TKT632=A312.*GX612 + A322.*GX622 + Ca2.*GX632;
	VBE31=-(A312.*AP12) - A322.*AP22 - AP32.*Ca2;
	JD1=1./TKT332;
	JU41=JD1.*TKT432;
	JU51=JD1.*TKT532;
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
	VR42=C2.*VR41 + A212.*VR51 + A312.*VR61;
	VR52=-(S2.*VR41) + A222.*VR51 + A322.*VR61;
	VR62=-(Sa2.*VR51) + Ca2.*VR61;
	GU2=JU12.*VR12 + JU22.*VR22 + JU32.*VR32 + JU42.*VR42 + JU52.*VR52 + JU62.*VR62;
	QDP2=-GU2 + GW2.*JD2;
	WP32=QDP2 + VR32;
	VR13=C3.*VR12 + A213.*VR22 + A313.*WP32 + WQ13;
	VR23=-(S3.*VR12) + A223.*VR22 + A323.*WP32 + WQ23;
	VR33=-(Sa3.*VR22) + Ca3.*WP32;
	VR43=LW13 + JPR123.*VR22 + C3.*VR42 + A213.*VR52 + A313.*VR62 + JPR133.*WP32;
	VR53=LW23 + JPR223.*VR22 - S3.*VR42 + A223.*VR52 + A323.*VR62 + JPR233.*WP32;
	GU3=JU13.*VR13 + JU23.*VR23 + JU33.*VR33 + JU43.*VR43 + JU53.*VR53;
	QDP3=-GU3 + GW3.*JD3;
	QDP4= 0;
    
    accelerations(1) = QDP1;
    accelerations(2) = QDP2;
    accelerations(3) = QDP3;
    
% *=*
% Number of operations : 339 '+' or '-', 488 '*' or '/'

