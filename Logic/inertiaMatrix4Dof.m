% (********************************************)
% (** SYMORO+ : SYmbolic MOdelling of RObots **)
% (**========================================**)
% (**      IRCCyN-ECN - 1, rue de la Noe     **)
% (**      B.P.92101                         **)
% (**      44321 Nantes cedex 3, FRANCE      **)
% (**      www.irccyn.ec-nantes.fr           **)
% (********************************************)


%    Name of file : C:\Users\dell.dell-PC\Desktop\arm_ioc\v3\arm_4dofv3.inm


%      Geometric parameters   

% j     ant   mu    sigma gamma b     alpha d     theta r
% 1     0     1     0     0     0     a1    0     t1    0
% 2     1     1     0     0     0     a2    0     t2    0
% 3     2     1     0     0     0     a3    0     t3    r3
% 4     3     1     0     0     0     a4    0     t4    0
% 5     4     0     2     0     0     a5    0     t5    r5



%              Inertial parameters

% j     XX    XY    XZ    YY    YZ    ZZ    MX    MY    MZ    M     Ia

% 1     0     0     0     0     0     0     0     0     0     0     0

% 2     XX2   XY2   XZ2   YY2   YZ2   ZZ2   MX2   MY2   MZ2   M2    0

% 3     0     0     0     0     0     0     0     0     0     0     0

% 4     XX4   XY4   XZ4   YY4   YZ4   ZZ4   MX4   MY4   MZ4   M4    0

% 5     0     0     0     0     0     0     0     0     0     0     0

%                       Inertia matrix
% Equations:

% Declaration of the function
function A = inertiaMatrix4Dof(armModel)

    % IMPORTANT: Consider transposing the matrix as shown in Liz's code

    % Arm state
    state = armModel.getState();
       
    % Arm current state
    [~, t2, t3, t4] = feval(@(x) x{:}, num2cell(state(1,:)));
    
    joint2 = armModel.joints(2);
    joint3 = armModel.joints(3);
    joint4 = armModel.joints(4);
    joint5 = armModel.joints(5);
    
    % Arm MDH parameters
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
        
% Function description:

    S2=sin(t2);
	C2=cos(t2);
	Sa2=sin(a2);
	Ca2=cos(a2);
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
	AS25=Ca5.*MY5 - MZ5.*Sa5;
	AS35=Ca5.*MZ5 + MY5.*Sa5;
	AJ215=Ca5.*XY5 - Sa5.*XZ5;
	AJ225=Ca5.*YY5 - Sa5.*YZ5;
	AJ235=Ca5.*YZ5 - Sa5.*ZZ5;
	AJ315=Sa5.*XY5 + Ca5.*XZ5;
	AJ325=Sa5.*YY5 + Ca5.*YZ5;
	AJ335=Sa5.*YZ5 + Ca5.*ZZ5;
	AJA225=AJ225.*Ca5 - AJ235.*Sa5;
	AJA325=AJ325.*Ca5 - AJ335.*Sa5;
	AJA335=AJ335.*Ca5 + AJ325.*Sa5;
	PAS115=-(AS25.*LOO25) - AS35.*LOO35;
	PAS125=LOO25.*MX5;
	PAS135=LOO35.*MX5;
	PAS225=-(AS35.*LOO35);
	PAS235=AS25.*LOO35;
	PAS325=AS35.*LOO25;
	PAS335=-(AS25.*LOO25);
	XXP4=(LOO25.^2 + LOO35.^2).*M5 - 2.*PAS115 + XX5;
	XYP4=AJ215 - PAS125;
	XZP4=AJ315 - PAS135;
	YYP4=AJA225 + LOO35.^2.*M5 - 2.*PAS225;
	YZP4=AJA325 - LOO25.*LOO35.*M5 - PAS235 - PAS325;
	ZZP4=AJA335 + LOO25.^2.*M5 - 2.*PAS335;
	MYP4=AS25 + LOO25.*M5;
	MZP4=AS35 + LOO35.*M5;
	AS14=C4.*MX5 - MYP4.*S4;
	AS24=A214.*MX5 + A224.*MYP4 - MZP4.*Sa4;
	AS34=A314.*MX5 + A324.*MYP4 + Ca4.*MZP4;
	AJ114=-((AJ215 - PAS125).*S4) + C4.*XXP4;
	AJ124=C4.*XYP4 - S4.*YYP4;
	AJ134=C4.*XZP4 - S4.*YZP4;
	AJ214=A224.*(AJ215 - PAS125) - (AJ315 - PAS135).*Sa4 + A214.*XXP4;
	AJ224=-((AJA325 - LOO25.*LOO35.*M5 - PAS235 - PAS325).*Sa4) + A214.*XYP4 + A224.*YYP4;
	AJ234=A214.*XZP4 + A224.*YZP4 - Sa4.*ZZP4;
	AJ314=A324.*(AJ215 - PAS125) + Ca4.*(AJ315 - PAS135) + A314.*XXP4;
	AJ324=Ca4.*(AJA325 - LOO25.*LOO35.*M5 - PAS235 - PAS325) + A314.*XYP4 + A324.*YYP4;
	AJ334=A314.*XZP4 + A324.*YZP4 + Ca4.*ZZP4;
	AJA114=AJ114.*C4 - AJ124.*S4;
	AJA214=AJ214.*C4 - AJ224.*S4;
	AJA314=AJ314.*C4 - AJ324.*S4;
	AJA224=A214.*AJ214 + A224.*AJ224 - AJ234.*Sa4;
	AJA324=A214.*AJ314 + A224.*AJ324 - AJ334.*Sa4;
	AJA334=A314.*AJ314 + A324.*AJ324 + AJ334.*Ca4;
	XXP3=AJA114 + XX3;
	XYP3=AJA214 + XY3;
	XZP3=AJA314 + XZ3;
	YYP3=AJA224 + YY3;
	YZP3=AJA324 + YZ3;
	ZZP3=AJA334 + ZZ3;
	MXP3=AS14 + MX3;
	MYP3=AS24 + MY3;
	MZP3=AS34 + MZ3;
	MP3=M3 + M5;
	AS13=C3.*MXP3 - MYP3.*S3;
	AS23=A213.*MXP3 + A223.*MYP3 - MZP3.*Sa3;
	AS33=A313.*MXP3 + A323.*MYP3 + Ca3.*MZP3;
	AJ113=C3.*XXP3 - S3.*(AJA214 + XY3);
	AJ123=C3.*XYP3 - S3.*YYP3;
	AJ133=C3.*XZP3 - S3.*YZP3;
	AJ213=A213.*XXP3 + A223.*(AJA214 + XY3) - Sa3.*(AJA314 + XZ3);
	AJ223=A213.*XYP3 + A223.*YYP3 - Sa3.*(AJA324 + YZ3);
	AJ233=A213.*XZP3 + A223.*YZP3 - Sa3.*ZZP3;
	AJ313=A313.*XXP3 + A323.*(AJA214 + XY3) + Ca3.*(AJA314 + XZ3);
	AJ323=A313.*XYP3 + A323.*YYP3 + Ca3.*(AJA324 + YZ3);
	AJ333=A313.*XZP3 + A323.*YZP3 + Ca3.*ZZP3;
	AJA113=AJ113.*C3 - AJ123.*S3;
	AJA213=AJ213.*C3 - AJ223.*S3;
	AJA313=AJ313.*C3 - AJ323.*S3;
	AJA223=A213.*AJ213 + A223.*AJ223 - AJ233.*Sa3;
	AJA323=A213.*AJ313 + A223.*AJ323 - AJ333.*Sa3;
	AJA333=A313.*AJ313 + A323.*AJ323 + AJ333.*Ca3;
	PAS113=-(AS23.*LOO23) - AS33.*LOO33;
	PAS123=AS13.*LOO23;
	PAS133=AS13.*LOO33;
	PAS223=-(AS33.*LOO33);
	PAS233=AS23.*LOO33;
	PAS323=AS33.*LOO23;
	PAS333=-(AS23.*LOO23);
	XXP2=AJA113 + (LOO23.^2 + LOO33.^2).*MP3 - 2.*PAS113;
	XYP2=AJA213 - PAS123;
	XZP2=AJA313 - PAS133;
	YYP2=AJA223 + LOO33.^2.*MP3 - 2.*PAS223;
	YZP2=AJA323 - LOO23.*LOO33.*MP3 - PAS233 - PAS323;
	ZZP2=AJA333 + LOO23.^2.*MP3 - 2.*PAS333;
	AJ312=A322.*(AJA213 - PAS123) + Ca2.*(AJA313 - PAS133) + A312.*XXP2;
	AJ322=Ca2.*(AJA323 - LOO23.*LOO33.*MP3 - PAS233 - PAS323) + A312.*XYP2 + A322.*YYP2;
	AJ332=A312.*XZP2 + A322.*YZP2 + Ca2.*ZZP2;
	AJA332=A312.*AJ312 + A322.*AJ322 + AJ332.*Ca2;
	EC12=-(C3.*MYP3) - MXP3.*S3;
	EC22=A223.*MXP3 - A213.*MYP3;
	EC32=A323.*MXP3 - A313.*MYP3;
	NC12=AJ133 + EC32.*LOO23 - EC22.*LOO33;
	NC22=AJ233 + EC12.*LOO33;
	NC32=AJ333 - EC12.*LOO23;
	NC33=A312.*NC12 + A322.*NC22 + Ca2.*NC32;
	ED12=-(C4.*MYP4) - MX5.*S4;
	ED22=A224.*MX5 - A214.*MYP4;
	ED32=A324.*MX5 - A314.*MYP4;
	ED13=C3.*ED12 - ED22.*S3;
	ED23=A213.*ED12 + A223.*ED22 - ED32.*Sa3;
	ED33=A313.*ED12 + A323.*ED22 + Ca3.*ED32;
	ND13=AJ134.*C3 + ED33.*LOO23 - ED23.*LOO33 - AJ234.*S3;
	ND23=A213.*AJ134 + A223.*AJ234 + ED13.*LOO33 - AJ334.*Sa3;
	ND33=A313.*AJ134 + A323.*AJ234 + AJ334.*Ca3 - ED13.*LOO23;
	ND34=A312.*ND13 + A322.*ND23 + Ca2.*ND33;
	
    A = [0 0 0 0;
        0 0 0 0;
        0 0 0 0;
        0 0 0 0];
    
    A(1,1)=AJA332;
	A(2,1)=AJ332;
	A(3,1)=NC33;
	A(4,1)=ND34;
	A(2,2)=ZZP2;
	A(3,2)=NC32;
	A(4,2)=ND33;
	A(3,3)=ZZP3;
	A(4,3)=AJ334;
	A(4,4)=ZZP4;
	    
    for i = 1:4
        for j = 1:4
            A(i,j) = A(j,i);
        end
    end


% *=*
% Number of operations : 211 '+' or '-', 278 '*' or '/'

%  QDP= 
% {QDP1, QDP2, QDP3, QDP4}

