% (********************************************)
% (** SYMORO+ : SYmbolic MOdelling of RObots **)
% (**========================================**)
% (**      IRCCyN-ECN - 1, rue de la Noe     **)
% (**      B.P.92101                         **)
% (**      44321 Nantes cedex 3, FRANCE      **)
% (**      www.irccyn.ec-nantes.fr           **)
% (********************************************)


%    Name of file : E:\simpleArm\simpleArm.inm




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

%                       Inertia matrix
% Equations

% Declaration of the function
function A = inertiaMatrix3Dof(armModel)

    % IMPORTANT: Consider transposing the matrix as shown in Liz's code

    % Arm state
    state = armModel.getState();
       
    % Arm current state
    [~, t2, t3] = feval(@(x) x{:}, num2cell(state(1,:)));
    
    joint2 = armModel.joints(2);
    joint3 = armModel.joints(3);
    
    % Arm MDH parameters
    a2 = joint2.twist;
    a3 = joint3.twist;
    d3 = joint3.length;
        
    % First moment of inertia
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
	AS13=C3.*MX3 - MY3.*S3;
	AS23=A213.*MX3 + A223.*MY3 - MZ3.*Sa3;
	AS33=A313.*MX3 + A323.*MY3 + Ca3.*MZ3;
	AJ113=C3.*XX3 - S3.*XY3;
	AJ123=C3.*XY3 - S3.*YY3;
	AJ133=C3.*XZ3 - S3.*YZ3;
	AJ213=A213.*XX3 + A223.*XY3 - Sa3.*XZ3;
	AJ223=A213.*XY3 + A223.*YY3 - Sa3.*YZ3;
	AJ233=A213.*XZ3 + A223.*YZ3 - Sa3.*ZZ3;
	AJ313=A313.*XX3 + A323.*XY3 + Ca3.*XZ3;
	AJ323=A313.*XY3 + A323.*YY3 + Ca3.*YZ3;
	AJ333=A313.*XZ3 + A323.*YZ3 + Ca3.*ZZ3;
	AJA113=AJ113.*C3 - AJ123.*S3;
	AJA213=AJ213.*C3 - AJ223.*S3;
	AJA313=AJ313.*C3 - AJ323.*S3;
	AJA223=A213.*AJ213 + A223.*AJ223 - AJ233.*Sa3;
	AJA323=A213.*AJ313 + A223.*AJ323 - AJ333.*Sa3;
	AJA333=A313.*AJ313 + A323.*AJ323 + AJ333.*Ca3;
	PAS213=AS23.*d3;
	PAS223=-(AS13.*d3);
	PAS313=AS33.*d3;
	PAS333=-(AS13.*d3);
	XXP2=AJA113 + XX2;
	XYP2=AJA213 - PAS213 + XY2;
	XZP2=AJA313 - PAS313 + XZ2;
	YYP2=AJA223 + d3.^2.*M3 - 2.*PAS223 + YY2;
	YZP2=AJA323 + YZ2;
	ZZP2=AJA333 + d3.^2.*M3 - 2.*PAS333 + ZZ2;
	AJ312=A312.*XXP2 + A322.*(AJA213 - PAS213 + XY2) + Ca2.*(AJA313 - PAS313 + XZ2);
	AJ322=A312.*XYP2 + A322.*YYP2 + Ca2.*(AJA323 + YZ2);
	AJ332=A312.*XZP2 + A322.*YZP2 + Ca2.*ZZP2;
	AJA332=A312.*AJ312 + A322.*AJ322 + AJ332.*Ca2;
	EC22=A223.*MX3 - A213.*MY3;
	EC32=A323.*MX3 - A313.*MY3;
	NC22=AJ233 - d3.*EC32;
	NC32=AJ333 + d3.*EC22;
	NC33=A312.*AJ133 + A322.*NC22 + Ca2.*NC32;
	
    A = [0 0 0;
        0 0 0;
        0 0 0];
    
	A(1,1)=AJA332;
	A(2,1)=AJ332;
	A(3,1)=NC33;
	A(2,2)=ZZP2;
	A(3,2)=NC32;
	A(3,3)=ZZ3;
    
    for i = 1:3
        for j = 1:3
            A(i,j) = A(j,i);
        end
    end


% *=*
% Number of operations : 211 '+' or '-', 278 '*' or '/'

%  QDP= 
% {QDP1, QDP2, QDP3, QDP4}

