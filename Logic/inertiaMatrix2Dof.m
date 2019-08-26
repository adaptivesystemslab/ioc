% (********************************************)
% (** SYMORO+ : SYmbolic MOdelling of RObots **)
% (**========================================**)
% (**      IRCCyN-ECN - 1, rue de la Noe     **)
% (**      B.P.92101                         **)
% (**      44321 Nantes cedex 3, FRANCE      **)
% (**      www.irccyn.ec-nantes.fr           **)
% (********************************************)


%    Name of file : E:\2DofArm_Symoro\Arm2Dof.inm




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


function A = inertiaMatrix2Dof(armModel)
    
    % Arm current state
    state = armModel.getState();
    [~, t2] = feval(@(x) x{:}, num2cell(state(1,:)));
    
    % Arm MDH parameters
    joint1 = armModel.joints(1);
    joint2 = armModel.joints(2);
    a2 = joint2.twist;
    d2 = joint2.length;    
    
    % Definition of first moments of inertia
    [MX2, MY2 , ~] = feval(@(x) x{:}, num2cell((joint2.com * joint2.mass)));
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
    
    % Function description:
	S2=sin(t2);
	C2=cos(t2);
	Sa2=sin(a2);
	Ca2=cos(a2);
	A212=Ca2.*S2;
	A222=C2.*Ca2;
	A312=S2.*Sa2;
	A322=C2.*Sa2;
	AS12=C2.*MX2 - MY2.*S2;
	AJ312=A312.*XX2 + A322.*XY2 + Ca2.*XZ2;
	AJ322=A312.*XY2 + A322.*YY2 + Ca2.*YZ2;
	AJ332=A312.*XZ2 + A322.*YZ2 + Ca2.*ZZ2;
	AJA332=A312.*AJ312 + A322.*AJ322 + AJ332.*Ca2;
	PAS332=-(AS12.*d2);
	ZZP1=AJA332 + d2.^2.*M2 - 2.*PAS332 + ZZ1;
	EB22=A222.*MX2 - A212.*MY2;
	NB32=AJ332 + d2.*EB22;
    
    A = zeros(2,2);
    
	A(1,1)=ZZP1;
	A(2,1)=NB32;
	A(2,2)=ZZ2;
    
    for i = 1:2
        for j = 1:2
            A(i,j) = A(j,i);
        end
    end
    
end

