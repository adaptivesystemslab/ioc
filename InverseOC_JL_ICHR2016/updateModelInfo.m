function model = updateModelInfo(model, param)

for i = 1:length(model.transforms)
    linkName = model.transforms(i).name;
    linkNameSplit = strsplit(linkName, 'length');
    
    if length(linkNameSplit) == 2
        linkNameInd = str2num(linkNameSplit{2});
        
%         model.bodies(linkNameInd).fX = [];
        model.bodies(linkNameInd).com = param.link_com{linkNameInd};
        model.bodies(linkNameInd).I   = param.link_i{linkNameInd};
        model.bodies(linkNameInd).m   = param.link_mass(linkNameInd);
%         model.bodies(linkNameInd).a = [];
%         model.bodies(linkNameInd).c = [];
%         model.bodies(linkNameInd).f = [];
%         model.bodies(linkNameInd).i = [];
%         model.bodies(linkNameInd).iA = [];
%         model.bodies(linkNameInd).pA = [];
%         model.bodies(linkNameInd).t = applyTranslation1(param.link_dh_d(linkNameInd), param.link_dh_r(linkNameInd));
%         model.bodies(linkNameInd).v = [];
%         model.bodies(linkNameInd).x = [];

        model.transforms(i).t = applyTranslation1(param.link_dh_d(linkNameInd), param.link_dh_r(linkNameInd));
%         model.transforms(i).x = [];

%         newL = applyTranslation(param.L{linkNameInd+1});
       
%         newL = applyTranslation2( ...
%             param.DH_d(linkNameInd), param.DH_r(linkNameInd), param.DH_alpha(linkNameInd));
%         newL = applyTranslation3( ...
%             param.DH_d(linkNameInd), param.DH_r(linkNameInd), param.DH_alpha(linkNameInd));
        
        

    end
end

end

function T = applyTranslation(d)
    tToApply = eye(4);
	tToApply(1:3, 4) = [d 0 0];
	T = tToApply;
end

function T = applyTranslation1(d, r)
    tToApply = eye(4);
	tToApply(1:3, 4) = [d r 0];
	T = tToApply;
end

function T = applyTranslation2(d, r, alpha)
% 	tToApply = eye(4);
% 	tToApply(1:3, 4) = [0 param.DH_d param.DH_r];
    tToApply = dhMatrix(0, d, r, alpha);
	T = tToApply;
end

function Tr=applyTranslation3(d, r, Alpha)

b=0;
% d=u(4);
% r=u(5);
Gamma=0;
% Alpha=u(7);
Theta=0;

Ct=cos(Theta);
St=sin(Theta);

Ca=cos(Alpha);
Sa=sin(Alpha);

Cg=cos(Gamma);
Sg=sin(Gamma);


Tr=[Cg*Ct-Sg*Ca*St     -Cg*St-Sg*Ca*Ct      Sg*Sa       d*Cg+r*Sg*Sa;
      Sg*Ct+Cg*Ca*St    -Sg*St+Cg*Ca*Ct     -Cg*Sa      d*Sg-r*Cg*Sa;
      Sa*St                      Sa*Ct                       Ca            r*Ca+b
        0                           0                              0              1];
end