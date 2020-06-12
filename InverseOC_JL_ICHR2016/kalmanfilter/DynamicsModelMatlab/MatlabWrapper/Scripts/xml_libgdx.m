%This is a script to match XML model to the ARS libGDX model,
%We add a bunch of bodies such that they align with LIBGDX bodies 


tic;
for i=1:1000
    mdl.position = rand(87,1);
    mdl.forwardPosition();
end
toc/1000