function [observation, fhmmIntermediate, sysparam] = fgInitialize
   sysparam = sysparamSet_threeTemplate;

observation.subject = 0;
observation.obsTime = [];
observation.obsData = [];
observation.obsVelo = [];
observation.currDataLength = 0;
observation.currDataChecked = 0;
observation.posIndex = 1:5;
observation.posIndexAccel = 6:11;

fhmmIntermediate.W = 3;
fhmmIntermediate.zi_ang = [];
fhmmIntermediate.zi_velo = [];