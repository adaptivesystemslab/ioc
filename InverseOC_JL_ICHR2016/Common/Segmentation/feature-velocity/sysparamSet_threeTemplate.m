function sysparam = sysparamSet_threeTemplate(sysparam)
    
sysparam.verbose = 1;
    
    sysparam.filter.filterType = 'low';
    sysparam.filter.typeTemplate = 'doublePass'; % noFilter, singlePass, doublePass
    sysparam.filter.typeObservation = 'doublePass';
    sysparam.filter.bwFreq = 0.04; % 0.04 is original, 0.03 is better fit it seems
    sysparam.filter.bwOrder = 2;
    
%     sysparam.filter.filterTemplate = ...
%         filter_butterworth(sysparam.filter.filterType, sysparam.filter.typeTemplate, ...
%         sysparam.filter.bwFreq, sysparam.filter.bwOrder, 0);
%     sysparam.filter.filterObservation = ...
%         filter_butterworth(sysparam.filter.filterType, sysparam.filter.typeObservation, ...
%         sysparam.filter.bwFreq, sysparam.filter.bwOrder, 0);
    
    sysparam.filter.filterTemplate = [];
    sysparam.filter.filterObservation = [];
    

    sysparam.template.LLThresholdModifierPositive = 6; %  currently unused - set so it's always negative, so this would drop the LLThres
    sysparam.template.LLThresholdModifierNegative = 1/2; % currently unused

    sysparam.template.autoLLThreshold = 1; % if 0, the LL threshold is set to some arb number, otherwise, use LOOCV
    sysparam.template.autoLLThresholdShift = -45000*10; % was -15000 on 2013/02 threshold is offset, then multiplied (was at -6000 or -10000 before)
    sysparam.template.autoLLThresholdMultiplier = 1;
    % -8000 normally, -200000 on single template, -15000 (so far) for modified filter
%     sysparam.template.manualLLThreshold = -8000; % if above value is 0, what should the arb LL be?
    % -15000 normally, -250000 on single template, -20000 (so far) for modified filter, -2e4 for twopeak
    sysparam.template.LLThresholdAbs = -50000*10; % if a given motion has a LL below this value, it'll be rejected, regardless of how high the actual LL threshold is

    sysparam.feature.attenuateTime = 5; % every 'x' seconds, the 'max/min' will drop by...
    sysparam.feature.attenuatePercent = 0.8; % by this much
    sysparam.feature.attenuateRate = 1; % every this much seconds seconds
    
    sysparam.feature.angNormalize = 1; % if yes, then the first timestamp of the ang template will be normalized to the incoming data
    sysparam.feature.angRange = 0.5; % if the ang features are within this range, pass it

    % Template matching
    sysparam.segment.Wzvc = 3; % window width examined for ZVC
    sysparam.segment.NspareZVC = 0; % how much ZVCs will we look ahead for the terminating point?
    sysparam.segment.recycleTime = 0.1; % offset for 'recycletime' for next cycle. set to zero if don't want lookback

    sysparam.segment.overlapExemplar = 0; % if 1, it suggests that motions could simultanously be two different motions

    sysparam.feature.lengthMultiplier = 0.5; % time bewteen segments
    sysparam.feature.peaktopeakMultiplier = 0.2; % time between peaks
    sysparam.template.zvcThreshold = 0.01; % in cases where velocity crosses really close to zero but doesn't actually cross, how close should it be before a ZVC is declared? (0.1 is good)
    sysparam.feature.smallMagRejectMultiplier = 0.25; % was 0.25
    sysparam.feature.velocityMultiplierThreshold = 0.25; % 1/1.5=0.67, 1/1.8=0.56 is current best, 1.5 and 2 works well (not used?)
    sysparam.feature.lookaheadZVC = 1.00; %0.5; % was 1.25 earlier

    sysparam.feature.windowThresholdMin = 0.25; % 0.55 bounds on the size that the template windows can be when inserted into HMM
    sysparam.feature.windowThresholdMax = 4;

%     sysparam.rejectRepeat = 0;
end