function [shiftValues, shift_rec] = JAShiftValues(JA)
% gets shift values from hip flex, knee extend, and ankle dorsiflex
% signals (6 total), then take average of these shift values

targCurr = JA.targAlign;
if(numel(JA.jointNames) == 35) % using shldrPrism model
    alignJoints = [22,29,25,32,27,34]; %hip flex, knee extend, ankle dorsiflex
elseif(numel(JA.jointNames) == 14) % using 2D sagittal model
    alignJoints = [9,10,11,12,13,14];
else %using fixed shldr position model
    alignJoints = [20,27,23,30,25,32]; %hip flex, knee extend, ankle dorsiflex
end

shiftValues = zeros(numel(targCurr(1).jump),numel(targCurr));
shift_rec = zeros(numel(alignJoints),numel(targCurr(1).jump),numel(targCurr));


for i_targ = 1:numel(targCurr)
    idx_LandFrame = JA.LandAlignJump(i_targ);

    for j = 1:numel(alignJoints)
        sig = zeros(size(targCurr(i_targ).jump(1).data,1),numel(targCurr(i_targ).jump));
        for i = 1:numel(targCurr(i_targ).jump)
            sig(:,i) = targCurr(i_targ).jump(i).data(:,alignJoints(j));
        end

        sig_shift = sig;
        %NOTE: only perform shift based on frames between TO and landing
        frame_st = JA.TOFrame(idx_LandFrame,i_targ) - 100;
        frame_end = JA.LandFrame(idx_LandFrame,i_targ);% + 100;
        maxshift = 50; %max samples shifted by w.r.t. reference
        
        refSig = sig_shift(frame_st:frame_end,end);
        for i = 1:size(sig_shift,2)-1 %use 12th (last) jump as reference
            currSig = sig_shift(frame_st:frame_end,i);
            shift = shiftMatchJAs(currSig,refSig,maxshift);
            shift_rec(j,i,i_targ) = shift;
            
%             figure(4); clf; hold on;
%             plot(sig_shift(frame_st:frame_end,end),'k');
%             plot(tmp,'r');
%             legend('Reference','Sig. To Shift');
            
        end
    end
    
    % average shift values for each set of JAs
    for i = 1:size(shiftValues,1)
        shiftValues(i,i_targ) = round(mean(shift_rec(:,i,i_targ)));
    end
    
end


