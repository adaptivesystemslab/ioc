function [sigCorrected, shiftSig, shift2Sig, flipSig, shiftFlipSig, shiftFlip2Sig, OL_smooth] ...
            = cleanAndAlignData(sig,alignWindow,removeNoisyData)
        
numSig = size(sig,2);
k = numSig;
nsigma = 4;
sf = alignWindow(1); %startFrame
ef = alignWindow(2); %endFrame

shiftSig = zeros(1,numSig);         % shift +/- 180 degrees
shift2Sig = zeros(1,numSig);        % shift +/- 360 degrees
flipSig = zeros(1,numSig);          % negative of original
shiftFlipSig = zeros(1,numSig);     % shift +/- 180 degrees and make negative
shiftFlip2Sig = zeros(1,numSig);    % shift +/- 360 degrees and make negative


[~,OL] = hampel(mean(sig),numSig,6); % "OL" = outliers

% flip/shift each signal and check for new outliers, compare original
% outlier detection to new one; if modified signal matches better then take
% that signal as correct
for i = 1:numSig
    sig_shiftUp = sig(sf:ef,:);
    sig_shiftDn = sig(sf:ef,:);
    sig_shift2Up = sig(sf:ef,:);
    sig_shift2Dn = sig(sf:ef,:);
    sig_flip = sig(sf:ef,:);
    sig_shiftFlipUp = sig(sf:ef,:);
    sig_shiftFlipDn = sig(sf:ef,:);
    sig_shiftFlip2Up = sig(sf:ef,:);
    sig_shiftFlip2Dn = sig(sf:ef,:);
    
    % modify original signal
    sig_shiftUp(:,i) = sig_shiftUp(:,i) + pi;
    sig_shiftDn(:,i) = sig_shiftDn(:,i) - pi;
    sig_shift2Up(:,i) = sig_shift2Up(:,i) + 2*pi;
    sig_shift2Dn(:,i) = sig_shift2Dn(:,i) - 2*pi;
    sig_flip(:,i) = -sig_flip(:,i);
    sig_shiftFlipUp(:,i) = -(sig_shiftFlipUp(:,i) + pi);
    sig_shiftFlipDn(:,i) = -(sig_shiftFlipDn(:,i) - pi);
    sig_shiftFlip2Up(:,i) = -(sig_shiftFlip2Up(:,i) + 2*pi);
    sig_shiftFlip2Dn(:,i) = -(sig_shiftFlip2Dn(:,i) - 2*pi);
    
    % detect outliers with modified signal
    [~,OL_shiftUp] = hampel(mean(sig_shiftUp),k,nsigma);
    [~,OL_shiftDn] = hampel(mean(sig_shiftDn),k,nsigma);
    [~,OL_shift2Up] = hampel(mean(sig_shift2Up),k,nsigma);
    [~,OL_shift2Dn] = hampel(mean(sig_shift2Dn),k,nsigma);
    [~,OL_flip] = hampel(mean(sig_flip));%,k,nsigma);
    [~,OL_shiftFlipUp] = hampel(mean(sig_shiftFlipUp),k,nsigma);
    [~,OL_shiftFlipDn] = hampel(mean(sig_shiftFlipDn),k,nsigma);
    [~,OL_shiftFlip2Up] = hampel(mean(sig_shiftFlip2Up),k,nsigma);
    [~,OL_shiftFlip2Dn] = hampel(mean(sig_shiftFlip2Dn),k,nsigma);
    
    % detect changes in outlier detection
    OLDiff_shiftUp = OL - OL_shiftUp; % "OLDiff" = outlier difference
    OLDiff_shiftDn = OL - OL_shiftDn;
    OLDiff_shift2Up = OL - OL_shift2Up;
    OLDiff_shift2Dn = OL - OL_shift2Dn;
    OLDiff_flip = OL - OL_flip;
    OLDiff_shiftFlipUp = OL - OL_shiftFlipUp;
    OLDiff_shiftFlipDn = OL - OL_shiftFlipDn;
    OLDiff_shiftFlip2Up = OL - OL_shiftFlip2Up;
    OLDiff_shiftFlip2Dn = OL - OL_shiftFlip2Dn;
    
    
    % set appropriate flags if find better match with modified signal
    if(OLDiff_shiftUp(i)==1)        shiftSig(i) = 1;        end
    if(OLDiff_shiftDn(i)==1)        shiftSig(i) = -1;       end
    if(OLDiff_shift2Up(i)==1)       shift2Sig(i) = 1;       end
    if(OLDiff_shift2Dn(i)==1)       shift2Sig(i) = -1;      end
    if(OLDiff_flip(i)==1)           flipSig(i) = 1;         end
    if(OLDiff_shiftFlipUp(i)==1)    shiftFlipSig(i) = 1;    end
    if(OLDiff_shiftFlipDn(i)==1)    shiftFlipSig(i) = -1;   end
    if(OLDiff_shiftFlip2Up(i)==1)   shiftFlip2Sig(i) = 1;   end
    if(OLDiff_shiftFlip2Dn(i)==1)   shiftFlip2Sig(i) = -1;  end
end


% modify signals based on ___Sig flags
sigCorrected = sig;
for i = 1:numSig
    if(shift2Sig(i)~=0)
        sigCorrected(:,i) = sig(:,i) + 2*pi*shift2Sig(i);
    elseif(shiftSig(i)~=0)
        sigCorrected(:,i) = sig(:,i) + pi*shiftSig(i);
    elseif(shiftFlip2Sig(i)~=0)
        sigCorrected(:,i) = -(sig(:,i) + 2*pi*shiftFlip2Sig(i));
    elseif(shiftFlipSig(i)~=0)
        sigCorrected(:,i) = -(sig(:,i) + pi*shiftFlipSig(i));
    elseif(flipSig(i)~=0)
        sigCorrected(:,i) = -sig(:,i);
    end
end
% NOTE: ordered signal mods above from least to most likely to have false positives


% Removes signals with "flickering"/bad smoothness
if(removeNoisyData)
    sigSmooth = trapz( diff(diff(sigCorrected(sf:ef,:))).^2 );
    [~,OL_smooth] = hampel(sigSmooth,k,1000);
    
    for i = 1:numSig
        if(OL_smooth(i))
            sigCorrected(:,i) = zeros(size(sigCorrected,1),1);
        end
    end
else
    OL_smooth = zeros(1,numSig);
end


