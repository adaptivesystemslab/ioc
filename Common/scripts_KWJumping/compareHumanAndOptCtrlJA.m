% compare opt. ctrl. generated JAs to human data EKF JAs

partNum = '02';


% load JA data struct from certain participant
mydir = pwd;
idcs = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1); 
loadFilePath = [newdir '\results\JA_P' partNum '.mat'];
if(exist(loadFilePath,'file')~=2)
    disp(['JA_P' partNum '.mat NOT FOUND']);
else
    load(loadFilePath); % JA data struct, human data
    load('opt_ctrl_th_rec.mat'); % optimal control solution joint angles
    th_rec = rad2deg(th_rec);
    scale = 0.005/0.002; %human data sampling freq / opt ctrl timestep
    optCtrlJA = [-90-th_rec(1:scale:end,2), -th_rec(1:scale:end,3), -th_rec(1:scale:end,4), 180+th_rec(1:scale:end,5)];
    optCtrlLandFrame = find(optCtrlJA(:,1)==-90,1);
    
    % Comparing 4 Joint Angles: ankleDorsi, kneeExt, hipFlex (relative to
    %    spine), and shoulderElev
    humanJAeul = [];
    humanJAlie = [];
    
    % EUL MODEL JAs
    for i_targ = 1:numel(JA.targAlign_eul)
        humanDataLength = size(JA.targAlign_eul(i_targ).jump(1).data,1);
        % these arrays will have dim2 > 1, but only add data if not zeros
        sig_ankleDorsi = zeros(humanDataLength,1); 
        sig_kneeExt = zeros(humanDataLength,1);
        sig_hipFlex = zeros(humanDataLength,1);
        sig_shldrElev = zeros(humanDataLength,1);
        sig_count = 1;
        
        
        for i = 1:numel(JA.targAlign_eul(i_targ).jump)
            if(JA.missingData(i,i_targ)==0 && JA.noisyData_eul(i,i_targ)==0) % JAs are useable
                data = rad2deg(JA.targAlign_eul(i_targ).jump(i).data);
                
                % for JAs with left and right signals, take average of both
                sig_ankleDorsi(:,sig_count) = (data(:,25) + data(:,32))/2; 
                sig_kneeExt(:,sig_count) = (data(:,22) + data(:,30))/2;
                sig_hipFlex(:,sig_count) = (data(:,20) + data(:,27))/2 - data(:,7); %hipFlex w.r.t. spine
                sig_shldrElev(:,sig_count) = (data(:,10) + data(:,15))/2;
                sig_count = sig_count+1;
            end
        end
        
        humanJAeul.targ(i_targ).mean = [mean(sig_ankleDorsi,2), mean(sig_kneeExt,2), mean(sig_hipFlex,2), mean(sig_shldrElev,2)];
        humanJAeul.targ(i_targ).var = [var(sig_ankleDorsi,0,2), var(sig_kneeExt,0,2), var(sig_hipFlex,0,2), var(sig_shldrElev,0,2)];
    end
    
    
    % LIE MODEL JAs
    for i_targ = 1:numel(JA.targAlign_lie)
        humanDataLength = size(JA.targAlign_lie(i_targ).jump(1).data,1);
        % these arrays will have dim2 > 1, but only add data if not zeros
        sig_ankleDorsi = zeros(humanDataLength,1); 
        sig_kneeExt = zeros(humanDataLength,1);
        sig_hipFlex = zeros(humanDataLength,1);
        sig_shldrElev = zeros(humanDataLength,1);
        sig_count = 1;
        
        
        for i = 1:numel(JA.targAlign_lie(i_targ).jump)
            if(JA.missingData(i,i_targ)==0 && JA.noisyData_lie(i,i_targ)==0) % JAs are useable
                data = rad2deg(JA.targAlign_lie(i_targ).jump(i).data);
                
                % for JAs with left and right signals, take average of both
                sig_ankleDorsi(:,sig_count) = (data(:,25) + data(:,32))/2; 
                sig_kneeExt(:,sig_count) = (data(:,22) + data(:,30))/2;
                sig_hipFlex(:,sig_count) = (data(:,20) + data(:,27))/2 - data(:,7); %hipFlex w.r.t. spine
                sig_shldrElev(:,sig_count) = (data(:,10) + data(:,15))/2;
                sig_count = sig_count+1;
            end
        end
        
        humanJAlie.targ(i_targ).mean = [mean(sig_ankleDorsi,2), mean(sig_kneeExt,2), mean(sig_hipFlex,2), mean(sig_shldrElev,2)];
        humanJAlie.targ(i_targ).var = [var(sig_ankleDorsi,0,2), var(sig_kneeExt,0,2), var(sig_hipFlex,0,2), var(sig_shldrElev,0,2)];
    end
end



%% Plots
JAtitles = {'ankleDorsi','kneeExt','hipFlex','shldrElev'};

% Compare to EKF Data
for figNum = 1:3
    figure(figNum); clf;
    startFrame = JA.startFrame(JA.TOAlignJump(figNum),figNum);
    TOFrame = JA.TOFrame(JA.TOAlignJump(figNum),figNum);
    LandFrame = JA.LandFrame(JA.TOAlignJump(figNum),figNum); 
    endFrame = JA.endFrame(JA.TOAlignJump(figNum),figNum);
    
    x = 1:size(humanJAeul.targ(figNum).mean,1); x = x';
    
    for i = 1:4 %each joint
        subplot(4,1,i); hold on;
        
        y = humanJAeul.targ(figNum).mean(:,i);
        dy = humanJAeul.targ(figNum).var(:,i);
        h = fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.6 .6 1],'linestyle','none');
        set(h,'facealpha',.5)
        plot(humanJAeul.targ(figNum).mean(:,i),'b','LineWidth',2);
        plot([TOFrame,TOFrame],[-200,200],'b--');
        plot([LandFrame,LandFrame],[-200,200],'k--');
        
        
        optCtrlFrames = (LandFrame - optCtrlLandFrame + 1):LandFrame;
        plot(optCtrlFrames,optCtrlJA(1:optCtrlLandFrame,i),'r','LineWidth',2);
        
        TOFrame_sim = (LandFrame - optCtrlLandFrame + 1) + 95*(0.002/0.005);
        plot([TOFrame_sim,TOFrame_sim],[-200,200],'r--'); % manually aligned simulation TOFrame
    
        title(['Target ' num2str(figNum) ': ' JAtitles{i}]);
%         ylim([-180, 180]);
        axis([startFrame, endFrame, -180, 180]);
        grid on;
    end
end


% Compare to LGEKF Data
for figNumLie = 1:3
    figure(figNumLie+3); clf;
    startFrame = JA.startFrame(JA.TOAlignJump(figNumLie),figNumLie);
    TOFrame = JA.TOFrame(JA.TOAlignJump(figNumLie),figNumLie);
    LandFrame = JA.LandFrame(JA.TOAlignJump(figNumLie),figNumLie); 
    endFrame = JA.endFrame(JA.TOAlignJump(figNumLie),figNumLie);
    
    x = 1:size(humanJAlie.targ(figNumLie).mean,1); x = x';
    
    for i = 1:4 %each joint
        subplot(4,1,i); hold on;
        
        y = humanJAlie.targ(figNumLie).mean(:,i);
        dy = humanJAlie.targ(figNumLie).var(:,i);
        h = fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.6 .6 1],'linestyle','none');
        set(h,'facealpha',.5)
        plot(humanJAlie.targ(figNumLie).mean(:,i),'b','LineWidth',2);
        plot([TOFrame,TOFrame],[-200,200],'b--');
        plot([LandFrame,LandFrame],[-200,200],'k--');
        
        
        optCtrlFrames = (LandFrame - optCtrlLandFrame + 1):LandFrame;
        plot(optCtrlFrames,optCtrlJA(1:optCtrlLandFrame,i),'r','LineWidth',2);
        
        TOFrame_sim = (LandFrame - optCtrlLandFrame + 1) + 95*(0.002/0.005);
        plot([TOFrame_sim,TOFrame_sim],[-200,200],'r--'); % manually aligned simulation TOFrame
    
        title(['Target ' num2str(figNumLie) ': ' JAtitles{i}]);
%         ylim([-180, 180]);
        axis([startFrame, endFrame, -180, 180]);
        grid on;
    end
end

