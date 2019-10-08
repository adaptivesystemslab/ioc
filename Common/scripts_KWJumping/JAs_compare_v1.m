% load JA data from multiple participants, compare mean and variance of
% perfect (P, P*) jumps to other (B,SB,SF,F) jumps


partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};
% partsToRun = {'04'};

target_percents = zeros(numel(partsToRun),2);


for partNum = 1:numel(partsToRun)
    % load JA data struct from certain participant
    mydir = pwd;
    idcs = strfind(mydir,'\');
    newdir = mydir(1:idcs(end)-1); 
    loadFilePath = [newdir '\results\JA_P' partsToRun{partNum} '.mat'];
    if(exist(loadFilePath,'file')~=2)
        disp(['JA_P' partNum '.mat NOT FOUND']);
    else
        load(loadFilePath); % JA data struct, human data
        
        target_percents(partNum,1) = JA.targetLengths(6);
        target_percents(partNum,2) = (JA.targetLengths(4)/0.55)*0.85;
        
        
%         for i_targ = 1:numel(JA.targAlign_eul)
%             
% %             startFrame(i,i_targ) = JA.startFrame(JA.LandAlignJump(i_targ),i_targ);
% %             TOFrame(i,i_targ) = JA.TOFrame(JA.LandAlignJump(i_targ),i_targ);
% %             LandFrame(i,i_targ) = JA.LandFrame(JA.LandAlignJump(i_targ),i_targ); 
% %             endFrame(i,i_targ) = JA.endFrame(JA.LandAlignJump(i_targ),i_targ);
%             
%             dataLength = size(JA.targAlign_eul(i_targ).jump(1).data,1);
%             % these arrays will have dim2 > 1, but only add data if not zeros
%             sig_ankleDorsi = zeros(dataLength,1); 
%             sig_kneeExt = zeros(dataLength,1);
%             sig_hipFlex = zeros(dataLength,1);
%             sig_shldrElev = zeros(dataLength,1);
%             sig_count = 1;
%             
%             sig_ankleDorsi_P = zeros(dataLength,1); 
%             sig_kneeExt_P = zeros(dataLength,1);
%             sig_hipFlex_P = zeros(dataLength,1);
%             sig_shldrElev_P = zeros(dataLength,1);
%             sig_count_P = 1;
%             
% 
%             for i = 1:numel(JA.targAlign_eul(i_targ).jump)
%                 if(JA.missingData(i,i_targ)==0 && JA.noisyData_eul(i,i_targ)==0) % JAs are useable
%                     data = rad2deg(JA.targAlign_eul(i_targ).jump(i).data);
% 
%                     % for JAs with left and right signals, take average of both
%                     sig_ankleDorsi(:,sig_count) = (data(:,25) + data(:,32))/2; 
%                     sig_kneeExt(:,sig_count) = (data(:,22) + data(:,30))/2;
%                     sig_hipFlex(:,sig_count) = (data(:,20) + data(:,27))/2;
%                     sig_shldrElev(:,sig_count) = (data(:,10) + data(:,15))/2;
%                     sig_count = sig_count+1;
%                 end
%                 
%                 % find jumps that are 'P' or 'P*'
%                 if(strcmp(jumpGrades{i_targ,i},'P') || strcmp(jumpGrades{i_targ,i},'P*'))
%                     data = rad2deg(JA.targAlign_eul(i_targ).jump(i).data);
% 
%                     % for JAs with left and right signals, take average of both
%                     sig_ankleDorsi_P(:,sig_count_P) = (data(:,25) + data(:,32))/2; 
%                     sig_kneeExt_P(:,sig_count_P) = (data(:,22) + data(:,30))/2;
%                     sig_hipFlex_P(:,sig_count_P) = (data(:,20) + data(:,27))/2;
%                     sig_shldrElev_P(:,sig_count_P) = (data(:,10) + data(:,15))/2;
%                     sig_count_P = sig_count_P+1;
%                 end
%                     
%             end
%             
%             humanJAeul.targ(i_targ).mean = [mean(sig_ankleDorsi,2), mean(sig_kneeExt,2), mean(sig_hipFlex,2), mean(sig_shldrElev,2)];
%             humanJAeul.targ(i_targ).var = [var(sig_ankleDorsi,0,2), var(sig_kneeExt,0,2), var(sig_hipFlex,0,2), var(sig_shldrElev,0,2)];
%             
%             humanJAeul.targ(i_targ).mean_P = [mean(sig_ankleDorsi_P,2), mean(sig_kneeExt_P,2), mean(sig_hipFlex_P,2), mean(sig_shldrElev_P,2)];
%             humanJAeul.targ(i_targ).var_P = [var(sig_ankleDorsi_P,0,2), var(sig_kneeExt_P,0,2), var(sig_hipFlex_P,0,2), var(sig_shldrElev_P,0,2)];
%         end
%         
%         
%         
%         for figNum = 1:3
%             figure(figNum); clf;
%             x = 1:size(humanJAeul.targ(figNum).mean,1); x = x';
%             
%             startFrame = JA.startFrame(JA.LandAlignJump(figNum),figNum);
%             TOFrame = JA.TOFrame(JA.LandAlignJump(figNum),figNum);
%             LandFrame = JA.LandFrame(JA.LandAlignJump(figNum),figNum); 
%             endFrame = JA.endFrame(JA.LandAlignJump(figNum),figNum);
%             
%             JAtitles = {'ankleDorsi','kneeExt','hipFlex','shldrElev'};
%             
%             for i = 1:4 %each joint
% %                 subplot(4,1,i); hold on;
%                 % mean and variance of ALL jumps
%                 subplot(4,2,2*i-1); hold on;
% 
%                 y = humanJAeul.targ(figNum).mean(:,i);
%                 dy = humanJAeul.targ(figNum).var(:,i);
%                 h = fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.6 .6 1],'linestyle','none');
%                 set(h,'facealpha',.5)
%                 plot(humanJAeul.targ(figNum).mean(:,i),'b','LineWidth',2);
%                 plot([TOFrame,TOFrame],[-200,200],'g--');
%                 plot([LandFrame,LandFrame],[-200,200],'r--');
%                 
%                 title(['Target ' num2str(figNum) ': ' JAtitles{i} ': ALL']);
%                 axis([startFrame, endFrame, -180, 180]);
%                 grid on;
%                 
%                 subplot(4,2,2*i); hold on;
% 
%                 y = humanJAeul.targ(figNum).mean_P(:,i);
%                 dy = humanJAeul.targ(figNum).var_P(:,i);
%                 h = fill([x;flipud(x)],[y-dy;flipud(y+dy)],[1 .6 .6],'linestyle','none');
%                 set(h,'facealpha',.5)
%                 plot(humanJAeul.targ(figNum).mean_P(:,i),'r','LineWidth',2);
%                 plot([TOFrame,TOFrame],[-200,200],'g--');
%                 plot([LandFrame,LandFrame],[-200,200],'k--');
%                 
%                 title(['Target ' num2str(figNum) ': ' JAtitles{i} ': On Target']);
%                 axis([startFrame, endFrame, -180, 180]);
%                 grid on;
%         
%             end
%         end
        
        
    end
end

        

