% signal outlier testing

colorVec = [0,0.6,0; 0,0.9,0; 0.7,0.95,0; 1,0.9,0; 1,0.75,0; 1,0.55,0; 1,0.1,0.1; 0.95,0,0.6; ...
            0.4,0.1,0.7; 0,0.5,1; 0,0.9,0.95; 0,1,0.9];


A = [10 10 -5 -20 -40 -30 -10 25 45 60 40 10 20 15 0 -10 -30 -35 -20 10] + 20;
B = [A(2:end) 5];
JA = [A+4; 1.1*A; 0.9*B; 0.85*A; A; 1.15*B-5; ...
      0.93*A; 1.04*A+5; B; 1.08*B-8; 0.88*A; B-2]';

toFlip = [3];
toShiftUp = [2,7];
toShiftDn = [];

for i = 1:size(JA,2)
    if(ismember(i,toFlip))      JA(:,i) = -JA(:,i);         end
    if(ismember(i,toShiftUp))   JA(:,i) = JA(:,i) + 180;    end
    if(ismember(i,toShiftDn))   JA(:,i) = JA(:,i) - 180;    end
end



figure(1); clf; hold on;
for i = 1:size(JA,2)
    plot(JA(:,i),'LineWidth',2,'Color',colorVec(i,:));
end
grid on;


% test for outliers with original data
[~,outliers] = hampel(mean(JA));

% test if flipping data results in better match
flipSig = zeros(1,size(JA,2));
for i = 1:size(JA,2)
    JA_flip = JA;
    JA_flip(:,i) = -JA_flip(:,i);
%     figure(2); clf; hold on;
%     for j = 1:size(JA_flip,2)
%         plot(JA_flip(:,j),'LineWidth',2,'Color',colorVec(j,:));
%     end
%     grid on;
    
    [~,outliers_flip] = hampel(mean(JA_flip));
    outliers_diff = outliers - outliers_flip;
    if(outliers_diff(i)==1)
        flipSig(i) = 1;
    end
end


% test is shifting +/- 180 results in better match
shiftSig = zeros(1,size(JA,2));
for i = 1:size(JA,2)
    JA_shiftUp = JA;
    JA_shiftDn = JA;
    JA_shiftUp(:,i) = JA_shiftUp(:,i) + 180;
    JA_shiftDn(:,i) = JA_shiftDn(:,i) - 180;
%     figure(3); clf; hold on;
%     for j = 1:size(JA_shiftUp,2)
%         plot(JA_shiftUp(:,j),'LineWidth',2,'Color',colorVec(j,:));
%     end
%     grid on;
%     figure(4); clf; hold on;
%     for j = 1:size(JA_shiftDn,2)
%         plot(JA_shiftDn(:,j),'LineWidth',2,'Color',colorVec(j,:));
%     end
%     grid on;

    [~,outliers_shiftUp] = hampel(mean(JA_shiftUp));
    [~,outliers_shiftDn] = hampel(mean(JA_shiftDn));
    
    outliers_diffUp = outliers - outliers_shiftUp;
    outliers_diffDn = outliers - outliers_shiftDn;
    if(outliers_diffUp(i)==1)
        shiftSig(i) = -1; %the signal IS 180 deg down, needs to be shifted up
    end
    if(outliers_diffDn(i)==1)
        shiftSig(i) = 1; %the signal IS 180 deg up, needs to be shifted down
    end
end

flipSig
shiftSig