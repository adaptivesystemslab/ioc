% load JA data from multiple participants, train LSTM on joint angle trajectories 
clear;

% NN LSTM parameters
inputSize = 35; % number of joints
classType = 'jumpGrade3'; % jumpGrade2, jumpGrade3, toeDistToTarg (others?)
testSetType = 'random'; % random, singlePart (others?)
singlePartTestSetNumber = '08';

numHiddenUnits = 100; % hidden units of LSTM, EXPERIMENT WITH THIS

partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};
% partsToRun = {'04'};


load('JA_NN_data.mat'); % X_all, Info_all




%% Create Y_all based on desired classification categories
Y_all = zeros(size(Info_all,1),1);

% {'B','SB','P','P*','SF','F'} is numbering order of "jumpGradeNumber" in Info_all
switch classType
    case 'jumpGrade2'
        numClasses = 2; % too short, on target, too long
        for i = 1:size(Info_all,1)
            if(Info_all(i,4)==3 || Info_all(i,4)==4) % P(3) or P*(4)
                Y_all(i) = 2; % perfect jump
            else
                Y_all(i) = 1; % jump too short or too long
            end
        end
        
    case 'jumpGrade3'
        numClasses = 3; % too short, on target, too long
        for i = 1:size(Info_all,1)
            if(Info_all(i,4)<3) % B(1) or SB(2)
                Y_all(i) = 1; % jump too short
            elseif(Info_all(i,4)>4) % SF(5) or F(6)
                Y_all(i) = 3; % jump too long
            else
                Y_all(i) = 2; % perfect jump
            end
        end
        
    case 'toeDistToTarg'
        shortJumpThresh = -0.02;
        longJumpThresh = 0.09;
        
        numClasses = 3;
        for i = 1:size(Info_all,1)
            if(Info_all(i,5) < shortJumpThresh)
                Y_all(i) = 1; % jump too short
            elseif(Info_all(i,4) > longJumpThresh)
                Y_all(i) = 3; % jump too long
            else
                Y_all(i) = 2; % perfect jump
            end
        end
end
Y_all = categorical(Y_all);


%% Split into training and test sets
switch testSetType
    case 'random'
        percentTest = 0.05;
        numObs = 21*36; % numParts x numJumps
        ordering = randperm(numObs);
        testObs = ordering(1:round(numObs*percentTest));
        trainObs = ordering(round(numObs*percentTest)+1:end);
    case 'singlePart'
        % something
end

X_train = X_all(trainObs);
X_test = X_all(testObs);
Y_train = Y_all(trainObs);
Y_test = Y_all(testObs);



%% Train LSTM
numObservations = numel(X_train);
% for i=1:numObservations
%     sequence = X_train{i};
%     sequenceLengths(i) = size(sequence,2);
% end

% [sequenceLengths,idx] = sort(sequenceLengths);
% X_train = X_train(idx);
% Y_train = Y_train(idx);
% 
% figure(1);
% bar(sequenceLengths);
% ylim([0 30]);
% xlabel('Sequence');
% ylabel('Length');
% title('Sorted Data');

layers = [ ...
    sequenceInputLayer(inputSize)
    bilstmLayer(numHiddenUnits,'OutputMode','sequence')
    bilstmLayer(numHiddenUnits,'OutputMode','last')
%     fullyConnect000000000000000000edLayer(numHiddenUnits)
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer]

maxEpochs = 100;

options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'GradientThreshold',1, ...
    'MaxEpochs',maxEpochs, ...
    'SequenceLength','longest', ...
    'Verbose',0, ...
    'Plots','training-progress');


jumpClassNN = trainNetwork(X_train,Y_train,layers,options);



%% Test LSTM
Y_pred = classify(jumpClassNN,X_test, ...
    'SequenceLength','longest');

figure(2); hold on;
plot(Y_test,'k','LineWidth',2);
plot(Y_pred,'r');

