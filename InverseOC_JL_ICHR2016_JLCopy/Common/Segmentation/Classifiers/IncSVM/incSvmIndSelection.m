function [ind0loc, ind1loc, out] = incSvmIndSelection(obj, out_inc, decisionValue_inc, currentInput, extraInput, ind_range, retrainFlagSideClass)

ind0locArray = {};
ind1locArray = {};
out = [];

for ind_incsvmmode = 1:length(obj.incsvmMode)
    incMode = obj.incsvmMode{ind_incsvmmode};
    incParam = obj.incsvmParam{ind_incsvmmode};
    
    switch incMode
        case 'highProb'
            % high prob      
            out = out_inc;
            prob_inc(ind_range) = sigmoid_predict(obj, decisionValue_inc(ind_range), obj.SVMmodel_inc.ProbA, obj.SVMmodel_inc.ProbB);
            highProb = find(abs(prob_inc(ind_range) - 0.5) > incParam.probThreshold);
            
            ind0locArray{ind_incsvmmode} = ind_range(out_inc(ind_range(highProb)) == 0);
            ind1locArray{ind_incsvmmode} = ind_range(out_inc(ind_range(highProb)) == 1);
            
        case {'fsmImu', 'fsmEkf', 'fsmPp'}
            switch incMode
                case 'fsmImu'
                    % fsm (on imu)
                    auxInput = extraInput.auxDataImu(ind_range, :);
                    [out, fsmProb] = obj.fsmImu.classify(auxInput);
                    
                case 'fsmEkf'
                    % fsm (on ekf)
                    auxInput = extraInput.auxDataEkf(ind_range, :);
                    [out, fsmProb] = obj.fsmEkf.classify(auxInput);
                    
                case 'fsmPp'
                    % fsm (on pp)
                    auxInput = extraInput.auxDataPp(ind_range, :);
                    [out, fsmProb] = obj.fsmPp.classify(auxInput);
            end
            
            switch incParam.labelSource
                case 'fsmLabels'
                    % use the labels determined by the
                    % FSM. feed into fsm and get the labels
                    
                    % highly certain non-seg entries
                    fsm0 = find(fsmProb(:, 1) == 1);
                    fsm1 = find(fsmProb(:, 2) == 1);
                    
                case 'svmfsmLabels'
                    % use the point when both
                    % the incsvm and the fsm agrees
                    % feed into fsm and get the labels
                    
                    % highly certain non-seg entries
                    fsm0 = intersect(find(out_inc(ind_range) == 0), find(fsmProb(:, 1) == 1));
                    fsm1 = intersect(find(out_inc(ind_range) == 1), find(fsmProb(:, 2) == 1));
                    
                case 'svmLabels'
                    % highly certain non-seg entries
                    highProb = union(fsmProb(:, 1) == 1, fsmProb(:, 2) == 1);
                    
                    fsm0 = intersect(find(out_inc(ind_range) == 0), highProb);
                    fsm1 = intersect(find(out_inc(ind_range) == 1), highProb);
            end
            
            ind0locArray{ind_incsvmmode} = ind_range(fsm0);
            ind1locArray{ind_incsvmmode} = ind_range(fsm1);
            
        case 'fgModel'
            out = out_inc;
            
%             startInd = obj.fgStateMachine.startInd + extraInput.settings.dimStrack;
%             endInd = ind_range(end) + extraInput.settings.dimStrack; % offset, since we're using the main data and not the testdata
      
            [startInd, endInd, observationNew] = temporallyAlignedDataLoad(ind_range, extraInput);
            
            if ind_range(1) == 1
                [observation, fhmmIntermediate, sysparam] = fgInitialize;
            else
                observation = obj.fgStateMachine.observation;
                fhmmIntermediate = obj.fgStateMachine.fhmmIntermediate;
                sysparam = obj.fgStateMachine.sysparam;
            end
            
            [observation, fhmmIntermediate] = updateDataStruct(observation, observationNew, sysparam, fhmmIntermediate);
            
            [segmentInfo, observation, fhmmIntermediate] = ...
                segmentModifyZvcMethod_runModel(observation, obj.fgModel, sysparam, fhmmIntermediate);
            
            obj.fgStateMachine.observation = observation;
            obj.fgStateMachine.fhmmIntermediate = fhmmIntermediate;
            obj.fgStateMachine.sysparam = sysparam;
            
            if ~isequal(obj.fgStateMachine.previousSegmentInfo, segmentInfo)
                % new entries in the segmentTemp. copy them out
                
%                 newStartInd = setxor(segmentTemp.indStart, obj.fgStateMachine.previousSegmentInfo.startTimeInd);
%                 newEndInd = setxor(segmentTemp.indEnd, obj.fgStateMachine.previousSegmentInfo.endTimeInd);

%                 subObj.label_notSegment = extraInput.label_notSegment;
%                 subObj.label_segment = extraInput.label_segment;
%                 subObj.label_unknown = extraInput.label_unknown;
% 
%                 subObj.settings = extraInput.settings;
                settings = extraInput.settings.segmentPointWindowSettings;

                data.time = observation.obsTime;
                data.jointAngles = observation.obsData';
                data.jointVelo = observation.obsVelo';

                [segLabel, segName, segTrainingInclude, segTestingInclude] = ...
                    segmentPointWindowExpand(extraInput, data, segmentInfo, length(data.time), settings);

                % points to include
                incStartInd = obj.fgStateMachine.startInd;
                incEndInd = segmentInfo.endTimeInd(end);
                inRangeToUse = incStartInd:incEndInd;
                segLabelToUse = segTrainingInclude(inRangeToUse);
                
                obj.fgStateMachine.startInd = incEndInd; % next start of the array to feed data in
                obj.fgStateMachine.previousSegmentInfo = segmentInfo;

                ind0locArray{ind_incsvmmode} = inRangeToUse(segLabelToUse == 0);
                ind1locArray{ind_incsvmmode} = inRangeToUse(segLabelToUse == 1);
                
                if 0
                    figure;
                    plot(data.jointVelo);
                    hold on
                    plot(segTrainingInclude, 'rx');
                    plot(inRangeToUse, segLabelToUse, 'bo');
                end
            else
                ind0locArray{ind_incsvmmode} = [];
                ind1locArray{ind_incsvmmode} = [];
            end
            
        case 'groundTruth'
             switch incParam.labelSource
                case 'testData'
                    out = extraInput.testingLabel(ind_range);
                    
                    fsm0 = find(out == 0);
                    fsm1 = find(out == 1);
             end

            ind0locArray{ind_incsvmmode} = ind_range(fsm0);
            ind1locArray{ind_incsvmmode} = ind_range(fsm1);
            
        case 'clusterNormal'
            
        case 'clusterDR'      
            switch obj.clusterDR.dataFormatName
                case 'PC'
%                     inputCluster = extraInput.testingData(ind_range, :);
                    testingDataParse = clusterDataGatewayTesting(extraInput, obj.clusterDR, ind_range);
                    
                case 'PP'
                    % PP uses normal q/dq and thus we need to add the
                    % dimstack to line up the data 
                    
%                     inputCluster = extraInput.data(ind_range, :);
                    testingDataParse = clusterDataGatewayTesting(extraInput, obj.clusterDR, [ind_range]+extraInput.settings.dimStrack);
            end
            
            [out, prob] = obj.clusterDR.classify(testingDataParse.input, extraInput.testingLabel(ind_range));
            
            fsm0 = find(out == 0);
            fsm1 = find(out == 1);
            ind0locArray{ind_incsvmmode} = ind_range(fsm0);
            ind1locArray{ind_incsvmmode} = ind_range(fsm1);
            
           if incParam.updateMachine && retrainFlagSideClass
               % how many of the old data points to keep
               percentageToKeep = incParam.percentageOriginalDataToKeepSubsequentReplacement;
               
               switch obj.clusterDR.dataFormatName
                   case {'PC', 'CR'}
                       inputRetrain = extraInput.testingData(ind_range, :);
%                        input0 = currentInput(out == 0, :);
%                        input1 = currentInput(out == 1, :);
                       
                   case 'PP'
                       inputRetrain = extraInput.data(ind_range, :);
%                        input0 = baseInput(out == 0, :);
%                        input1 = baseInput(out == 1, :);
               end
               
%                input01 = [input0; input1];
%                output01 = [out(out == 0) out(out == 1)]';
               
               testingDataParse.input = inputRetrain; % pre-PCA data
               testingDataParse.output = out';
               testingDataParse.distance = prob';
               testingDataParse.outputGndTruth = extraInput.testingLabel(ind_range); % INDEXING
               
               obj.clusterDR.retrain(testingDataParse, percentageToKeep, incParam.percentageDistanceFromCoreToKeep);
           end
           
        case 'pickPlaceDR'
            [out, prob] = obj.pickPlaceDR.classify(currentInput);
            
            fsm0 = find(out == 0);
            fsm1 = find(out == 1);
            ind0locArray{ind_incsvmmode} = ind_range(fsm0);
            ind1locArray{ind_incsvmmode} = ind_range(fsm1);
            
            if incParam.updateMachine && retrainFlagSideClass
                input0 = currentInput(out == 0, :);
                input1 = currentInput(out == 1, :);
                input01 = [input0; input1];
                output01 = [out(out == 0) out(out == 1)]';
                obj.pickPlaceDR.retrain(input01, output01);
            end
            
        otherwise
            % set it to 'all good' for inclusion, since
            % the setxor would remove the poor ones
            ind0locArray{ind_incsvmmode} = ind_range;
            ind1locArray{ind_incsvmmode} = ind_range;
    end
end

ind0loc = ind_range;
ind1loc = ind_range;
for ind_setxor = 1:length(ind0locArray)
    % pick the entry that all activated voters select
    ind0loc = intersect(ind0loc, ind0locArray{ind_setxor});
    ind1loc = intersect(ind1loc, ind1locArray{ind_setxor});
end
end

function [startInd, endInd, observationNew] = temporallyAlignedDataLoad(ind_range, extraInput)
    startInd = ind_range(1) + extraInput.settings.dimStrack;
    endInd = ind_range(end) + extraInput.settings.dimStrack; % offset, since we're using the main data and not the testdata

    observationNew.obsTime = extraInput.time(startInd:endInd);
    observationNew.obsData = extraInput.data(startInd:endInd, 1:5);
    observationNew.obsVelo = extraInput.data(startInd:endInd, 6:10);
end