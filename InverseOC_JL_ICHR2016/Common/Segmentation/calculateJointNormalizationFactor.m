function [normVal_q, maxInd_q1, maxInd_q2, normVal_q_all, maxInd_q1_all ,maxInd_q2_all] = ...
    calculateJointNormalizationFactor(normalizeFlag, ekfQ, ekfTime, segTimeStart, segTimeEnd)

    lenTimeStart = length(segTimeStart);

    % calculate the joint angle magnitude normalization
    switch normalizeFlag
        case 1
            [normVal_q, maxInd_q1, maxInd_q2] = normalizeValueFromSegments(ekfQ, ekfTime, segTimeStart, segTimeEnd, 1:lenTimeStart);
            normVal_q_all = 0;
            maxInd_q1_all = 0;
            maxInd_q2_all = 0;
%             fprintf('Training - full normalization: %f \n', normVal_q);

        case 2
            switch obj.settings.mode
                case 'Training'
                    [normVal_q, maxInd_q1, maxInd_q2] = normalizeValueFromSegments(ekfQ, ekfTime, segTimeStart, segTimeEnd, 1:lenTimeStart);
                    normVal_q_all = 0;
                    maxInd_q1_all = 0;
                    maxInd_q2_all = 0;
%                     fprintf('Training - full normalization: %f \n', normVal_q);

                case 'Testing'
                    [normVal_q_all, maxInd_q1_all, maxInd_q2_all] = normalizeValueFromSegments(ekfQ, ekfTime, segTimeStart, segTimeEnd, 1:lenTimeStart);
                    [normVal_q, maxInd_q1, maxInd_q2]             = normalizeValueFromSegments(ekfQ, ekfTime, segTimeStart, segTimeEnd, 1);
%                     fprintf('Testing         - full normalization: %f \n', normVal_q);
%                     fprintf('Testing if full - full normalization: %f \n', normVal_q_all);
            end

        otherwise
            normVal_q = 1;
            maxInd_q1 = 0;
            maxInd_q2 = 0;
            normVal_q_all = 0;
            maxInd_q1_all = 0;
            maxInd_q2_all = 0;
    end