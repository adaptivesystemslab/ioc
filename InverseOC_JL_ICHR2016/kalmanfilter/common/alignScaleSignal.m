function [dt, timeArray, featureSet_imu_aligned, featureSet_mocap_aligned, dataInstance_imu_aligned, dataInstance_mocap_aligned, outputLog, h1] = ...
    alignScaleSignal(featureSet_mocap, featureSet_imu, dataInstance_mocap, dataInstance_imu, currFileEntry)

    cropEdgesFlag = 1;

    allJointStrMocap = featureSet_mocap.joint_labels;
    allJointStrImu = featureSet_imu.joint_labels;
    allSensorStrMocap = featureSet_mocap.measurement_labels;
    allSensorStrImu = featureSet_imu.measurement_labels;

    switch 'scale'
        case 'interpolate'
            [dt, timeArray, q_mocap, q_imu, x_mocap, x_imu] = alignInterpolate(featureSet_mocap, featureSet_imu, allJointStrMocap, allJointStrImu);
            h1 = [];
            featureSet_imu_aligned.q = q_imu;
            featureSet_mocap_aligned.q = q_mocapl;
            
            featureSet_imu_aligned.x = x_imu;
            featureSet_mocap_aligned.x = x_mocapl;
            
        case 'scale'
            if length(currFileEntry.exerciseName) < 12
                jointToAlign1 = 'joint_rknee_0';
                jointToAlign2 = 'joint_rknee_0';
            else
                switch currFileEntry.exerciseName(1:12)
                    case {'HEFO_STD_NON', 'KHEF_STD_NON'}
                        jointToAlign1 = 'joint_rhip_0';
                        jointToAlign2 = 'joint_rhip_0';
                        
                    case 'HAAO_STD_NON'
                        jointToAlign1 = 'joint_rhip_0';
                        jointToAlign2 = 'joint_rhip_1';
                        
                    case {'KEFO_SIT_NON', 'KFEO_STD_NON', 'SQUA_STD_NON', ...
                            'LUNG_STD_NON'}
                        jointToAlign1 = 'joint_rknee_0';
                        jointToAlign2 = 'joint_rknee_0';
                        
                    case {'LCAL_STD_NON', 'LFLA_STD_NON', 'ACAL_STD_NON'}
                        jointToAlign1 = 'joint_rknee_0';
                        jointToAlign2 = 'joint_rknee_0';
                        
                    case {'STEP_STD_NON', 'GAIT_STD_NON', 'OVER_STD_NON'}
                        jointToAlign1 = 'joint_rknee_0';
                        jointToAlign2 = 'joint_rknee_0';
                        
                    otherwise
                        jointToAlign1 = 'joint_rknee_0';
                        jointToAlign2 = 'joint_rknee_0';
                end
            end
            

            fieldsToAlign = {'q', 'dq', 'ddq', 'x', 'measurement_input', 'measurement_output'};
            
            % determine the time scaling
            [dt, timeArray1, featureSet_imu_aligned1, featureSet_mocap_aligned1, dataInstance_mocap_aligned1, dataInstance_imu_aligned1, h1, ...
                imu_first, imu_last, mocap_first, mocap_last, currJointInd_imu, currJointInd_mocap, startingOffset, endingOffset, output_N] = ...
                alignScale(featureSet_mocap, featureSet_imu, dataInstance_mocap, dataInstance_imu, jointToAlign1, allJointStrMocap, allJointStrImu, fieldsToAlign);
            if ishandle(h1)
                close(h1);
            end
            N1 = output_N/(mocap_last-mocap_first);
            N2 = output_N/(imu_last-imu_first);
            
            % then do a second alignment
            [dt, timeArray2, featureSet_imu_aligned2, featureSet_mocap_aligned2, dataInstance_mocap_aligned2, dataInstance_imu_aligned2, h1, D] = ...
                alignXcorr(featureSet_mocap_aligned1, featureSet_imu_aligned1, dataInstance_mocap_aligned1, dataInstance_imu_aligned1, ...
                jointToAlign2, allJointStrMocap, allJointStrImu, fieldsToAlign);
            if ishandle(h1)
                close(h1);
            end
            
            timeArray = timeArray2;
            featureSet_imu_aligned = featureSet_imu_aligned2;
            featureSet_mocap_aligned = featureSet_mocap_aligned2;
            dataInstance_imu_aligned = dataInstance_imu_aligned2;
            dataInstance_mocap_aligned = dataInstance_mocap_aligned2;
            
            if cropEdgesFlag
                % then remove the calibrating squat from the motions
                [timeArray3, indArray, featureSet_imu_aligned3, featureSet_mocap_aligned3, dataInstance_imu_aligned3, dataInstance_mocap_aligned3] = ...
                    alignCropZVC(featureSet_mocap_aligned2, featureSet_imu_aligned2, dataInstance_mocap_aligned2, dataInstance_imu_aligned2, jointToAlign1, allJointStrMocap, allJointStrImu, fieldsToAlign, ...
                    mocap_first, mocap_last);
                if ishandle(h1)
                    close(h1);
                end
                timeArray = timeArray3;
                featureSet_imu_aligned = featureSet_imu_aligned3;
                featureSet_mocap_aligned = featureSet_mocap_aligned3;
                dataInstance_imu_aligned = dataInstance_imu_aligned3;
                dataInstance_mocap_aligned = dataInstance_mocap_aligned3;
            end
            
            h1 = [];
                     
            h1 = figure;
            time_imu = featureSet_imu.time - featureSet_imu.time(1);
            time_mocap = featureSet_mocap.time - featureSet_mocap.time(1);
            signal_imu = filter_dualpassBW(featureSet_imu.q);
            signal_mocap = filter_dualpassBW(featureSet_mocap.q);
            ax(1) = subplot(411); hold on;
            plot(time_imu, signal_imu(:, currJointInd_imu), 'DisplayName', 'IMU');
            plot(time_imu(imu_first), signal_imu(imu_first, currJointInd_imu), 'kx');
            plot(time_imu(imu_last), signal_imu(imu_last, currJointInd_imu), 'kx');
            plot(time_mocap, signal_mocap(:, currJointInd_mocap), 'DisplayName', 'Mocap');
            plot(time_mocap(mocap_first), signal_mocap(mocap_first, currJointInd_mocap), 'kx');
            plot(time_mocap(mocap_last), signal_mocap(mocap_last, currJointInd_mocap), 'kx');
            title('orig filt(q)');
            
            %         figure;
            time_imu = featureSet_imu_aligned1.time;
            time_mocap = featureSet_mocap_aligned1.time;
            signal_imu = featureSet_imu_aligned1.q;
            signal_mocap = featureSet_mocap_aligned1.q;
            ax(2) = subplot(412);  hold on;
            plot(time_imu, signal_imu(:, currJointInd_imu), 'DisplayName', 'IMU');
            plot(time_imu(startingOffset), signal_imu(startingOffset, currJointInd_imu), 'kx');
            plot(time_imu(size(signal_imu, 1)-endingOffset), signal_imu(size(signal_imu, 1)-endingOffset, currJointInd_imu), 'kx');
            plot(time_mocap, signal_mocap(:, currJointInd_mocap), 'DisplayName', 'Mocap');
            plot(time_mocap(startingOffset), signal_mocap(startingOffset, currJointInd_mocap), 'kx');
            plot(time_mocap(size(signal_imu, 1)-endingOffset), signal_mocap(size(signal_imu, 1)-endingOffset, currJointInd_mocap), 'kx');
            title(['output q (N = mocap ' num2str(N1) ',imu ' num2str(N2) ')']);
            
            time_imu = featureSet_imu_aligned2.time;
            time_mocap = featureSet_mocap_aligned2.time;
            signal_imu = featureSet_imu_aligned2.q;
            signal_mocap = featureSet_mocap_aligned2.q;
            ax(3) = subplot(413);  hold on;
            plot(time_imu, signal_imu(:, currJointInd_imu), 'DisplayName', 'IMU');
            plot(time_mocap, signal_mocap(:, currJointInd_mocap), 'DisplayName', 'Mocap');
            title(['output q (D = ' num2str(D) ')']);
            
            if cropEdgesFlag
                time_imu = featureSet_imu_aligned3.time;
                time_mocap = featureSet_mocap_aligned3.time;
                signal_imu = featureSet_imu_aligned3.q;
                signal_mocap = featureSet_mocap_aligned3.q;
                ax(4) = subplot(414);  hold on;
                plot(time_imu, signal_imu(:, currJointInd_imu), 'DisplayName', 'IMU');
                plot(time_mocap, signal_mocap(:, currJointInd_mocap), 'DisplayName', 'Mocap');
                title(['output q (D = ' num2str(D) ')']);
            end
            
            linkaxes(ax, 'xy');
            
            outputLog.time_imu_align1 = featureSet_imu.time - featureSet_imu.time(1);
            outputLog.time_mocap_align1 = featureSet_mocap.time - featureSet_mocap.time(1);
            outputLog.signal_imu_align1 = filter_dualpassBW(featureSet_imu.q(:, currJointInd_imu));
            outputLog.signal_mocap_align1 = filter_dualpassBW(featureSet_mocap.q(:, currJointInd_imu));
            outputLog.timefirst_imu_align1 = imu_first;
            outputLog.timelast_imu_align1 = imu_last;
            outputLog.timefirst_mocap_align1 = mocap_first;
            outputLog.timelast_mocap_align1 = mocap_last;
            
            outputLog.time_imu_align2 = featureSet_imu_aligned1.time;
            outputLog.time_mocap_align2 = featureSet_mocap_aligned1.time;
            outputLog.signal_imu_align2 = featureSet_imu_aligned1.q(:, currJointInd_imu);
            outputLog.signal_mocap_align2 = featureSet_mocap_aligned1.q(:, currJointInd_imu);
            outputLog.timefirst_align2 = startingOffset;
            outputLog.timelast_align2 = endingOffset;
            
            if ~cropEdgesFlag
                outputLog.time_imu_align3 = featureSet_imu_aligned2.time;
                outputLog.time_mocap_align3 = featureSet_mocap_aligned2.time;
                outputLog.signal_imu_align3 = featureSet_imu_aligned2.q(:, currJointInd_imu);
                outputLog.signal_mocap_align3 = featureSet_mocap_aligned2.q(:, currJointInd_imu);
            else
                outputLog.time_imu_align3 = featureSet_imu_aligned3.time;
                outputLog.time_mocap_align3 = featureSet_mocap_aligned3.time;
                outputLog.signal_imu_align3 = featureSet_imu_aligned3.q(:, currJointInd_imu);
                outputLog.signal_mocap_align3 = featureSet_mocap_aligned3.q(:, currJointInd_imu);
            end
            
        case 'alignXcorr'
            [dt, timeArray, featureSet_imu_aligned, featureSet_mocap_aligned, h1] = ...
                alignXcorr(featureSet_mocap, featureSet_imu, 'joint_rknee_0', allJointStrMocap, allJointStrImu);
    end