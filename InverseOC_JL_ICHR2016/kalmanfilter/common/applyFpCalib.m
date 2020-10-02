function [F1, F2, M1, M2] = applyFpCalib(cortex_f1_adc, cortex_f1_adc_offset, cortex_f2_adc, cortex_f2_adc_offset)
    % defining the calibration matrix, as is from the AMTI calibration files
    % FP1 = 0541
    % FP2 = 0561
    
    % setting up conversion factors
    analogToVoltFactor = (2^16)/20; % analog units to voltage
    unitConv_force = 4.4482216282509; % lb to N
    unitConv_moment = 0.11298482933333; % lb in to N m
    
    % fp1 and fp2 calibration matrix
    FP1_calib_mtx = [0.1857718	-0.1510818	0.2560428	-0.2383812	0.3552629	-11.8969181	11.8100539	0.3578005
        -0.7145054	-0.3994613	0.1534758	0.452656	11.5410966	-0.1482715	0.1471889	11.6235325
        -24.0784063	-24.5258619	-24.1982515	-24.2939251	-0.3326175	-0.0151676	0.1636783	0.1749859
        173.1750103	175.2073443	-171.8125152	-174.1176904	0	0	0	0
        180.4522012	-181.6292807	180.1535382	-180.2768301	0	0	0	0
        -1.7264629	-1.7585462	-1.735056	-1.7419159	122.264218	-122.6549011	-121.7593482	-123.1375287];
    
    FP2_calib_mtx = [0.3119342	-0.3027807	0.4253523	-0.0834534	0.3247007	-11.6861478	11.6962775	0.3193138
        -0.7610162	-0.3878758	0.6419344	0.6092338	11.7931326	0.0639227	-0.0639781	11.5974797
        -23.8892571	-23.6714047	-24.3150895	-23.5307501	-0.3019417	-0.0247871	-0.0251555	-0.0965701
        173.2423527	170.311134	-170.8479385	-167.1206904	0	0	0	0
        179.1098509	-173.5546663	181.7210017	-174.9241195	0	0	0	0
        -1.5781653	-1.5637736	-1.6062965	-1.5544817	122.4718391	-120.9667541	-121.0716099	-120.4399804];
        
    % apply calibration matrix
    [F1_volt, F1_lb] = calibrateFPdata(cortex_f1_adc, FP1_calib_mtx, analogToVoltFactor, cortex_f1_adc_offset);
    [F2_volt, F2_lb] = calibrateFPdata(cortex_f2_adc, FP2_calib_mtx, analogToVoltFactor, cortex_f2_adc_offset);
    
    % apply unit conversion
    F1 = F1_lb(:, 1:3)* unitConv_force;
    F2 = F2_lb(:, 1:3) * unitConv_force;
    
    M1 = F1_lb(:, 4:6) * unitConv_moment;
    M2 = F2_lb(:, 4:6) * unitConv_moment; 
end

function [F1_volt, F1_lb] = calibrateFPdata(F1_adc, FP1_calib, analogToVoltFactor, F1_adc_unloaded)
    % converting from bits to voltage, and removing offset
    F1_volt = F1_adc / analogToVoltFactor; % ADC to voltage
    
    if ~isempty(F1_adc_unloaded)
        % external zero-offset calibration file exists
        F1_volt_unloaded = F1_adc_unloaded / analogToVoltFactor; % ADC to voltage
        F1_volt = F1_volt - repmat(mean(F1_volt_unloaded), size(F1_volt, 1), 1); % remove offset
    else
        % no exiternal calibration file exist, using the starting timestamp
         F1_volt = F1_volt - repmat(F1_volt(1, :), size(F1_volt, 1), 1); % remove offset
    end
   
    % apply calibration matrix
    F1_calib = [];
    for i = 1:size(F1_adc, 1)
        rawDataRow = F1_volt(i, :);
        calibDataRow = FP1_calib*rawDataRow';
        F1_calib = [F1_calib; calibDataRow'];
    end

    F1_lb = F1_calib;
end