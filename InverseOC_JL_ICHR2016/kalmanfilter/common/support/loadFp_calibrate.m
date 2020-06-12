function fp = loadFp_calibrate(filepath_data, filepath_offset)
    % defining the calibration matrix, as is from the AMTI calibration files
    % FP1 = 0541
    % FP2 = 0561
    
    % if filepath_offset is set, it will use this file as the unloaded
    % forceplate ('unloaded_fp1.anc' typically). if it is not set, it will
    % use the starting value of filepath_data to zero the data
    
    if ~exist('filepath_data', 'var')
        % sample code
        filepath_data = '..\data_fp\amtitest1.anc';
        filepath_offset = '';
    end
    
    % loading the test data
    if exist(filepath_data, 'file')
        [cortex_time, cortex_f1_adc, cortex_f2_adc] = loadFPDataCortex(filepath_data);
    else
        % if it doesn't exist...return nulls
        fp.time = [];
        fp.F1 = [];
        fp.F2 = [];
        fp.M1 = [];
        fp.M2 = [];
        return
    end
    
    % loading unloaded fp data if exist
    if exist(filepath_offset, 'file')
          [cortex_time_offset, cortex_f1_adc_offset, cortex_f2_adc_offset] = loadFPDataCortex(filepath_offset);
    else
          cortex_time_offset = [];
          cortex_f1_adc_offset = [];
          cortex_f2_adc_offset = [];
    end
        
     [F1, F2, M1, M2] = applyFpCalib(cortex_f1_adc, cortex_f1_adc_offset, cortex_f2_adc, cortex_f2_adc_offset);
     
     fp.time = cortex_time;
     fp.F1 = F1;
     fp.F2 = F2;
     fp.M1 = M1;
     fp.M2 = M2;
end

function [time, fp1, fp2] = loadFPDataCortex(fileLoad)
    fId = fopen(fileLoad);
    for i = 1:11
        % remove header data that Cortex inserts
        fgetl(fId);
    end
    data = textscan(fId,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter','\t');
    fclose(fId);
    
    time = data{1};
    F1XAB = data{2};
    F1YBD = data{3};
    F1YAC = data{4};
    F1XDC = data{5};
    F1ZA = data{6};
    F1ZB = data{7};
    F1ZC = data{8};
    F1ZD = data{9};
    F2XAB = data{10};
    F2YBD = data{11};
    F2YAC = data{12};
    F2XDC = data{13};
    F2ZA = data{14};
    F2ZB = data{15};
    F2ZC = data{16};
    F2ZD = data{17};

    % restack the array in the same layout as the AMTI calib mtx
    fp1 = [F1ZC F1ZD F1ZA F1ZB F1YAC F1XDC F1XAB F1YBD];
    fp2 = [F2ZC F2ZD F2ZA F2ZB F2YAC F2XDC F2XAB F2YBD];
end