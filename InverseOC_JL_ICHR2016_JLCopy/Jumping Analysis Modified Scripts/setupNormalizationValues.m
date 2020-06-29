function param = setupNormalizationValues(param, traj_load, segmentInfo, config_feature, config_cf)

    % use the input data
    feature_full.t = param.t_full;
    feature_full.q = traj_load.q;
    feature_full.dq = traj_load.dq;
    feature_full.ddq = traj_load.ddq;
    feature_full.dddq = traj_load.dddq;

    % generate the feature set with no normalization enabled
    param_nonorm = feature_norm('none', [], param); 
    param_nonorm = feature_cf('none', [], param_nonorm, []);
    param_nonorm.windowMode = 'narrow';
    
   [~, indStart] = findClosestValue(segmentInfo.timeStart(1), feature_full.t);
   [~, indEnd] = findClosestValue(segmentInfo.timeEnd(1), feature_full.t);

    feature_win = feature_windowed(feature_full, indStart, indEnd, param_nonorm);
    
    param_nonorm.cost_function_names = {'ddq', 'dddq', 'ddx', 'dddx', 'tau', 'dtau', 'ddtau', ...
        'potential_en', 'kinetic_en', 'kinetic_en_sqrt', 'dq_tau', 'cop', 'dcop', ...
        'COM', 'dCOM', 'ddCOM', 'dddCOM', 'dCOM_mag', 'ddCOM_mag', 'dddCOM_mag', ...
        'COM_height', 'dCOM_height', 'ddCOM_height', 'dddCOM_height', ...
        'dToeZ', 'ddToeZ', 'ToeTargX', 'dToeTargX', 'COMToeX', 'dCOMToeX', 'COMTargX', 'dCOMTargX'};

    [~, ~, J_cost_array] = calc_direct_cost(ones(size(param_nonorm.cost_function_names)), feature_win, param_nonorm);
    
%     if 1
%         % generate value for feature_cf
%         param_test = feature_cf('rang', J_cost_array, param_nonorm, param_nonorm.cost_function_names);
%         pullVal = param_test.coeff_cf.array
%     end
    
    param = feature_norm(config_feature, feature_win, param);
    param = feature_cf(config_cf, J_cost_array, param, param_nonorm.cost_function_names);
end

function rangeout = rangeMod(data)
    rangeout = range(data);
    
    for i = 1:length(rangeout)
        if rangeout(i) == 0
            rangeout(i) = 1;
        end
    end
end

function param = feature_cf(config, J_cost_array, param, cost_function_names)
    switch config
        case 'rang'
            % sylla style
            param.coeff_cf.ddq   = 1 ./ J_cost_array(strcmpi('ddq', cost_function_names));
            param.coeff_cf.dddq  = 1 ./ J_cost_array(strcmpi('dddq', cost_function_names));
            param.coeff_cf.ddx   = 1 ./ J_cost_array(strcmpi('ddx', cost_function_names));
            param.coeff_cf.dddx  = 1 ./ J_cost_array(strcmpi('dddx', cost_function_names));
            
            param.coeff_cf.tau   = 1 ./ J_cost_array(strcmpi('tau', cost_function_names));
            param.coeff_cf.dtau  = 1 ./ J_cost_array(strcmpi('dtau', cost_function_names));
            param.coeff_cf.ddtau = 1 ./ J_cost_array(strcmpi('ddtau', cost_function_names));
            param.coeff_cf.ep    = 1 ./ J_cost_array(strcmpi('potential_en', cost_function_names));
            param.coeff_cf.ek    = 1 ./ J_cost_array(strcmpi('kinetic_en', cost_function_names));
            param.coeff_cf.geo   = 1 ./ J_cost_array(strcmpi('kinetic_en_sqrt', cost_function_names));
            param.coeff_cf.en    = 1 ./ J_cost_array(strcmpi('dq_tau', cost_function_names));
            param.coeff_cf.cop    = 1 ./ J_cost_array(strcmpi('cop', cost_function_names));
            param.coeff_cf.dcop   = 1 ./ J_cost_array(strcmpi('dcop', cost_function_names));
            
            param.coeff_cf.COM   = 1 ./ J_cost_array(strcmpi('COM', cost_function_names));
            param.coeff_cf.dCOM   = 1 ./ J_cost_array(strcmpi('dCOM', cost_function_names));
            param.coeff_cf.ddCOM   = 1 ./ J_cost_array(strcmpi('ddCOM', cost_function_names));
            param.coeff_cf.dddCOM   = 1 ./ J_cost_array(strcmpi('dddCOM', cost_function_names));
            param.coeff_cf.dCOM_mag   = 1 ./ J_cost_array(strcmpi('dCOM_mag', cost_function_names));
            param.coeff_cf.ddCOM_mag   = 1 ./ J_cost_array(strcmpi('ddCOM_mag', cost_function_names));
            param.coeff_cf.dddCOM_mag   = 1 ./ J_cost_array(strcmpi('dddCOM_mag', cost_function_names));
            param.coeff_cf.COM_height   = 1 ./ J_cost_array(strcmpi('COM_height', cost_function_names));
            param.coeff_cf.dCOM_height   = 1 ./ J_cost_array(strcmpi('dCOM_height', cost_function_names));
            param.coeff_cf.ddCOM_height   = 1 ./ J_cost_array(strcmpi('ddCOM_height', cost_function_names));
            param.coeff_cf.dddCOM_height   = 1 ./ J_cost_array(strcmpi('dddCOM_height', cost_function_names));
            
            param.coeff_cf.dToeZ   = 1 ./ J_cost_array(strcmpi('dToeZ', cost_function_names));
            param.coeff_cf.ddToeZ   = 1 ./ J_cost_array(strcmpi('ddToeZ', cost_function_names));
            
            param.coeff_cf.ToeTargX   = 1 ./ J_cost_array(strcmpi('ToeTargX', cost_function_names));
            param.coeff_cf.dToeTargX   = 1 ./ J_cost_array(strcmpi('dToeTargX', cost_function_names));
            param.coeff_cf.COMToeX   = 1 ./ J_cost_array(strcmpi('COMToeX', cost_function_names));
            param.coeff_cf.dCOMToeX   = 1 ./ J_cost_array(strcmpi('dCOMToeX', cost_function_names));
            param.coeff_cf.COMTargX   = 1 ./ J_cost_array(strcmpi('COMTargX', cost_function_names));
            param.coeff_cf.dCOMTargX   = 1 ./ J_cost_array(strcmpi('dCOMTargX', cost_function_names));

        case 'none'
            % sylla style
            param.coeff_cf.ddq   = 1;
            param.coeff_cf.dddq  = 1;
            param.coeff_cf.ddx   = 1;
            param.coeff_cf.dddx  = 1;
            
            param.coeff_cf.tau   = 1;
            param.coeff_cf.dtau  = 1;
            param.coeff_cf.ddtau = 1;
            param.coeff_cf.ep    = 1;
            param.coeff_cf.ek    = 1;
            param.coeff_cf.geo   = 1;
            param.coeff_cf.en    = 1;
            param.coeff_cf.cop    = 1;
            param.coeff_cf.dcop   = 1;
            
            param.coeff_cf.COM   = 1;
            param.coeff_cf.dCOM   = 1;
            param.coeff_cf.ddCOM   = 1;
            param.coeff_cf.dddCOM   = 1;
            param.coeff_cf.dCOM_mag   = 1;
            param.coeff_cf.ddCOM_mag   = 1;
            param.coeff_cf.dddCOM_mag   = 1;
            param.coeff_cf.COM_height   = 1;
            param.coeff_cf.dCOM_height   = 1;
            param.coeff_cf.ddCOM_height   = 1;
            param.coeff_cf.dddCOM_height   = 1;
            
            param.coeff_cf.dToeZ   = 1;
            param.coeff_cf.ddToeZ   = 1;
            
            param.coeff_cf.ToeTargX   = 1;
            param.coeff_cf.dToeTargX   = 1;
            param.coeff_cf.COMToeX   = 1;
            param.coeff_cf.dCOMToeX   = 1;
            param.coeff_cf.COMTargX   = 1;
            param.coeff_cf.dCOMTargX   = 1;
            
        case 'rset'
            % to obtain values for this, use 'rang' to and average over all
            % the subjects
            switch param.dataset
                case {'squats_tuat_2011'}
%                     convArray = [0.105125133245040,0.00240354826516813,0.551764742792730,0.00755871910030054,9.84947320988321e-06,4.18297151919837e-06,2.94437170598407e-08,0.00228608574516544,0.0193066399340372,0.106863794781893,0.00433856866000618,17.5619860909456,4.71352242480247;];
                    convArray = [0.109811900666667,0.00252979850000000,0.747077835000000,0.0111124616666667,0.000350317833333333,0.000217953333333333,2.93216666666667e-06,0.0190849576666667,0.137866198000000,0.305067012833333,0.0257237905000000,140.171846162167,82.4768502316667;];
                    normProfile = 'squats_tuat_2011';
                    
                case {'sim', 'squats_tuat_2015'}
%                     convArray = [0.109811900616738,0.00252979861696423,0.690534132247458,0.0103504869674637,4.05850558960241e-05,1.93561064030306e-05,2.02777628683924e-07,0.00446049483991481,0.0219897985847883,0.111569461876832,0.00761042445413306,109.244959168577,22.1185024541362;];
                    convArray = [0.0209193477500000,0.000610628637500000,0.391311233500000,0.00917582787500000,0.000555224862500000,4.87641200000000e-05,1.04980645000000e-06,0.0286230961250000,0.0615042142500000,0.227054622250000,0.0259565217500000,36.2111514456250,25.9765244890000;];
                    normProfile = 'squats_tuat_2015';
                       
                case {'healthy1ankle'}
%                     convArray = [0.0718071301521757,0.00173238702742724,0.381796073178321,0.00864465904688942,8.01134768988754e-05,1.40086764607366e-05,1.90335711866653e-07,0.00749498968821210,0.0223835908729775,0.109384977028953,0.0131423993992226,120.063327275424,12.2626834339111;];
                    convArray = [0.0718071301750000,0.00173238705000000,0.381796073175000,0.00864465910000000,0.000162312695000000,4.88828002500000e-05,6.42907222500000e-07,0.0113850658000000,0.0459814886000000,0.156623211700000,0.0209052461000000,57.1828467592500,24.2696260045000;];
                    normProfile = 'healthy1ankle';

                case {'healthy1hip'}
%                     convArray = [0.189608878671599,0.00300635555801744,0.158954475971594,0.00251107669665943,0.000253641032073400,0.000408490880293129,3.23341504849206e-06,0.0136817133157400,0.421313169332929,0.517479173197018,0.0669053251439090,12.6599493430575,45.9417691625808;];
                    convArray = [0.189588507298246,0.00300635624561404,0.158954475982456,0.00251107666666667,0.000253640543859649,0.000408492385964912,3.23156140350877e-06,0.0136820203333333,0.421313169298246,0.517479173192982,0.0669053251578947,442.802011495263,208.424977970649;];
                    normProfile = 'healthy1hip';
                    
                case {'jump2D'}
%                     convArray = [0.189588507298246,0.00300635624561404,0.158954475982456,0.00251107666666667,0.000253640543859649,0.000408492385964912,3.23156140350877e-06,0.0136820203333333,0.421313169298246,0.517479173192982,0.0669053251578947,442.802011495263,208.424977970649;];
                    convArray = ones(1,13);
                    normProfile = 'jump2D';
            end
            
            param.coeff_cf = array2param(convArray);
            param.normProfile = normProfile;
    end
    
    param.coeff_cf.array = param2array(param);
end

function array = param2array(param)
    array = [param.coeff_cf.ddq param.coeff_cf.dddq param.coeff_cf.ddx param.coeff_cf.dddx ...
        param.coeff_cf.tau param.coeff_cf.dtau param.coeff_cf.ddtau param.coeff_cf.ep param.coeff_cf.ek ...
        param.coeff_cf.geo param.coeff_cf.en param.coeff_cf.cop param.coeff_cf.dcop];
end

function coeff_cf = array2param(array)
    coeff_cf.ddq = array(1);
    coeff_cf.dddq = array(2);
    coeff_cf.ddx = array(3);
    coeff_cf.dddx = array(4);
    coeff_cf.tau = array(5);
    coeff_cf.dtau = array(6);
    coeff_cf.ddtau = array(7);
    coeff_cf.ep = array(8);
    coeff_cf.ek = array(9);
    coeff_cf.geo = array(10);
    coeff_cf.en = array(11);
    coeff_cf.cop = array(12);
    coeff_cf.dcop= array(13);
    
    
end

function param = feature_norm(config, feature_win, param)
    switch config(1:4) % feature-style normalization
        case 'none'
            param.coeff_feature.q     = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.dq    = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.ddq   = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.dddq  = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.x     = ones(3, 1);
            param.coeff_feature.dx    = ones(3, 1);
            param.coeff_feature.ddx   = ones(3, 1);
            param.coeff_feature.dddx  = ones(3, 1);
            
            param.coeff_feature.tau   = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.dtau  = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.ddtau = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.ep    = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.ek    = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.geo   = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.en    = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.cop    = ones(3, 1); %
            param.coeff_feature.dcop   = ones(3, 1); %
            param.coeff_feature.ddcop  = ones(3, 1); %
            param.coeff_feature.com    = ones(3, 1); %
            param.coeff_feature.dcom   = ones(3, 1); %
            param.coeff_feature.ddcom  = ones(3, 1); %
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            param.coeff_feature.COM   = ones(3, 1);
            param.coeff_feature.dCOM   = ones(3, 1);
            param.coeff_feature.ddCOM   = ones(3, 1);
            param.coeff_feature.dddCOM   = ones(3, 1);
            param.coeff_feature.dCOM_mag   = 1;
            param.coeff_feature.ddCOM_mag   = 1;
            param.coeff_feature.dddCOM_mag   = 1;
            param.coeff_feature.COM_height   = 1;
            param.coeff_feature.dCOM_height   = 1;
            param.coeff_feature.ddCOM_height   = 1;
            param.coeff_feature.dddCOM_height   = 1;
            
            param.coeff_feature.ToeZ   = 1;
            param.coeff_feature.dToeZ  = 1;
            
            param.coeff_feature.ToeTargX   = 1;
            param.coeff_feature.dToeTargX   = 1;
            param.coeff_feature.COMToeX   = 1;
            param.coeff_feature.dCOMToeX   = 1;
            param.coeff_feature.COMTargX   = 1;
            param.coeff_feature.dCOMTargX   = 1;
    end
    
    if length(config) == 4
        masterMultiplier = 1;
    else
        switch config(5:8)
            case 'Half'
                masterMultiplier = 0.5;

            case 'Doub'
                masterMultiplier = 2;
                
            otherwise
                masterMultiplier = 1;
        end
    end
    
    param.coeff_feature.q     = param.coeff_feature.q*masterMultiplier;
    param.coeff_feature.dq    = param.coeff_feature.dq*masterMultiplier;
    param.coeff_feature.ddq   = param.coeff_feature.ddq*masterMultiplier;
    param.coeff_feature.dddq  = param.coeff_feature.dddq*masterMultiplier;
    param.coeff_feature.x     = param.coeff_feature.x*masterMultiplier;
    param.coeff_feature.dx    = param.coeff_feature.dx*masterMultiplier;
    param.coeff_feature.ddx   = param.coeff_feature.ddx*masterMultiplier;
    param.coeff_feature.dddx  = param.coeff_feature.dddx*masterMultiplier;
    param.coeff_feature.tau   = param.coeff_feature.tau*masterMultiplier;
    param.coeff_feature.dtau  = param.coeff_feature.dtau*masterMultiplier;
    param.coeff_feature.ddtau = param.coeff_feature.ddtau*masterMultiplier;
    param.coeff_feature.ep    = param.coeff_feature.ep*masterMultiplier;
    param.coeff_feature.ek    = param.coeff_feature.ek*masterMultiplier;
    param.coeff_feature.geo   = param.coeff_feature.geo*masterMultiplier;
    param.coeff_feature.en    = param.coeff_feature.en*masterMultiplier;
    param.coeff_feature.cop    = param.coeff_feature.cop*masterMultiplier;
    param.coeff_feature.dcop   = param.coeff_feature.dcop*masterMultiplier;
    param.coeff_feature.ddcop  = param.coeff_feature.ddcop*masterMultiplier;
    param.coeff_feature.com    = param.coeff_feature.com*masterMultiplier;
    param.coeff_feature.dcom   = param.coeff_feature.dcom*masterMultiplier;
    param.coeff_feature.ddcom  = param.coeff_feature.ddcom*masterMultiplier;
end