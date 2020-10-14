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
    
   [~, indStart] = findClosestValue(segmentInfo.timeStart(1), feature_full.t); % normalizes over first frame of first window...
   [~, indEnd] = findClosestValue(segmentInfo.timeEnd(end), feature_full.t); % ...to last frame of last window 

    feature_win = feature_windowed(feature_full, indStart, indEnd, param_nonorm);
    
%     param_nonorm.cost_function_names = {'ddq', 'dddq', 'ddx', 'dddx', 'tau', 'dtau', 'ddtau', ...
%         'potential_en', 'kinetic_en', 'kinetic_en_sqrt', 'dq_tau', 'cop', 'dcop', ...
%         'COM', 'dCOM', 'ddCOM', 'dddCOM', 'dCOM_mag', 'ddCOM_mag', 'dddCOM_mag', ...
%         'COM_height', 'dCOM_height', 'ddCOM_height', 'dddCOM_height', ...
%         'dToeZ', 'ddToeZ', 'ToeTargX', 'dToeTargX', 'COMToeX', 'dCOMToeX', 'COMTargX', 'dCOMTargX', ...
%         'ddq_leg', 'dddq_leg', 'tau_leg',  'dq_tau_leg', 'kinetic_en_sqrt_leg', ...
%         'ddq_arm', 'dddq_arm', 'tau_arm', 'dq_tau_arm', 'kinetic_en_sqrt_arm', ...
%         'ddq_tor', 'dddq_tor', 'tau_tor', 'dq_tau_tor', 'kinetic_en_sqrt_tor', ...
%         'ang_mom', 'dang_mom', 'ddang_mom'};
   
   param_nonorm.cost_function_names = {'half_joint_task', 'joint_length', 'joint_limit', ...
       'manip_rot', 'manip_trans', 'manipulability', 'orientation_length', 'task_length'};
   
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
            
            param.coeff_cf.half_joint_task   = 1 ./ J_cost_array(strcmpi('half_joint_task', cost_function_names));
            param.coeff_cf.joint_length   = 1 ./ J_cost_array(strcmpi('joint_length', cost_function_names));
            param.coeff_cf.joint_limit   = 1 ./ J_cost_array(strcmpi('joint_limit', cost_function_names));
            param.coeff_cf.manip_rot   = 1 ./ J_cost_array(strcmpi('manip_rot', cost_function_names));
            param.coeff_cf.manip_trans   = 1 ./ J_cost_array(strcmpi('manip_trans', cost_function_names));
            param.coeff_cf.manipulability   = 1 ./ J_cost_array(strcmpi('manipulability', cost_function_names));
            param.coeff_cf.orientation_length   = 1 ./ J_cost_array(strcmpi('orientation_length', cost_function_names));
            param.coeff_cf.task_length   = 1 ./ J_cost_array(strcmpi('task_length', cost_function_names));
            
        case 'none'
            % sylla style
            param.coeff_cf.ddq   = 1;
            param.coeff_cf.dddq  = 1;
            param.coeff_cf.ddx   = 1;
            param.coeff_cf.dddx  = 1;
            
            param.coeff_cf.half_joint_task   = 1;
            param.coeff_cf.joint_length   = 1;
            param.coeff_cf.joint_limit   = 1;
            param.coeff_cf.manip_rot   = 1;
            param.coeff_cf.manip_trans   = 1;
            param.coeff_cf.manipulability   = 1;
            param.coeff_cf.orientation_length   = 1;
            param.coeff_cf.task_length   = 1;
            
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
                    
                case {'expressive_ioc'}
%                     convArray = [0.0151567656717017,0.000141203809158701,0.408202998389100,0.00529903065978939,0.0438651754029304,0.0582138403879714,0.000463160753182444,7.10256892089768,1.73091816733913,0.552881966878168,4.83250045750189e-07,3165.19715217695,1.18925808373141,2.17528253105143,0.133299135074705,0.0302740929280632,3.81342449106885,3.81342449106885,25.4926049459811];
%                     convArray = [0.0151567656717017,0.000141203809158701,0.408202998389100,0.00529903065978939,0.0438651754029304,...
%                         0.0582138403879714,0.000463160753182444,7.10256892089768,1.73091816733913,0.552881966878168,4.83250045750189e-07,...
%                         3165.19715217695,1.18925808373141,2.17528253105143,0.133299135074705,0.0302740929280632,3.21993761330724,...
%                         2.62767410765346,25.4926049459811,9.58393434815072, 15.0455091167740, 0.0714144429792718, 0.000711912139210026,...
%                         1.14278503267674e-09];
%                     convArray = [0.0438651754029304,0.552881966878168,0.552881966878168,4.83250045750189e-07,3165.19715217695,1.18925808373141,...
%                         2.17528253105143,0.133299135074705,3.21993761330724,2.62767410765346,25.4926049459811,9.58393434815072,15.0455091167740,...
%                         0.0714144429792718,0.000711912139210026];
                    convArray = [0.0151567656717017,0.000141203809158701,0.408202998389099,0.00529903065978938,0.0518887469154429,0.0651536180356405,0.000522336280642349,...
                        7.50966260934148,1.78059348677811,0.588658903502858,4.85812173758764e-07,3165.19715217695,1.19775859687547,2.18272357556875,0.135474296170114,...
                        0.0303632098137464,16.9974736278158,29.9724760071953,1.40728003100660,2.17755952718146,27.0101817617975,9.58393434815072,15.0455091167740,...
                        0.0714144429792716,0.000766085082581547, 1.14278503271832e-09,0.166420432680378];
                    normProfile = 'expressive_ioc';
            end
            
            param.coeff_cf = array2param(convArray);
            param.normProfile = normProfile;
    end
    
    param.coeff_cf.array = param2array(param);
end

function array = param2array(param)
%     array = [param.coeff_cf.ddq param.coeff_cf.dddq param.coeff_cf.ddx param.coeff_cf.dddx ...
%         param.coeff_cf.tau param.coeff_cf.dtau param.coeff_cf.ddtau ...
%         param.coeff_cf.ek param.coeff_cf.geo param.coeff_cf.en, ...
%         param.coeff_cf.quantity_motion_dx, param.coeff_cf.volume_bounding_box, param.coeff_cf.weight_effort, ...
%         param.coeff_cf.time_effort, param.coeff_cf.space_effort, param.coeff_cf.flow_effort, ...
%         param.coeff_cf.x_anchor_1, param.coeff_cf.x_anchor_2, param.coeff_cf.rot_anchor_1, param.coeff_cf.rot_anchor_2, param.coeff_cf.com, ...
%         param.coeff_cf.x_displace,  param.coeff_cf.dx_displace,  param.coeff_cf.cartCurv,  param.coeff_cf.dcom,  param.coeff_cf.shapeDir,...
%         param.coeff_cf.x_cartCurv];

%     array = [param.coeff_cf.tau param.coeff_cf.en param.coeff_cf.geo, ...         
%         param.coeff_cf.quantity_motion_dx, param.coeff_cf.volume_bounding_box, param.coeff_cf.weight_effort, ...         
%         param.coeff_cf.time_effort, param.coeff_cf.space_effort, ...         
%         param.coeff_cf.x_anchor_1, param.coeff_cf.x_anchor_2, param.coeff_cf.com, ...         
%         param.coeff_cf.x_displace,  param.coeff_cf.dx_displace,  param.coeff_cf.cartCurv,  param.coeff_cf.dcom];

    array = [...
        param.coeff_cf.half_joint_task 
        param.coeff_cf.joint_length 
        param.coeff_cf.joint_limit
        param.coeff_cf.manip_rot
        param.coeff_cf.manip_trans
        param.coeff_cf.manipulability
        param.coeff_cf.orientation_length
        param.coeff_cf.task_length];
end

function coeff_cf = array2param(array)
%     coeff_cf.ddq = array(1);
%     coeff_cf.dddq = array(2);
%     coeff_cf.ddx = array(3);
%     coeff_cf.dddx = array(4);
%     coeff_cf.tau = array(5);
%     coeff_cf.dtau = array(6);
%     coeff_cf.ddtau = array(7);
%     coeff_cf.ek = array(8);
%     coeff_cf.geo = array(9);
%     coeff_cf.en = array(10);
%     coeff_cf.quantity_motion_dx = array(11);
%     coeff_cf.volume_bounding_box = array(12);
%     coeff_cf.weight_effort = array(13);
%     coeff_cf.time_effort = array(14);
%     coeff_cf.space_effort = array(15);
%     coeff_cf.flow_effort = array(16);
%     coeff_cf.x_anchor_1 = array(17);
%     coeff_cf.x_anchor_2 = array(18);
%     coeff_cf.rot_anchor_1 = array(19);
%     coeff_cf.rot_anchor_2 = array(20);
%     coeff_cf.com = array(21);
%     coeff_cf.x_displace = array(22);
%     coeff_cf.dx_displace = array(23);
%     coeff_cf.cartCurv = array(24);
%     coeff_cf.dcom = array(25);
%     coeff_cf.shapeDir = array(26);
%     coeff_cf.x_cartCurv = array(27);
% %       coeff_cf.tau = array(1);   
% %       coeff_cf.en = array(2);     
% %       coeff_cf.geo = array(3);     
% %       coeff_cf.quantity_motion_dx = array(4);     
% %       coeff_cf.volume_bounding_box = array(5);     
% %       coeff_cf.weight_effort = array(6);     
% %       coeff_cf.time_effort = array(7);     
% %       coeff_cf.space_effort = array(8);        
% %       coeff_cf.x_anchor_1 = array(9);     
% %       coeff_cf.x_anchor_2 = array(10);     
% %       coeff_cf.com = array(11);         
% %       coeff_cf.x_displace = array(12);     
% %       coeff_cf.dx_displace = array(13);     
% %       coeff_cf.cartCurv = array(14);     
% %       coeff_cf.dcom = array(15);
end

function param = feature_norm(config, feature_win, param)
    switch config(1:4) % feature-style normalization
        case 'none'
            % Set to size of variable at a single frame
            param.coeff_feature.q     = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.dq    = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.ddq   = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.dddq  = ones(length(param.dofsFromFull), 1);
            param.coeff_feature.x     = ones(3, 1);
            param.coeff_feature.dx    = ones(3, 1);
            param.coeff_feature.ddx   = ones(3, 1);
            param.coeff_feature.dddx  = ones(3, 1);

            param.coeff_feature.half_joint_task   = 1;
            param.coeff_feature.joint_length   = 1;
            param.coeff_feature.joint_limit   = 1;
            param.coeff_feature.manip_rot   = 1;
            param.coeff_feature.manip_trans   = 1;
            param.coeff_feature.manipulability   = 1;
            param.coeff_feature.orientation_length   = 1;
            param.coeff_feature.task_length   = 1;
    end
end