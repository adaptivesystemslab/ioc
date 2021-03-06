% load a specific .mat file and using the fixed window cf functions,
% recover the ind window values
addpath(genpath('Common'));
addpath(genpath('kalmanfilter'));

% README so these pathing would all need to be updated. 
matPath = 'C:\Users\jf2lin\Downloads\TransferXL-00j4BNd2YW3mTp\Subject01_SingleArm60Time_4071_7071_end.mat'; % the IOC mat file
filepathModel = 'D:\aslab_github\ioc\InverseOC_JL_ICHR2016\kalmanfilter\instance_expressiveioc\model/ioc_v4_rightarm_fixedbase.xml'; % the EKF path to the RL model
filepathMat = 'D:\aslab\data_IK\expressive_ioc\2019_04_11_rightarm3\matEkfIk\Subject01_SingleArm60Time_OnlyRightArm_01_mocap_mocap_X00_floating_ekfId1_ekfIk.mat'; % the original IK mat file

% README some hardcoded values from setup_main
%      case 'Subject01'
pickPose = [0.9405    0.3285   -0.0873    0.1286;
    -0.1124    0.0582   -0.9920    0.3439;
    -0.3207    0.9427    0.0917   -0.2761;
    0         0         0    1.0000];
placementPose = [0.8497    0.2090    0.4842    0.0965;
    0.4769    0.0872   -0.8746    0.3725;
    -0.2250    0.9740   -0.0256   -0.1028;
    0         0         0    1.0000];

if ~exist('saveVar', 'var')
    load(matPath);
end

modelInstance = rlModelInstance_expressiveioc_rightArm('');
modelInstance.loadModelFromModelSpecsNoSensor(filepathModel, filepathMat);
      
param = saveVar.param;
param.cost_function_names = saveVar.cost_function_names;
param.dt_spline = mean(diff(saveVar.feature_win_save{1}.t));
param.model = modelInstance.model;
param.dofsFromFull = 1:size(saveVar.feature_win_save{1}.q, 1);
param.dofsFull = 1:size(saveVar.feature_win_save{1}.q, 1);

param.endEffectorName = 'frame_rhand_end';
param.end_eff_frames = {'frame_relbow_4', 'frame_rwrist_4', 'frame_rhand_end'};
param.end_eff_framesInds = {1:3, 4:6, 7:9};

param = feature_cf(param); % rset normalization

param.x_anchor_1 = pickPose(1:3, 4);
param.rot_anchor_1 = pickPose(1:3, 1:3);
param.x_anchor_2 = placementPose(1:3, 4);
param.rot_anchor_2 = placementPose(1:3, 1:3);

massWeights = [param.model.bodies(2).m param.model.bodies(3).m param.model.bodies(4).m];
massWeights = massWeights / sum(massWeights);

param.weights.quantity_motion_dx = massWeights;
param.weights.weight_effort = massWeights;
param.weights.time_effort = massWeights;
param.weights.space_effort = massWeights;
param.weights.flow_effort = massWeights;
param.weights.shapeDir = massWeights;

lenWin = length(saveVar.feature_win_save);
for ind_win = 1:lenWin
   [cf_name{ind_win}, len{ind_win}, feat{ind_win}, J_all{ind_win}, ...
       J_norm{ind_win}, J_ind{ind_win}, c_use{ind_win}] = loadStuff(saveVar, ind_win, param);
end

% generate the time varying weights
weight_avg = [];
counter_avg = [];
t = saveVar.feature_full.t;
for ind_t = 1:length(t)
    currWeight = [];
    avgCount = 0;
    curr_t = t(ind_t);
    for ind_win = 1:lenWin
        t_curr_array = saveVar.t_recon_plot_array{ind_win};
        if find(abs(t_curr_array - curr_t) < 0.005)
            avgCount = avgCount + 1;
            currWeight = [currWeight; c_use{ind_win}];
        end
    end
    
    meanedWeight = mean(currWeight, 1);
    weight_avg(ind_t, :) = meanedWeight;
    counter_avg(ind_t) = avgCount;
end

function [cf_name, len, feat, J_all, J_norm, J_ind, c_use] = loadStuff(saveVar, ind_win, param)
    % README so the original IOC features are generataed by taking the
    % original path, resample it as a spline trajectory, then calculate the
    % features from that. I had a lot of trouble duplicating the proper
    % parts of the splining because it is a minor nighmare mess of code,
    % but I've gotten it working with just the straight up original
    % trajectory. The difference should not be that much since our windows
    % are small, 60 timesteps under a 100 Hz sampling rate, so it's less
    % than a second of movement. but as a result, the numbers will not
    % match up perfectly
    
    % README2 I am also giving you the "recovered" weights (which is what I
    % feed back into DOC), which has been normalized so that the weights
    % all sum to 1. This is different then the "generated" weights, where
    % the pivot weight is set to 1 (so by def, the sum of weights is more
    % than 1). So there will be some minor scaling issues with J as well
    
    feature_opt = saveVar.feature_win_save{ind_win}; 
    minRmseInd = saveVar.minRmseIndArray(ind_win);
    c_use = saveVar.output_inverse{ind_win}(minRmseInd).c_recovered;
    len = length(feature_opt.t);
    
    % compute the cost function values corresponding to this combination.
    % that is, how does the cost function change when we change the
    feature_use = calc_features(feature_opt.q, [], param);
    [J_all, J_contrib, J_ind, J_cost_array] = calc_direct_cost(c_use, feature_use, param);
    
    for i = 1:length(param.cost_function_names)
        switch param.cost_function_names{i}
            case 'kinetic_en_sqrt'
                J_norm(i) = param.coeff_cf.geo;
                
            case 'dq_tau'
                J_norm(i) = param.coeff_cf.en;
                
            otherwise
                J_norm(i) = param.coeff_cf.(param.cost_function_names{i});
        end
        
        switch param.cost_function_names{i}
            case 'kinetic_en_sqrt'
                feat{i} = feature_use.geo;
                
            case 'dq_tau'
                feat{i} = feature_use.en;
                
            case 'quantity_motion_dx'
                dx = feature_use.dx_all;
                J = 0;
                for i = 1:size(dx, 1)
                    J = J + param.weights.quantity_motion_dx(i) * ...
                        sum(abs(dx(i, :)));
                end
                feat{i} = J / sum(param.weights.quantity_motion_dx);
                
            case 'space_effort'
                x = feature_use.x_all;
                eff = [0 0 0]';
                for i = 2:size(x, 2)
                    eff = eff + abs(x(:, i) - x(:, i-1));
                end
                feat{i} = param.weights.space_effort(:) .* eff(:);
                
            case 'time_effort'
                ddx = feature_use.ddx_all;
                eff = sum(abs(ddx), 2) / size(ddx, 2);
                feat{i} = param.weights.time_effort(:) .* eff(:);
                
            case 'volume_bounding_box'
                x = feature_use.x_max;
                vol = x(1, :) .* x(2, :) .* x(3, :);
                feat{i} = cost_fct_calc_sq(vol, len);
                
            case 'weight_effort'
                dx = feature_use.dx_all;
                eff = dx.^2;
                for i = 1:size(dx, 1)
                    eff(i, :) = eff(i, :)*param.weights.weight_effort(i);
                end
                effSum = sum(eff, 1); % sum across joints
                feat{i} = max(effSum);
                
            otherwise
                feat{i} = feature_use.(param.cost_function_names{i});
        end
    end
    
    % PAMELA THE VARIABLES YOU WANT ARE... (AT EACH WINDOW)
    % J_all = sum w * J_ind/J_norm
    % J_ind = sum(sum(feat .^ 2)) / len
    cf_name = param.cost_function_names; % order of the cfs
%     len % len
%     feat % the individual features used to calc the basis cf...sort of. some have been squared and some have not, so it'll be inconsistent unless you want to copy them all out
%     J_all % the total J val
%     J_norm % the normalization coeff for each basis cf
%     J_ind % individual cf terms
%     c_use % the weights, 
end

function J = cost_fct_calc_sq(feat, len)
	J = sum(sum(feat .^ 2)) / len;
% 	J = sum(sum(abs(feat))) / len;
end

function param = feature_cf(param)
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
            
            param.coeff_cf = array2param(convArray);
            param.normProfile = normProfile;
    end
    
    param.coeff_cf.array = param2array(param);
end

function array = param2array(param)

    array = [param.coeff_cf.ddq param.coeff_cf.dddq param.coeff_cf.ddx param.coeff_cf.dddx ...
        param.coeff_cf.tau param.coeff_cf.dtau param.coeff_cf.ddtau ...
        param.coeff_cf.ek param.coeff_cf.geo param.coeff_cf.en, ...
        param.coeff_cf.quantity_motion_dx, param.coeff_cf.volume_bounding_box, param.coeff_cf.weight_effort, ...
        param.coeff_cf.time_effort, param.coeff_cf.space_effort, param.coeff_cf.flow_effort, ...
        param.coeff_cf.x_anchor_1, param.coeff_cf.x_anchor_2, param.coeff_cf.rot_anchor_1, param.coeff_cf.rot_anchor_2, param.coeff_cf.com, ...
        param.coeff_cf.x_displace,  param.coeff_cf.dx_displace,  param.coeff_cf.cartCurv,  param.coeff_cf.dcom,  param.coeff_cf.shapeDir,...
        param.coeff_cf.x_cartCurv];

%     array = [param.coeff_cf.tau param.coeff_cf.en param.coeff_cf.geo, ...         
%         param.coeff_cf.quantity_motion_dx, param.coeff_cf.volume_bounding_box, param.coeff_cf.weight_effort, ...         
%         param.coeff_cf.time_effort, param.coeff_cf.space_effort, ...         
%         param.coeff_cf.x_anchor_1, param.coeff_cf.x_anchor_2, param.coeff_cf.com, ...         
%         param.coeff_cf.x_displace,  param.coeff_cf.dx_displace,  param.coeff_cf.cartCurv,  param.coeff_cf.dcom];

    
end

function coeff_cf = array2param(array)
    coeff_cf.ddq = array(1);
    coeff_cf.dddq = array(2);
    coeff_cf.ddx = array(3);
    coeff_cf.dddx = array(4);
    coeff_cf.tau = array(5);
    coeff_cf.dtau = array(6);
    coeff_cf.ddtau = array(7);
    coeff_cf.ek = array(8);
    coeff_cf.geo = array(9);
    coeff_cf.en = array(10);
    coeff_cf.quantity_motion_dx = array(11);
    coeff_cf.volume_bounding_box = array(12);
    coeff_cf.weight_effort = array(13);
    coeff_cf.time_effort = array(14);
    coeff_cf.space_effort = array(15);
    coeff_cf.flow_effort = array(16);
    coeff_cf.x_anchor_1 = array(17);
    coeff_cf.x_anchor_2 = array(18);
    coeff_cf.rot_anchor_1 = array(19);
    coeff_cf.rot_anchor_2 = array(20);
    coeff_cf.com = array(21);
    coeff_cf.x_displace = array(22);
    coeff_cf.dx_displace = array(23);
    coeff_cf.cartCurv = array(24);
    coeff_cf.dcom = array(25);
    coeff_cf.shapeDir = array(26);
    coeff_cf.x_cartCurv = array(27);
%       coeff_cf.tau = array(1);   
%       coeff_cf.en = array(2);     
%       coeff_cf.geo = array(3);     
%       coeff_cf.quantity_motion_dx = array(4);     
%       coeff_cf.volume_bounding_box = array(5);     
%       coeff_cf.weight_effort = array(6);     
%       coeff_cf.time_effort = array(7);     
%       coeff_cf.space_effort = array(8);        
%       coeff_cf.x_anchor_1 = array(9);     
%       coeff_cf.x_anchor_2 = array(10);     
%       coeff_cf.com = array(11);         
%       coeff_cf.x_displace = array(12);     
%       coeff_cf.dx_displace = array(13);     
%       coeff_cf.cartCurv = array(14);     
%       coeff_cf.dcom = array(15);
end