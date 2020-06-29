if ispc
    addpath(genpath('D:\aslab_gitlab\kalmanfilter\General_FKEKF\DynamicsModelMatlab\MatlabWrapper'));
    addpath(genpath('D:\aslab_gitlab\kalmanfilter\ik_framework\common'));
    addpath(genpath('D:\aslab_gitlab\kalmanfilter\ik_framework\instance_expressiveioc'));
    xmlPath = 'D:\aslab_gitlab\kalmanfilter\ik_framework\instance_expressiveioc\model';
    matPath = 'D:\aslab\data_IK\PamelasData\2019_04_11_fullbody2\matEkfIk\';
    trcPath = 'D:\aslab\data\PamelasData\';
elseif isunix
    addpath(genpath('../../General_FKEKF/DynamicsModelMatlab/MatlabWrapper/'));
    addpath(genpath('../../ik_framework/common'));
    addpath(genpath('../../ik_framework/instance_expressiveioc'));
    xmlPath = '../../ik_framework/instance_expressiveioc/model/';
    %matPath = '../../../data/2019_04_11_fullBody/matEkfIk/';
    matPath = '../../../data/2019_04_11_singleArm/matEkfIk/';
    trcPath = '../../../data/Vicon/';
    disp('Linux');
else
    disp('Platform not supported') 
end

% path to the IK mat file and initialize the varables
subjectId = 0;

% addpath(genpath('D:\aslab_gitlab\kalmanfilter\ik_framework\instance_expressiveioc'));
% modelInstance = rlModelInstance_expressiveioc(subjectId);
% xmlFilePath = 'D:\aslab_gitlab\kalmanfilter\ik_framework\instance_expressiveioc\model\ioc_v4_upperbody.xml';
% matFilePath = 'D:\aslab\data_IK\PamelasData\2019_04_11_fullbody2\matEkfIk\Subject01_SingleArm60Time_Full_body_01_mocap_mocap_X00_floating_ekfId1_ekfIk.mat';
% trcFilePath = 'D:\aslab\data\PamelasData\Subject01\Full_body\Subject01_SingleArm60Time_Full_body.trc';

% addpath(genpath('D:\aslab_gitlab\kalmanfilter\ik_framework\instance_expressiveioc'));
% modelInstance = rlModelInstance_expressiveioc_rightArm(subjectId);
% xmlFilePath = 'D:\aslab_gitlab\kalmanfilter\ik_framework\instance_expressiveioc\model\ioc_v4_rightarm.xml';
% matFilePath = 'D:\aslab\data_IK\PamelasData\2019_04_11_rightarm3\matEkfIk\Subject01_SingleArm60Time_OnlyRightArm_01_mocap_mocap_X00_floating_ekfId1_ekfIk.mat';
% trcFilePath = 'D:\aslab\data\PamelasData\Subject01\Only_Right_arm\Subject01_SingleArm60Time_OnlyRightArm.trc';

% addpath(genpath('D:\aslab_gitlab\kalmanfilter\ik_framework\instance_iit'));
% modelInstance = rlModelInstance_iit(subjectId);
% xmlFilePath = 'D:\aslab_gitlab\kalmanfilter\ik_framework\instance_iit\model\iit_v10.xml';
% matFilePath = 'D:\aslab\data_IK\FullBody_IIT_2017\2019_04_01_arm01\matEkfIk\PICK_STD_RES_Subject15_mocap_mocap_X00_floating_ekfId1_ekfIk.mat';
% trcFilePath = 'D:\aslab\data\Fullbody_IIT_2017\Subject15\mocap_fp\exercise5.trc';

addpath(genpath('D:\aslab_gitlab\kalmanfilter\ik_framework\instance_iit'));
modelInstance = rlModelInstance_iit(subjectId);
xmlFilePath = 'D:\aslab_gitlab\kalmanfilter\ik_framework\instance_iit\model\iit_v10_fixedbase.xml';
matFilePath = 'D:\aslab\data_IK\FullBody_IIT_2017\2019_04_10_squat02\matEkfIk\SQUA_STD_FAT_Subject01_mocap_mocap_X00_floating_ekfId1_ekfIk.mat';
trcFilePath = 'D:\aslab\data\Fullbody_IIT_2017\Subject01\mocap_fp\exercise7.trc';

% modelInstance = rlModelInstance_expressiveioc_rightArm(subjectId);
% xmlFilePath = strcat(xmlPath, 'ioc_v4_rightarm.xml');
% %xmlFilePath = strcat(xmlPath, 'ioc_v4_upperbody.xml');
% matFilePath = strcat(matPath, 'Subject01_SingleArm60Time_OnlyRightArm_01_mocap_mocap_X00_floating_ekfId1_ekfIk.mat');
% trcFilePath = strcat(trcPath, 'Subject01/RightArm/Subject01_SingleArm60Time_OnlyRightArm.trc');

load(matFilePath);

startTime = 0;
time = saveVar.time;
% frameData = saveVar.frameData;
q = saveVar.jointAngle.array;
jointNameMat = saveVar.jointLabels';

% [startVal, startInd] = findClosestValue(startTime, time);

% % % % % Example: Visualize joint centers without base translation
% % % % m_names = frameData.name;
% % % % 
% % % % vis = rlVisualizer('vis',640,480);
% % % % vis.update();
% % % % 
% % % % for i=1:2000
% % % %     frameBase = find(ismember(m_names, 'body_base'));
% % % %     tBase = reshape(frameData(frameBase).position(i, :), 4, 4);
% % % %     for j=2:numel(m_names)
% % % %         frameInd = find(ismember(m_names, m_names{j}));
% % % %         t = reshape(frameData(frameInd).position(i, :), 4, 4);
% % % %         pos = t(1:3,4) - tBase(1:3,4);
% % % %         vis.addMarker(m_names{j}, pos + [0 0 1]');
% % % %         pause(0.001);
% % % %     end
% % % %     vis.update
% % % % end

% % % Example: Visualize markers to make sure initial rotatio is applied
% % % vis = rlVisualizer('vis',640,480);
% % % vis.addModel(mdl);
% % % mdl.position = q(1, :);
% % % mdl.forwardPosition();
% % % vis.update();
% % 
% % % trc = loadDataFromTrc_expressiveioc(trcFilePath, algorithmParam);
% % % 
% % % m_names = fieldnames(trc.data);
% % 
% % % %% Visualize only markers
% % % for i=1:1000
% % %    for j=1:numel(m_names)
% % %        vis.addMarker(m_names{j},trc.data.(m_names{j})(i,:));
% % %    end
% % %    vis.update
% % % end

% % % % %% example: regenerate the transform matrix
% % % % frameStr = 'frame_rhand_end';
% % % % frameInd = find(ismember({frameData.name}, frameStr));
% % % % 
% % % % for i = 1:length(time)
% % % %     t = reshape(frameData(frameInd).position(i, :), 4, 4);
% % % % end
% % % % 

% %% example: load the model and trc data
modelInstance.loadModelFromModelSpecsNoSensor(xmlFilePath, matFilePath);
mdl = modelInstance.model;
jointNameMdl = {mdl.joints.name}';

mdl.forwardPosition();
mdl.base = 'frame_6dof_root';

algorithmParam.ekfRun = [];
algorithmParam.addUnknownMarkersToEkf = 1;
algorithmParam.missingMarkerValue = 0;
dataInstance = rlDataInstance_trc(modelInstance);
dataInstance.loadData(trcFilePath, algorithmParam);
dataInstance.dataProcessing(algorithmParam);

% visualize the playback
vis = rlVisualizer('vis',640,480);
vis.addModel(mdl);
vis.update();

for i = startInd:length(time)
    fprintf('Index %u/%u, Time %f/%f\n', i, length(time), time(i), max(time));
    mdl.position = q(i, :);
    mdl.forwardPosition();
    
    applyMarkersToVisualization(vis, mdl, dataInstance, i);
    
    vis.update();
end