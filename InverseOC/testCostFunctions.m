% test cost functions
clear all;

setPaths();

dt = 0.01;
y = [0 0 0 0;
    60 -60 30 -600;
    0 0 0 0];
x = [1
    20
     40];

configFile = jsondecode(fileread('../Data/IOC_testCostFunction3.json'));
trialInfo = configFile.Files(1);

% cfs = featuresEnums.parseEnum(trialInfo);
rlModel = getModel(trialInfo.model, trialInfo.modelType);

dy = zeros(size(y));
ddy = zeros(size(y));

% generate quintic spline
[q, dq, ddq, t] = quinticSpline(deg2rad(y), dy, ddy, x, dt);
dq(x, :) = 0;
ddq(x, :) = 0;

% now use RL to generate tau
tau = zeros(size(q));

for i = 1:size(q, 1)
    currQ = q(i, :);
    currDq = dq(i, :);
    currDdq = ddq(i, :);
    
    rlModel.updateState(currQ(1:end), currDq(1:end));
    tau(i, :) = rlModel.inverseDynamics(currDdq);
    ddq2(i, :) = rlModel.forwardDynamics(tau(i, :));
end

% calculate cost functions
ioc = IOCInstanceNew(rlModel, dt);
ioc.init(trialInfo);

state = encodeState(q, dq);
control = tau;
% ioc.setFeatureNormalization(state, control);
f = ioc.calcFeatures(state, control);
p = ioc.calcDynamics(state, control);

% calculate the gradient for the full length
[df_dx, df_du, dp_dx, dp_du] = ioc.getDerivativesNewObservation(state, control);

if 0 
    figure;
    ax(1) = subplot(221); plot(q); title('q');
    ax(2) = subplot(222); plot(dq); title('dq');
    ax(3) = subplot(223); plot(ddq); title('ddq');
    ax(4) = subplot(224); plot(ddq); title('ddq2');
%     ax(4) = subplot(224); plot(tau); title('tau');
    linkaxes(ax, 'x');
end
