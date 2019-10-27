function [cost,currStateTraj]=costEval(currControlTraj,initialState,iocObj,weights)

if size(weights,2)>1; weights=weights'; end
if size(initialState,1)>1; initialState=initialState'; end

% Use current control trajectory to generate the state trajectory
timeHorizon=length(currControlTraj);
dimState=length(initialState);
currStateTraj=zeros(timeHorizon,dimState);
currState=initialState;
for t=1:timeHorizon
    currControl=currControlTraj(t,:);
    nextState=iocObj.calcDynamics(currState,currControl);
    currStateTraj(t,:)=nextState;
    currState=nextState;
end

% Evaluate the cost
featuresVal=zeros(timeHorizon,length(weights));
minWinLen=3;
for t=1:minWinLen-1
    currU=currControlTraj(t:t+minWinLen,:);
    currU=flipud(currU);
    currX=currStateTraj(t:t+minWinLen,:);
    currX=flipud(currX);
    cffull=iocObj.calcFeatures(currX, currU);
    featuresVal(t,:)=cffull(end, :);
end

currU=currControlTraj;
currX=currStateTraj;
cf2full=iocObj.calcFeatures(currX, currU);
featuresVal(minWinLen:end,:)=cf2full(minWinLen:end,:);


stageCost=featuresVal*weights;
cost=sum(stageCost);


end
