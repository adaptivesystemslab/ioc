function gradControl=trajGradient(currControlTraj,initialState, iocObj,weights)

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

% compute the recovery matrix
currStateTraj=[initialState;currStateTraj];
minWinLen=3;
H1=[];
H2=[];
F=[];
for t=1:timeHorizon
    if t<minWinLen % this is used to solve `no previous data' problem in computing jerk information
        currU=currControlTraj(t:t+minWinLen-1,:);
        currU=flipud(currU);
        currX=currStateTraj(t+1:t+minWinLen,:);
        currX=flipud(currX);
        nextU=currControlTraj(t+1:t+minWinLen,:);
        nextU=flipud(nextU);
        prevX=currStateTraj(t:t+minWinLen-1,:);
        prevX=flipud(prevX);
    elseif t==timeHorizon
        currU=currControlTraj(t-minWinLen+1:t,:);
        currX=currStateTraj(t+2-minWinLen:t+1,:);
        nextU=currControlTraj(t+1-minWinLen:t,:); % actually, this can be set anything
        prevX=currStateTraj(t-minWinLen+1:t,:);
    else
        currU=currControlTraj(t-minWinLen+1:t,:);
        currX=currStateTraj(t+2-minWinLen:t+1,:);
        nextU=currControlTraj(t+2-minWinLen:t+1,:);
        prevX=currStateTraj(t-minWinLen+1:t,:);
    end

    % Compute the deriviatives(jacobilaoge an matrix) for the system
    Jx=iocObj.getDynamicsStateDerivatives(currX,nextU);
    Ju=iocObj.getDynamicsControlDerivatives(prevX,currU);
    % Compute the deriviatives(jacobina matrix) for the features
    Px=iocObj.getFeaturesStateDerivatives(currX,currU);
    Pu=iocObj.getFeaturesControlDerivatives(currX,currU);
    
    % Compute the recovery matrix 
    Jx=Jx';
    Ju=Ju';
    Px=Px';
    Pu=Pu';
    
    % assemble the recovery matrix
    [H1,H2]=assembleH1H2(Jx, Ju, Px, Pu, H1,H2);
    
    
%     % assemble another matrix to prevent the final point changing....
%     prevJx=iocObj.getDynamicsStateDerivatives(prevX,currU);
%     F=assembleF(prevJx,Ju,F);
    
end

% Discard H2 because it contain nothing in DOC
G=H1*weights;
% invK=-inv(F'*F);
% VecG=G+F*invK*F'*G;+
VecG=G;
gradControl=reshape(VecG,size(currControlTraj,2),timeHorizon);
gradControl=gradControl';
end


function F=assembleF(prevJx,Ju,F)
if isempty(F)
    F=Ju;
else
    F=[F*prevJx;Ju];
end
end
