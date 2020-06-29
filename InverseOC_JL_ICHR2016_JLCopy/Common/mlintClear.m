% Sometimes MLint would just fail for no good reason. This will
% clear the MLintFailureFile to reset that            
fprintf('Clearing MLintFailureFiles \n');

tmpcurrdir = pwd;
cd(prefdir);
delete('MLintFailureFiles');
cd(tmpcurrdir);