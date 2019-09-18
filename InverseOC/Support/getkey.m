function rsp=getkey(t2wait)


% improve portability of your code acorss operating systems 
KbName('UnifyKeyNames');
% specify key names of interest in the study
activeKeys = [KbName('LeftArrow') KbName('RightArrow') KbName('s')];
% set value for maximum time to wait for response (in seconds)
% t2wait = 2; 
% if the wait for presses is in a loop, 
% then the following two commands should come before the loop starts
% restrict the keys for keyboard input to the keys we want
RestrictKeysForKbCheck(activeKeys);
% suppress echo to the command line for keypresses
ListenChar(2);
% get the time stamp at the start of waiting for key input 
% so we can evaluate timeout and reaction time
% tStart can also be the timestamp for the onset of the stimuli, 
% for example the VBLTimestamp returned by the 'Flip'
tStart = GetSecs;
% repeat until a valid key is pressed or we time out
timedout = false;
% initialise fields for rsp variable 
% that would contain details about the response if given
rsp.RT = NaN; rsp.keyCode = []; rsp.keyName = [];
while ~timedout
    % check if a key is pressed
    % only keys specified in activeKeys are considered valid
    [ keyIsDown, keyTime, keyCode ] = KbCheck; 
      if(keyIsDown), break; end
      if( (keyTime - tStart) > t2wait), timedout = true; end
end
% store code for key pressed and reaction time
if(~timedout)
  rsp.RT      = keyTime - tStart;
  rsp.keyCode = keyCode;
  rsp.keyName = KbName(rsp.keyCode);
end
% if the wait for presses is in a loop, 
% then the following two commands should come after the loop finishes
% reset the keyboard input checking for all keys
RestrictKeysForKbCheck;
% re-enable echo to the command line for key presses
% if code crashes before reaching this point 
% CTRL-C will reenable keyboard input
ListenChar(1)