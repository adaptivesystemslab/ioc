function data = process_data(bagFile)
% if nargin == 0
    bagFile = 'D:\data\data_pr2\RRTstarLengthJoint_trial00_button0_.bag';
% end

bag = rosbag(bagFile); %change the name for different files
rosbag('info', bagFile);

bSel = select(bag, 'Topic', '/joint_states');
joint_states = readMessages(bSel, 'DataFormat','struct');
num_index = length(joint_states)

%% processing the time stamps
for i = 1:num_index
js_stamp(i) = joint_states{i}.Header.Stamp;
end

init_time = double(js_stamp(1).Sec) + double(js_stamp(1).Nsec)/1.0e-9;
TIME = zeros(1,num_index);

for i = 1:num_index
TIME(i) = double(js_stamp(i).Sec) + double(js_stamp(i).Nsec)/1.0e-9 - init_time;
end

%% for processing torqe if the bag file contains efforts

% EFFORT = zeros(45, num_index - 1);
% joint_states{1}.Effort
% for i = 1:num_index - 1
% EFFORT(:,i) = joint_states{i}.Effort;
% end
% 
% TORQUE = zeros(7,num_index - 1);
% TORQUE(1,:) = EFFORT(19,:);
% TORQUE(2,:) = EFFORT(20,:);
% TORQUE(3,:) = EFFORT(18,:);
% TORQUE(4,:) = EFFORT(22,:);
% TORQUE(5,:) = EFFORT(21,:);
% TORQUE(6,:) = EFFORT(23,:);
% TORQUE(7,:) = EFFORT(24,:);

%% process joint positions
POS_FULL = zeros(45,num_index);

for i = 1:num_index
POS_FULL(:,i) = joint_states{i}.Position;
end

POSITION = zeros(7,num_index);
POSITION(1,:) = POS_FULL(19,:);
POSITION(2,:) = POS_FULL(20,:);
POSITION(3,:) = POS_FULL(18,:);
POSITION(4,:) = POS_FULL(22,:);
POSITION(5,:) = POS_FULL(21,:);
POSITION(6,:) = POS_FULL(23,:);
POSITION(7,:) = POS_FULL(24,:);

%% for processing velocity if the bag file contains velocities
% VEL_FULL = zeros(45,num_index - 1);
% 
% for i = 1:num_index - 1
% VEL_FULL(:,i) = joint_states{i}.Velocity;
% end
% 
% VELOCITY = zeros(7,num_index);
% VELOCITY(1,:) = VEL_FULL(19,:);
% VELOCITY(2,:) = VEL_FULL(20,:);
% VELOCITY(3,:) = VEL_FULL(18,:);
% VELOCITY(4,:) = VEL_FULL(22,:);
% VELOCITY(5,:) = VEL_FULL(21,:);
% VELOCITY(6,:) = VEL_FULL(23,:);
% VELOCITY(7,:) = VEL_FULL(24,:);

%% get velocities from joint positions

VELOCITY = zeros(7,num_index - 1);

% calculate velocity using finite difference between points t and t+1
for i = 1:7
for j = 1:num_index - 1
VELOCITY(i,j) = ((POSITION(i,j+1) - POSITION(i,j))) ./ ((TIME(j+1) - TIME(j)));
end
end

%%
ACCELERATION = zeros(7,num_index - 2);

for i = 1:7
for j = 1:num_index - 2
ACCELERATION(i,j) = ((VELOCITY(i,j+1) - VELOCITY(i,j))) ./ ((TIME(j+1) - TIME(j)));
end
end

% make POSITION and VELOCITY the same length as ACCELERATION, 7*num_index-2
POSITION(:,num_index) = [];
POSITION(:,num_index - 1) = [];
VELOCITY(:,num_index - 1) = [];
% VELOCITY(:,num_index) = [];

% compute inverse
% TORQUE = TORQUE';

data.t = TIME;
data.q = POSITION';
data.dq = VELOCITY';
data.ddq = ACCELERATION';

% plot(TORQUE)

% note: Pam said something about adding a low pass filter before computing
% ddq from dq then another filter after, but I haven't done that (wasn't
% sure what low pass filter to use)
end