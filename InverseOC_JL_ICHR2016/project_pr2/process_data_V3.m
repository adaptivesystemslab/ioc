function data = process_data_V3(bagFile)

if nargin == 0
    bagFile = ('C:\Users\Maram\Data\button5.bag');
end
[pathstr,name,ext] = fileparts(bagFile);
% make sure it's a bag file
if ~strcmp(ext, '.bag') 
warning('Is not a bag file: %s\n', bagFile)
return
end
% define mat file path
matf = fullfile(pathstr, [name '.mat']);

bag = rosbag(bagFile);
% extract three topics
% EE_topic = select(bag, 'Topic','/joint_states/link_pos/r_gripper_tool_frame');
% EE_data = readMessages(EE_topic, 'DataFormat', 'struct');
% 
% manip_topic = select(bag, 'Topic','/joint_states/manipulability_measures');
% manip_data = readMessages(manip_topic, 'DataFormat', 'struct');

joint_topic = select(bag, 'Topic','/joint_states');
joint_states = readMessages(joint_topic, 'DataFormat', 'struct');

%% processing the time stamps

% define the start and end time 
% the time at which the robot was in move to adjust all the data to include
% only this frame 

MOTION_THRESHOLD = 1e-5;
num_index = length(joint_states);
for i = 1:num_index
    js_stamp(i) = joint_states{i}.Header.Stamp;
end

init_time = double(js_stamp(1).Sec) + double(js_stamp(1).Nsec)/1e9;
last_time = double(js_stamp(end).Sec) + double(js_stamp(end).Nsec)/1e9;

start_pos = joint_states{1,1}.Position;

% once there is a motion, the time for this specific position is saved to
% t_start_offset
for msg = 1:num_index
   if(norm(start_pos - joint_states{msg,1}.Position) > MOTION_THRESHOLD)
       init_time_real  = double(joint_states{msg,1}.Header.Stamp.Sec)+ ...
           double(joint_states{msg,1}.Header.Stamp.Nsec)/1.0e9;
       init_idx_real = msg;
       break;
   end
end

end_pos = joint_states{end,1}.Position;
for msg = num_index:-1:1
   if(norm(joint_states{msg,1}.Position - end_pos) > MOTION_THRESHOLD)
       time_x = double(joint_states{msg,1}.Header.Stamp.Sec)+ ...
           double(joint_states{msg,1}.Header.Stamp.Nsec)/1.0e9;
       last_time_real = time_x;
       last_idx_real = msg;
       break;
   end
end


idx = 1;
for i = init_idx_real:last_idx_real
    TIME(idx) = double(js_stamp(i).Sec) + double(js_stamp(i).Nsec)/1.0e9 - init_time_real;
    idx = idx +1;
end

%% process joint positions
num_index_real = length(TIME);
POS_FULL = zeros(45,num_index_real);

idx = 1;
for i = init_idx_real:last_idx_real
    POS_FULL(:,idx) = joint_states{i}.Position;
    idx = idx +1;
end
% extract the position for the right arm joints only
POSITION = zeros(7,num_index_real);
POSITION(1,:) = POS_FULL(19,:);
POSITION(2,:) = POS_FULL(20,:);
POSITION(3,:) = POS_FULL(18,:);
POSITION(4,:) = POS_FULL(22,:);
POSITION(5,:) = POS_FULL(21,:);
POSITION(6,:) = POS_FULL(23,:);
POSITION(7,:) = POS_FULL(24,:);

% POSITION = filter_bw_lpf(POSITION');
% POSITION = POSITION';
%% get velocities from joint positions

VELOCITY = zeros(7,num_index_real - 1);

%calculate velocity using finite difference between points t and t+1
for i = 1:7
    for j = 1:num_index_real - 1
        VELOCITY(i,j) = ((POSITION(i,j+1) - POSITION(i,j))) ./ ((TIME(j+1) - TIME(j)));
    end
end
%% get acceleration from joint velocities

ACCELERATION = zeros(7,num_index_real - 2);

for i = 1:7
    for j = 1:num_index_real - 2
        ACCELERATION(i,j) = ((VELOCITY(i,j+1) - VELOCITY(i,j))) ./ ((TIME(j+1) - TIME(j)));
    end
end

% make POSITION and VELOCITY the same length as ACCELERATION, 7*num_index-2
POSITION(:,num_index_real) = [];
POSITION(:,num_index_real - 1) = [];
VELOCITY(:,num_index_real - 1) = [];

data.t = TIME;
data.q = POSITION;
data.dq = VELOCITY;
data.ddq = ACCELERATION;

save(matf, 'TIME','POSITION', 'VELOCITY', 'ACCELERATION')
% save(matf, 'data')

end