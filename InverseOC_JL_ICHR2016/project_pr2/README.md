# IOC RL

Load the Robotics Library (RL) model, process ROSBAG trajectories data, compute cost for various objectives

`load_model.m` loads the RL model and runs `process_data.m`.

`process_data.m` is for processing the [task and joint trajectories](https://drive.google.com/drive/u/0/folders/1GMlAUqx0Xcv6XFsNGpDpAw9AUsJXiUOW)

`Process_data_V3.m` read the joints data from bag file and cut the parts that the robot was steady in the beginning and the end. The Time, joint position, velocity, and acceleration are saved to mat file in the same directory as the bag file with the same name.

`robot.rlmdl.xml` is the RL model of the PR2. It is the updated version that contains the finger links. [Here](https://github.com/slin14/PR2-DH) is the full RL model. Note: the `robot.rlmdl.xml` file there does not contain the finger links.

The cost functions are in based on the [cost functions that Jesse developed](https://gitlab.com/jesse.li2002/moveit/-/tree/skill_transfer/moveit_planners/ompl/ompl_interface/include/moveit/ompl_interface/detail). I just re-wrote them for RL in MATLAB.
