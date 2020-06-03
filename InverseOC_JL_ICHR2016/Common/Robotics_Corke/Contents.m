% Robotics Toolbox.
% Version 8  December 2008
%
% What's new.
%   Readme      - New features and enhancements in this version.
%
% Homogeneous transformations
%   angvec2r    - angle/vector to rotation matrix 3x3
%   angvec2tr   - angle/vector to transform 4x4
%   eul2r       - Euler angle to rotation matrix 3x3
%   eul2tr      - Euler angle to transform 4x4
%   oa2r        - orientation and approach vector to rotation matrix 3x3
%   oa2tr       - orientation and approach vector to transform 4x4
%   r2t         - rotation submatrix to transform
%   rotx        - transform for rotation about X-axis 3x3
%   roty        - transform for rotation about Y-axis 3x3
%   rotz        - transform for rotation about Z-axis 3x3
%   rpy2r       - roll/pitch/yaw angles to rotation matrix 3x3
%   rpy2tr      - roll/pitch/yaw angles to transform 4x4
%   t2r         - transform to rotation submatrix
%   tr2angvec   - transform to angle/vector form
%   tr2eul      - transform to Euler angles 
%   tr2rpy      - transform to roll/pitch/yaw angles
%   transl      - set or extract the translational component of a transform 4x4
%   trnorm      - normalize a transform 
%   trplot      - plot a transform as a coordinate frame
%   trotx       - transform for rotation about X-axis 4x4
%   troty       - transform for rotation about Y-axis 4x4
%   trotz       - transform for rotation about Z-axis 4x4
%   
% Quaternion methods:
%   /           - divide quaternion by quaternion or scalar
%   *           - multiply quaternion by a quaternion or vector
%   inv         - invert a quaternion 
%   norm        - norm of a quaternion 
%   plot        - display a quaternion as a 3D rotation
%   unit        - unitize a quaternion 
%   qinterp     - interpolate a quaternion
%
% Kinematics
%   diff2tr     - differential motion vector to transform 
%   fkine       - compute forward kinematics 
%   ikine       - compute inverse kinematics 
%   ikine560    - compute inverse kinematics for Puma 560 like arm
%   jacob0      - compute Jacobian in base coordinate frame
%   jacobn      - compute Jacobian in end-effector coordinate frame
%   tr2diff     - transform to differential motion vector 
%   tr2jac      - transform to Jacobian 
%   
% Dynamics
%   accel       - compute forward dynamics
%   cinertia    - compute Cartesian manipulator inertia matrix 
%   coriolis    - compute centripetal/coriolis torque 
%   fdyn        - forward dynamics
%   ftrans      - transform force/moment 
%   gravload    - compute gravity loading 
%   inertia     - compute manipulator inertia matrix 
%   itorque     - compute inertia torque 
%   rne         - inverse dynamics 
%   
% Trajectory generation
%   ctraj       - Cartesian trajectory 
%   jtraj       - joint space trajectory 
%   trinterp    - interpolate transform s
%   
% Graphics
%   drivebot    - drive a graphical  robot 
%   plot        - plot/animate robot 
%   
% Other
%   ishomog     - true if argument is a 4x4 matrix
%   isrot       - true if argument is a 3x3 matrix
%   isvec       - true if argument is a 3-vector
%   maniplty    - compute manipulability 
%   unit        - unitize a vector
%
% Robot methods:
%   *           - compound two robots
%   friction    - return joint friction torques
%   nofriction  - return a robot object with no friction
%   perturb     - return a robot object with perturbed parameters
%   plot        - plot/animate a robot
%
% Creation of robot models.
%   Fanuc10L    - Fanuc 10L (DH, kine)
%   link        - construct a robot link object 
%   MotomanHP6  - Motoman HP6 (DH, kine)
%   puma560     - Puma 560 data (DH, kine, dyn)
%   puma560akb  - Puma 560 data (MDH, kine, dyn)
%   robot       - construct a robot object 
%   stanford    - Stanford arm data (DH, kine, dyn)
%   S4ABB2p8    - ABB S4 2.8 (DH, kine)
%   twolink     - simple 2-link example (DH, kine)
%
% Demonstrations.
%   demos      - toolbox demonstration
%   
% Copyright (C) 2008, by Peter I. Corke
