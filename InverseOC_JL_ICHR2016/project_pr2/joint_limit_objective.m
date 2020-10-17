function cost = joint_limit_objective(mdl, q, dq, ddq, index1)
% Purpose: compute the cost for the joint_limit_objective
% Input: a time index
% Prereq: run process_data.m so that the following are loaded
% - POSITION - joint positions of the trajectory

% These are the PR2 joint limits
shoulder_pan_LOWER = -2.28539816;
shoulder_lift_LOWER = -0.5236;
upper_arm_roll_LOWER = -3.9;
elbow_flex_LOWER = -2.3213;
wrist_flex_LOWER = -2.18;

shoulder_pan_UPPER = 0.71460184;
shoulder_lift_UPPER = 1.3963;
upper_arm_roll_UPPER = 0.8;
elbow_flex_UPPER = 0;
wrist_flex_UPPER = 0;

ang = q(:, index1);

% Store the positions for each joint that has limits
shoulder_pan = ang(1);
shoulder_lift = ang(2);
upper_arm_roll = ang(3);
elbow_flex = ang(4);
wrist_flex = ang(6);

% Compute the cost for the joint_limit_objective
cost =        4.0 * (shoulder_pan - shoulder_pan_LOWER)     * (shoulder_pan_UPPER - shoulder_pan)     / (shoulder_pan_UPPER - shoulder_pan_LOWER)^2;
cost = cost * 4.0 * (shoulder_lift - shoulder_lift_LOWER)   * (shoulder_lift_UPPER - shoulder_lift)   / (shoulder_lift_UPPER - shoulder_lift_LOWER)^2;
cost = cost * 4.0 * (upper_arm_roll - upper_arm_roll_LOWER) * (upper_arm_roll_UPPER - upper_arm_roll) / (upper_arm_roll_UPPER - upper_arm_roll_LOWER)^2;
cost = cost * 4.0 * (elbow_flex - elbow_flex_LOWER)         * (elbow_flex_UPPER - elbow_flex)         / (elbow_flex_UPPER - elbow_flex_LOWER)^2;
cost = cost * 4.0 * (wrist_flex - wrist_flex_LOWER)         * (wrist_flex_UPPER - wrist_flex)         / (wrist_flex_UPPER - wrist_flex_LOWER)^2;

% old formula
% prod = 4.0 / (shoulder_pan_UPPER - shoulder_pan_LOWER) / (shoulder_pan_UPPER - shoulder_pan_LOWER) * (shoulder_pan - shoulder_pan_LOWER) * (shoulder_pan_UPPER - shoulder_pan);
% prod = prod * 4.0 / (shoulder_lift_UPPER - shoulder_lift_LOWER) / (shoulder_lift_UPPER - shoulder_lift_LOWER) * (shoulder_lift - shoulder_lift_LOWER) * (shoulder_lift_UPPER - shoulder_lift);
% prod = prod * 4.0 / (upper_arm_roll_UPPER - upper_arm_roll_LOWER) / (upper_arm_roll_UPPER - upper_arm_roll_LOWER) * (upper_arm_roll - upper_arm_roll_LOWER) * (upper_arm_roll_UPPER - upper_arm_roll);
% prod = prod * 4.0 / (elbow_flex_UPPER - elbow_flex_LOWER) / (elbow_flex_UPPER - elbow_flex_LOWER) * (elbow_flex - elbow_flex_LOWER) * (elbow_flex_UPPER - elbow_flex);
% prod = prod * 4.0 / (wrist_flex_UPPER - wrist_flex_LOWER) / (wrist_flex_UPPER - wrist_flex_LOWER) * (wrist_flex - wrist_flex_LOWER) * (wrist_flex_UPPER - wrist_flex);
% cost = (1 - prod);