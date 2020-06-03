% knots:   1  11  21  31  41  51  61  71  81  91 101
%        101 111 121 131 141 151 161 171 181 191 201
% -> 11 knots in 101 points, so 10 knots in settings
% double length, so 20 knots in settings

variableFactors.win_length_sim = 201;
variableFactors.doc_sim_win_length_sim = variableFactors.win_length_sim/100; % in [s]
variableFactors.spline_length_sim = variableFactors.win_length_sim;
variableFactors.knots_sim = 20;

if 0
length_traj = variableFactors.win_length_sim;
skip_count = floor(length_traj/variableFactors.knots_sim);
intermed_ind = unique([1:skip_count:length_traj-1 length_traj]); % this ensures there is a knot point at the end, and on constraints
end