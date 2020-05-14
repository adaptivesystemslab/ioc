% knots:   1  21  41  61  81 101
%        101 121 141 161 181 201
% -> 6 knots in 101 points, so 5 knots in settings
% double length, so 10 knots in settings

variableFactors.win_length_sim = 201;
variableFactors.doc_sim_win_length_sim = variableFactors.win_length_sim/100; % in [s]
variableFactors.spline_length_sim = variableFactors.win_length_sim;
variableFactors.knots_sim = 10;