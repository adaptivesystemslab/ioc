% from: http://www.mathworks.com/matlabcentral/answers/36999-how-do-i-regression-fit-a-sinwave-to-a-dataset

%% Generate some data

X = 2* pi*rand(100,1);
X = sortrows(X);
Y = 9 + 7*sin(2*X + 4*pi) + randn(100,1);

scatter(X,Y)

%%  Generate a fit

% Note that we need to pass three sets of input arguments to NonLinearModel
% # The X and Y data
% # A string describing our model
% # Starting conditions for the optimization solvers

% Generate some good starting conditions for the solvers

scatter(X, Y)
hold on

B0 = mean(Y);  % Vertical shift
B1 = (max(Y) - min(Y))/2; % Amplitude
B2 = 2; % Phase (Number of peaks)
B3 = 0; % Phase shift (eyeball the Curve)

myFit = NonLinearModel.fit(X,Y, 'y ~ b0 + b1*sin(b2*x1 + b3)', [B0, B1, B2, B3])

% Note that all the coefficient estimates are very good except for b3 where
% any even integer is equally valid

%% look at the complete set of methods

methods(myFit)

%% Generate a plot

hold on
plot(X, myFit.Fitted)
hold off

%% Generate a fit using an alternative syntax

myFit2 = NonLinearModel.fit(X,Y, @(b,x)(b(1) + b(2)*sin(b(3)*x + b(4))), [B0, B1, B2, B3])  
