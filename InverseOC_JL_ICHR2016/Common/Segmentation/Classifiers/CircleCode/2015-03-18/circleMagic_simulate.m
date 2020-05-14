function circleMagic
%     % generate simulatation data
    wx = [3 0.5];
    Ax = [1 11010];
    phix = [0 1];
    betax = [-pi/2 0];
    X = 0:pi/8:4*pi/2;
    Y = zeros(size(X));
    
    for ind = 1:length(wx)
        Y = Y + wx(ind) * sin(Ax(ind) * X + phix(ind)) + betax(ind) + rand(1, length(X))/3;
    end   
%     plot(Y, 'b');
    
%     % now fit the data

% return



    
%% Generate some data

% X = 2* pi*rand(100,1);
% X = sortrows(X);
% Y = 9 + 7*sin(2*X + 4*pi) + randn(100,1);

scatter(X,Y)

%%  Generate a fit

% Note that we need to pass three sets of input arguments to NonLinearModel
% # The X and Y data
% # A string describing our model
% # Starting conditions for the optimization solvers

% Generate some good starting conditions for the solvers

scatter(X, Y)
hold on

B01 = mean(Y);  % Vertical shift
B11 = (max(Y) - min(Y))/2; % Amplitude
B21 = 1; % Phase (Number of peaks)
B31 = 0; % Phase shift (eyeball the Curve)

B02 = mean(Y);  % Vertical shift
B12 = (max(Y) - min(Y))/2; % Amplitude
B22 = 1; % Phase (Number of peaks)
B32 = 0; % Phase shift (eyeball the Curve)

myFit = NonLinearModel.fit(X,Y, 'y ~ b01 + b11*sin(b21*x1 + b31) + b02 + b12*sin(b22*x1 + b32)', [B01, B11, B21, B31, B02, B12, B22, B32])

% Note that all the coefficient estimates are very good except for b3 where
% any even integer is equally valid

%% look at the complete set of methods

methods(myFit)

%% Generate a plot

hold on
plot(X, myFit.Fitted)
hold off

%% Generate a fit using an alternative syntax

% myFit2 = NonLinearModel.fit(X,Y, @(b,x)(b(1) + b(2)*sin(b(3)*x + b(4))), [B01, B11, B21, B31])  

end