% function circlefitting
% clearvars
load('temp', 'trainingData');
x = trainingData.trainingDataDR;
label = trainingData.trainingLabel;


c = 0*ones(size(x)); % weither or not we feed in a full array for c does not change P
% c = x;
[c(:, 1), c(:, 2)] = mapUnitCircle(x(:, 1), x(:, 2));

P = 0.01*eye(size(x, 2), size(x, 2));
% xP = c
Paug = x\c;
% P = Paug;
P(1:2, :) = Paug(1:2, :);

% for i = 1:size(P, 2)
% %     for j = 1:size(P, 1)
%     j = i;
%     P(i,j) = 1/P(i,j);
% %     end
% end

xp = x*P;

figure;
 
scatter3(x(label == 0, 1), x(label == 0, 2), x(label == 0, 3), 'r.'); hold on
scatter3(x(label == 1, 1), x(label == 1, 2), x(label == 1, 3), 'b.');

scatter3(xp(label == 0, 1), xp(label == 0, 2), xp(label == 0, 3), 'm.'); hold on
scatter3(xp(label == 1, 1), xp(label == 1, 2), xp(label == 1, 3), 'g.');

scatter3(c(label == 0, 1), c(label == 0, 2), c(label == 0, 3), 'k.'); hold on
scatter3(c(label == 1, 1), c(label == 1, 2), c(label == 1, 3), 'c.');
grid on

% P0 = zeros(size(x, 2), 1);
% A = [];
% b = [];
% for i = 1:size(x, 2) % iterate through dim
%     Acol = zeros(size(x));
%     Acol(:, i) = x(:, i);
%     
%     switch i
%         case {1, 2}
%             bcol = c(:, i);
%             
%         otherwise
%             bcol = zeros(size(c, 1), 1);
%     end
%     
%     A = [A; Acol];
%     b = [b; bcol];
% end
% 
% P = fmincon(@myfun,P0,A,b);

% x0 = x(:, 1);
% y0 = x(:, 2);
% [x1, y1] = mapUnitCircle(x0, y0);

% figure; 
% plot(x0, y0, '.');
% hold on
% plot(x0(1:45), y0(1:45), 'r.');
% 
% plot(x1, y1, 'g.');
% plot(x1(1:45), y1(1:45), 'm.');
% end
% 
% function f = myfun(P)
% % we want to min this
% % Pmat = diag(P);
% % f = Pmat*x;
% f = sum(P);
% end