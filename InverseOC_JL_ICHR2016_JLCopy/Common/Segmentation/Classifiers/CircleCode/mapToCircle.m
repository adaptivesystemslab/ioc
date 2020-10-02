function [xp, P, h] = mapToCircle(x, label, fig)

if ~exist('fig', 'var')
    fig = 0;
end

% assume it's zero mean already coming ing

% apply a zero mean formula first and see
% x(:, 1) = x(:, 1) - mean(x(:, 1));
% x(:, 2) = x(:, 2) - mean(x(:, 2));

% ellipsePartStruct = fit_ellipse(x(:, 1), x(:, 2));
% x(:, 1) = x(:, 1) - ellipsePartStruct.X0_in;
% x(:, 2) = x(:, 2) - ellipsePartStruct.Y0_in;

c = 0*ones(size(x)); % weither or not we feed in a full array for c does not change P
% c = x;
[c(:, 1), c(:, 2)] = mapUnitCircle(x(:, 1), x(:, 2));

Peye = 0.01*eye(size(x, 2), size(x, 2)); % the 0.01 is to shrink the magnitude of off-axis terms % xP = c
P = x\c; % perform least square
% P = Paug;
P(3:end, 3:end) = Peye(3:end, 3:end);

% for i = 1:size(P, 2)
% %     for j = 1:size(P, 1)
%     j = i;
%     P(i,j) = 1/P(i,j);
% %     end
% end

xp = x*P;

if fig
h = figure('position', [    1419         169        1765         907]);
 
scatter3(x(label == 0, 1), x(label == 0, 2), x(label == 0, 3), 'r.', 'DisplayName', 'orig non-seg'); hold on
scatter3(x(label == 1, 1), x(label == 1, 2), x(label == 1, 3), 'b.', 'DisplayName', 'orig seg');

scatter3(xp(label == 0, 1), xp(label == 0, 2), xp(label == 0, 3), 'm.', 'DisplayName', 'mapped non-seg'); hold on
scatter3(xp(label == 1, 1), xp(label == 1, 2), xp(label == 1, 3), 'g.', 'DisplayName', 'mapped seg');

scatter3(c(label == 0, 1), c(label == 0, 2), c(label == 0, 3), 'k.', 'DisplayName', 'circle non-seg'); hold on
scatter3(c(label == 1, 1), c(label == 1, 2), c(label == 1, 3), 'c.', 'DisplayName', 'circle seg');
grid on
view([-35 65]);
legend
else 
    h = [];
end
