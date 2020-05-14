function plotSupportVectors(input, output, SVMmodel)
% pull out all the points
p1 = input(output == 1, :);
p0 = input(output == 0, :);

% now pull out all the svs
svIndices = SVMmodel.sv_indices;
svInput = input(svIndices, :);
svOutput = output(svIndices);
svp1 = svInput(svOutput == 1, :);
svp0 = svInput(svOutput == 0, :);

figure;
hold on
%             plot(p1(:, 1), p1(:, 2), 'b.');
%             plot(p0(:, 1), p0(:, 2), 'r.');
%
%             plot(svp1(:, 1), svp1(:, 2), 'bo');
%             plot(svp0(:, 1), svp0(:, 2), 'ro');

skipind = 8;
p1indtoplot = 1:skipind:size(p1, 1);
p0indtoplot = 1:skipind:size(p0, 1);
svp1indtoplot = 1:skipind:size(svp1, 1);
svp0indtoplot = 1:skipind:size(svp0, 1);

%             scatter3(p1(p1indtoplot, 1), p1(p1indtoplot, 2), p1(p1indtoplot, 3), ...
%                 'CData', [0 0 1], 'Marker', '.', 'DisplayName', 'p1 training');
%
%             scatter3(p0(p0indtoplot, 1), p0(p0indtoplot, 2), p0(p0indtoplot, 3), ...
%                 'CData', [1 0 0], 'Marker', '.', 'DisplayName', 'p0 training');

scatter3(svp1(svp1indtoplot, 1), svp1(svp1indtoplot, 2), svp1(svp1indtoplot, 3), ...
    'CData', [0 0 1], 'Marker', 'o', 'DisplayName', 'p1 sv');

scatter3(svp0(svp0indtoplot, 1), svp0(svp0indtoplot, 2), svp0(svp0indtoplot, 3), ...
    'CData', [1 0 0], 'Marker', 'x', 'DisplayName', 'p0 sv');

xlabel('X');
ylabel('Y');
zlabel('Z');

view([-250 10]);

grid on
end