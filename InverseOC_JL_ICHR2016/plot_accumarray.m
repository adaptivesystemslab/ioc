% testing

function [X, Y] = plot_accumarray(x_input, y_input)
%  data = [
%      8,10
%      8,11
%      8,11
%      9,10
%      9,11
%      9,11
%      9,11];

for i = 1:size(y_input, 2)
    x_use{i} = (1:size(y_input, 1))';
    y_use{i} = y_input(:, i);
end
    
x = vertcat(x_use{:});
y = vertcat(y_use{:});

 
ageValues = unique(y);          %# Vector of unique age values
barData = accumarray(x,y,[],@(x) {hist(x,1:length(ageValues))});
X = find(~cellfun('isempty',barData));  %# Find x values for bars
Y = vertcat(barData{:});                %# Matrix of y values for bars
hBar = bar(X,Y,'stacked');              %# Create a stacked histogram
set(hBar,{'FaceColor'},{'g';'r'});      %# Set bar colors
legend(cellstr(num2str(ageValues)),'Location','NorthWest');  %# Add a legend