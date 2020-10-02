% make correlation
function lala

xyzData;

colourList = {'b', 'c', 'g', 'm', 'r', 'k'};
symsList = {'.', 'x', 'o', '^', 'v', '+'};
xName = {'RmseMean', 'RmseStd', 'ResnormMean', 'ResnormStd', 'SegLenMean', 'SegLenStd'}; 

runningTableCounter = 1;
currTable{1} = 'ManVar';
currTable{2} = 'Basis';
currTable{3} = 'Exercise';
currTable{4} = 'R';
currTable{5} = 'ItemsConsidered';
runningTable = currTable;

uniqueDep = unique(dep);
uniqueDep{end+1} = 'all';
uniqueMotion = unique(motion);
uniqueMotion{end+1} = 'all';

% h = figure;

for k = length(uniqueMotion):-1:1
for i = length(uniqueDep):-1:1
    currX = [];
    currY = [];
    
    currDep = uniqueDep{i};
    currMotion = uniqueMotion{k};
    
    for j = 1:length(y)
        if strcmpi(currDep, dep{j}) || strcmpi(currDep, 'all')
            
        else
            continue
        end
        
        if strcmpi(currMotion, motion{j}) || strcmpi(currMotion, 'all')
            
        else
            continue
        end
        
        currX = [currX; x(j, :)];
        currY = [currY; y(j)];
    end
    
    if length(currY) < 5
        continue
    end
    
    for j = 1:size(currX, 2)
        [rmse_b, rmse_y, rmse_Rsq] = calcRegression(currX(:, j), currY);
%       h = generateAndPlot3(h, rmse_b, rmse_y, rmse_Rsq, [xName{j} '-' currDep], colourList{i}, symsList{j});
        
        currTable{1} = xName{j};
        currTable{2} = currDep;
        currTable{3} = currMotion;
        currTable{4} = rmse_Rsq;
        currTable{5} = length(currY);        
        
        runningTableCounter = runningTableCounter + 1;
        runningTable(runningTableCounter, :) = currTable;
    end
end
end

title('squat h1');
ylabel('Seg Acc');
legend show

runningTable
end

function h = generateAndPlot3(h, rmse_b, rmse_y, rmse_Rsq, name, color, sym)
%     [rmse_b, rmse_y, rmse_Rsq] = calcRegression(x, y);

    if ~isempty(h)
        figure(h);
        hold on
        plot(x, y, [color sym], 'LineWidth', 5, 'DisplayName', ['R=' num2str(rmse_Rsq) ', ' name ]);
        plot(x, rmse_y, color, 'DisplayName', ' ');
    end
end