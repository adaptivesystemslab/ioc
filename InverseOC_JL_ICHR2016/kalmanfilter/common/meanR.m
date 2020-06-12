function t = meanR(tCell)
    t = zeros(size(tCell{1}));
    ind = 0;
    for i = 1:length(tCell)
        if isempty(tCell{i}) || ~isempty(find(isnan(tCell{i}), 1)) % if the cell is empty
            continue;
        end
        
        if isempty(t) % if tcell{1} is empty
            t = zeros(size(tCell{i}));
        end
        
        ind = ind + 1;
        t = t + tCell{i};
    end
    
    t = t / ind;
end