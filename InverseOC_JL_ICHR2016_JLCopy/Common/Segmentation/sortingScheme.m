function [nameMapping, nameInd] = sortingScheme(tName, pPoints)

    nameMapping = {};
    for ind_p1name = 1:length(pPoints)
        currNameTemp = tName(pPoints(ind_p1name));
        currNameTemp = unique(currNameTemp{1});
        for ind_inner = 1:length(currNameTemp)

            if ind_inner > 1
                blah = 1;
            end

            currName = currNameTemp{ind_inner};

            if sum(strcmpi(currName, nameMapping))
                % if this entry already exist
                indToProduce = find(strcmpi(currName, nameMapping));
                nameInd{indToProduce} = [nameInd{indToProduce} ind_p1name];
            else
                nameMapping = [nameMapping; currName];
                indToProduce = length(nameMapping);
                nameInd{indToProduce} = ind_p1name;
            end
        end
    end
end
