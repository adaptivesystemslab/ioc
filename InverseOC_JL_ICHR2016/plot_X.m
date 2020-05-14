function plotX(x, param)
    figure;
    for j = 1:5:size(x, 2)
%         clf
        i = 3;
        prevX = x((1:3) + (i-1)*3, j);
        for i= 4:6%1:param.NF
            currX = x((1:3) + (i-1)*3, j);
            hold on
            if mod(i, 2)
                scatter3(currX(1), currX(2), currX(3), 'o')
            else
                scatter3(currX(1), currX(2), currX(3), 'x')
            end
            plot3([prevX(1) currX(1)], [prevX(2) currX(2)], [prevX(3) currX(3)]);
            prevX = currX;
        end
        xlim([-0.5 0.1]+1);
        ylim([-0.1 1.6]);
        gg = 1;
    end
end