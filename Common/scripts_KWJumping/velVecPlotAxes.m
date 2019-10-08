function velVecPlotAxes(velVecStart,indivPlot)

% inset "origin axes" for x_dot and z_dot

if(indivPlot)
    plot3([velVecStart(1), velVecStart(1)+0.2], [velVecStart(2), velVecStart(2)], [velVecStart(3), velVecStart(3)],'k');
    plot3([velVecStart(1), velVecStart(1)], [velVecStart(2), velVecStart(2)], [velVecStart(3), velVecStart(3)+0.2],'k');
    text(velVecStart(1)+0.155, velVecStart(2), velVecStart(3)-0.035,'X','FontSize',12);
    text(velVecStart(1)-0.04, velVecStart(2), velVecStart(3)+0.17,'Z','FontSize',12);
    plot3(velVecStart(1)+0.17, velVecStart(2), velVecStart(3)-0.02,'k.');
    plot3(velVecStart(1)-0.025, velVecStart(2), velVecStart(3)+0.19,'k.');
    
else
    plot3([velVecStart(1), velVecStart(1)+0.2], [velVecStart(2), velVecStart(2)], [velVecStart(3), velVecStart(3)],'k');
    plot3([velVecStart(1), velVecStart(1)], [velVecStart(2), velVecStart(2)], [velVecStart(3), velVecStart(3)+0.2],'k');
    text(velVecStart(1)+0.15, velVecStart(2), velVecStart(3)-0.055,'X');
    text(velVecStart(1)-0.05, velVecStart(2), velVecStart(3)+0.15,'Z');
    plot3(velVecStart(1)+0.17, velVecStart(2), velVecStart(3)-0.02,'k.');
    plot3(velVecStart(1)-0.025, velVecStart(2), velVecStart(3)+0.19,'k.');
end

