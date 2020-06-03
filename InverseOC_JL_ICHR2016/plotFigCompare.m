function [h, r1, r2] = plotFigCompare(h, rmse_fct, feature_optsim1, feature_optsim2, param_docsim1, param_docsim2, secondWindowOffset, colour1, colour2)
indOffset = (1:50) + secondWindowOffset;

subplot(131);
hold on
plot(1:100, feature_optsim1.q', '-', 'DisplayName', 'Full');
plot(indOffset, feature_optsim2.q', [colour2 '.'], 'DisplayName', 'Partial');
for jj = 1:length(param_docsim1.const_x)
    plot(param_docsim1.const_x(jj), param_docsim1.const_y(:, jj), 'x', 'MarkerSize', 20);
end
for jj = 1:length(param_docsim2.const_x)
    plot(param_docsim2.const_x(jj) + secondWindowOffset, param_docsim2.const_y(:, jj), [colour2 'o'], 'MarkerSize', 20);
end
r1 = rmse_fct(feature_optsim1.q(:, indOffset), feature_optsim2.q(:, 1:50), 50);
title(['q RMSE: ' num2str(r1)]);

subplot(132);
hold on
plot(1:100, feature_optsim1.dq', '-');
for jj = 1:length(param_docsim1.const_dx)
    plot(param_docsim1.const_dx(jj), param_docsim1.const_dy(:, jj), 'x', 'MarkerSize', 20);
end
plot(indOffset, feature_optsim2.dq', [colour2 '.']);
for jj = 1:length(param_docsim2.const_dx)
    plot(param_docsim2.const_dx(jj) + secondWindowOffset, param_docsim2.const_dy(:, jj), [colour2 'o'], 'MarkerSize', 20);
end
r2 = rmse_fct(feature_optsim1.dq(:, indOffset), feature_optsim2.dq(:, 1:50), 50);
title(['dq RMSE: ' num2str(r2)]);

subplot(133);
hold on
plot(1:100, feature_optsim1.ddq', '-', 'DisplayName', 'Full');
for jj = 1:length(param_docsim1.const_ddx)
    plot(param_docsim1.const_ddx(jj), param_docsim1.const_ddy(:, jj), 'x', 'MarkerSize', 20, 'DisplayName', 'Full');
end
plot(indOffset, feature_optsim2.ddq', [colour2 '.'], 'DisplayName', 'Half');
for jj = 1:length(param_docsim2.const_ddx)
    plot(param_docsim2.const_ddx(jj) + secondWindowOffset, param_docsim2.const_ddy(:, jj), [colour2 'o'], 'MarkerSize', 20, 'DisplayName', 'Half');
end
r3 = rmse_fct(feature_optsim1.ddq(:, indOffset), feature_optsim2.ddq(:, 1:50), 50);
title(['ddq RMSE: ' num2str(r3)]);
legend show
end