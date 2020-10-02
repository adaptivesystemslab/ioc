function [h, r1_1, r1_2, r2_1, r2_2] = plotFigCompare2(h, rmse_fct, feature_optsim1, feature_optsim2, feature_optsim3, ...
    param_docsim1, param_docsim2, param_docsim3, secondWindowOffset, thirdWindowOffset)

len1 = length(feature_optsim1.q);
len2 = length(feature_optsim2.q);
len3 = length(feature_optsim3.q);

indSet1 = 1:len1;
indOffsetSecond = (1:len2) + secondWindowOffset;
indOffsetThird = (1:len3) + thirdWindowOffset;

colour2 = 'r';
colour3 = 'g';

r1_1 = rmse_fct(feature_optsim1.q(:, indOffsetSecond), feature_optsim2.q(:, 1:len2), len2);
r1_2 = rmse_fct(feature_optsim1.q(:, indOffsetThird), feature_optsim3.q(:, 1:len3), len3);
r2_1 = rmse_fct(feature_optsim1.dq(:, indOffsetSecond), feature_optsim2.dq(:, 1:len2), len2);
r2_2 = rmse_fct(feature_optsim1.dq(:, indOffsetThird), feature_optsim3.dq(:, 1:len3), len3);

subplot(121);
hold on
plot(indSet1, feature_optsim1.q', '-', 'DisplayName', 'Full');
for jj = 1:length(param_docsim1.const_x)
    plot(param_docsim1.const_x(jj), param_docsim1.const_y(:, jj), 'o', 'MarkerSize', 20);
end

plot(indOffsetSecond, feature_optsim2.q', [colour2 '.'], 'DisplayName', 'Partial1');
for jj = 1:length(param_docsim2.const_x)
    plot(param_docsim2.const_x(jj) + secondWindowOffset, param_docsim2.const_y(:, jj), [colour2 'x'], 'MarkerSize', 20);
end

plot(indOffsetThird, feature_optsim3.q', [colour3 '.'], 'DisplayName', 'Partial2');
for jj = 1:length(param_docsim3.const_x)
    plot(param_docsim3.const_x(jj) + thirdWindowOffset, param_docsim3.const_y(:, jj), [colour3 '+'], 'MarkerSize', 20);
end
title(['q RMSE: ' num2str(r1_1) ', ' num2str(r1_2)]);

subplot(122);
hold on
plot(indSet1, feature_optsim1.dq', '-');
for jj = 1:length(param_docsim1.const_dx)
    plot(param_docsim1.const_dx(jj), param_docsim1.const_dy(:, jj), 'o', 'MarkerSize', 20);
end

plot(indOffsetSecond, feature_optsim2.dq', [colour2 '.']);
for jj = 1:length(param_docsim2.const_dx)
    plot(param_docsim2.const_dx(jj) + secondWindowOffset, param_docsim2.const_dy(:, jj), [colour2 'x'], 'MarkerSize', 20);
end

plot(indOffsetThird, feature_optsim3.dq', [colour3 '.']);
for jj = 1:length(param_docsim3.const_dx)
    plot(param_docsim3.const_dx(jj) + thirdWindowOffset, param_docsim3.const_dy(:, jj), [colour3 '+'], 'MarkerSize', 20);
end
title(['dq RMSE: ' num2str(r2_1) ', ' num2str(r2_2)]);

% subplot(133);
% hold on
% plot(1:100, feature_optsim1.ddq', '-', 'DisplayName', 'Full1');
% for jj = 1:length(param_docsim1.const_ddx)
%     plot(param_docsim1.const_ddx(jj), param_docsim1.const_ddy(:, jj), 'x', 'MarkerSize', 20, 'DisplayName', 'Full');
% end
% plot(indOffsetSecond, feature_optsim2.ddq', [colour2 '.'], 'DisplayName', 'Half2');
% for jj = 1:length(param_docsim2.const_ddx)
%     plot(param_docsim2.const_ddx(jj) + secondWindowOffset, param_docsim2.const_ddy(:, jj), [colour2 'o'], 'MarkerSize', 20, 'DisplayName', 'Half');
% end
% plot(indOffsetThird, feature_optsim3.ddq', [colour3 '.'], 'DisplayName', 'Half3');
% for jj = 1:length(param_docsim3.const_ddx)
%     plot(param_docsim3.const_ddx(jj) + thirdWindowOffset, param_docsim3.const_ddy(:, jj), [colour3 'o'], 'MarkerSize', 20, 'DisplayName', 'Half');
% end
% r3_1 = rmse_fct(feature_optsim1.ddq(:, indOffsetSecond), feature_optsim2.ddq(:, 1:50), 50);
% r3_2 = rmse_fct(feature_optsim1.ddq(:, indOffsetThird), feature_optsim3.ddq(:, 1:50), 50);
% title(['ddq RMSE: ' num2str(r3_1) ', ' num2str(r3_2)]);
legend show
end