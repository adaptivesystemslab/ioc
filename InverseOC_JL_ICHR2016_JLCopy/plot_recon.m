figure;

plot(q_opt'); hold on
plot(q_recon{1}', 'x', 'DisplayName', num2str(rmse(1)));
plot(q_recon{2}', 'o', 'DisplayName', num2str(rmse(2)));
plot(q_recon{3}', '+', 'DisplayName', num2str(rmse(3)));