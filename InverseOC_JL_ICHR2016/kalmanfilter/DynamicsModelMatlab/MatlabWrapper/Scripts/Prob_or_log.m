%Here we test if using probabilities or log(P) is ok for matching markers
%to sensors on model

mes_size = 2;
A = rand(mes_size,mes_size); % generate a random n x n matrix
% construct a symmetric matrix
A = 0.5*(A+A');
A = A + mes_size*eye(mes_size);

[mes_x,mes_y] = meshgrid(linspace(-5,5,100));
mes = [mes_x(:), mes_y(:)];
mes_mean = zeros(1,2);
mes_mean_all = repmat(mes_mean,size(mes,1),1);
A_all = repmat(A,1,1,size(mes,1));
tic;
p = mvnpdf(mes,mes_mean_all,A_all);
p_norm = p./(1/sqrt((2*pi)^size(A,1)*det(A)));
time_p = toc;
disp(['Probs time: ' num2str(time_p)]);
surf(mes_x,mes_y,reshape(p_norm,size(mes_x,1),size(mes_x,2)));

%Now lets look at the log 
n = size(A,1);
er = mes-repmat(mes_mean,size(mes,1),1);
%log_norm = n*log(2*pi)+log(det(A))+arrayfun(@(n) er(n,:)*(A\er(n,:)'), 1:size(er,1));
tic;
log_norm = arrayfun(@(n) er(n,:)*(A\er(n,:)'), 1:size(er,1));
time_p = toc;
disp(['Log time: ' num2str(time_p)]);
surf(mes_x,mes_y,reshape(-1/2*log_norm,size(mes_x,1),size(mes_x,2)));


