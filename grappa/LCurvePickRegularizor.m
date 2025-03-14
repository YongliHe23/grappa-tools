A=calib(src);
b=calib(trg);

s=norm(A*A'); %biggest singular value, ~8e9 

lambdas = logspace(-9, -1, 100).*s; % Range of lambda values
residuals = zeros(size(lambdas));
solutions = zeros(size(lambdas));

pinv_tikhonov=@(A,lambda) A'/(A*A'+lambda*eye(size(A,1)));

for i = 1:length(lambdas)
    i
    A_pinv = pinv_tikhonov(A, lambdas(i));
    w_lambda = b*A_pinv; % Regularized solution
    residuals(i) = norm(w_lambda*A - b);
    solutions(i) = norm(w_lambda);
end

loglog(residuals, solutions, '-o');
xlabel('Residual Norm ||wA - b||');
ylabel('Solution Norm ||w||');
title('L-Curve for Regularization Parameter Selection');
grid on;