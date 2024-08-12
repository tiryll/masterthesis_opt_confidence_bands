
function results = fit_var(X, p, d)

K = size(X, 2);
Y = X((p+1):end, :);
T = size(Y, 1);
Z = zeros(T, p*K);

for ii = 1:p
    Z(:, ((ii-1)*K+1):(ii*K)) = X((p+1-ii):(end-ii), :);
end

if d == 1
    Z = [ones(T, 1), Z];
end

Y = Y';
Z = Z';

A = Y*Z' / (Z * Z');
U = Y - A * Z;

Sigma_U_OLS = (1/(T-p*K-d)) * (U * U');
Sigma_U_MLE = ((T-K*p-d)/T) * Sigma_U_OLS;
Gamma = (1/T) * (Z * Z');

yhat = (A * Z)';
llike = -(T*K)/2 * (1 + log(2*pi)) - T/2 * log(det(Sigma_U_MLE));

results.A = A;
results.T = T;
results.p = p;
results.d = d;
results.Sigma_U_OLS = Sigma_U_OLS;
results.Sigma_U_MLE = Sigma_U_MLE;
results.Gamma = Gamma;
results.yhat = yhat;
results.llike = llike;
results.U = U;

end