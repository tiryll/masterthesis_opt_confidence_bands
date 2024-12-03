
function results = fit_var_vec(X, p, d)

K = size(X, 2);
Y = X((p+1):end, :,:);
T = size(Y, 1);
mc_runs=size(X,3);
Z = zeros(T, p*K,size(X,3));

for ii = 1:p
    Z(:, ((ii-1)*K+1):(ii*K),:) = X((p+1-ii):(end-ii), :,:);
end

if d == 1
    Z = cat(2,ones(T, 1,mc_runs),Z);
end

Y = permute(Y,[2,1,3]);
Z = permute(Z,[2,1,3]);

A = pagemtimes(pagemtimes(Y,'none',Z,'transpose'),pageinv(pagemtimes(Z,'none',Z,'transpose')));
U = Y - pagemtimes(A,Z);

Sigma_U_OLS = (1/(T-p*K-d)).*pagemtimes(U,'none',U,'transpose');
Sigma_U_MLE = ((T-K*p-d)/T) .* Sigma_U_OLS;
Gamma = (1/T) .* pagemtimes(Z,'none',Z,'transpose');

yhat = permute(pagemtimes(A,Z),[2,1,3]);

results.A = A;
results.T = T;
results.p = p;
results.d = d;
results.Sigma_U_OLS = Sigma_U_OLS;
results.Sigma_U_MLE = Sigma_U_MLE;
results.Gamma = Gamma;
results.yhat = yhat;
results.U = U;

end
